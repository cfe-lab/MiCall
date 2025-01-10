from typing import Iterable, Iterator, Optional, FrozenSet, Tuple, Sequence, TextIO, MutableMapping, Literal, Union
from dataclasses import dataclass
from fractions import Fraction
from Bio import SeqIO, Seq
from Bio.SeqRecord import SeqRecord
import logging
from mappy import Aligner
from functools import cached_property

from micall.utils.contig_stitcher_context import StitcherContext
from micall.utils.consensus_aligner import Alignment
from micall.utils.overlap_stitcher import align_queries, \
    calculate_concordance, sort_concordance_indexes, calc_overlap_pvalue
from micall.utils.contig_stitcher_contigs import Contig
from micall.utils.find_maximum_overlap import find_maximum_overlap, OverlapFinder


logger = logging.getLogger(__name__)

ACCEPTABLE_STITCHING_PROB = Fraction(1, 20)


@dataclass(frozen=True)
class ContigWithAligner(Contig):
    @cached_property
    def aligner(self) -> Aligner:
        return Aligner(seq=self.seq)

    @staticmethod
    def make(contig: Contig) -> 'ContigWithAligner':
        return ContigWithAligner(name=contig.name, seq=contig.seq)

    @staticmethod
    def empty() -> 'ContigWithAligner':
        return ContigWithAligner.make(Contig.empty())

    def map_overlap(self, overlap: str) -> Iterator[Alignment]:
        for x in self.aligner.map(overlap):
            if x.is_primary:
                yield x


@dataclass(frozen=True)
class ContigsPath:
    # Contig representing all combined contigs in the path.
    whole: ContigWithAligner

    # Id's of contigs that comprise this path.
    parts_ids: FrozenSet[int]

    # Lower is better. This is an estimated probability that
    # all the components in this path came together by accident.
    probability: Fraction
    pessimisstic_probability: Fraction

    def score(self) -> Tuple[Fraction, Fraction, int]:
        return (1-self.pessimisstic_probability, 1-self.probability, len(self.parts_ids))

    def has_contig(self, contig: Contig) -> bool:
        return contig.id in self.parts_ids

    @property
    def is_empty(self) -> bool:
        return len(self.parts_ids) == 0

    @staticmethod
    def empty() -> 'ContigsPath':
        return ContigsPath(ContigWithAligner.empty(), frozenset(),
                           probability=Fraction(1),
                           pessimisstic_probability=ACCEPTABLE_STITCHING_PROB)


@dataclass(frozen=True)
class Overlap:
    # A negative integer.
    # Represents the shift value for left contig relative to right
    # contig to achieve the maximum overlap.
    shift: int
    size: int


ContigId = int
GET_OVERLAP_CACHE: MutableMapping[
    Tuple[ContigId, ContigId],
    Optional[Overlap]] = {}


def get_overlap(finder: OverlapFinder, left: ContigWithAligner, right: ContigWithAligner) -> Optional[Overlap]:
    if len(left.seq) == 0 or len(right.seq) == 0:
        return None

    key = (left.id, right.id)
    cached: Union[Overlap, None, Literal[0]] = GET_OVERLAP_CACHE.get(key, 0)
    if cached != 0:
        return cached

    (shift, value) = find_maximum_overlap(left.seq, right.seq, finder=finder)
    if shift == 0:
        ret = None
        GET_OVERLAP_CACHE[key] = ret
        return ret

    size = min((abs(shift), len(left.seq), len(right.seq)))
    ret = Overlap(shift=shift, size=size)
    GET_OVERLAP_CACHE[key] = ret
    return ret


def combine_probability(current: Fraction, new: Fraction) -> Fraction:
    return current * new


TRY_COMBINE_CACHE: MutableMapping[
    Tuple[ContigId, ContigId],
    Optional[Tuple[ContigWithAligner, Fraction]]] = {}


def try_combine_contigs(finder: OverlapFinder,
                        current_prob: Fraction,
                        max_acceptable_prob: Fraction,
                        a: ContigWithAligner, b: ContigWithAligner,
                        ) -> Optional[Tuple[ContigWithAligner, Fraction]]:

    if len(b.seq) == 0:
        return (a, Fraction(1))

    if len(a.seq) == 0:
        return (b, Fraction(1))

    overlap = get_overlap(finder, a, b)
    if overlap is None:
        return None

    if abs(overlap.shift) < len(a.seq):
        shift = overlap.shift
        left = a
        right = b
    else:
        shift = (len(b.seq) + len(a.seq) - abs(overlap.shift)) * -1
        left = b
        right = a

    optimistic_number_of_matches = overlap.size
    optimistic_result_probability = calc_overlap_pvalue(L=overlap.size, M=optimistic_number_of_matches)
    if combine_probability(optimistic_result_probability, current_prob) > max_acceptable_prob:
        return None

    left_initial_overlap = left.seq[len(left.seq) - abs(shift):(len(left.seq) - abs(shift) + len(right.seq))]
    right_initial_overlap = right.seq[:abs(shift)]

    key = (left.id, right.id)
    existing: Union[Tuple[ContigWithAligner, Fraction], None, Literal[0]] \
        = TRY_COMBINE_CACHE.get(key, 0)

    if existing != 0:
        return existing

    elif len(left_initial_overlap) < len(right_initial_overlap):
        left_overlap_alignments = left.map_overlap(right_initial_overlap)
        left_cutoff = min((al.r_st for al in left_overlap_alignments), default=None)
        if left_cutoff is None:
            TRY_COMBINE_CACHE[key] = None
            return None

        right_overlap_alignments = right.map_overlap(left_initial_overlap)
        right_cutoff = max((al.r_en for al in right_overlap_alignments), default=None)
        if right_cutoff is None:
            TRY_COMBINE_CACHE[key] = None
            return None
    else:
        right_overlap_alignments = right.map_overlap(left_initial_overlap)
        right_cutoff = max((al.r_en for al in right_overlap_alignments), default=None)
        if right_cutoff is None:
            TRY_COMBINE_CACHE[key] = None
            return None

        left_overlap_alignments = left.map_overlap(right_initial_overlap)
        left_cutoff = min((al.r_st for al in left_overlap_alignments), default=None)
        if left_cutoff is None:
            TRY_COMBINE_CACHE[key] = None
            return None

    left_overlap = left.seq[left_cutoff:(left_cutoff + len(right.seq))]
    left_remainder = left.seq[:left_cutoff]
    right_overlap = right.seq[:right_cutoff]
    right_remainder = right.seq[right_cutoff:]
    aligned_left, aligned_right = align_queries(left_overlap, right_overlap)

    number_of_matches = sum(1 for x, y
                            in zip(aligned_left, aligned_right)
                            if x == y and x != '-')
    result_probability = calc_overlap_pvalue(L=len(left_overlap), M=number_of_matches)
    if combine_probability(current_prob, result_probability) > max_acceptable_prob:
        return None

    is_covered = len(right.seq) < abs(shift)
    if is_covered:
        logger.debug("Between %s and %s, returning the covering contig, %s.",
                     a.unique_name, b.unique_name, left.unique_name)
        TRY_COMBINE_CACHE[key] = (left, Fraction(1))
        return (left, Fraction(1))

    else:
        concordance = calculate_concordance(aligned_left, aligned_right)
        max_concordance_index = next(iter(sort_concordance_indexes(concordance)))
        left_overlap_chunk = ''.join(x for x in aligned_left[:max_concordance_index] if x != '-')
        right_overlap_chunk = ''.join(x for x in aligned_right[max_concordance_index:] if x != '-')

        result_seq = left_remainder + left_overlap_chunk + right_overlap_chunk + right_remainder
        result_contig = ContigWithAligner(None, result_seq)

        logger.debug("Joined %s and %s together in a contig %s with lengh %s.",
                     a.unique_name, b.unique_name,
                     result_contig.unique_name, len(result_contig.seq))

        TRY_COMBINE_CACHE[key] = (result_contig, result_probability)
        return (result_contig, result_probability)


def extend_by_1(finder: OverlapFinder,
                max_acceptable_prob: Fraction,
                path: ContigsPath,
                candidate: ContigWithAligner,
                ) -> Iterator[ContigsPath]:
    if path.has_contig(candidate):
        return

    combination = try_combine_contigs(finder, path.probability, max_acceptable_prob, path.whole, candidate)
    if combination is None:
        return

    (combined, prob) = combination
    probability = combine_probability(path.probability, prob)
    pessimisstic_probability = min(path.pessimisstic_probability, prob)
    new_elements = path.parts_ids.union([candidate.id])
    new_path = ContigsPath(combined, new_elements, probability, pessimisstic_probability)
    yield new_path


def calc_extension(finder: OverlapFinder,
                   max_acceptable_prob: Fraction,
                   contigs: Sequence[ContigWithAligner],
                   path: ContigsPath,
                   ) -> Iterator[ContigsPath]:

    for contig in contigs:
        yield from extend_by_1(finder, max_acceptable_prob, path, contig)


def calc_multiple_extensions(finder: OverlapFinder,
                             max_acceptable_prob: Fraction,
                             paths: Iterable[ContigsPath],
                             contigs: Sequence[ContigWithAligner],
                             ) -> Iterator[ContigsPath]:
    for path in paths:
        yield from calc_extension(finder, max_acceptable_prob, contigs, path)


def filter_extensions(existing: MutableMapping[str, ContigsPath],
                      extensions: Iterable[ContigsPath],
                      ) -> Iterator[ContigsPath]:
    ret = {}

    for path in extensions:
        key = path.whole.seq
        alternative = existing.get(key)
        if alternative is None or path.score() > alternative.score():
            existing[key] = path
            ret[key] = path

    yield from ret.values()


def calculate_all_paths(contigs: Sequence[ContigWithAligner]) -> Iterator[ContigsPath]:
    max_acceptable_prob = ACCEPTABLE_STITCHING_PROB
    existing: MutableMapping[str, ContigsPath] = {}
    finder = OverlapFinder.make('ACTG')
    extensions = calc_extension(finder, max_acceptable_prob, contigs, ContigsPath.empty())
    paths = tuple(filter_extensions(existing, extensions))
    yield from paths

    logger.debug("Calculating all paths...")
    cycle = 1
    while paths:
        logger.debug("Cycle %s started with %s paths.", cycle, len(paths))

        extensions = calc_multiple_extensions(finder, max_acceptable_prob, paths, contigs)
        paths = tuple(filter_extensions(existing, extensions))

        if paths:
            longest = max(paths, key=lambda path: len(path.whole.seq))
            size = len(longest.whole.seq)
            parts = len(longest.parts_ids)
            logger.debug("Cycle %s finished with %s new paths, %s [%s parts] being the longest.",
                         cycle, len(paths), size, parts)

        MAX_ALTERNATIVES = 30
        if len(paths) > MAX_ALTERNATIVES:
            # This is necessary so that the time complexity does not explode.
            # Note that if MAX_ALTERNATIVES = 1, then this algorithm becomes
            # the standard greedy algorithm that simply chooses
            # the most promising alternative at each point.
            logger.debug("Dropping %s paths that have the lowest scores.",
                         len(paths) - MAX_ALTERNATIVES)
            paths = tuple(sorted(paths, key=ContigsPath.score)[-MAX_ALTERNATIVES:])

        cycle += 1
        yield from paths
        if paths:
            max_acceptable_prob = max(x.probability for x in paths)


def find_most_probable_path(contigs: Sequence[ContigWithAligner]) -> ContigsPath:
    paths = calculate_all_paths(contigs)
    return max(paths, key=ContigsPath.score)


def stitch_consensus(contigs: Iterable[ContigWithAligner]) -> Iterator[ContigWithAligner]:
    def contig_size(contig: Contig) -> int:
        return -len(contig.seq)

    remaining = tuple(sorted(contigs, key=contig_size))
    while remaining:
        most_probable = find_most_probable_path(remaining)
        yield most_probable.whole
        remaining = tuple(contig for contig in remaining
                          if not most_probable.has_contig(contig))


def write_contigs(output_fasta: TextIO, contigs: Iterable[ContigWithAligner]):
    records = (SeqRecord(Seq.Seq(contig.seq),
                         description='',
                         id=contig.unique_name,
                         name=contig.unique_name)
               for contig in contigs)
    SeqIO.write(records, output_fasta, "fasta")


def read_contigs(input_fasta: TextIO) -> Iterable[ContigWithAligner]:
    for record in SeqIO.parse(input_fasta, "fasta"):
        yield ContigWithAligner(name=record.name, seq=str(record.seq))


def referenceless_contig_stitcher(input_fasta: TextIO,
                                  output_fasta: Optional[TextIO],
                                  ) -> int:
    with StitcherContext.fresh():
        contigs = tuple(read_contigs(input_fasta))
        logger.debug("Loaded %s contigs.", len(contigs))

        if output_fasta is not None:
            contigs = tuple(stitch_consensus(contigs))
            logger.debug("Outputting %s contigs", len(contigs))

        if output_fasta is not None:
            write_contigs(output_fasta, contigs)

        return len(contigs)
