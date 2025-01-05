from typing import Iterable, Iterator, Optional, FrozenSet, Tuple, Sequence, TextIO, MutableMapping
from dataclasses import dataclass
from fractions import Fraction
from Bio import SeqIO, Seq
from Bio.SeqRecord import SeqRecord
import logging
from mappy import Aligner

from micall.utils.contig_stitcher_context import StitcherContext
from micall.utils.consensus_aligner import Alignment
from micall.utils.overlap_stitcher import align_queries, \
    calculate_concordance, sort_concordance_indexes, calc_overlap_pvalue
from micall.utils.contig_stitcher_contigs import Contig
from micall.utils.find_maximum_overlap import find_maximum_overlap


logger = logging.getLogger(__name__)

ACCEPTABLE_STITCHING_PROB = Fraction(1, 1000000)


@dataclass(frozen=True)
class ContigsPath:
    # Contig representing all combined contigs in the path.
    whole: Contig

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
        return ContigsPath(Contig.empty(), frozenset(), Fraction(1), Fraction(1))


@dataclass(frozen=True)
class Overlap:

    # A negative integer.
    # Represents the shift value for left contig relative to right
    # contig to achieve the maximum overlap.
    shift: int


def get_overlap(left: Contig, right: Contig) -> Optional[Overlap]:
    if len(left.seq) == 0 or len(right.seq) == 0:
        return None

    shift = find_maximum_overlap(left.seq, right.seq)
    if shift == 0:
        return None

    return Overlap(shift)


def map_overlap_onto_candidate(overlap: str, candidate: str) -> Iterator[Alignment]:
    # TODO: Move this implementation into consensus_aligner maybe.
    aligner = Aligner(seq=candidate, bw=500, bw_long=500, preset='map-ont')
    for x in aligner.map(overlap):
        if x.is_primary:
            yield x


def try_combine_contigs(a: Contig, b: Contig,
                        ) -> Optional[Tuple[Contig, Fraction]]:
    # TODO: Memoize this function.
    #       Two-layer caching seems most optimal:
    #       first by key=contig.id, then by key=contig.seq.

    if len(b.seq) == 0:
        return (a, Fraction(1))

    if len(a.seq) == 0:
        return (b, Fraction(1))

    overlap = get_overlap(a, b)
    if overlap is None:
        return None

    logger.debug("Trying to put together %s with %s (with lengts %s and %s).",
                 a.unique_name, b.unique_name, len(a.seq), len(b.seq))

    if abs(overlap.shift) < len(a.seq):
        shift = overlap.shift
        left = a
        right = b
    else:
        shift = (len(b.seq) + len(a.seq) - abs(overlap.shift)) * -1
        left = b
        right = a

    logger.debug("Overlap of size %s detected between %s and %s.",
                 min([abs(shift), len(a.seq), len(b.seq)]),
                 a.unique_name, b.unique_name)

    is_covered = len(right.seq) < abs(shift)
    if is_covered:
        logger.debug("Contig %s is covered by %s.", right.unique_name, left.unique_name)
    else:
        logger.debug("Contig %s comes before %s.", left.unique_name, right.unique_name)

    left_initial_overlap = left.seq[len(left.seq) - abs(shift):(len(left.seq) - abs(shift) + len(right.seq))]
    right_initial_overlap = right.seq[:abs(shift)]

    left_overlap_alignments = map_overlap_onto_candidate(str(right_initial_overlap), str(left.seq))
    right_overlap_alignments = map_overlap_onto_candidate(str(left_initial_overlap), str(right.seq))
    left_cutoff = min((al.r_st for al in left_overlap_alignments), default=None)
    right_cutoff = max((al.r_en for al in right_overlap_alignments), default=None)
    if left_cutoff is None:
        logger.debug("Overlap alignment between %s and %s failed.", a.unique_name, b.unique_name)
        return None
    if right_cutoff is None:
        logger.debug("Overlap alignment between %s and %s failed.", a.unique_name, b.unique_name)
        return None

    left_overlap = left.seq[left_cutoff:(left_cutoff + len(right.seq))]
    left_remainder = left.seq[:left_cutoff]
    right_overlap = right.seq[:right_cutoff]
    right_remainder = right.seq[right_cutoff:]

    logger.debug("Cut off sizes between %s and %s are %s and %s.",
                 a.unique_name, b.unique_name,
                 len(left_overlap), len(right_overlap))

    aligned_left, aligned_right = align_queries(str(left_overlap), str(right_overlap))

    number_of_matches = sum(1 for x, y
                            in zip(aligned_left, aligned_right)
                            if x == y and x != '-')
    result_probability = calc_overlap_pvalue(L=len(left_overlap), M=number_of_matches)
    logger.debug("Overlap probability between %s and %s is %s.",
                 a.unique_name, b.unique_name, f"{float(1 - result_probability):.5f}")

    if result_probability > ACCEPTABLE_STITCHING_PROB:
        # FIXME: Adjust the threold to something more based.
        # FIXME: Print a warning that this happened.
        logger.debug("Overlap probability between %s and %s is too small.", a.unique_name, b.unique_name)
        return None

    if is_covered:
        logger.debug("Between %s and %s, returning the covering contig, %s.",
                     a.unique_name, b.unique_name, left.unique_name)
        return (left, Fraction(1))

    else:
        concordance = calculate_concordance(aligned_left, aligned_right)
        max_concordance_index = next(iter(sort_concordance_indexes(concordance)))
        left_overlap_chunk = ''.join(x for x in aligned_left[:max_concordance_index] if x != '-')
        right_overlap_chunk = ''.join(x for x in aligned_right[max_concordance_index:] if x != '-')

        result_seq = left_remainder + left_overlap_chunk + right_overlap_chunk + right_remainder
        result_contig = Contig(None, result_seq)

        logger.debug("Joined %s and %s together in a contig %s with lengh %s.",
                     a.unique_name, b.unique_name,
                     result_contig.unique_name, len(result_contig.seq))

        return (result_contig, result_probability)


def extend_by_1(path: ContigsPath, candidate: Contig) -> Iterator[ContigsPath]:
    if path.has_contig(candidate):
        return

    combination = try_combine_contigs(path.whole, candidate)
    if combination is None:
        return

    (combined, prob) = combination
    probability = path.probability * prob
    if prob == 1:
        pessimisstic_probability = path.pessimisstic_probability
    else:
        pessimisstic_probability = max(path.pessimisstic_probability, prob)
    new_elements = path.parts_ids.union([candidate.id])
    new_path = ContigsPath(combined, new_elements, probability, pessimisstic_probability)
    yield new_path


def calc_extension(contigs: Sequence[Contig],
                   path: ContigsPath,
                   ) -> Iterator[ContigsPath]:

    for contig in contigs:
        yield from extend_by_1(path, contig)


def calc_multiple_extensions(paths: Iterable[ContigsPath],
                             contigs: Sequence[Contig],
                             ) -> Iterator[ContigsPath]:
    for path in paths:
        yield from calc_extension(contigs, path)


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


def calculate_all_paths(contigs: Sequence[Contig]) -> Iterator[ContigsPath]:
    existing: MutableMapping[str, ContigsPath] = {}

    extensions = calc_extension(contigs, ContigsPath.empty())
    paths = tuple(filter_extensions(existing, extensions))
    yield from paths

    logger.debug("Calculating all paths...")
    cycle = 1
    while paths:
        logger.debug("Cycle %s started with %s paths.", cycle, len(paths))

        extensions = calc_multiple_extensions(paths, contigs)
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


def find_most_probable_path(contigs: Sequence[Contig]) -> ContigsPath:
    paths = calculate_all_paths(contigs)
    return max(paths, key=ContigsPath.score)


def stitch_consensus(contigs: Iterable[Contig]) -> Iterable[Contig]:
    remaining = tuple(contigs)
    while remaining:
        most_probable = find_most_probable_path(remaining)
        yield most_probable.whole
        remaining = tuple(contig for contig in remaining
                          if not most_probable.has_contig(contig))


def write_contigs(output_fasta: TextIO, contigs: Iterable[Contig]):
    records = (SeqRecord(Seq.Seq(contig.seq),
                         description='',
                         id=contig.unique_name,
                         name=contig.unique_name)
               for contig in contigs)
    SeqIO.write(records, output_fasta, "fasta")


def read_contigs(input_fasta: TextIO) -> Iterable[Contig]:
    for record in SeqIO.parse(input_fasta, "fasta"):
        yield Contig(name=record.name, seq=record.seq)


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
