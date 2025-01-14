from typing import Iterable, Iterator, Optional, FrozenSet, Tuple, \
    Sequence, TextIO, MutableMapping, Literal, Union
from dataclasses import dataclass
from Bio import SeqIO, Seq
from Bio.SeqRecord import SeqRecord
import logging
from mappy import Aligner
from functools import cached_property
from sortedcontainers import SortedList

from micall.utils.contig_stitcher_context import StitcherContext
from micall.utils.consensus_aligner import Alignment
from micall.utils.overlap_stitcher import align_queries, \
    calculate_concordance, sort_concordance_indexes, calc_overlap_pvalue
from micall.utils.contig_stitcher_contigs import Contig
from micall.utils.find_maximum_overlap import find_maximum_overlap, OverlapFinder


logger = logging.getLogger(__name__)

Score = int
SCORE_NOTHING = 0
SCORE_EPSILON = 1

ACCEPTABLE_STITCHING_SCORE: Score = calc_overlap_pvalue(L=15, M=15)
MAX_ALTERNATIVES = 30


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

    # Higher is better. This is an estimated probability that
    # all the components in this path came together by accident.
    probability: int

    def score(self) -> int:
        return self.probability

    def __lt__(self, other: 'ContigsPath') -> bool:
        return self.score() < other.score()

    def __le__(self, other: 'ContigsPath') -> bool:
        return self.score() <= other.score()

    def __gt__(self, other: 'ContigsPath') -> bool:
        return self.score() > other.score()

    def __ge__(self, other: 'ContigsPath') -> bool:
        return self.score() >= other.score()

    def has_contig(self, contig: Contig) -> bool:
        return contig.id in self.parts_ids

    @staticmethod
    def empty() -> 'ContigsPath':
        return ContigsPath(ContigWithAligner.empty(),
                           frozenset(),
                           SCORE_NOTHING,
                           )


@dataclass
class Pool:
    paths: SortedList[ContigsPath]
    existing: MutableMapping[str, ContigsPath]
    size: int

    @staticmethod
    def empty() -> 'Pool':
        return Pool(SortedList(), {}, 0)

    @property
    def min_acceptable_score(self) -> Score:
        if self.size > 0:
            return self.paths[0].probability
        else:
            return ACCEPTABLE_STITCHING_SCORE

    def add(self, path: ContigsPath) -> bool:
        key = path.whole.seq
        alternative = self.existing.get(key)
        if alternative is not None and alternative.score() >= path.score():
            return False

        if self.size >= MAX_ALTERNATIVES:
            to_delete = self.paths[0]
            if to_delete.score() >= path.score():
                return False

            del self.paths[0]
            self.size -= 1

        self.size += 1
        self.paths.add(path)
        self.existing[key] = path
        return True


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


def combine_probability(current: Score, new: Score) -> Score:
    return current + new


TRY_COMBINE_CACHE: MutableMapping[
    Tuple[ContigId, ContigId],
    Optional[Tuple[ContigWithAligner, Score]]] = {}


def try_combine_contigs(finder: OverlapFinder,
                        current_prob: Score,
                        pool: Pool,
                        a: ContigWithAligner, b: ContigWithAligner,
                        ) -> Optional[Tuple[ContigWithAligner, Score]]:

    if len(b.seq) == 0:
        return (a, SCORE_NOTHING)

    if len(a.seq) == 0:
        return (b, SCORE_NOTHING)

    maximum_overlap_size = min(len(a.seq), len(b.seq))
    maximum_number_of_matches = maximum_overlap_size
    maximum_result_probability = calc_overlap_pvalue(L=maximum_overlap_size, M=maximum_number_of_matches)
    if combine_probability(maximum_result_probability, current_prob) < pool.min_acceptable_score:
        return None

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
    if combine_probability(optimistic_result_probability, current_prob) < pool.min_acceptable_score:
        return None

    left_initial_overlap = left.seq[len(left.seq) - abs(shift):(len(left.seq) - abs(shift) + len(right.seq))]
    right_initial_overlap = right.seq[:abs(shift)]

    key = (left.id, right.id)
    existing: Union[Tuple[ContigWithAligner, Score], None, Literal[False]] \
        = TRY_COMBINE_CACHE.get(key, False)

    if existing is not False:
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
    if combine_probability(current_prob, result_probability) < pool.min_acceptable_score:
        return None

    is_covered = len(right.seq) < abs(shift)
    if is_covered:
        ret = (left, SCORE_EPSILON)
        TRY_COMBINE_CACHE[key] = ret
        return ret

    else:
        concordance = calculate_concordance(aligned_left, aligned_right)
        max_concordance_index = next(iter(sort_concordance_indexes(concordance)))
        left_overlap_chunk = ''.join(x for x in aligned_left[:max_concordance_index] if x != '-')
        right_overlap_chunk = ''.join(x for x in aligned_right[max_concordance_index:] if x != '-')

        result_seq = left_remainder + left_overlap_chunk + right_overlap_chunk + right_remainder
        result_contig = ContigWithAligner(None, result_seq)

        ret = (result_contig, result_probability)
        TRY_COMBINE_CACHE[key] = ret
        return ret


def extend_by_1(finder: OverlapFinder,
                pool: Pool,
                path: ContigsPath,
                candidate: ContigWithAligner,
                ) -> bool:
    if path.has_contig(candidate):
        return False

    combination = try_combine_contigs(finder, path.probability, pool, path.whole, candidate)
    if combination is None:
        return False

    (combined, prob) = combination
    probability = combine_probability(path.probability, prob)
    new_elements = path.parts_ids.union([candidate.id])
    new_path = ContigsPath(combined, new_elements, probability)
    return pool.add(new_path)


def calc_extension(finder: OverlapFinder,
                   pool: Pool,
                   contigs: Sequence[ContigWithAligner],
                   path: ContigsPath,
                   ) -> bool:
    ret = False
    for contig in contigs:
        ret = extend_by_1(finder, pool, path, contig) or ret
    return ret


def calc_multiple_extensions(finder: OverlapFinder,
                             pool: Pool,
                             paths: Iterable[ContigsPath],
                             contigs: Sequence[ContigWithAligner],
                             ) -> bool:
    ret = False
    for path in paths:
        ret = calc_extension(finder, pool, contigs, path) or ret
    return ret


def calculate_all_paths(contigs: Sequence[ContigWithAligner]) -> Iterable[ContigsPath]:
    pool = Pool.empty()
    finder = OverlapFinder.make('ACTG')
    extending = calc_extension(finder, pool, contigs, ContigsPath.empty())

    logger.debug("Calculating all paths...")
    cycle = 1
    while extending:
        logger.debug("Cycle %s started with %s paths.", cycle, pool.size)

        extending = calc_multiple_extensions(finder, pool, pool.paths, contigs)

        if pool.paths:
            longest = pool.paths[-1]
            size = len(longest.whole.seq)
            parts = len(longest.parts_ids)
            logger.debug("Cycle %s finished with %s new paths, %s [%s parts] being the fittest.",
                         cycle, pool.size, size, parts)

        cycle += 1

    return pool.paths


def find_most_probable_path(contigs: Sequence[ContigWithAligner]) -> ContigsPath:
    paths = calculate_all_paths(contigs)
    return max(paths, key=lambda path: path.score())


def contig_size_fun(contig: Contig) -> int:
    return -len(contig.seq)


def stitch_consensus_overlaps(contigs: Iterable[ContigWithAligner]) -> Iterator[ContigWithAligner]:
    remaining = tuple(sorted(contigs, key=contig_size_fun))
    while remaining:
        most_probable = find_most_probable_path(remaining)
        yield most_probable.whole
        remaining = tuple(contig for contig in remaining
                          if not most_probable.has_contig(contig))


def try_combine_1(finder: OverlapFinder,
                  contigs: Iterable[ContigWithAligner],
                  ) -> Optional[Tuple[ContigWithAligner,
                                      ContigWithAligner,
                                      ContigWithAligner]]:
    for first in contigs:
        for second in contigs:
            if first.id >= second.id:
                continue

            pool = Pool.empty()
            result = try_combine_contigs(finder=finder,
                                         current_prob=SCORE_NOTHING,
                                         pool=pool,
                                         a=first, b=second,
                                         )
            if result is not None:
                (combined, prob) = result
                return (first, second, combined)

    return None


def o2_loop(contigs: Iterable[ContigWithAligner],
            ) -> Iterator[ContigWithAligner]:
    finder = OverlapFinder.make('ACTG')
    buf = {contig.id: contig for contig in contigs}

    while True:
        combination = try_combine_1(finder, buf.values())
        if combination is None:
            break

        (first, second, combined) = combination
        del buf[first.id]
        del buf[second.id]
        buf[combined.id] = combined

    yield from sorted(buf.values(), key=contig_size_fun)


def stitch_consensus(contigs: Iterable[ContigWithAligner],
                     ) -> Iterator[ContigWithAligner]:
    contigs = stitch_consensus_overlaps(contigs)
    contigs = o2_loop(contigs)
    yield from contigs


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
