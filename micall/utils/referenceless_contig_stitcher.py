from typing import Iterable, Iterator, Optional, Tuple, \
    Sequence, TextIO, MutableMapping, Literal, Union
from dataclasses import dataclass
from Bio import SeqIO, Seq
from Bio.SeqRecord import SeqRecord
import logging
from sortedcontainers import SortedList
import itertools
from functools import cache

from micall.utils.contig_stitcher_context import ReferencelessStitcherContext
from micall.utils.overlap_stitcher import align_queries, \
    calculate_concordance, sort_concordance_indexes, calculate_overlap_score
import micall.utils.referenceless_contig_stitcher_events as events
from micall.utils.referenceless_contig_with_aligner import ContigWithAligner
from micall.utils.referenceless_contig_path import ContigsPath
from micall.utils.referenceless_score import Score, SCORE_EPSILON, SCORE_NOTHING
from micall.utils.contig_stitcher_contigs import Contig


logger = logging.getLogger(__name__)


@cache
def calculate_referenceless_overlap_score(L: int, M: int) -> Score:
    return 1024 * calculate_overlap_score(L=L, M=M)


ACCEPTABLE_STITCHING_SCORE: Score = calculate_referenceless_overlap_score(L=15, M=14)
MAX_ALTERNATIVES = 30


def log(e: events.EventType) -> None:
    ReferencelessStitcherContext.get().emit(e)
    logger.debug("%s", e)


@dataclass
class Pool:
    # TODO: Factor out a SortedRing structure out of here.
    paths: SortedList[ContigsPath]
    existing: MutableMapping[str, ContigsPath]
    size: int
    capacity: int
    smallest_score: Score

    @staticmethod
    def empty() -> 'Pool':
        return Pool(SortedList(), {}, 0, 999999, ACCEPTABLE_STITCHING_SCORE)

    @property
    def min_acceptable_score(self) -> Score:
        return self.smallest_score

    def resize(self, new_capacity: int) -> None:
        if new_capacity < self.size:
            stop = self.size - new_capacity
            del self.paths[:stop]
            self.size = new_capacity

        self.capacity = new_capacity

        if self.size > 0:
            smallest_path = self.paths[0]
            self.smallest_score = max(ACCEPTABLE_STITCHING_SCORE,
                                      smallest_path.score)
        else:
            self.smallest_score = ACCEPTABLE_STITCHING_SCORE

    def add(self, path: ContigsPath) -> bool:
        key = path.whole.seq
        alternative = self.existing.get(key)
        if alternative is not None and alternative.get_score() >= path.get_score():
            return False

        if self.size > 0:
            smallest_path = self.paths[0]
            if self.size >= self.capacity:
                if smallest_path.get_score() >= path.get_score():
                    return False

                del self.paths[0]
                self.size -= 1

            if path.get_score() < smallest_path.get_score():
                smallest_path = path

            self.smallest_score = max(ACCEPTABLE_STITCHING_SCORE,
                                      smallest_path.score)
        else:
            self.smallest_score = max(ACCEPTABLE_STITCHING_SCORE,
                                      path.get_score())

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


def get_overlap(left: ContigWithAligner, right: ContigWithAligner) -> Optional[Overlap]:
    if len(left.seq) == 0 or len(right.seq) == 0:
        return None

    key = (left.id, right.id)
    cached: Union[Overlap, None, Literal[0]] = GET_OVERLAP_CACHE.get(key, 0)
    if cached != 0:
        return cached

    (shift, value) = left.find_maximum_overlap(right)
    if shift == 0:
        ret = None
        GET_OVERLAP_CACHE[key] = ret
        return ret

    if abs(shift) <= len(left.seq):
        size = min(abs(shift), len(right.seq))
    else:
        size = min(len(left.seq), len(left.seq) + len(right.seq) + shift)

    assert size > 0, f"{shift}, {len(left.seq)}, {len(right.seq)}"
    assert size <= len(left.seq), f"{shift}, {size}, {len(left.seq)}, {len(right.seq)}"
    assert size <= len(right.seq), f"{shift}, {size}, {len(left.seq)}, {len(right.seq)}"

    ret = Overlap(shift=shift, size=size)
    GET_OVERLAP_CACHE[key] = ret
    return ret


def combine_scores(current: Score, new: Score) -> Score:
    return current + new


def get_minimum_base_score(current: Score, minimum: Score) -> Score:
    """
    Returns minimum score that when combined with current is still acceptable.
    """

    return minimum - current


ALIGN_CACHE: MutableMapping[Tuple[str, str], Tuple[str, str]] = {}


def align_overlaps(left_overlap: str, right_overlap: str) -> Tuple[str, str]:
    key = (left_overlap, right_overlap)
    existing = ALIGN_CACHE.get(key)
    if existing is not None:
        return existing

    result = align_queries(left_overlap, right_overlap)
    ALIGN_CACHE[key] = result
    return result


CutoffsCacheResult = Optional[Tuple[int, int]]
CutoffsCache = MutableMapping[Tuple[ContigId, ContigId], CutoffsCacheResult]
CUTOFFS_CACHE: CutoffsCache = {}


def find_overlap_cutoffs(minimum_score: Score,
                         left: ContigWithAligner,
                         right: ContigWithAligner,
                         shift: int,
                         left_initial_overlap: str,
                         right_initial_overlap: str,
                         ) -> CutoffsCacheResult:

    key = (left.id, right.id)

    existing: Union[CutoffsCacheResult, Literal[False]] \
        = CUTOFFS_CACHE.get(key, False)

    if existing is None:
        return None

    elif existing is not False:
        (left_cutoff, right_cutoff) = existing

    elif len(left.seq) == len(left_initial_overlap):
        overlap_alignments = tuple(right.map_overlap(minimum_score, "cover", left_initial_overlap))
        right_cutoff = max((end for start, end in overlap_alignments), default=-1)
        if right_cutoff < 0:
            ret = (len(right.seq) - abs(shift) + 1, len(right.seq) - abs(shift) + len(left_initial_overlap) + 1)
            CUTOFFS_CACHE[key] = ret
            return ret

        left_cutoff = min((start for start, end in overlap_alignments), default=-1)

    elif len(right.seq) == len(right_initial_overlap):
        overlap_alignments = tuple(left.map_overlap(minimum_score, "cover", right_initial_overlap))
        left_cutoff = min((start for start, end in overlap_alignments), default=-1)
        if left_cutoff < 0:
            ret = (len(left.seq) - abs(shift) + 1, len(left.seq) - abs(shift) + len(right_initial_overlap) + 1)
            CUTOFFS_CACHE[key] = ret
            return ret

        right_cutoff = max((end for start, end in overlap_alignments), default=-1)

    elif len(left.seq) < len(right.seq):
        left_overlap_alignments = left.map_overlap(minimum_score, "left", right_initial_overlap)
        left_cutoff = min((start for start, end in left_overlap_alignments), default=-1)
        if left_cutoff < 0:
            CUTOFFS_CACHE[key] = None
            return None

        right_overlap_alignments = right.map_overlap(minimum_score, "right", left_initial_overlap)
        right_cutoff = max((end for start, end in right_overlap_alignments), default=-1)
        if right_cutoff < 0:
            CUTOFFS_CACHE[key] = None
            return None
    else:
        right_overlap_alignments = right.map_overlap(minimum_score, "right", left_initial_overlap)
        right_cutoff = max((end for start, end in right_overlap_alignments), default=-1)
        if right_cutoff < 0:
            CUTOFFS_CACHE[key] = None
            return None

        left_overlap_alignments = left.map_overlap(minimum_score, "left", right_initial_overlap)
        left_cutoff = min((start for start, end in left_overlap_alignments), default=-1)
        if left_cutoff < 0:
            CUTOFFS_CACHE[key] = None
            return None

    ret = (left_cutoff, right_cutoff)
    CUTOFFS_CACHE[key] = ret
    return ret


def try_combine_contigs(is_debug2: bool,
                        current_score: Score,
                        pool: Pool,
                        a: ContigWithAligner, b: ContigWithAligner,
                        ) -> Optional[Tuple[ContigWithAligner, Score]]:

    minimum_base_score = get_minimum_base_score(current_score, pool.min_acceptable_score)

    maximum_overlap_size = min(len(a.seq), len(b.seq)) - 1
    maximum_number_of_matches = maximum_overlap_size
    maximum_result_score = calculate_referenceless_overlap_score(L=maximum_overlap_size+1, M=maximum_number_of_matches)
    if maximum_result_score < minimum_base_score:
        return None

    overlap = get_overlap(a, b)
    if overlap is None:
        return None

    optimistic_overlap_size = overlap.size
    optimistic_number_of_matches = optimistic_overlap_size
    optimistic_result_score = calculate_referenceless_overlap_score(L=optimistic_overlap_size+1, M=optimistic_number_of_matches)
    if optimistic_result_score < minimum_base_score:
        return None

    if abs(overlap.shift) < len(a.seq):
        shift = overlap.shift
        left = a
        right = b
    else:
        shift = (len(b.seq) + len(a.seq) - abs(overlap.shift)) * -1
        left = b
        right = a

    left_initial_overlap = left.seq[len(left.seq) - abs(shift):(len(left.seq) - abs(shift) + len(right.seq))]
    right_initial_overlap = right.seq[:abs(shift)]

    assert len(right_initial_overlap) == overlap.size, f"{len(right_initial_overlap)} == {overlap.size}"
    assert len(left_initial_overlap) == overlap.size, f"{len(left_initial_overlap)} == {overlap.size}"
    assert calculate_referenceless_overlap_score(L=len(left_initial_overlap)+1, M=len(left_initial_overlap)) >= minimum_base_score

    cutoffs = find_overlap_cutoffs(minimum_base_score,
                                   left, right, shift,
                                   left_initial_overlap,
                                   right_initial_overlap,
                                   )

    if is_debug2:
        log(events.CalculatedCutoffs(left.unique_name, right.unique_name,
                                     len(left_initial_overlap),
                                     cutoffs))

    if cutoffs is None:
        return None

    (left_cutoff, right_cutoff) = cutoffs

    left_is_covered = len(left.seq) <= overlap.size
    right_is_covered = len(right.seq) <= overlap.size
    is_covered = left_is_covered or right_is_covered
    if is_covered:
        if left_is_covered:
            covered = left
            bigger = right
        else:
            covered = right
            bigger = left

        covered_overlap = covered.seq
        bigger_overlap = bigger.seq[left_cutoff:right_cutoff]
        aligned_1, aligned_2 = align_overlaps(covered_overlap, bigger_overlap)
        result_length = max(len(bigger_overlap), len(covered_overlap)) + 2

    else:
        left_overlap = left.seq[left_cutoff:(left_cutoff + len(right.seq))]
        left_remainder = left.seq[:left_cutoff]
        right_overlap = right.seq[:right_cutoff]
        right_remainder = right.seq[right_cutoff:]
        aligned_1, aligned_2 = align_overlaps(left_overlap, right_overlap)
        result_length = max(len(left_overlap), len(right_overlap)) + 1

    assert result_length > 0, "The overlap cannot be empty."
    number_of_matches = sum(1 for x, y
                            in zip(aligned_1, aligned_2)
                            if x == y and x != '-')

    # Note that result_length is not necessarily == len(left_overlap_chunk) + len(right_overlap_chunk).
    # The addition would give a more precise value for the result_score, but it's much more expensive to calculate.
    result_score = calculate_referenceless_overlap_score(L=result_length, M=number_of_matches)
    if is_debug2:
        log(events.DeterminedOverlap(left.unique_name, right.unique_name,
                                     result_length - is_covered - 1,
                                     number_of_matches, result_score))

    if result_score < minimum_base_score:
        return None

    if is_covered:
        if is_debug2:
            log(events.Covered(left.unique_name, right.unique_name))
        return (bigger, SCORE_EPSILON)

    else:
        concordance = calculate_concordance(aligned_1, aligned_2)
        max_concordance_index = next(iter(sort_concordance_indexes(concordance)))
        left_overlap_chunk = ''.join(x for x in aligned_1[:max_concordance_index] if x != '-')
        right_overlap_chunk = ''.join(x for x in aligned_2[max_concordance_index:] if x != '-')

        result_seq = left_remainder + left_overlap_chunk + right_overlap_chunk + right_remainder
        result_contig = ContigWithAligner(None, result_seq)

        if is_debug2:
            log(events.CombinedContings(left_contig=left.unique_name,
                                        right_contig=right.unique_name,
                                        result_contig=result_contig.unique_name,
                                        overlap_size=len(left_overlap_chunk) + len(right_overlap_chunk),
                                        ))
        return (result_contig, result_score)


def extend_by_1(is_debug2: bool,
                pool: Pool,
                path: ContigsPath,
                candidate: ContigWithAligner,
                ) -> bool:
    if path.has_contig(candidate):
        return False

    combination = try_combine_contigs(is_debug2, path.score, pool, path.whole, candidate)
    if combination is None:
        return False

    (combined, additional_score) = combination
    score = combine_scores(path.score, additional_score)
    new_elements = path.parts_ids.union([candidate.id])
    new_path = ContigsPath(combined, new_elements, score)
    return pool.add(new_path)


def calc_extension(is_debug2: bool,
                   pool: Pool,
                   contigs: Sequence[ContigWithAligner],
                   path: ContigsPath,
                   ) -> bool:
    ret = False
    for contig in contigs:
        ret = extend_by_1(is_debug2, pool, path, contig) or ret
    return ret


def calc_multiple_extensions(is_debug2: bool,
                             pool: Pool,
                             paths: Iterable[ContigsPath],
                             contigs: Sequence[ContigWithAligner],
                             ) -> bool:
    ret = False
    for path in paths:
        ret = calc_extension(is_debug2, pool, contigs, path) or ret
    return ret


def calculate_all_paths(paths: Sequence[ContigsPath],
                        contigs: Sequence[ContigWithAligner],
                        ) -> Iterable[ContigsPath]:
    is_debug2 = ReferencelessStitcherContext.get().is_debug2

    pool = Pool.empty()
    pool.resize(MAX_ALTERNATIVES)

    for path in sorted(paths):
        if pool.size >= pool.capacity:
            break
        pool.add(path)

    log(events.CalculatingAll())
    for cycle in itertools.count(1):
        log(events.CycleStart(cycle, pool.size))

        if not calc_multiple_extensions(is_debug2, pool, pool.paths, contigs):
            break

        log(events.CycleEnd(cycle, pool.size, pool))

    return pool.paths


def find_most_probable_path(seeds: Sequence[ContigsPath],
                            contigs: Sequence[ContigWithAligner],
                            ) -> ContigsPath:
    paths = calculate_all_paths(seeds, contigs)
    return max(paths, key=ContigsPath.get_score)


def contig_size_fun(contig: Contig) -> int:
    return -len(contig.seq)


def stitch_consensus_overlaps(contigs: Iterable[ContigWithAligner]) -> Iterator[ContigWithAligner]:
    remaining = tuple(sorted(contigs, key=contig_size_fun))

    log(events.InitializingSeeds())
    seeds = tuple(sorted(map(ContigsPath.singleton, remaining),
                         key=lambda path: len(path.whole.seq),
                         reverse=True))
    log(events.Starting(len(seeds)))

    while remaining:
        most_probable = find_most_probable_path(seeds, remaining)
        log(events.Constructed(most_probable))
        yield most_probable.whole
        remaining = tuple(contig for contig in remaining
                          if not most_probable.has_contig(contig))
        seeds = tuple(path for path in seeds
                      if most_probable.parts_ids.isdisjoint(path.parts_ids))
        log(events.Remove(len(most_probable.parts_ids), len(remaining)))
        if len(most_probable.parts_ids) == 1:
            log(events.GiveUp())
            yield from remaining
            return


def try_combine_1(contigs: Iterable[ContigWithAligner],
                  ) -> Optional[Tuple[ContigWithAligner,
                                      ContigWithAligner,
                                      ContigWithAligner]]:
    is_debug2 = ReferencelessStitcherContext.get().is_debug2

    for first in contigs:
        for second in contigs:
            if first.id >= second.id:
                continue

            pool = Pool.empty()
            result = try_combine_contigs(is_debug2=is_debug2,
                                         current_score=SCORE_NOTHING,
                                         pool=pool,
                                         a=first, b=second,
                                         )
            if result is not None:
                (combined, additional_score) = result
                return (first, second, combined)

    return None


def o2_loop(contigs: Iterable[ContigWithAligner],
            ) -> Iterator[ContigWithAligner]:
    buf = {contig.id: contig for contig in contigs}

    while True:
        combination = try_combine_1(buf.values())
        if combination is None:
            break

        (first, second, combined) = combination
        del buf[first.id]
        del buf[second.id]
        buf[combined.id] = combined

    yield from sorted(buf.values(), key=contig_size_fun)


def stitch_consensus(contigs: Iterable[ContigWithAligner],
                     ) -> Iterator[ContigWithAligner]:
    contigs = tuple(stitch_consensus_overlaps(contigs))
    log(events.InitiallyProduced(len(contigs)))
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


def referenceless_contig_stitcher_with_ctx(
        input_fasta: TextIO,
        output_fasta: Optional[TextIO],
) -> int:
    contigs = tuple(read_contigs(input_fasta))
    log(events.Loaded(len(contigs)))

    if output_fasta is not None:
        contigs = tuple(stitch_consensus(contigs))
        log(events.Outputting(len(contigs)))

    if output_fasta is not None:
        write_contigs(output_fasta, contigs)

    return len(contigs)


def referenceless_contig_stitcher(input_fasta: TextIO,
                                  output_fasta: Optional[TextIO],
                                  ) -> int:
    with ReferencelessStitcherContext.fresh() as ctx:
        if logger.level == logging.DEBUG - 1:
            ctx.is_debug2 = True
        return referenceless_contig_stitcher_with_ctx(input_fasta, output_fasta)
