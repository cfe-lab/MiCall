import itertools
import logging
from dataclasses import dataclass
from functools import cache
from typing import Iterable, Iterator, Optional, Tuple, Sequence, TextIO, MutableMapping, Literal, Union, Set

from Bio import Seq, SeqIO
from Bio.SeqRecord import SeqRecord

from micall.utils.contig_stitcher_context import ReferencelessStitcherContext
from micall.utils.contig_stitcher_contigs import Contig
from micall.utils.overlap_stitcher import (
    align_queries,
    calculate_concordance,
    sort_concordance_indexes,
    calculate_overlap_score,
)
from micall.utils.sorted_ring import SortedRing
import micall.utils.referenceless_contig_stitcher_events as events
from micall.utils.referenceless_contig_path import ContigsPath
from micall.utils.referenceless_contig_with_aligner import ContigWithAligner
from micall.utils.referenceless_score import Score, SCORE_EPSILON, SCORE_NOTHING


logger = logging.getLogger(__name__)


@cache
def calculate_referenceless_overlap_score(L: int, M: int) -> Score:
    """
    Transform overlap scores for optimal contig path selection in referenceless stitching.

    This function takes the raw statistical overlap score and applies transformations to:
    1. Emphasize high-quality overlaps by scaling
    2. Ensure clear separation from SCORE_EPSILON

    Amplification
    -----------------------
    The (score+9)^4 transformation addresses a pathfinding problem: when combining multiple
    overlaps into contig paths, excellent overlaps should dominate the total path score
    rather than being averaged out by mediocre ones. This creates strong preference for
    paths with consistently high-quality connections.

    Separation from SCORE_EPSILON
    ------------------------
    The large offset (coming from +9 factor) ensures that no genuine
    overlap score can accidentally equal SCORE_EPSILON. This
    separation is critical because:
    - SCORE_EPSILON triggers different algorithmic behavior (ex: covered contigs)
    - All other scoring is purely relative (only <, >, >= comparisons matter)
    - The specific magnitude is irrelevant, only the ordering and separation

    Monotonicity Preservation
    -------------------------
    Both transformations preserve the monotonic ordering of the original scores,
    ensuring that better statistical overlaps always result in higher transformed scores.

    Parameters
    ----------
    L : int
        Length of the overlap region between two contigs (must be >= 0)
    M : int
        Number of matching nucleotides within the overlap region (0 <= M <= L)

    Returns
    -------
    Score
        A transformed score where:
        - Ordering reflects overlap quality (higher = better)
        - Positive values indicate genuine overlaps
        - Negative values indicate spurious overlaps
        - Zero only when L=0 (no overlap)
        - Guaranteed to be far from SCORE_EPSILON
    """

    if L == 0:
        return SCORE_NOTHING

    base = calculate_overlap_score(L=L, M=M)
    sign = 1 if base >= 0 else -1
    return sign * (abs(base) + 9) ** 4


ACCEPTABLE_STITCHING_SCORE: Score = calculate_referenceless_overlap_score(L=99, M=98)
MAX_ALTERNATIVES = 999
MIN_ALTERNATIVES = 1


def intrapolate_number_of_alternatives(n_candidates: int) -> int:
    """
    Compute how many alternative paths to explore in one stitching cycle.

    To avoid exponential blow-up in comparisons when extending
    contig-paths, we limit the total work per cycle:

        T(cycle) ≈ n_candidates × N_alt  ≤  CONSTANT

    Solving for the per-cycle alternative count:

        N_alt ≈ CONSTANT / max(1, n_candidates - 2)

    We pick CONSTANT = 999 experimentally.  The "-2" in the
    denominator compensates for the fact that at very small
    n_candidates (≤ 2) we still want at least one alternative,
    and we avoid division by zero.

    Finally, the result is clamped into [MIN_ALTERNATIVES, MAX_ALTERNATIVES]
    = [1, 999] and rounded to the nearest integer.

    Args:
        n_candidates: Number of contigs available for extension
                      in the current cycle.

    Returns:
        The (rounded and clamped) maximum number of alternative
        contig-paths to consider for this cycle.
    """

    ret = MAX_ALTERNATIVES / max(1, n_candidates - 2)
    clamped = max(MIN_ALTERNATIVES, min(MAX_ALTERNATIVES, ret))
    rounded = round(clamped)
    return rounded


def log(e: events.EventType) -> None:
    """
    Emit an event to the current stitching context and log it at debug level.
    """
    ReferencelessStitcherContext.get().emit(e)
    logger.debug("%s", e)


@dataclass
class Pool:
    ring: SortedRing[ContigsPath]
    set: Set[str]
    existing: MutableMapping[str, ContigsPath]
    smallest_score: Score

    @staticmethod
    def empty(capacity: int) -> 'Pool':
        # initial capacity is large; smallest_score starts at threshold
        ring: SortedRing[ContigsPath] = SortedRing(capacity=capacity)
        return Pool(ring, set(), {}, ACCEPTABLE_STITCHING_SCORE)

    @property
    def min_acceptable_score(self) -> Score:
        return self.smallest_score

    def add(self, path: ContigsPath) -> bool:
        key = path.whole.seq
        alternative = self.existing.get(key)
        if alternative is not None and alternative.get_score() >= path.get_score():
            return False

        # insert the new path and record it
        self.existing[key] = path

        deleted_seq = None
        old_size = len(self.ring)

        if path.whole.seq in self.set:
            to_delete_index = -1
            for i, candidate in enumerate(self.ring):
                if candidate.whole.seq == path.whole.seq:
                    to_delete_index = i
                    break

            assert to_delete_index >= 0
            deleted_seq = path.whole.seq
            del self.ring[to_delete_index]

        else:
            if len(self.ring) > 0:
                deleted_seq = self.ring[0].whole.seq

        if self.ring.insert(path):
            new_size = len(self.ring)

            if new_size == old_size and deleted_seq is not None:
                self.set.remove(deleted_seq)

            self.set.add(path.whole.seq)
            self.smallest_score = max(self.ring[0].get_score(), ACCEPTABLE_STITCHING_SCORE)
            return True

        return False


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
    """
    Combine two scores by summing them.
    """
    return current + new


def get_minimum_base_score(current: Score, minimum: Score) -> Score:
    """
    Calculate the minimum additional score required so that
    current + additional >= minimum.
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

    else:
        left_overlap = left.seq[left_cutoff:(left_cutoff + len(right.seq))]
        left_remainder = left.seq[:left_cutoff]
        right_overlap = right.seq[:right_cutoff]
        right_remainder = right.seq[right_cutoff:]
        aligned_1, aligned_2 = align_overlaps(left_overlap, right_overlap)

    result_length = len(aligned_1)
    assert result_length > 0, "The overlap cannot be empty."
    number_of_matches = sum(1 for x, y
                            in zip(aligned_1, aligned_2)
                            if x == y and x != '-')

    # Note that result_length is not necessarily == len(left_overlap_chunk) + len(right_overlap_chunk).
    # The addition would give a more precise value for the result_score, but it's much more expensive to calculate.
    result_score = calculate_referenceless_overlap_score(L=result_length, M=number_of_matches)
    if is_debug2:
        denominator = max(minimum_base_score, ACCEPTABLE_STITCHING_SCORE)
        relative_score = result_score / denominator if denominator != 0 else float("inf")
        log(events.DeterminedOverlap(left.unique_name, right.unique_name,
                                     result_length - is_covered - 1,
                                     number_of_matches, relative_score))

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
                ) -> Optional[ContigsPath]:
    """
    Attempt to extend a contig path by one candidate contig.
    If a valid combination is found, return it.
    Otherwise return None.
    """

    if path.has_contig(candidate):
        return None

    combination = try_combine_contigs(is_debug2, path.score, pool, path.whole, candidate)
    if combination is None:
        return None

    combined, additional_score = combination
    score = combine_scores(path.score, additional_score)
    is_covered = additional_score == SCORE_EPSILON
    if is_covered:
        new_elements = path.contigs_ids
    else:
        new_elements = path.contigs_ids.union([candidate.id])
    new_contained_elements = path.contains_contigs_ids.union([candidate.id])
    new_path = ContigsPath(combined, new_elements, new_contained_elements, score)
    return new_path


def calc_extension(is_debug2: bool,
                   pool: Pool,
                   contigs: Sequence[ContigWithAligner],
                   path: ContigsPath,
                   ) -> bool:
    """
    Try to extend a single path with each contig in contigs.
    Return True if any extension was added to the pool.
    """

    ret = False
    for contig in contigs:
        new = extend_by_1(is_debug2, pool, path, contig)
        if new is not None:
            ret = pool.add(new) or ret

    return ret


def calc_multiple_extensions(is_debug2: bool,
                             pool: Pool,
                             paths: Iterable[ContigsPath],
                             contigs: Sequence[ContigWithAligner],
                             ) -> bool:
    """
    Attempt to extend multiple paths with a set of contigs.
    Returns True if any new extensions were added to the pool.
    """
    ret = False
    for path in paths:
        ret = calc_extension(is_debug2, pool, contigs, path) or ret
    return ret


def calculate_all_paths(paths: Sequence[ContigsPath],
                        contigs: Sequence[ContigWithAligner],
                        ) -> Iterable[ContigsPath]:
    """
    Iteratively extend seed paths with contigs to generate candidate contig paths.
    Returns an iterable of best paths up to a max number of alternatives.
    """
    is_debug2 = ReferencelessStitcherContext.get().is_debug2

    max_alternatives = intrapolate_number_of_alternatives(len(contigs))
    pool = Pool.empty(max_alternatives)
    for path in sorted(paths):
        # stop if we reached capacity
        if len(pool.ring) >= pool.ring.capacity:
            break
        pool.add(path)

    log(events.CalculatingAll())
    for cycle in itertools.count(1):
        # log start with current pool size
        log(events.CycleStart(cycle, len(pool.ring)))

        if not calc_multiple_extensions(is_debug2, pool, pool.ring, contigs):
            break

        # log end with updated pool size
        log(events.CycleEnd(cycle, len(pool.ring), pool))

    return pool.ring


def find_most_probable_path(seeds: Sequence[ContigsPath],
                            contigs: Sequence[ContigWithAligner],
                            ) -> ContigsPath:
    """
    Select the single most probable contig path from seeds extended by contigs.
    """
    paths = calculate_all_paths(seeds, contigs)
    return max(paths, key=ContigsPath.get_score)


def contig_size_fun(contig: Contig) -> int:
    """
    Sorting key to order contigs by descending sequence length.
    """
    return -len(contig.seq)


def stitch_consensus_overlaps(contigs: Iterable[ContigWithAligner]) -> Iterator[ContigWithAligner]:
    """
    Produce a consensus stitching of contigs based on overlaps.
    Iteratively selects and stitches the most probable contig paths.
    """
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
                      if most_probable.contigs_ids.isdisjoint(path.contigs_ids))
        log(events.Remove(len(most_probable.contigs_ids), len(remaining)))
        if len(most_probable.contains_contigs_ids) == 1:
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

            pool = Pool.empty(1)
            result = try_combine_contigs(
                is_debug2=is_debug2,
                current_score=SCORE_NOTHING,
                pool=pool,
                a=first,
                b=second,
            )
            if result is not None:
                combined, additional_score = result
                return first, second, combined

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
    """
    Execute the full referenceless contig stitching algorithm.
    First perform consensus overlap stitching, then O2 greedy merging.
    """
    contigs = tuple(stitch_consensus_overlaps(contigs))
    log(events.InitiallyProduced(len(contigs)))
    contigs = o2_loop(contigs)
    yield from contigs


def write_contigs(output_fasta: TextIO, contigs: Iterable[ContigWithAligner]):
    """
    Write an iterable of contigs to the given output FASTA file handle.
    """
    records = (SeqRecord(Seq.Seq(contig.seq),
                         description='',
                         id=contig.unique_name,
                         name=contig.unique_name)
               for contig in contigs)
    SeqIO.write(records, output_fasta, "fasta")


def read_contigs(input_fasta: TextIO) -> Iterable[ContigWithAligner]:
    """
    Read contigs from a FASTA file handle, yielding ContigWithAligner objects.
    """
    for record in SeqIO.parse(input_fasta, "fasta"):
        yield ContigWithAligner(name=record.name, seq=str(record.seq))


def referenceless_contig_stitcher_with_ctx(
        input_fasta: TextIO,
        output_fasta: Optional[TextIO],
) -> int:
    """
    Main entrypoint for referenceless contig stitching with context.
    Reads input contigs, performs stitching if output is specified, and writes results.
    Returns the number of contigs (stitched or original).
    """
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
    """
    Wrapper that initializes a fresh stitching context and calls the core stitching function.
    """
    with ReferencelessStitcherContext.fresh() as ctx:
        if logger.level == logging.DEBUG - 1:
            ctx.is_debug2 = True
        return referenceless_contig_stitcher_with_ctx(input_fasta, output_fasta)
