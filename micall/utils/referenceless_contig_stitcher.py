import itertools
import logging
from dataclasses import dataclass
from functools import cache
from typing import (
    Iterable,
    Iterator,
    Optional,
    Tuple,
    Sequence,
    TextIO,
    MutableMapping,
    Set,
)

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
from .referenceless_contig_stitcher_pool import Pool
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
    The score^exp transformation addresses a pathfinding problem: when combining multiple
    overlaps into contig paths, excellent overlaps should dominate the total path score
    rather than being averaged out by mediocre ones. This creates strong preference for
    paths with consistently high-quality connections.

    Separation from SCORE_EPSILON
    ------------------------
    The large offset (coming from 999 factor) ensures that no genuine
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

    base = calculate_overlap_score(L=L, M=M) - ACCEPTABLE_BASE_STITCHING_SCORE()
    sign = 1 if base >= 0 else -1
    magnitude = 999 + (999 * base) ** 2
    return sign * magnitude


MIN_MATCHES = 99


@cache
def ACCEPTABLE_BASE_STITCHING_SCORE():
    return calculate_overlap_score(L=MIN_MATCHES + 1, M=MIN_MATCHES)


@cache
def ACCEPTABLE_STITCHING_SCORE():
    return calculate_referenceless_overlap_score(L=MIN_MATCHES + 1, M=MIN_MATCHES)


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


@dataclass(frozen=True)
class Overlap:
    """Represents a maximal-overlap placement between two contigs.

    Attributes:
        shift: Relative shift to place `left` vs `right` at their best overlap.
            The value is negative when `right` starts before the end of `left`.
            A value of 0 indicates no overlap (we never construct Overlap with 0).
        size: Length (in bases) of the overlapping window induced by `shift`.
    """

    shift: int
    size: int


ContigId = int
GET_OVERLAP_CACHE: MutableMapping[Tuple[ContigId, ContigId], Overlap] = {}
GET_OVERLAP_NEGATIVE: Set[Tuple[ContigId, ContigId]] = set()


def compute_overlap_size(left_len: int, right_len: int, shift: int) -> int:
    """Compute the overlap size given contig lengths and a shift.

    The shift is defined as returned by `ContigWithAligner.find_maximum_overlap`.
    """
    if abs(shift) <= left_len:
        size = min(abs(shift), right_len)
    else:
        size = min(left_len, left_len + right_len + shift)
    assert size > 0, f"{shift}, {left_len}, {right_len}"
    assert size <= left_len, f"{shift}, {size}, {left_len}, {right_len}"
    assert size <= right_len, f"{shift}, {size}, {left_len}, {right_len}"
    return size


def get_overlap(left: ContigWithAligner, right: ContigWithAligner) -> Optional[Overlap]:
    """Return the best overlap placement between two contigs, if any.

    Uses an internal cache to avoid repeated expensive overlap calculations.

    Returns:
        Overlap(shift, size) if an overlap was found; otherwise None.
    """
    if len(left.seq) == 0 or len(right.seq) == 0:
        return None

    key = (left.id, right.id)
    if key in GET_OVERLAP_NEGATIVE:
        return None
    cached = GET_OVERLAP_CACHE.get(key)
    if cached is not None:
        return cached

    shift, _ = left.find_maximum_overlap(right)
    if shift == 0:
        GET_OVERLAP_NEGATIVE.add(key)
        return None

    size = compute_overlap_size(len(left.seq), len(right.seq), shift)
    ret = Overlap(shift=shift, size=size)
    GET_OVERLAP_CACHE[key] = ret
    return ret


def combine_scores(current: Score, new: Score) -> Score:
    return current + new


def get_minimum_base_score(current: Score, minimum: Score) -> Score:
    """
    Calculate the minimum additional score required so that
    combine_scores(current, additional) >= minimum.
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


# Cutoff cache types and stores
CutoffsCacheResult = Optional[Tuple[int, int]]
CutoffsCache = MutableMapping[Tuple[ContigId, ContigId], CutoffsCacheResult]
CUTOFFS_CACHE: MutableMapping[Tuple[ContigId, ContigId], Tuple[int, int]] = {}
CUTOFFS_NEGATIVE: Set[Tuple[ContigId, ContigId]] = set()


def precheck_and_prepare_overlap(
    a: ContigWithAligner,
    b: ContigWithAligner,
    minimum_base_score: Score,
) -> Optional[Tuple[ContigWithAligner, ContigWithAligner, int, str, str, Overlap]]:
    """
    Run early-bound checks and prepare normalized overlap windows.

    Returns (left, right, shift, left_initial_overlap, right_initial_overlap, overlap)
    or None if any precondition fails. No events are emitted here.
    """

    # Quick upper bound check: if even a perfect overlap cannot reach the minimum, abort early.
    if max_possible_overlap_score(len(a.seq), len(b.seq)) < minimum_base_score:
        return None

    # Identify best overlap placement.
    overlap = get_overlap(a, b)
    if overlap is None:
        return None

    # Optimistic check using the discovered overlap window size.
    if optimistic_overlap_score(overlap.size) < minimum_base_score:
        return None

    # Normalize orientation and derive initial overlap windows.
    left, right, shift = normalize_orientation(a, b, overlap)
    left_initial_overlap, right_initial_overlap = initial_overlap_windows(
        left, right, shift
    )

    # Maintain safety assertions from original implementation.
    assert len(right_initial_overlap) == overlap.size, (
        f"{len(right_initial_overlap)} == {overlap.size}"
    )
    assert len(left_initial_overlap) == overlap.size, (
        f"{len(left_initial_overlap)} == {overlap.size}"
    )
    assert (
        calculate_referenceless_overlap_score(
            L=len(left_initial_overlap) + 1, M=len(left_initial_overlap)
        )
        >= minimum_base_score
    )

    return left, right, shift, left_initial_overlap, right_initial_overlap, overlap


def coverage_flags(
    left: ContigWithAligner,
    right: ContigWithAligner,
    overlap_size: int,
) -> Tuple[bool, Optional[ContigWithAligner], Optional[ContigWithAligner]]:
    """Compute coverage flags and identify covered/bigger contigs if applicable."""
    left_is_covered = len(left.seq) <= overlap_size
    right_is_covered = len(right.seq) <= overlap_size
    is_covered = left_is_covered or right_is_covered
    if is_covered:
        covered = left if left_is_covered else right
        bigger = right if left_is_covered else left
        return True, covered, bigger
    return False, None, None


def compute_alignment_and_score(
    is_covered: bool,
    covered: Optional[ContigWithAligner],
    bigger: Optional[ContigWithAligner],
    left: ContigWithAligner,
    right: ContigWithAligner,
    left_cutoff: int,
    right_cutoff: int,
) -> Tuple[str, str, str, str, int, int, Score]:
    """Align according to covered/non-covered case and compute overlap score.

    Returns (aligned_1, aligned_2, left_remainder, right_remainder,
             result_length, number_of_matches, result_score).
    """
    if is_covered:
        assert covered is not None and bigger is not None
        aligned_1, aligned_2 = align_for_merge_covered(
            covered, bigger, left_cutoff, right_cutoff
        )
        left_remainder = right_remainder = ""
    else:
        aligned_1, aligned_2, left_remainder, right_remainder = (
            align_for_merge_noncovered(left, right, left_cutoff, right_cutoff)
        )

    result_length, number_of_matches, result_score = score_alignment(
        aligned_1, aligned_2
    )
    return (
        aligned_1,
        aligned_2,
        left_remainder,
        right_remainder,
        result_length,
        number_of_matches,
        result_score,
    )


# Cutoff computation helpers
def cutoffs_left_covered(
    minimum_score: Score,
    left: ContigWithAligner,
    right: ContigWithAligner,
    shift: int,
    left_initial_overlap: str,
) -> Tuple[int, int]:
    overlap_alignments = tuple(
        right.map_overlap(minimum_score, "cover", left_initial_overlap)
    )
    right_cutoff = max((end for start, end in overlap_alignments), default=-1)
    if right_cutoff < 0:
        return (
            len(right.seq) - abs(shift) + 1,
            len(right.seq) - abs(shift) + len(left_initial_overlap) + 1,
        )
    left_cutoff = min((start for start, end in overlap_alignments), default=-1)
    return left_cutoff, right_cutoff


def cutoffs_right_covered(
    minimum_score: Score,
    left: ContigWithAligner,
    right: ContigWithAligner,
    shift: int,
    right_initial_overlap: str,
) -> Tuple[int, int]:
    overlap_alignments = tuple(
        left.map_overlap(minimum_score, "cover", right_initial_overlap)
    )
    left_cutoff = min((start for start, end in overlap_alignments), default=-1)
    if left_cutoff < 0:
        return (
            len(left.seq) - abs(shift) + 1,
            len(left.seq) - abs(shift) + len(right_initial_overlap) + 1,
        )
    right_cutoff = max((end for start, end in overlap_alignments), default=-1)
    return left_cutoff, right_cutoff


def cutoffs_left_shorter(
    minimum_score: Score,
    left: ContigWithAligner,
    right: ContigWithAligner,
    left_initial_overlap: str,
    right_initial_overlap: str,
) -> Optional[Tuple[int, int]]:
    left_overlap_alignments = left.map_overlap(
        minimum_score, "left", right_initial_overlap
    )
    left_cutoff = min((start for start, end in left_overlap_alignments), default=-1)
    if left_cutoff < 0:
        return None
    right_overlap_alignments = right.map_overlap(
        minimum_score, "right", left_initial_overlap
    )
    right_cutoff = max((end for start, end in right_overlap_alignments), default=-1)
    if right_cutoff < 0:
        return None
    return left_cutoff, right_cutoff


def cutoffs_right_shorter_or_equal(
    minimum_score: Score,
    left: ContigWithAligner,
    right: ContigWithAligner,
    left_initial_overlap: str,
    right_initial_overlap: str,
) -> Optional[Tuple[int, int]]:
    right_overlap_alignments = right.map_overlap(
        minimum_score, "right", left_initial_overlap
    )
    right_cutoff = max((end for start, end in right_overlap_alignments), default=-1)
    if right_cutoff < 0:
        return None
    left_overlap_alignments = left.map_overlap(
        minimum_score, "left", right_initial_overlap
    )
    left_cutoff = min((start for start, end in left_overlap_alignments), default=-1)
    if left_cutoff < 0:
        return None
    return left_cutoff, right_cutoff


def compute_overlap_cutoffs(
    minimum_score: Score,
    left: ContigWithAligner,
    right: ContigWithAligner,
    shift: int,
    left_initial_overlap: str,
    right_initial_overlap: str,
) -> Optional[Tuple[int, int]]:
    if len(left.seq) == len(left_initial_overlap):
        return cutoffs_left_covered(
            minimum_score, left, right, shift, left_initial_overlap
        )
    if len(right.seq) == len(right_initial_overlap):
        return cutoffs_right_covered(
            minimum_score, left, right, shift, right_initial_overlap
        )
    if len(left.seq) < len(right.seq):
        return cutoffs_left_shorter(
            minimum_score, left, right, left_initial_overlap, right_initial_overlap
        )
    return cutoffs_right_shorter_or_equal(
        minimum_score, left, right, left_initial_overlap, right_initial_overlap
    )


def find_overlap_cutoffs(
    minimum_score: Score,
    left: ContigWithAligner,
    right: ContigWithAligner,
    shift: int,
    left_initial_overlap: str,
    right_initial_overlap: str,
) -> CutoffsCacheResult:
    """Locate cutoffs delimiting the high-confidence overlap region.

    The function determines indices into the original left/right contig sequences
    that bound the overlap to be aligned for scoring and merging. Cutoffs are
    chosen via `map_overlap` using the provided `minimum_score`.

    Caching:
        Results are cached per (left.id, right.id). When no valid overlap region
        is possible for the given pair (under `minimum_score`), a negative entry
        is recorded and None is returned on subsequent calls.

    Returns:
        Tuple (left_cutoff, right_cutoff) on success, or None if no valid
        overlap region satisfies the minimum score.
    """
    key = (left.id, right.id)

    if key in CUTOFFS_NEGATIVE:
        return None
    existing = CUTOFFS_CACHE.get(key)
    if existing is not None:
        return existing

    value = compute_overlap_cutoffs(
        minimum_score, left, right, shift, left_initial_overlap, right_initial_overlap
    )
    if value is None:
        CUTOFFS_NEGATIVE.add(key)
    else:
        CUTOFFS_CACHE[key] = value
    return value


def normalize_orientation(
    a: ContigWithAligner, b: ContigWithAligner, overlap: Overlap
) -> Tuple[ContigWithAligner, ContigWithAligner, int]:
    """Return (left, right, shift) ensuring `shift` corresponds to left->right."""
    if abs(overlap.shift) < len(a.seq):
        return a, b, overlap.shift
    shift = (len(b.seq) + len(a.seq) - abs(overlap.shift)) * -1
    return b, a, shift


def initial_overlap_windows(
    left: ContigWithAligner, right: ContigWithAligner, shift: int
) -> Tuple[str, str]:
    """Extract initial overlap windows on left/right given a placement shift."""
    left_initial = left.seq[
        len(left.seq) - abs(shift) : (len(left.seq) - abs(shift) + len(right.seq))
    ]
    right_initial = right.seq[: abs(shift)]
    return left_initial, right_initial


def max_possible_overlap_score(len_a: int, len_b: int) -> Score:
    """Upper bound score for a perfect overlap between two contigs.

    Uses L = min(len_a, len_b) - 1 and M = L.
    """
    maximum_overlap_size = min(len_a, len_b) - 1
    maximum_number_of_matches = maximum_overlap_size
    return calculate_referenceless_overlap_score(
        L=maximum_overlap_size + 1,
        M=maximum_number_of_matches,
    )


def optimistic_overlap_score(overlap_size: int) -> Score:
    """Optimistic score using discovered overlap window size (assume M = L)."""
    return calculate_referenceless_overlap_score(L=overlap_size + 1, M=overlap_size)


def score_alignment(aligned_1: str, aligned_2: str) -> Tuple[int, int, Score]:
    """Compute (result_length, number_of_matches, result_score) from an alignment."""
    result_length = len(aligned_1)
    assert result_length > 0, "The overlap cannot be empty."
    number_of_matches = sum(
        1 for x, y in zip(aligned_1, aligned_2) if x == y and x != "-"
    )
    result_score = calculate_referenceless_overlap_score(
        L=result_length, M=number_of_matches
    )
    return result_length, number_of_matches, result_score


def align_for_merge_covered(
    covered: ContigWithAligner,
    bigger: ContigWithAligner,
    left_cutoff: int,
    right_cutoff: int,
) -> Tuple[str, str]:
    """Align sequences for the covered-contig case."""
    covered_overlap = covered.seq
    bigger_overlap = bigger.seq[left_cutoff:right_cutoff]
    return align_overlaps(covered_overlap, bigger_overlap)


def align_for_merge_noncovered(
    left: ContigWithAligner,
    right: ContigWithAligner,
    left_cutoff: int,
    right_cutoff: int,
) -> Tuple[str, str, str, str]:
    """Align overlap windows and return alignment plus remainders for merging."""
    left_overlap = left.seq[left_cutoff : (left_cutoff + len(right.seq))]
    left_remainder = left.seq[:left_cutoff]
    right_overlap = right.seq[:right_cutoff]
    right_remainder = right.seq[right_cutoff:]
    aligned_1, aligned_2 = align_overlaps(left_overlap, right_overlap)
    return aligned_1, aligned_2, left_remainder, right_remainder


def merge_by_concordance(
    aligned_1: str, aligned_2: str, left_remainder: str, right_remainder: str
) -> Tuple[str, int]:
    """Produce merged sequence by splitting at the best concordance index.

    Returns tuple of (result_seq, overlap_size_for_event).
    """
    concordance = calculate_concordance(aligned_1, aligned_2)
    max_concordance_index = next(iter(sort_concordance_indexes(concordance)))
    left_overlap_chunk = "".join(
        x for x in aligned_1[:max_concordance_index] if x != "-"
    )
    right_overlap_chunk = "".join(
        x for x in aligned_2[max_concordance_index:] if x != "-"
    )

    result_seq = (
        left_remainder + left_overlap_chunk + right_overlap_chunk + right_remainder
    )
    overlap_size = len(left_overlap_chunk) + len(right_overlap_chunk)
    return result_seq, overlap_size


def try_combine_contigs(
    is_debug2: bool,
    current_score: Score,
    pool: Pool,
    a: ContigWithAligner,
    b: ContigWithAligner,
) -> Optional[Tuple[ContigWithAligner, Score]]:
    """Attempt to combine two contigs respecting the pool's score threshold.

    Fast-fails using upper bounds on achievable overlap, normalizes orientation,
    determines precise cutoffs for the overlap region, aligns the overlaps, and
    either returns a merged contig and its score or indicates that merging is
    not beneficial under the current threshold.
    """

    # Compute the minimum additional score needed for acceptance.
    minimum_base_score = get_minimum_base_score(
        current_score, pool.min_acceptable_score
    )

    # Prechecks and preparation (no events emitted here).
    prepared = precheck_and_prepare_overlap(a, b, minimum_base_score)
    if prepared is None:
        return None
    left, right, shift, left_initial_overlap, right_initial_overlap, overlap = prepared

    # Determine refined cutoffs within the initial windows.
    cutoffs = find_overlap_cutoffs(
        minimum_base_score,
        left,
        right,
        shift,
        left_initial_overlap,
        right_initial_overlap,
    )

    if is_debug2:
        log(
            events.CalculatedCutoffs(
                left.unique_name, right.unique_name, len(left_initial_overlap), cutoffs
            )
        )

    if cutoffs is None:
        return None

    left_cutoff, right_cutoff = cutoffs

    # Covered-contig short-circuit: if one contig is fully covered by the overlap.
    is_covered, covered, bigger = coverage_flags(left, right, overlap.size)

    (
        aligned_1,
        aligned_2,
        left_remainder,
        right_remainder,
        result_length,
        number_of_matches,
        result_score,
    ) = compute_alignment_and_score(
        is_covered, covered, bigger, left, right, left_cutoff, right_cutoff
    )

    if is_debug2:
        denominator = max(minimum_base_score, ACCEPTABLE_STITCHING_SCORE())
        relative_score = (
            result_score / denominator if denominator != 0 else float("inf")
        )
        log(
            events.DeterminedOverlap(
                left.unique_name,
                right.unique_name,
                result_length,
                number_of_matches,
                relative_score,
            )
        )

    if result_score < minimum_base_score:
        return None

    if is_covered:
        assert bigger is not None
        if is_debug2:
            log(events.Covered(left.unique_name, right.unique_name))
        return (bigger, SCORE_EPSILON)

    # Merge by concordance position to produce the combined contig.
    result_seq, overlap_size = merge_by_concordance(
        aligned_1, aligned_2, left_remainder, right_remainder
    )
    result_contig = ContigWithAligner(None, result_seq)

    if is_debug2:
        log(
            events.CombinedContings(
                left_contig=left.unique_name,
                right_contig=right.unique_name,
                result_contig=result_contig.unique_name,
                overlap_size=overlap_size,
            )
        )
    return (result_contig, result_score)


def extend_by_1(
    is_debug2: bool,
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

    combination = try_combine_contigs(
        is_debug2, path.score, pool, path.whole, candidate
    )
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


def calc_extension(
    is_debug2: bool,
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


def calc_multiple_extensions(
    is_debug2: bool,
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


def calculate_all_paths(
    paths: Sequence[ContigsPath],
    contigs: Sequence[ContigWithAligner],
) -> Iterable[ContigsPath]:
    """
    Iteratively extend seed paths with contigs to generate candidate contig paths.
    Returns an iterable of best paths up to a max number of alternatives.
    """
    is_debug2 = ReferencelessStitcherContext.get().is_debug2

    max_alternatives = intrapolate_number_of_alternatives(len(contigs))
    pool = Pool.empty(max_alternatives, ACCEPTABLE_STITCHING_SCORE())
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


def find_most_probable_path(
    seeds: Sequence[ContigsPath],
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


def stitch_consensus_overlaps(
    contigs: Iterable[ContigWithAligner],
) -> Iterator[ContigWithAligner]:
    """
    Produce a consensus stitching of contigs based on overlaps.
    Iteratively selects and stitches the most probable contig paths.
    """
    remaining = tuple(sorted(contigs, key=contig_size_fun))

    log(events.InitializingSeeds())
    seeds = tuple(
        sorted(
            map(ContigsPath.singleton, remaining),
            key=lambda path: len(path.whole.seq),
            reverse=True,
        )
    )
    log(events.Starting(len(seeds)))

    while remaining:
        most_probable = find_most_probable_path(seeds, remaining)
        log(events.Constructed(most_probable))
        yield most_probable.whole
        remaining = tuple(
            contig for contig in remaining if not most_probable.has_contig(contig)
        )
        seeds = tuple(
            path
            for path in seeds
            if most_probable.contigs_ids.isdisjoint(path.contigs_ids)
        )
        log(events.Remove(len(most_probable.contigs_ids), len(remaining)))
        if len(most_probable.contains_contigs_ids) == 1:
            log(events.GiveUp())
            yield from remaining
            return


def try_combine_1(
    contigs: Iterable[ContigWithAligner],
) -> Optional[Tuple[ContigWithAligner, ContigWithAligner, ContigWithAligner]]:
    is_debug2 = ReferencelessStitcherContext.get().is_debug2

    for first in contigs:
        for second in contigs:
            if first.id >= second.id:
                continue

            pool = Pool.empty(1, ACCEPTABLE_STITCHING_SCORE())
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


def o2_loop(
    contigs: Iterable[ContigWithAligner],
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


def stitch_consensus(
    contigs: Iterable[ContigWithAligner],
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
    records = (
        SeqRecord(
            Seq.Seq(contig.seq),
            description="",
            id=contig.unique_name,
            name=contig.unique_name,
        )
        for contig in contigs
    )
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


def referenceless_contig_stitcher(
    input_fasta: TextIO,
    output_fasta: Optional[TextIO],
) -> int:
    """
    Wrapper that initializes a fresh stitching context and calls the core stitching function.
    """
    with ReferencelessStitcherContext.fresh() as ctx:
        if logger.level == logging.DEBUG - 1:
            ctx.is_debug2 = True
        return referenceless_contig_stitcher_with_ctx(input_fasta, output_fasta)
