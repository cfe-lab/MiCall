import itertools
import logging
from functools import cache
from typing import (
    Iterable,
    Iterator,
    Optional,
    Tuple,
    Sequence,
    TextIO,
    MutableMapping,
    AbstractSet,
)

from Bio import Seq, SeqIO
from Bio.SeqRecord import SeqRecord

from micall.utils.contig_stitcher_context import ReferencelessStitcherContext
from micall.utils.contig_stitcher_contigs import Contig, ContigId
from micall.utils.overlap_stitcher import (
    align_queries,
    calculate_concordance,
    sort_concordance_indexes,
    calculate_overlap_score,
)
from .referenceless_contig_stitcher_pool import Pool
import micall.utils.referenceless_contig_stitcher_events as events
from micall.utils.referenceless_contig_path import ContigsPath
from micall.utils.referenceless_contig_with_aligner import ContigWithAligner, find_maximum_overlap, map_overlap
from micall.utils.referenceless_score import Score, SCORE_EPSILON, SCORE_NOTHING
from micall.utils.referenceless_contig_stitcher_overlap import Overlap


# Kmer size for filtering non-overlapping contigs
KMER_SIZE = 50

# Minimum number of matches to consider an overlap acceptable
MIN_MATCHES = 99


logger = logging.getLogger(__name__)


def cache_kmers(cache: MutableMapping[str, AbstractSet[str]], contig_sequence: str) -> AbstractSet[str]:
    """
    Extract all k-mers from a contig sequence.

    This is used as a fast filter to quickly reject contig pairs that
    cannot possibly overlap. If two contigs don't share any k-mers,
    they are considered to not overlap.

    Args:
        cache: A mutable mapping used to cache previously computed k-mer sets.
        contig_sequence: The DNA sequence of the contig

    Returns:
        A set of all k-mers of size _KMER_SIZE found in the sequence.
        If the sequence is shorter than _KMER_SIZE, returns an empty set.
    """

    existing = cache.get(contig_sequence)
    if existing is not None:
        return existing

    kmer_size = KMER_SIZE
    if len(contig_sequence) < kmer_size:
        return set()

    kmers = set()
    for i in range(len(contig_sequence) - kmer_size + 1):
        kmer = contig_sequence[i:i + kmer_size]
        kmers.add(kmer)

    return kmers


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


@cache
def ACCEPTABLE_BASE_STITCHING_SCORE():
    return calculate_overlap_score(L=MIN_MATCHES + 1, M=MIN_MATCHES)


@cache
def ACCEPTABLE_STITCHING_SCORE():
    return calculate_referenceless_overlap_score(L=MIN_MATCHES + 1, M=MIN_MATCHES)


MAX_ALTERNATIVES = 999


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

    Finally, the result is clamped into [1, MAX_ALTERNATIVES]
    and rounded to the nearest integer.

    Args:
        n_candidates: Number of contigs available for extension
                      in the current cycle.

    Returns:
        The (rounded and clamped) maximum number of alternative
        contig-paths to consider for this cycle.
    """

    ret = MAX_ALTERNATIVES / max(1, n_candidates - 2)
    clamped = max(1, min(MAX_ALTERNATIVES, ret))
    rounded = round(clamped)
    return rounded


def log(e: events.EventType) -> None:
    """
    Emit an event to the current stitching context and log it at debug level.
    """
    ReferencelessStitcherContext.get().emit(e)
    logger.debug("%s", e)


def compute_overlap_size(left_len: int, right_len: int, shift: int) -> int:
    """Compute the overlap size given contig lengths and a shift.

    The shift is defined as returned by `find_maximum_overlap`.
    """
    if abs(shift) <= left_len:
        size = min(abs(shift), right_len)
    else:
        size = min(left_len, left_len + right_len + shift)
    assert size > 0, f"{shift}, {left_len}, {right_len}"
    assert size <= left_len, f"{shift}, {size}, {left_len}, {right_len}"
    assert size <= right_len, f"{shift}, {size}, {left_len}, {right_len}"
    return size


def get_overlap(
    left: ContigWithAligner,
    right: ContigWithAligner,
    get_overlap_cache: MutableMapping[Tuple[ContigId, ContigId], Optional[Overlap]],
    kmers_cache: MutableMapping[str, AbstractSet[str]],
) -> Optional[Overlap]:
    """Return the best overlap placement between two contigs, if any.

    Uses a provided cache to avoid repeated expensive overlap calculations.

    Returns:
        Overlap(shift, size) if an overlap was found; otherwise None.
    """
    if len(left.seq) == 0 or len(right.seq) == 0:
        return None

    key = (left.id, right.id)
    existing = get_overlap_cache.get(key, -1)
    if existing != -1:
        ret: Optional[Overlap] = existing # type: ignore[assignment]
        return ret

    # Fast k-mer filter: if two contigs don't share any k-mers, they cannot overlap
    # Only apply the filter if both contigs are long enough (at least kmer_size)
    # This prevents false negatives when the overlap region is smaller than kmer_size
    min_length_for_filter = KMER_SIZE
    if len(left.seq) >= min_length_for_filter and len(right.seq) >= min_length_for_filter:
        left_kmers = cache_kmers(kmers_cache, left.seq)
        right_kmers = cache_kmers(kmers_cache, right.seq)
        # We only filter if both have kmers AND they don't share any
        if left_kmers and right_kmers and not (left_kmers & right_kmers):
            # No shared k-mers - so no significant overlap
            get_overlap_cache[key] = None
            return None

    shift, _ = find_maximum_overlap(left, right)
    if shift == 0:
        get_overlap_cache[key] = None
        return None

    size = compute_overlap_size(len(left.seq), len(right.seq), shift)
    ret = Overlap(shift=shift, size=size)
    get_overlap_cache[key] = ret
    return ret


def combine_scores(current: Score, new: Score) -> Score:
    return current + new


def get_minimum_base_score(current: Score, minimum: Score) -> Score:
    """
    Calculate the minimum additional score required so that
    combine_scores(current, additional) >= minimum.
    """
    return minimum - current


def align_overlaps(
    left_overlap: str,
    right_overlap: str,
    align_cache: MutableMapping[Tuple[str, str], Tuple[str, str]],
) -> Tuple[str, str]:
    key = (left_overlap, right_overlap)
    existing = align_cache.get(key)
    if existing is not None:
        return existing

    result = align_queries(left_overlap, right_overlap)
    align_cache[key] = result
    return result


# Cutoff cache types
CutoffsCacheResult = Optional[Tuple[int, int]]
CutoffsCache = MutableMapping[Tuple[ContigId, ContigId], CutoffsCacheResult]


def precheck_and_prepare_overlap(
    a: ContigWithAligner,
    b: ContigWithAligner,
    minimum_base_score: Score,
    get_overlap_cache: MutableMapping[Tuple[ContigId, ContigId], Optional[Overlap]],
    kmers_cache: MutableMapping[str, AbstractSet[str]],
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
    overlap = get_overlap(a, b, get_overlap_cache, kmers_cache)
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

    # Sanity checks.
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


def calculate_covered(
    left: ContigWithAligner,
    right: ContigWithAligner,
    overlap_size: int,
) -> Optional[Tuple[ContigWithAligner, ContigWithAligner]]:
    """Compute coverage flags and identify covered/bigger contigs if applicable."""
    left_is_covered = len(left.seq) <= overlap_size
    right_is_covered = len(right.seq) <= overlap_size
    is_covered = left_is_covered or right_is_covered
    if is_covered:
        covered = left if left_is_covered else right
        bigger = right if left_is_covered else left
        return covered, bigger
    return None


def compute_alignment_and_score(
    covered: Optional[ContigWithAligner],
    bigger: Optional[ContigWithAligner],
    left: ContigWithAligner,
    right: ContigWithAligner,
    left_cutoff: int,
    right_cutoff: int,
    align_cache: MutableMapping[Tuple[str, str], Tuple[str, str]],
) -> Tuple[str, str, str, str, int, int, Score]:
    """Align according to covered/non-covered case and compute overlap score.

    Returns (aligned_1, aligned_2, left_remainder, right_remainder,
             result_length, number_of_matches, result_score).
    """

    is_covered = covered is not None
    if is_covered:
        assert covered is not None and bigger is not None
        aligned_1, aligned_2 = align_for_merge_covered(
            covered, bigger, left_cutoff, right_cutoff, align_cache
        )
        left_remainder = right_remainder = ""
    else:
        aligned_1, aligned_2, left_remainder, right_remainder = (
            align_for_merge_noncovered(left, right, left_cutoff, right_cutoff, align_cache)
        )

    # Apply +1 length bonus for non-covered overlaps and +2 for covered overlaps
    result_length, number_of_matches, result_score = score_alignment(
        aligned_1, aligned_2, is_covered,
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
    minimum_acceptable: Score,
    left: ContigWithAligner,
    right: ContigWithAligner,
    shift: int,
    left_initial_overlap: str,
) -> Tuple[int, int]:
    overlap_alignments = tuple(
        map_overlap(right, minimum_acceptable, "cover", left_initial_overlap)
    )
    right_cutoff = max((end for start, end in overlap_alignments), default=-1)
    if right_cutoff < 0:
        return (
            len(right.seq) - abs(shift),
            len(right.seq) - abs(shift) + len(left_initial_overlap),
        )
    left_cutoff = min((start for start, end in overlap_alignments), default=-1)
    return left_cutoff, right_cutoff


def cutoffs_right_covered(
    minimum_acceptable: Score,
    left: ContigWithAligner,
    right: ContigWithAligner,
    shift: int,
    right_initial_overlap: str,
) -> Tuple[int, int]:
    overlap_alignments = tuple(
        map_overlap(left, minimum_acceptable, "cover", right_initial_overlap)
    )
    left_cutoff = min((start for start, end in overlap_alignments), default=-1)
    if left_cutoff < 0:
        return (
            len(left.seq) - abs(shift),
            len(left.seq) - abs(shift) + len(right_initial_overlap),
        )
    right_cutoff = max((end for start, end in overlap_alignments), default=-1)
    return left_cutoff, right_cutoff


def cutoffs_left_shorter(
    minimum_acceptable: Score,
    left: ContigWithAligner,
    right: ContigWithAligner,
    left_initial_overlap: str,
    right_initial_overlap: str,
) -> Optional[Tuple[int, int]]:
    left_overlap_alignments = map_overlap(
        left, minimum_acceptable, "left", right_initial_overlap
    )
    left_cutoff = min((start for start, end in left_overlap_alignments), default=-1)
    if left_cutoff < 0:
        return None
    right_overlap_alignments = map_overlap(
        right, minimum_acceptable, "right", left_initial_overlap
    )
    right_cutoff = max((end for start, end in right_overlap_alignments), default=-1)
    if right_cutoff < 0:
        return None
    return left_cutoff, right_cutoff


def cutoffs_right_shorter_or_equal(
    minimum_acceptable: Score,
    left: ContigWithAligner,
    right: ContigWithAligner,
    left_initial_overlap: str,
    right_initial_overlap: str,
) -> Optional[Tuple[int, int]]:
    right_overlap_alignments = map_overlap(
        right, minimum_acceptable, "right", left_initial_overlap
    )
    right_cutoff = max((end for start, end in right_overlap_alignments), default=-1)
    if right_cutoff < 0:
        return None
    left_overlap_alignments = map_overlap(
        left, minimum_acceptable, "left", right_initial_overlap
    )
    left_cutoff = min((start for start, end in left_overlap_alignments), default=-1)
    if left_cutoff < 0:
        return None
    return left_cutoff, right_cutoff


def compute_overlap_cutoffs(
    minimum_acceptable: Score,
    left: ContigWithAligner,
    right: ContigWithAligner,
    shift: int,
    left_initial_overlap: str,
    right_initial_overlap: str,
) -> Optional[Tuple[int, int]]:
    if len(left.seq) == len(left_initial_overlap):
        return cutoffs_left_covered(
            minimum_acceptable, left, right, shift, left_initial_overlap
        )
    if len(right.seq) == len(right_initial_overlap):
        return cutoffs_right_covered(
            minimum_acceptable, left, right, shift, right_initial_overlap
        )
    if len(left.seq) < len(right.seq):
        return cutoffs_left_shorter(
            minimum_acceptable, left, right, left_initial_overlap, right_initial_overlap
        )
    return cutoffs_right_shorter_or_equal(
        minimum_acceptable, left, right, left_initial_overlap, right_initial_overlap
    )


def find_overlap_cutoffs(
    minimum_acceptable: Score,
    left: ContigWithAligner,
    right: ContigWithAligner,
    shift: int,
    left_initial_overlap: str,
    right_initial_overlap: str,
    cutoffs_cache: MutableMapping[Tuple[ContigId, ContigId], Optional[Tuple[int, int]]],
) -> CutoffsCacheResult:
    """Locate cutoffs delimiting the high-confidence overlap region.

    The function determines indices into the original left/right contig sequences
    that bound the overlap to be aligned for scoring and merging. Cutoffs are
    chosen via `map_overlap`.

    Caching:
        Results are cached per (left.id, right.id). When no valid overlap region
        is possible for the given pair, a None entry is
        recorded and None is returned on subsequent calls.

    Returns:
        Tuple (left_cutoff, right_cutoff) on success, or None if no valid
        overlap region satisfies the minimum acceptable score or something lower than it.
    """

    # Note:
    # It is fine to omit `minimum_acceptable` from the cache key because
    # the cutoffs are monotonic with respect to `minimum_acceptable`.
    # Increasing `minimum_acceptable` can only reduce the valid overlap region,
    # never expand it. Therefore, the cutoffs computed for a lower
    # `minimum_acceptable` are always valid for a higher one.
    # The `minimum_acceptable` is monotonic.
    key = (left.id, right.id)
    existing = cutoffs_cache.get(key, -1)
    if existing != -1:
        ret: CutoffsCacheResult = existing  # type: ignore[assignment]
        return ret

    value = compute_overlap_cutoffs(
        minimum_acceptable, left, right, shift, left_initial_overlap, right_initial_overlap
    )
    cutoffs_cache[key] = value
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


def score_alignment(aligned_1: str, aligned_2: str, is_covered: bool) -> Tuple[int, int, Score]:
    """Compute (result_length, number_of_matches, result_score) from an alignment.

    The statistical length L used for scoring is result_length + length_bonus.
    Use length_bonus = 1 for non-covering overlaps and 2 for covering overlaps.
    """
    result_length = len(aligned_1)
    assert result_length > 0, "The overlap cannot be empty."
    number_of_matches = sum(
        1 for x, y in zip(aligned_1, aligned_2) if x == y and x != "-"
    )

    # We add a length bonus because the statistical model assumes
    # that the overlap is flanked by non-matching positions on both sides.
    length_bonus = 2 if is_covered else 1
    result_score = calculate_referenceless_overlap_score(
        L=result_length + length_bonus, M=number_of_matches
    )

    return result_length, number_of_matches, result_score


def align_for_merge_covered(
    covered: ContigWithAligner,
    bigger: ContigWithAligner,
    left_cutoff: int,
    right_cutoff: int,
    align_cache: MutableMapping[Tuple[str, str], Tuple[str, str]],
) -> Tuple[str, str]:
    """Align sequences for the covered-contig case."""
    covered_overlap = covered.seq
    bigger_overlap = bigger.seq[left_cutoff:right_cutoff]
    return align_overlaps(covered_overlap, bigger_overlap, align_cache)


def align_for_merge_noncovered(
    left: ContigWithAligner,
    right: ContigWithAligner,
    left_cutoff: int,
    right_cutoff: int,
    align_cache: MutableMapping[Tuple[str, str], Tuple[str, str]],
) -> Tuple[str, str, str, str]:
    """Align overlap windows and return alignment plus remainders for merging."""
    left_overlap = left.seq[left_cutoff : left_cutoff + right_cutoff]
    left_remainder = left.seq[:left_cutoff]
    right_overlap = right.seq[:right_cutoff]
    right_remainder = right.seq[right_cutoff:]
    aligned_1, aligned_2 = align_overlaps(left_overlap, right_overlap, align_cache)
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
    get_overlap_cache: MutableMapping[Tuple[ContigId, ContigId], Optional[Overlap]],
    kmers_cache: MutableMapping[str, AbstractSet[str]],
    cutoffs_cache: MutableMapping[Tuple[ContigId, ContigId], Optional[Tuple[int, int]]],
    align_cache: MutableMapping[Tuple[str, str], Tuple[str, str]],
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

    prepared = precheck_and_prepare_overlap(a, b, minimum_base_score, get_overlap_cache, kmers_cache)
    if prepared is None:
        return None
    left, right, shift, left_initial_overlap, right_initial_overlap, overlap = prepared

    cutoffs = find_overlap_cutoffs(
        pool.min_acceptable_score,
        left,
        right,
        shift,
        left_initial_overlap,
        right_initial_overlap,
        cutoffs_cache,
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

    coverage = calculate_covered(left, right, overlap.size)
    if coverage is None:
        covered = bigger = None
    else:
        covered, bigger = coverage

    (
        aligned_1,
        aligned_2,
        left_remainder,
        right_remainder,
        result_length,
        number_of_matches,
        result_score,
    ) = compute_alignment_and_score(
        covered, bigger, left, right, left_cutoff, right_cutoff, align_cache
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

    if covered is not None:
        # Check if the coverage is mutual (i.e., both contigs are fully covered).
        is_both_covered = left_cutoff == 0 and right_cutoff >= len(right.seq)
        if is_both_covered:
            # In case of mutual coverage, keep both contigs.
            if left.seq == right.seq:
                # But if they are identical, keep only one to avoid redundancy.
                if is_debug2:
                    log(events.Covered(left.unique_name, right.unique_name))
                return (left, SCORE_EPSILON)

            return None

        if is_debug2:
            log(events.Covered(left.unique_name, right.unique_name))

        assert bigger is not None
        return (bigger, SCORE_EPSILON)

    result_seq, overlap_size = merge_by_concordance(
        aligned_1, aligned_2, left_remainder, right_remainder
    )
    # When merging contigs, sum their read counts if both are available
    if left.reads_count is not None and right.reads_count is not None:
        combined_reads_count = left.reads_count + right.reads_count
    else:
        combined_reads_count = None
    result_contig = ContigWithAligner(None, result_seq, combined_reads_count)

    if is_debug2:
        log(
            events.CombinedContigs(
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
    get_overlap_cache: MutableMapping[Tuple[ContigId, ContigId], Optional[Overlap]],
    kmers_cache: MutableMapping[str, AbstractSet[str]],
    cutoffs_cache: MutableMapping[Tuple[ContigId, ContigId], Optional[Tuple[int, int]]],
    align_cache: MutableMapping[Tuple[str, str], Tuple[str, str]],
) -> Optional[ContigsPath]:
    """
    Attempt to extend a contig path by one candidate contig.
    If a valid combination is found, return it.
    Otherwise return None.
    """

    if path.has_contig(candidate):
        return None

    combination = try_combine_contigs(
        is_debug2, path.score, pool, path.whole, candidate,
        get_overlap_cache, kmers_cache, cutoffs_cache, align_cache
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
    get_overlap_cache: MutableMapping[Tuple[ContigId, ContigId], Optional[Overlap]],
    kmers_cache: MutableMapping[str, AbstractSet[str]],
    cutoffs_cache: MutableMapping[Tuple[ContigId, ContigId], Optional[Tuple[int, int]]],
    align_cache: MutableMapping[Tuple[str, str], Tuple[str, str]],
) -> bool:
    """
    Try to extend a single path with each contig in contigs.
    Return True if any extension was added to the pool.
    """

    ret = False
    for contig in contigs:
        new = extend_by_1(
            is_debug2, pool, path, contig,
            get_overlap_cache, kmers_cache, cutoffs_cache, align_cache
        )
        if new is not None:
            ret = pool.add(new) or ret

    return ret


def calc_multiple_extensions(
    is_debug2: bool,
    pool: Pool,
    paths: Iterable[ContigsPath],
    contigs: Sequence[ContigWithAligner],
    get_overlap_cache: MutableMapping[Tuple[ContigId, ContigId], Optional[Overlap]],
    kmers_cache: MutableMapping[str, AbstractSet[str]],
    cutoffs_cache: MutableMapping[Tuple[ContigId, ContigId], Optional[Tuple[int, int]]],
    align_cache: MutableMapping[Tuple[str, str], Tuple[str, str]],
) -> bool:
    """
    Attempt to extend multiple paths with a set of contigs.
    Returns True if any new extensions were added to the pool.
    """
    ret = False
    for path in paths:
        ret = calc_extension(
            is_debug2, pool, contigs, path,
            get_overlap_cache, kmers_cache, cutoffs_cache, align_cache,
        ) or ret
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
    ctx = ReferencelessStitcherContext.get()
    get_overlap_cache = ctx.get_overlap_cache
    kmers_cache = ctx.kmers_cache
    cutoffs_cache = ctx.cutoffs_cache
    align_cache = ctx.align_cache

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

        if not calc_multiple_extensions(
            is_debug2, pool, pool.ring, contigs,
            get_overlap_cache, kmers_cache, cutoffs_cache, align_cache,
        ):
            break

        # log end with updated pool size
        log(events.CycleEnd(cycle, len(pool.ring), pool))

    return pool.ring


def try_combine_1(
    contigs: Iterable[ContigWithAligner],
) -> Optional[Tuple[ContigWithAligner, ContigWithAligner, ContigWithAligner]]:
    is_debug2 = ReferencelessStitcherContext.get().is_debug2
    ctx = ReferencelessStitcherContext.get()
    get_overlap_cache = ctx.get_overlap_cache
    kmers_cache = ctx.kmers_cache
    cutoffs_cache = ctx.cutoffs_cache
    align_cache = ctx.align_cache

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
                get_overlap_cache=get_overlap_cache,
                kmers_cache=kmers_cache,
                cutoffs_cache=cutoffs_cache,
                align_cache=align_cache,
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
        yield ContigWithAligner(name=record.name, seq=str(record.seq), reads_count=None)


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
        if logger.getEffectiveLevel() == logging.DEBUG - 1:
            ctx.is_debug2 = True
        return referenceless_contig_stitcher_with_ctx(input_fasta, output_fasta)


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
