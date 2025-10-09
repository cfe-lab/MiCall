from fractions import Fraction
from typing import Sequence, Iterator, Tuple, TypeVar
from operator import itemgetter
import numpy as np
from micall.utils.referenceless_score import Score
from Bio.Align import PairwiseAligner


# set up the aligner
ALIGNER = PairwiseAligner()
ALIGNER.mode = "global"
# By default end gaps are free (score 0).
# Make them penalized by picking -1 per gap position at ends.
ALIGNER.end_gap_score = -1


def align_queries(seq1: str, seq2: str) -> Tuple[str, str]:
    """
    Globally align two query sequences against each other
    and return the resulting aligned sequences in MSA format.
    """

    # do the alignment, grab top hit
    aln = next(iter(ALIGNER.align(seq1, seq2)))

    # aln.aligned is a tuple of two lists of (start,end) pairs:
    #    aln.aligned[0] describes the aligned blocks on seq1
    #    aln.aligned[1] describes the aligned blocks on seq2
    blocks1, blocks2 = aln.aligned

    g1 = []
    g2 = []
    i = j = 0

    for (a0, a1), (b0, b1) in zip(blocks1, blocks2):
        # handle any gap in seq1 before this block
        if b0 > j:
            g1.append("-" * (b0 - j))
            g2.append(seq2[j:b0])
        # handle any gap in seq2 before this block
        if a0 > i:
            g1.append(seq1[i:a0])
            g2.append("-" * (a0 - i))

        # handle the matching/mismatching block
        g1.append(seq1[a0:a1])
        g2.append(seq2[b0:b1])

        i = a1
        j = b1

    # tail-end gaps
    if i < len(seq1):
        g1.append(seq1[i:])
        g2.append("-" * (len(seq1) - i))
    if j < len(seq2):
        g1.append("-" * (len(seq2) - j))
        g2.append(seq2[j:])

    aseq1 = "".join(g1)
    aseq2 = "".join(g2)
    return aseq1, aseq2


def normalize_array(array: Sequence[float]) -> Sequence[Fraction]:
    if len(array) == 0:
        return ()

    low = min(array)
    high = max(array)
    if low == high:
        if low > 0:
            value = 1
        else:
            value = 0
        return (Fraction(value),) * len(array)

    diff = Fraction(high - low)
    return tuple(Fraction(x - low) / diff for x in array)


def calculate_concordance_norm(left: Sequence[object], right: Sequence[object],
                               ) -> Sequence[Fraction]:
    absolute = calculate_concordance(left, right)
    return normalize_array(absolute)


def consecutive_true_counts(arr: np.ndarray) -> np.ndarray:
    """
    Given a boolean array arr, returns an integer array
    where out[i] is the count of consecutive True values
    up to (and including) index i, resetting on False.
    """

    c = np.add.accumulate(arr)                  # running count of 1s
    z = np.where(arr == 0, c, 0)                # counts at zeros
    m = np.maximum.accumulate(z)                # running max of that
    return c - m

    """
    Example trace:

    arr   = [0, 1, 1, 1, 0, 1, 1, 0]
    c     = [0, 1, 2, 3, 3, 4, 5, 5]     # running count of 1s
    z     = [0, 0, 0, 0, 3, 0, 0, 5]     # counts at zeros
    m     = [0, 0, 0, 0, 3, 3, 3, 5]     # running max of that
    c - m = [0, 1, 2, 3, 0, 1, 2, 0]

    """


def exp_accumulate_array_positive(array: np.ndarray) -> np.ndarray:
    forward = consecutive_true_counts(array)
    reverse = np.flip(consecutive_true_counts(np.flip(array)))
    return forward ** 0.5 + reverse ** 0.5


def exp_accumulate_array(array: np.ndarray) -> np.ndarray:
    positive = exp_accumulate_array_positive(array)
    negative = exp_accumulate_array_positive(1 - array)
    return positive - negative


def calculate_concordance(left: Sequence[object], right: Sequence[object],
                          ) -> Sequence[float]:
    """
    Calculate concordance for two given sequences using a sliding
    average.

    The function compares the two strings character by character,
    simultaneously from both left to right and right to left,
    calculating a score that represents a moving average of matches at
    each position. If characters match at a given position, a score of
    1 is added; otherwise, a score of 0 is added. This "sliding
    average" score is then processed again, but in reverse direction.
    The score is not normalized here.

    :param left: string representing first sequence
    :param right: string representing second sequence
    :return: list representing concordance score for each position
    """

    if len(left) != len(right):
        raise ValueError("Can only calculate concordance for same sized sequences")

    array = np.fromiter((x == y for x, y in zip(left, right)),
                        count=len(left), dtype=np.bool)
    return exp_accumulate_array(array)  # type: ignore


T = TypeVar("T")


def disambiguate_concordance(concordance: Sequence[T],
                             ) -> Iterator[Tuple[T, int]]:
    for i, x in enumerate(concordance):
        if i < len(concordance) / 2:
            global_rank = i
        else:
            global_rank = len(concordance) - i - 1
        yield x, global_rank


def sort_concordance_indexes(concordance: Sequence[object]) -> Iterator[int]:
    concordance_d = disambiguate_concordance(concordance)
    for i, v in sorted(enumerate(concordance_d),
                       key=itemgetter(1),
                       reverse=True,
                       ):
        yield i


def exp_dropoff_array_iter(array: np.ndarray,
                           direction: int,
                           factor: int = 2,
                           ) -> None:
    if direction > 0:
        iterator = range(len(array))
    else:
        iterator = range(len(array) - 1, -1, -1)

    previous = 0.0
    for i in iterator:
        current = array[i]
        dropoff = previous / factor
        if current < dropoff:
            array[i] = dropoff
            previous = dropoff
        else:
            previous = current


def exp_dropoff_array(array: np.ndarray, factor: int = 2) -> None:
    exp_dropoff_array_iter(array=array, direction=1, factor=factor)
    exp_dropoff_array_iter(array=array, direction=-1, factor=factor)


def calculate_overlap_score(L: int, M: int) -> Score:
    """
    Computes a monotonic scoring metric for an overlap between two
    sequences over a four-letter alphabet. It measures how much the
    observed match count M in an overlap of length L deviates from its
    expectation under a four-letter model, with a scaling exponent
    chosen to reflect the correlated nature of real genomic sequences.

    Derivation
    ----------
    1.  Uniform four-letter alphabet (match probability p = 1/4):
          - Expected matches:
                E[M] = L * p = L * 1/4
          - Independent-match variance:
                Var[M] = L * p * (1 - p) = L * 1/4 * (1 - 1/4) = L * (1/4 - 1/16)
          - Independent-match standard deviation:
                SD[M] = sqrt(Var[M]) = sqrt(L * (1/4 - 1/16))
          - Classic z-score for M:
                z = (M - E[M]) / SD[M] = (M - L * 1/4) / sqrt(L * (1/4 - 1/16))

    2.  Correlated-match model:
          - Real DNA (repeats, conserved motifs, low-complexity regions)
            exhibits long-range correlations so that
                SD[M] ∝ L^a
            with empirical estimates
                a ≈ 0.8.

    3.  Generalized overlap score:
          - Replace the sqrt(L) scaling with ^a, giving
                score = (M - L * 1/4) / ((L * (1/4 - 1/16))^a)

    4.  Optimizations:
          - Multiplying by L^-a avoids an explicit (slower) division:
                score = (M - L * 1/4) * ((L * (1/4 - 1/16))^-a)
          Constant factors rescale every score equally and do not affect ordering.
          - With that we multiply the first factor by 4 to eliminate division operation:
                score = (4 * M - L) * ((L * (1/4 - 1/16))^-a)
          - Expanding exponentiation we get another constant factor:
                score = (4 * M - L) * (L^-a * (1/4 - 1/16)^-a)
          - We eliminate the constant part:
                score = (4 * M - L) * (L^-a)

    This score preserves the monotonic "rarity" ordering of overlaps
    (higher ⇒ more unexpected), while penalizing long overlaps more
    strongly than the independent-match model.

    Parameters
    ----------
    L : int
        Length of the overlap (must be > 0).
    M : int
        Number of matching characters (0 <= M <= L).

    Returns
    -------
    Score
        A monotonic overlap score based on z-score.
    """

    alpha = -0.90
    return (4 * M - L) * (L ** alpha)


def find_max_overlap_length(M: int, X: Score, L_low: int = -1, L_high: int = -1) -> int:
    """
    Find the maximum integer L for which calc_overlap_pvalue(L, M) is
    greater than a given threshold X using binary search.

    This function assumes that the calc_overlap_pvalue function is
    monotonic (or nearly so) in L over the range of interest.  It
    performs a binary search within the bounds [L_low, L_high] to
    efficiently determine the largest L for which the score exceeds X.

    Args:
        M (int): The parameter M used in calc_overlap_pvalue.
        X (Score): The threshold value; the function finds the maximum
        L such that calc_overlap_pvalue(L, M) > X.
        L_low (int, optional): The lower bound of the search interval
        for L. Defaults to M.
        L_high (int, optional): The upper bound of the search interval
        for L. Defaults to ~M * M.

    Returns:
        int: The largest integer L for which calc_overlap_pvalue(L, M) > X.
    """

    if L_low < 0:
        L_low = M

    score = 0.0
    old_score = score
    ret = L_low
    while True:
        if L_high >= 0 and ret >= L_high:
            return L_high

        score = calculate_overlap_score(L=ret, M=M)
        if score < X:
            return ret

        if score == old_score:
            return ret - 1

        old_score = score
        ret += 1
