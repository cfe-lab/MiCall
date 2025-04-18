from fractions import Fraction
from typing import Sequence, Iterator, Tuple, TypeVar
from operator import itemgetter
from gotoh import align_it
import numpy as np
from micall.utils.referenceless_score import Score


def align_queries(seq1: str, seq2: str) -> Tuple[str, str]:
    """
    Globally align two query sequences against each other
    and return the resulting aligned sequences in MSA format.
    """

    gap_open_penalty = 15
    gap_extend_penalty = 3
    use_terminal_gap_penalty = 1
    aseq1, aseq2, score = \
        align_it(
            seq1, seq2,
            gap_open_penalty,
            gap_extend_penalty,
            use_terminal_gap_penalty)

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

    diff = high - low
    return tuple(Fraction(x - low) / Fraction(diff) for x in array)


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

    n = len(arr)
    out = np.empty(n, dtype=np.int64)
    count = 0

    for i in range(n):
        count = (count + 1) * arr[i]
        out[i] = count

    return out


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
    Computes a monotonic scoring metric for an overlap event between
    two sequences using a four-letter alphabet. This function returns
    a z-score that quantifies how far the observed number of matching
    characters deviates from the expected value under the assumption
    that each character match occurs with a probability of 1/4.

    Given:
        - M: the observed number of matching characters in the overlap,
        - L: the length of the overlap.

    Under the uniform model for a four-letter alphabet:
        - The probability of a match per position is 1/4.
        - Expected number of matches: L / 4.
        - Variance: L * (1/4) * (3/4) = 3L / 16.
        - Standard deviation: sqrt(3L) / 4.

    The z-score is computed using the formula:
        z = (4*M - L) / sqrt(3*L)

    This z-score provides a monotonic measure for ranking overlap
    events: a higher z-score indicates a larger deviation from the
    expected match count and thus a rarer event.

    :param M: Number of matching characters.
    :param L: Length of the overlap (must be greater than 0).
    :return: A z-score as a Score that serves as a ranking metric for the event.
    :raises ValueError: If L is not greater than 0.
    """

    # Compute z-score for a four-letter alphabet where P(match)=1/4:
    # Expected matches = L / 4 and standard deviation = sqrt(3L) / 4.
    return (4 * M - L) / ((3 * L) ** 0.5)


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
