from fractions import Fraction
from typing import Sequence, Iterator, Tuple, TypeVar
from operator import itemgetter
from gotoh import align_it
import numpy as np
import math


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


def calc_overlap_pvalue(L: int, M: int) -> float:
    """
    Compute the probability (p-value) of observing at least M matches
    out of L under a binomial model where each position has
    probability `match_prob` of matching.
    Lower return value suggest that the matches are not coincidental.

    NOTE: This is a very simplistic model.

    :param L: Total length of the overlap region (int)
    :param M: Number of matching positions observed (int)
    :return: p-value (float), the probability of seeing at least M
    matches by chance
    """

    L += 1
    baseline = L / 4
    extra = M - baseline

    return 9 + extra * (-1 * math.log(1/L))
