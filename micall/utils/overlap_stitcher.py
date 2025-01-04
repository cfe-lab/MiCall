from fractions import Fraction
from typing import Sequence, Iterator, Tuple
from operator import itemgetter
from functools import lru_cache
from gotoh import align_it
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


def calculate_concordance(left: Sequence[object], right: Sequence[object],
                          ) -> Sequence[Fraction]:
    """
    Calculate concordance for two given sequences using a sliding
    average.

    The function compares the two strings character by character,
    simultaneously from both left to right and right to left,
    calculating a score that represents a moving average of matches at
    each position. If characters match at a given position, a score of
    1 is added; otherwise, a score of 0 is added. The score is then
    averaged with the previous scores using a weighted sliding average
    where the current score has a weight of 1/3 and the accumulated
    score has a weight of 2/3.  This sliding average score is halved
    and then processed again, but in reverse direction.

    :param left: string representing first sequence
    :param right: string representing second sequence
    :return: list representing concordance ratio for each position
    """

    if len(left) != len(right):
        raise ValueError("Can only calculate concordance for same sized sequences")

    result = [Fraction(0)] * len(left)

    def slide(start, end):
        scores_sum = Fraction(0)
        inputs = list(zip(left, right))
        increment = 1 if start <= end else -1

        for i in range(start, end, increment):
            (a, b) = inputs[i]
            current = Fraction(1) if a == b else Fraction(0)
            scores_sum = (scores_sum * 2 / 3 + current * 1 / 3)
            result[i] += scores_sum / 2

    # Slide forward, then in reverse, adding the scores at each position.
    slide(0, len(left))
    slide(len(left) - 1, -1)

    return result


def disambiguate_concordance(concordance: Sequence[Fraction],
                             ) -> Iterator[Tuple[Fraction, int]]:
    for i, x in enumerate(concordance):
        if i < len(concordance) / 2:
            global_rank = i
        else:
            global_rank = len(concordance) - i - 1
        yield x, global_rank


def sort_concordance_indexes(concordance: Sequence[Fraction]) -> Iterator[int]:
    concordance_d = disambiguate_concordance(concordance)
    for i, v in sorted(enumerate(concordance_d),
                       key=itemgetter(1),
                       reverse=True,
                       ):
        yield i


@lru_cache
def calc_overlap_pvalue(L: int, M: int) -> Fraction:
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

    # TODO: Implement an evidence-based p-value calculation.

    pval = Fraction(0)
    match_prob = Fraction(1, 4)

    # Accounting for the fact that the contigs differ at the ends
    # (left end is different and right end is different).
    L += 2

    # Summation of
    # Binomial(L, x) * match_prob^x * (1-match_prob)^(L-x)
    #   from x = M to L
    for x in range(M, L + 1):
        pval += (math.comb(L, x) *
                 (match_prob ** x) *
                 ((1 - match_prob) ** (L - x)))

    return pval
