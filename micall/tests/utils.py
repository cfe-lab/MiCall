from dataclasses import dataclass
from math import floor, ceil
from contextlib import contextmanager
import random

from micall.utils.consensus_aligner import CigarActions


@dataclass
class MockAlignment:
    is_rev: bool
    mapq: int
    cigar: list
    cigar_str: str
    q_st: int
    q_en: int
    r_st: int
    r_en: int


class MockAligner:
    """
    Mock for the mappy's aligner class.
    Only reports exact matches.
    """

    def __init__(self, seq, *args, **kwargs):
        self.seq = seq
        self.max_matches = 5
        self.min_length = 3


    def map(self, seq):
        max_matches = self.max_matches
        returned = set()
        for length in range(len(seq), self.min_length - 1, -1):
            for start in range(len(seq) - length + 1):
                end = start + length
                substring = seq[start:end]
                if substring not in self.seq:
                    continue

                mapq = 60
                is_rev = False # Doesn't handle reverse complements in this mock.
                r_st = self.seq.index(substring)
                r_en = r_st + len(substring)
                q_st = start
                q_en = end
                cigar = [[q_en - q_st, CigarActions.MATCH]]
                cigar_str = f'{(q_en - q_st)}M'
                al = MockAlignment(is_rev, mapq, cigar, cigar_str, q_st, q_en, r_st, r_en)
                if (q_st, q_en, r_st, r_en) not in returned:
                    returned.add((q_st, q_en, r_st, r_en))
                    yield MockAlignment(is_rev, mapq, cigar, cigar_str, q_st, q_en, r_st, r_en)

                    max_matches -= 1
                    if max_matches < 1:
                        return


@contextmanager
def fixed_random_seed(seed):
    original_state = random.getstate()
    random.seed(seed)
    try:
        yield
    finally:
        random.setstate(original_state)
