from contextlib import contextmanager
from dataclasses import dataclass
import random
from aligntools import CigarActions
from typing import Tuple, List, Iterator

from micall.utils.alignment import Alignment


def find_all_occurrences(s, substring):
    start = 0
    while True:
        start = s.find(substring, start)
        if start == -1:  # no more occurrences found
            return
        yield start
        start += len(substring)


@dataclass(frozen=True)
class MockAlignment(Alignment):
    is_primary: bool


class MockAligner:
    """
    Mock for the mappy's aligner class.
    Only reports exact matches.
    """

    def __init__(self, seq, *args, **kwargs):
        self.seq = seq
        self.max_matches = 5
        self.min_length = 3

    def map(self, seq) -> Iterator[Alignment]:
        max_matches = self.max_matches
        returned = set()
        for length in range(len(seq), self.min_length - 1, -1):
            for start in range(len(seq) - length + 1):
                end = start + length
                substring = seq[start:end]
                for r_st in find_all_occurrences(self.seq, substring):
                    mapq = 60
                    strand = 1  # Doesn't handle reverse complements in this mock.
                    r_en = r_st + len(substring)
                    q_st = start
                    q_en = end
                    cigar = [(q_en - q_st, CigarActions.MATCH)]
                    cigar_str = f'{(q_en - q_st)}M'
                    if (q_st, q_en, r_st, r_en) not in returned:
                        returned.add((q_st, q_en, r_st, r_en))
                        yield MockAlignment(ctg='N/A', ctg_len=0,
                                            strand=strand, mapq=mapq,
                                            cigar=cigar, cigar_str=cigar_str,
                                            q_st=q_st, q_en=q_en,
                                            r_st=r_st, r_en=r_en,
                                            is_primary=True)

                        max_matches -= 1
                        if max_matches < 1:
                            return


def mock_align_consensus(reference_seq: str, consensus: str) -> Tuple[List[Alignment], str]:
    alignments = list(MockAligner(reference_seq).map(consensus))
    algorithm = 'mock'
    return (alignments, algorithm)


@contextmanager
def fixed_random_seed(seed):
    original_state = random.getstate()
    random.seed(seed)
    try:
        yield
    finally:
        random.setstate(original_state)
