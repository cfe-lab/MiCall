
import numpy as np

from abc import ABC, abstractmethod
from typing import Iterator, Tuple, Mapping, Literal, NoReturn
from dataclasses import dataclass
from mappy import Aligner as OriginalMappyAligner
from functools import cached_property

from micall.utils.referenceless_score import Score
from micall.utils.contig_stitcher_contigs import Contig
from micall.utils.overlap_stitcher import \
    exp_dropoff_array, find_max_overlap_length
from micall.utils.find_maximum_overlap import \
    get_overlap_results, choose_convolution_method


OverlapRelation = Literal["left", "right", "cover"]


class LocalAligner(ABC):
    @abstractmethod
    def map(self, query: str) -> Iterator[Tuple[int, int]]: ...


class MappyAligner(LocalAligner):
    def __init__(self, seq: str) -> None:
        self.aligner = OriginalMappyAligner(seq=seq)

    def map(self, query: str) -> Iterator[Tuple[int, int]]:
        for x in self.aligner.map(query):
            if x.is_primary:
                yield (x.r_st, x.r_en)


PADDING = 40


class ForwardAligner(LocalAligner):
    def __init__(self, seq: str) -> None:
        if seq.startswith('A'):
            self.seq = 'C' * PADDING + seq
        else:
            self.seq = 'A' * PADDING + seq
        self.aligner = OriginalMappyAligner(seq=self.seq)

    def map(self, query: str) -> Iterator[Tuple[int, int]]:
        if self.seq.startswith('A'):
            query = 'A' * PADDING + query
        else:
            query = 'C' * PADDING + query

        for x in self.aligner.map(query):
            if x.is_primary:
                end = x.r_en - PADDING
                if end > 0:
                    yield (-99999999999, end)


class ReversedAligner(LocalAligner):
    def __init__(self, seq: str) -> None:
        if seq.endswith('A'):
            self.seq = seq + 'C' * PADDING
        else:
            self.seq = seq + 'A' * PADDING
        self.aligner = OriginalMappyAligner(seq=self.seq)

    def map(self, query: str) -> Iterator[Tuple[int, int]]:
        if self.seq.endswith('A'):
            query = query + 'A' * PADDING
        else:
            query = query + 'C' * PADDING

        for x in self.aligner.map(query):
            if x.is_primary:
                start = x.r_st
                if start < len(self.seq) - PADDING:
                    yield (start, -99999999999999999)


@dataclass(frozen=True)
class ContigWithAligner(Contig):
    @cached_property
    def mappy_aligner(self) -> LocalAligner:
        return MappyAligner(seq=self.seq)

    @cached_property
    def forward_aligner(self) -> LocalAligner:
        return ForwardAligner(seq=self.seq)

    @cached_property
    def reversed_aligner(self) -> LocalAligner:
        return ReversedAligner(seq=self.seq)

    @staticmethod
    def make(contig: Contig) -> 'ContigWithAligner':
        return ContigWithAligner(name=contig.name, seq=contig.seq)

    @staticmethod
    def empty() -> 'ContigWithAligner':
        return ContigWithAligner.make(Contig.empty())

    @cached_property
    def nucleotide_seq(self) -> np.ndarray:
        ret = np.frombuffer(self.seq.encode('utf-8'), dtype='S1')
        return ret

    @cached_property
    def alphabet(self) -> Tuple[str, ...]:
        return tuple(sorted(set(self.seq)))

    @cached_property
    def alignment_seqs(self) -> Mapping[str, np.ndarray]:
        def to_array(letter: str) -> np.ndarray:
            value = letter.encode('utf-8')
            ret = np.zeros(len(self.nucleotide_seq))
            ret[self.nucleotide_seq == value] = 1
            exp_dropoff_array(ret, factor=8)
            return ret

        return {x: to_array(x) for x in self.alphabet}

def find_maximum_overlap(
    left: ContigWithAligner, right: ContigWithAligner
) -> Tuple[int, float]:

    total = np.zeros(len(left.seq) + len(right.seq) - 1)
    method = choose_convolution_method(len(left.seq), len(right.seq))
    keys = sorted(set(list(left.alphabet) + list(right.alphabet)))
    for key in keys:
        if key not in left.alignment_seqs or \
            key not in right.alignment_seqs:
            continue

        x = left.alignment_seqs[key]
        y = np.flip(right.alignment_seqs[key])
        total += method(x, y, mode='full')

    return get_overlap_results(total, len(left.seq), len(right.seq))


def map_overlap(self: ContigWithAligner,
                minimum_score: Score,
                relation: OverlapRelation,
                overlap: str,
                ) -> Iterator[Tuple[int, int]]:

    if relation == "left":
        aligner = self.reversed_aligner
    elif relation == "right":
        aligner = self.forward_aligner
    else:
        aligner = self.mappy_aligner
    shift = 0

    if relation != "cover":
        optimistic_number_of_matches = len(overlap)
        max_length = find_max_overlap_length(M=optimistic_number_of_matches,
                                                X=minimum_score,
                                                L_high=len(self.seq),
                                                )

        assert max_length > 0
        assert max_length >= len(overlap)
        assert max_length <= len(self.seq)

        if max_length < len(self.seq):
            if relation == "left":
                seq = self.seq[-max_length:]
                shift = len(self.seq) - max_length
                aligner = ReversedAligner(seq=seq)
            elif relation == "right":
                seq = self.seq[:max_length]
                shift = 0
                aligner = ForwardAligner(seq=seq)
            else:
                _x: NoReturn = relation

    for start, end in aligner.map(overlap):
        yield (start + shift, end + shift)
