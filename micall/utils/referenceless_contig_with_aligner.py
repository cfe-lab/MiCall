
import numpy as np

from typing import Iterator, Tuple, Mapping, Literal, NoReturn
from dataclasses import dataclass
from mappy import Aligner
from functools import cached_property

from micall.utils.consensus_aligner import Alignment
from micall.utils.referenceless_score import Score
from micall.utils.contig_stitcher_contigs import Contig
from micall.utils.overlap_stitcher import \
    exp_dropoff_array, find_max_overlap_length
from micall.utils.find_maximum_overlap import \
    get_overlap_results, choose_convolution_method


OverlapRelation = Literal["left", "right", "cover"]


@dataclass(frozen=True)
class ContigWithAligner(Contig):
    @cached_property
    def aligner(self) -> Aligner:
        return Aligner(seq=self.seq)

    @staticmethod
    def make(contig: Contig) -> 'ContigWithAligner':
        return ContigWithAligner(name=contig.name, seq=contig.seq)

    @staticmethod
    def empty() -> 'ContigWithAligner':
        return ContigWithAligner.make(Contig.empty())

    def map_overlap(self,
                    minimum_score: Score,
                    relation: OverlapRelation,
                    overlap: str,
                    ) -> Iterator[Tuple[Alignment, int]]:

        optimistic_number_of_matches = len(overlap)
        max_length = find_max_overlap_length(M=optimistic_number_of_matches,
                                             X=minimum_score,
                                             L_high=len(self.seq),
                                             )

        assert max_length > 0
        assert max_length >= len(overlap)
        assert max_length <= len(self.seq)

        if max_length < len(self.seq) and relation != "cover":
            if relation == "left":
                seq = self.seq[-max_length:]
                shift = len(self.seq) - max_length
            elif relation == "right":
                seq = self.seq[:max_length]
                shift = 0
            else:
                _x: NoReturn = relation
            aligner = Aligner(seq=seq)
        else:
            aligner = self.aligner
            shift = 0

        for x in aligner.map(overlap):
            if x.is_primary:
                yield (x, shift)

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

    def find_maximum_overlap(self, other: 'ContigWithAligner',
                             ) -> Tuple[int, float]:

        total = np.zeros(len(self.seq) + len(other.seq) - 1)
        method = choose_convolution_method(len(self.seq), len(other.seq))
        keys = sorted(set(list(self.alphabet) + list(other.alphabet)))
        for key in keys:
            if key not in self.alignment_seqs or \
               key not in other.alignment_seqs:
                continue

            x = self.alignment_seqs[key]
            y = np.flip(other.alignment_seqs[key])
            total += method(x, y, mode='full')

        return get_overlap_results(total, len(self.seq), len(other.seq))
