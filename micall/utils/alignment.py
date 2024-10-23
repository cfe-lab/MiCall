from typing import Tuple, List, Sequence, Optional
from dataclasses import dataclass
from operator import attrgetter

from aligntools import CigarActions, Cigar, CigarHit
from gotoh import align_it
from mappy import Aligner
import mappy


@dataclass(frozen=True)
class Alignment:
    """
    Our representation of mappy's Alignment object.
    """

    ctg: str
    ctg_len: int
    r_st: int
    r_en: int
    strand: int
    q_st: int
    q_en: int
    mapq: int
    cigar: Sequence[Tuple[int, CigarActions]]
    cigar_str: str

    @staticmethod
    def coerce(obj: object) -> 'Alignment':
        if isinstance(obj, Alignment):
            return obj
        elif isinstance(obj, mappy.Alignment):
            cigar: List[Tuple[int, CigarActions]] = []
            for (size, action) in obj.cigar:
                cigar.append((size, CigarActions(action)))

            return Alignment(ctg=obj.ctg,
                             ctg_len=obj.ctg_len,
                             r_st=obj.r_st, r_en=obj.r_en,
                             strand=obj.strand,
                             q_st=obj.q_st, q_en=obj.q_en,
                             mapq=obj.mapq,
                             cigar=cigar,
                             cigar_str=obj.cigar_str,
                             )
        else:
            raise TypeError(f"Cannot coerce from {obj!r}.")

    def to_cigar_hit(self) -> CigarHit:
        return CigarHit(Cigar(self.cigar),
                        r_st=self.r_st, r_ei=self.r_en - 1,
                        q_st=self.q_st, q_ei=self.q_en - 1)

    @staticmethod
    def from_cigar_hit(hit: CigarHit, ctg: str, ctg_len: int, strand: int, mapq: int) -> 'Alignment':
        return Alignment(ctg=ctg,
                         ctg_len=ctg_len,
                         r_st=hit.r_st, r_en=hit.r_ei + 1,
                         strand=strand,
                         q_st=hit.q_st, q_en=hit.q_ei + 1,
                         mapq=mapq,
                         cigar=list(hit.cigar._data),
                         cigar_str=str(hit.cigar),
                         )


def align_gotoh(coordinate_seq: str, consensus: str) -> Optional[Alignment]:
    gap_open_penalty = 15
    gap_extend_penalty = 3
    use_terminal_gap_penalty = 1
    assert '&' not in consensus, "Consensus contains forbidden character '&'"
    consensus = ''.join('&' if x == '-' else x for x in consensus)
    aligned_coordinate, aligned_consensus, score = align_it(
        coordinate_seq,
        consensus,
        gap_open_penalty,
        gap_extend_penalty,
        use_terminal_gap_penalty)

    if min(len(coordinate_seq), len(consensus)) < score:
        cigar = Cigar.from_msa(aligned_coordinate, aligned_consensus)
        hit = CigarHit(cigar,
                       q_st=0, q_ei=len(consensus)-1,
                       r_st=0, r_ei=len(coordinate_seq)-1)
        hit = hit.lstrip_query().lstrip_reference().rstrip_query().rstrip_reference()
        return Alignment.from_cigar_hit(
            hit,
            ctg='N/A',
            ctg_len=len(coordinate_seq),
            strand=1,
            mapq=0)
    else:
        return None


def align_consensus(coordinate_seq: str, consensus: str) -> Tuple[List[Alignment], str]:
    aligner = Aligner(seq=coordinate_seq, preset='map-ont')
    mappy_alignments: List[mappy.Alignment] = list(aligner.map(consensus))
    if mappy_alignments or 10_000 < len(consensus):
        algorithm = 'minimap2'
        alignments = [Alignment.coerce(alignment)
                      for alignment in mappy_alignments
                      if alignment.is_primary]
    else:
        algorithm = 'gotoh'
        gotoh_alignment = align_gotoh(coordinate_seq, consensus)
        if gotoh_alignment:
            alignments = [gotoh_alignment]
        else:
            alignments = []

    alignments.sort(key=attrgetter('q_st'))
    return (alignments, algorithm)
