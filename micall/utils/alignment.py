from typing import Tuple, List, Sequence
from dataclasses import dataclass

from aligntools import CigarActions, Cigar, CigarHit
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
    def from_cigar_hit(hit: CigarHit, ctg='', ctg_len=0, strand=1, mapq=0) -> 'Alignment':
        return Alignment(ctg=ctg,
                         ctg_len=ctg_len,
                         r_st=hit.r_st, r_en=hit.r_ei + 1,
                         strand=strand,
                         q_st=hit.q_st, q_en=hit.q_ei + 1,
                         mapq=mapq,
                         cigar=hit.cigar._data,
                         cigar_str=str(hit.cigar),
                         )
