from typing import Iterable, Tuple, List

from aligntools import CigarActions, Cigar, CigarHit
import mappy


class Alignment:
    """
    Our representation of mappy's Alignment object.
    """

    def __init__(self,
                 ctg='',
                 ctg_len=0,
                 r_st=0,
                 r_en=0,
                 strand=1,
                 q_st=0,
                 q_en=0,
                 mapq=0,
                 cigar: Iterable[Tuple[int, CigarActions]] = tuple(),
                 cigar_str=None):

        cigar = list(cigar)
        if not cigar:
            cigar = [(max(q_en-q_st, r_en-r_st), CigarActions.MATCH)]
        if cigar_str is None:
            cigar_str = str(Cigar(cigar))

        self.ctg = ctg
        self.ctg_len = ctg_len
        self.r_st = r_st
        self.r_en = r_en
        self.strand = strand
        self.q_st = q_st
        self.q_en = q_en
        self.mapq = mapq
        self.cigar = cigar
        self.cigar_str = cigar_str

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

    def __eq__(self, other: object):
        # Filter out private attributes (those starting with an underscore)
        self_public_attrs = {k: v for k, v in self.__dict__.items() if not k.startswith('_')}
        other_public_attrs = {k: v for k, v in other.__dict__.items() if not k.startswith('_')}
        return self_public_attrs == other_public_attrs

    def __repr__(self):
        return (f'Alignment({self.ctg!r}, {self.ctg_len}, '
                f'{self.r_st}, {self.r_en}, {self.strand}, '
                f'{self.q_st}, {self.q_en})')
