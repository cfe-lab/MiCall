from dataclasses import dataclass
from typing import Optional, Literal

from micall.utils.cigar_tools import CigarHit


@dataclass(frozen=True)
class Contig:
    name: str
    seq: str


@dataclass(frozen=True)
class GenotypedContig(Contig):
    ref_name: str
    group_ref: str
    ref_seq: Optional[str]          # The sequence of self.group_ref. None in cases where the reference organism is unknown.
    match_fraction: float           # Approximated overall concordance between `seq` and `ref_seq`. It is calculated by BLAST as qcovhsp/100, where qcovhsp means Query Coverage Per HSP.


@dataclass(frozen=True)
class AlignedContig(GenotypedContig):
    alignment: CigarHit
    strand: Literal["forward", "reverse"]

    @staticmethod
    def make(query: GenotypedContig, alignment: CigarHit, strand: Literal["forward", "reverse"]):
        return AlignedContig(
            alignment=alignment,
            strand=strand,
            seq=query.seq,
            name=query.name,
            ref_name=query.ref_name,
            group_ref=query.group_ref,
            ref_seq=query.ref_seq,
            match_fraction=query.match_fraction)
