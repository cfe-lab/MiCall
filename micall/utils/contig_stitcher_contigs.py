from dataclasses import dataclass
from typing import Optional, Literal
from functools import cached_property
from aligntools import CigarHit


ID_STATE = 0

def generate_new_id() -> int:
    global ID_STATE
    ID_STATE += 1
    return ID_STATE


@dataclass(frozen=True)
class Contig:
    name: Optional[str]
    seq: str

    @cached_property
    def id(self) -> int:
        return generate_new_id()

    @cached_property
    def unique_name(self) -> str:
        index = self.register()
        unqualified = repr(self.name) if self.name is not None else ""
        if index == 1 and self.name:
            return unqualified
        else:
            return unqualified + f'({index})'

    def register(self) -> int:
        from micall.utils.contig_stitcher_context import context
        ctx = context.get()
        return ctx.register(key=self.id, value=self.name)


@dataclass(frozen=True)
class GenotypedContig(Contig):
    ref_name: str
    group_ref: Optional[str]

    # The sequence of self.group_ref. None in cases where the reference organism is unknown.
    ref_seq: Optional[str]

    # Approximated overall concordance between `seq` and `ref_seq`.
    # It is calculated by BLAST as qcovhsp/100, where qcovhsp means Query Coverage Per HSP.
    match_fraction: float


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
