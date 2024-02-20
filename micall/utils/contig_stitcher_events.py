from typing import Union, List, Tuple, Literal
from dataclasses import dataclass
from fractions import Fraction

from micall.utils.cigar_tools import CigarHit
from micall.utils.contig_stitcher_contigs import GenotypedContig, AlignedContig


@dataclass(frozen=True)
class Cut:
    original: AlignedContig
    left: AlignedContig
    right: AlignedContig
    cut_point: float

    def __str__(self) -> str:
        return f"Created contigs {self.left.name!r} at {self.left.alignment} and {self.right.name!r} at {self.right.alignment} by cutting {self.original.name!r} at {self.original.alignment} at cut point = {round(self.cut_point, 1)}."


@dataclass(frozen=True)
class LStrip:
    original: AlignedContig
    result: AlignedContig

    def __str__(self) -> str:
        return f"Doing lstrip of {self.original.name!r} at {self.original.alignment} (len {len(self.original.seq)}) resulted in {self.result.name!r} at {self.result.alignment} (len {len(self.result.seq)})."


@dataclass(frozen=True)
class RStrip:
    original: AlignedContig
    result: AlignedContig

    def __str__(self) -> str:
        return f"Doing rstrip of {self.original.name!r} at {self.original.alignment} (len {len(self.original.seq)}) resulted in {self.result.name!r} at {self.result.alignment} (len {len(self.result.seq)})."


@dataclass(frozen=True)
class Munge:
    left: AlignedContig
    right: AlignedContig
    result: AlignedContig

    def __str__(self) -> str:
        return f"Munged contigs {self.left.name!r} at {self.left.alignment} with {self.right.name!r} at {self.right.alignment} resulting in {self.result.name!r} at {self.result.alignment}."


@dataclass(frozen=True)
class Combine:
    contigs: List[AlignedContig]
    result: AlignedContig

    def __str__(self) -> str:
        return f"Created a frankenstein {self.result.name!r} at {self.result.alignment} (len {len(self.result.seq)}) from {[f'{x.name!r} at {x.alignment} (len {len(x.seq)})' for x in self.contigs]}."


@dataclass(frozen=True)
class NoRef:
    contig: GenotypedContig

    def __str__(self) -> str:
        return f"Contig {self.contig.name!r} not aligned - no reference."


@dataclass(frozen=True)
class InitialHit:
    contig: GenotypedContig
    index: int
    hit: CigarHit
    strand: Literal["forward", "reverse"]

    def __str__(self) -> str:
        return f"Part {self.index} of contig {self.contig.name!r} aligned at {self.hit}{'' if self.strand == 'forward' else ' (rev)'}."


@dataclass(frozen=True)
class ZeroHits:
    contig: GenotypedContig

    def __str__(self) -> str:
        return f"Contig {self.contig.name!r} not aligned - backend's choice."


@dataclass(frozen=True)
class StrandConflict:
    contig: GenotypedContig

    def __str__(self) -> str:
        return f"Discarding contig {self.contig.name!r} because it aligned both in forward and reverse sense."


@dataclass(frozen=True)
class ReverseComplement:
    contig: GenotypedContig
    result: GenotypedContig

    def __str__(self) -> str:
        return f"Reverse complemented contig {self.contig.name!r}."


@dataclass(frozen=True)
class HitNumber:
    contig: GenotypedContig
    initial: List[Tuple[CigarHit, Literal["reverse", "forward"]]]
    connected: List[CigarHit]

    def __str__(self) -> str:
        return f"Contig {self.contig.name!r} produced {len(self.initial)} aligner hits. After connecting them, the number became {len(self.connected)}."


@dataclass(frozen=True)
class ConnectedHit:
    contig: GenotypedContig
    part: AlignedContig
    index: int

    def __str__(self) -> str:
        return f"Part {self.index} of contig {self.contig.name!r} re-aligned as {self.part.name!r} at {self.part.alignment}{'' if self.part.strand == 'forward' else ' (rev)'}."


@dataclass(frozen=True)
class InitialStrip:
    contig: AlignedContig
    q_st: int
    q_ei: int

    def __str__(self) -> str:
        return f"Trimming (strip) contig {self.contig.name!r} from {self.q_st} to {self.q_ei}."


@dataclass(frozen=True)
class StitchCut:
    left: AlignedContig
    right: AlignedContig
    left_overlap: AlignedContig
    right_overlap: AlignedContig
    left_remainder: AlignedContig
    right_remainder: AlignedContig

    def __str__(self) -> str:
        return f"Stitching {self.left.name!r} at {self.left.alignment} (len {len(self.left.seq)}) with {self.right.name!r} at {self.right.alignment} (len {len(self.right.seq)}). The left_overlap {self.left_overlap.name!r} is at {self.left_overlap.alignment} (len {len(self.left_overlap.seq)}) and the right_overlap {self.right_overlap.name!r} is at {self.right_overlap.alignment} (len {len(self.right_overlap.seq)})."


@dataclass(frozen=True)
class Overlap:
    left: AlignedContig
    right: AlignedContig
    left_overlap: AlignedContig
    right_overlap: AlignedContig
    left_remainder: AlignedContig
    right_remainder: AlignedContig
    left_take: AlignedContig
    right_take: AlignedContig
    concordance: List[Fraction]
    average: Fraction
    cut_point: int
    cut_point_scaled: Fraction

    def __str__(self) -> str:
        average_concordance = round(self.average * 100)
        cut_point_location_scaled = round(self.cut_point_scaled * 100)
        concordance_str = ', '.join(str(int(round(x * 100)) / 100) for x in self.concordance)
        return f"Created overlap contigs {self.left_take.name!r} at {self.left_overlap.alignment} and {self.right_take.name!r} at {self.right_take.alignment} based on parts of {self.left.name!r} and {self.right.name!r}, with avg. concordance {average_concordance}%, cut point at {cut_point_location_scaled}%, and full concordance [{concordance_str}]."


@dataclass(frozen=True)
class NoOverlap:
    contig: AlignedContig

    def __str__(self) -> str:
        return f"Nothing overlaps with {self.contig.name!r}."


@dataclass(frozen=True)
class Stitch:
    left: AlignedContig
    right: AlignedContig
    result: AlignedContig

    def __str__(self) -> str:
        return f"Stitching {self.left.name!r} with {self.right.name!r} results in {self.result.name!r} at {self.result.alignment} (len {len(self.result.seq)})."


@dataclass(frozen=True)
class Drop:
    contig: AlignedContig
    covering: List[AlignedContig]

    def __str__(self) -> str:
        return f"Dropped contig {self.contig.name!r} as it is completely covered by these contigs: {', '.join(repr(x.name) for x in self.covering)}."


@dataclass(frozen=True)
class IgnoreGap:
    contig: AlignedContig
    gap: CigarHit

    def __str__(self) -> str:
        return f"Ignored insignificant gap of {self.contig.name!r}, {self.gap}."


@dataclass(frozen=True)
class SplitGap:
    contig: AlignedContig
    gap: CigarHit
    left: AlignedContig
    right: AlignedContig

    def __str__(self) -> str:
        return f"Split contig {self.contig.name!r} at {self.contig.alignment} around its gap at [{self.gap.q_st}, {self.gap.q_ei}]->[{self.gap.r_st}, {self.gap.r_ei}]. Left part: {self.left.name!r} at {self.left.alignment}, right part: {self.right.name!r} at {self.right.alignment}."


@dataclass(frozen=True)
class Intro:
    contig: GenotypedContig

    def __str__(self) -> str:
        return f"Introduced contig {self.contig.name!r} (seq = {self.contig.seq}) of ref {self.contig.ref_name!r}, group_ref {self.contig.group_ref} (seq = {self.contig.ref_seq}), and length {len(self.contig.seq)}."


@dataclass(frozen=True)
class FinalCombine:
    contigs: List[AlignedContig]
    result: AlignedContig

    def __str__(self) -> str:
        return f"Combining these contigs for final output for {self.result.group_ref}: {['%r at %s (len %s)' % (x.name, x.alignment, len(x.seq)) for x in self.contigs]}."


AlignmentEvent = Union[NoRef, InitialHit, ZeroHits, StrandConflict, ReverseComplement, HitNumber, ConnectedHit]
ModifyEvent = Union[LStrip, RStrip]
EventType = Union[Cut, ModifyEvent, Munge, Combine, AlignmentEvent, InitialStrip, StitchCut, Overlap, NoOverlap, Stitch, Drop, IgnoreGap, SplitGap, Intro, FinalCombine]
