from typing import Union, Sequence, Tuple, Literal, List
from dataclasses import dataclass
from fractions import Fraction
from aligntools import CigarHit

from micall.utils.contig_stitcher_contigs import GenotypedContig, AlignedContig


class Warning:
    pass


@dataclass(frozen=True)
class Cut:
    original: AlignedContig
    left: AlignedContig
    right: AlignedContig
    cut_point: float

    def __str__(self) -> str:
        return (
            f"Created contigs {self.left.unique_name} at {self.left.alignment} and "
            f"{self.right.unique_name} at {self.right.alignment} by cutting "
            f"{self.original.unique_name} at {self.original.alignment} at cut point = "
            f"{round(self.cut_point, 1)}."
        )


@dataclass(frozen=True)
class LStrip:
    original: AlignedContig
    result: AlignedContig

    def __str__(self) -> str:
        return (
            f"Doing lstrip of {self.original.unique_name} at {self.original.alignment} (len "
            f"{len(self.original.seq)}) resulted in {self.result.unique_name} at "
            f"{self.result.alignment} (len {len(self.result.seq)})."
        )


@dataclass(frozen=True)
class RStrip:
    original: AlignedContig
    result: AlignedContig

    def __str__(self) -> str:
        return (
            f"Doing rstrip of {self.original.unique_name} at {self.original.alignment} (len "
            f"{len(self.original.seq)}) resulted in {self.result.unique_name} at "
            f"{self.result.alignment} (len {len(self.result.seq)})."
        )


@dataclass(frozen=True)
class Munge:
    left: AlignedContig
    right: AlignedContig
    result: AlignedContig

    def __str__(self) -> str:
        return (
            f"Munged contigs {self.left.unique_name} at {self.left.alignment} with "
            f"{self.right.unique_name} at {self.right.alignment} resulting in "
            f"{self.result.unique_name} at {self.result.alignment}."
        )


@dataclass(frozen=True)
class Combine:
    contigs: Sequence[AlignedContig]
    result: AlignedContig

    def __str__(self) -> str:
        contigs_str = ', '.join(
            [f"{x.unique_name} at {x.alignment} (len {len(x.seq)})" for x in self.contigs])
        return (
            f"Created a frankenstein {self.result.unique_name} at {self.result.alignment} "
            f"(len {len(self.result.seq)}) from [{contigs_str}]."
        )


@dataclass(frozen=True)
class NoRef:
    contig: GenotypedContig

    def __str__(self) -> str:
        return f"Contig {self.contig.unique_name} not aligned - no reference."


@dataclass(frozen=True)
class InitialHit:
    contig: GenotypedContig
    index: int
    hit: CigarHit
    strand: Literal["forward", "reverse"]

    def __str__(self) -> str:
        strand_info = '' if self.strand == 'forward' else ' (rev)'
        return (
            f"Part {self.index} of contig {self.contig.unique_name} aligned at {self.hit}"
            f"{strand_info}."
        )


@dataclass(frozen=True)
class ZeroHits:
    contig: GenotypedContig

    def __str__(self) -> str:
        return f"Contig {self.contig.unique_name} not aligned - backend's choice."


@dataclass(frozen=True)
class StrandConflict:
    contig: GenotypedContig

    def __str__(self) -> str:
        return (
            f"Discarding contig {self.contig.unique_name} because it aligned both in forward "
            "and reverse sense."
        )


@dataclass(frozen=True)
class ReverseComplement:
    contig: GenotypedContig
    result: GenotypedContig

    def __str__(self) -> str:
        return f"Reverse complemented contig {self.contig.unique_name}."


@dataclass(frozen=True)
class HitNumber:
    contig: GenotypedContig
    initial: Sequence[Tuple[CigarHit, Literal["reverse", "forward"]]]
    connected: Sequence[CigarHit]

    def __str__(self) -> str:
        return (
            f"Contig {self.contig.unique_name} produced {len(self.initial)} aligner hits. "
            f"After connecting them, the number became {len(self.connected)}."
        )


@dataclass(frozen=True)
class ConnectedHit:
    contig: GenotypedContig
    part: AlignedContig
    index: int

    def __str__(self) -> str:
        part_strand_info = '' if self.part.strand == 'forward' else ' (rev)'
        return (
            f"Part {self.index} of contig {self.contig.unique_name} re-aligned as "
            f"{self.part.unique_name} at {self.part.alignment}{part_strand_info}."
        )


@dataclass(frozen=True)
class InitialStrip:
    contig: AlignedContig
    q_st: int
    q_ei: int

    def __str__(self) -> str:
        return (
            f"Trimming (strip) contig {self.contig.unique_name} from {self.q_st} to "
            f"{self.q_ei}."
        )


@dataclass(frozen=True)
class StitchCut:
    left: AlignedContig
    right: AlignedContig
    left_overlap: AlignedContig
    right_overlap: AlignedContig
    left_remainder: AlignedContig
    right_remainder: AlignedContig

    def __str__(self) -> str:
        return (
            f"Stitching {self.left.unique_name} at {self.left.alignment} (len {len(self.left.seq)}) "
            f"with {self.right.unique_name} at {self.right.alignment} (len {len(self.right.seq)}). "
            f"The left_overlap {self.left_overlap.unique_name} is at {self.left_overlap.alignment} "
            f"(len {len(self.left_overlap.seq)}) and the right_overlap {self.right_overlap.unique_name} is "
            f"at {self.right_overlap.alignment} (len {len(self.right_overlap.seq)})."
        )


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
    concordance: Sequence[Fraction]
    cut_point: int
    cut_point_scaled: float

    def __str__(self) -> str:
        cut_point_location_scaled = round(self.cut_point_scaled * 100)
        concordance_str = ', '.join(str(int(round(x * 100)) / 100) for x in self.concordance)
        return (
            f"Created overlap contigs {self.left_take.unique_name} at {self.left_overlap.alignment} and "
            f"{self.right_take.unique_name} at {self.right_take.alignment} based on parts of "
            f"{self.left.unique_name} and {self.right.unique_name}, with "
            f"cut point at {cut_point_location_scaled}%, and full concordance [{concordance_str}]."
        )


@dataclass(frozen=True)
class NoOverlap:
    contig: AlignedContig

    def __str__(self) -> str:
        return f"Nothing overlaps with {self.contig.unique_name}."


@dataclass(frozen=True)
class Stitch:
    left: AlignedContig
    right: AlignedContig
    result: AlignedContig

    def __str__(self) -> str:
        return (
            f"Stitching {self.left.unique_name} with {self.right.unique_name} results in "
            f"{self.result.unique_name} at {self.result.alignment} (len {len(self.result.seq)})."
        )


@dataclass(frozen=True)
class Drop:
    contig: AlignedContig
    covering: Sequence[AlignedContig]

    def __str__(self) -> str:
        covering_contig_names = ', '.join(repr(x.unique_name) for x in self.covering)
        return (
            f"Dropped contig {self.contig.unique_name} as it is completely covered by these contigs: "
            f"{covering_contig_names}."
        )


@dataclass(frozen=True)
class IgnoreGap:
    contig: AlignedContig
    gap: CigarHit

    def __str__(self) -> str:
        return f"Ignored insignificant gap of {self.contig.unique_name}, {self.gap}."


@dataclass(frozen=True)
class SplitGap:
    contig: AlignedContig
    gap: CigarHit
    left: AlignedContig
    right: AlignedContig

    def __str__(self) -> str:
        return (
            f"Split contig {self.contig.unique_name} at {self.contig.alignment} around its gap at "
            f"[{self.gap.q_st}, {self.gap.q_ei}]->[{self.gap.r_st}, {self.gap.r_ei}]. Left part: "
            f"{self.left.unique_name} at {self.left.alignment}, right part: {self.right.unique_name} at "
            f"{self.right.alignment}."
        )


@dataclass(frozen=True)
class Intro:
    contig: GenotypedContig

    def __str__(self) -> str:
        return (
            f"Introduced contig {self.contig.unique_name} (seq = {self.contig.seq}) of ref "
            f"{self.contig.ref_name!r}, group_ref {self.contig.group_ref} (seq = {self.contig.ref_seq}), "
            f"and length {len(self.contig.seq)}."
        )


@dataclass(frozen=True)
class FinalCombine:
    contigs: Sequence[AlignedContig]
    result: AlignedContig

    def __str__(self) -> str:
        contigs_str = [f"{x.unique_name} at {x.alignment} (len {len(x.seq)})" for x in self.contigs]
        contigs_format = ', '.join(contigs_str)
        return (
            f"Combining these contigs for final output for {self.result.group_ref}: "
            f"[{contigs_format}]."
        )

@dataclass(frozen=True)
class IgnoreCoverage(Warning):
    current: AlignedContig
    overlaping_contigs: List[AlignedContig]

    def __str__(self) -> str:
        overlapping_names = ', '.join([c.unique_name for c in self.overlaping_contigs])
        return (
            f"Ignoring coverage comparison for contig {self.current.unique_name} at "
            f"{self.current.alignment} with overlapping contigs [{overlapping_names}] "
            "due to unknown read counts."
        )


AlignmentEvent = Union[NoRef, InitialHit, ZeroHits, StrandConflict, ReverseComplement,
                       HitNumber, ConnectedHit]
ModifyEvent = Union[LStrip, RStrip]
EventType = Union[Cut, ModifyEvent, Munge, Combine, AlignmentEvent, InitialStrip, StitchCut,
                  Overlap, NoOverlap, Stitch, Drop, IgnoreGap, SplitGap, Intro, FinalCombine, IgnoreCoverage]
