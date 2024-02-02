from typing import Union, List, Tuple, Literal
from dataclasses import dataclass
from fractions import Fraction
from micall.utils.cigar_tools import Cigar, connect_cigar_hits, CigarHit


@dataclass(frozen=True)
class Cut:
    original: 'Contig'
    left: 'Contig'
    right: 'Contig'


@dataclass(frozen=True)
class LStrip:
    original: 'AlignedContig'
    result: 'AlignedContig'


@dataclass(frozen=True)
class RStrip:
    original: 'AlignedContig'
    result: 'AlignedContig'


@dataclass(frozen=True)
class Munge:
    left: 'AlignedContig'
    right: 'AlignedContig'
    result: 'AlignedContig'


@dataclass(frozen=True)
class Combine:
    contigs: List['AlignedContig']
    result: 'AlignedContig'


@dataclass(frozen=True)
class NoRef:
    contig: 'GenotypedContig'


@dataclass(frozen=True)
class InitialHit:
    contig: 'GenotypedContig'
    hit: CigarHit
    strand: Literal["forward", "reverse"]


@dataclass(frozen=True)
class ZeroHits:
    contig: 'GenotypedContig'


@dataclass(frozen=True)
class StrandConflict:
    contig: 'GenotypedContig'


@dataclass
class ReverseComplement:
    contig: 'GenotypedContig'
    result: 'GenotypedContig'


@dataclass(frozen=True)
class HitNumber:
    contig: 'GenotypedContig'
    initial: List[Tuple[CigarHit, Literal["reverse", "forward"]]]
    connected: List[CigarHit]


@dataclass(frozen=True)
class ConnectedHit:
    contig: 'GenotypedContig'
    part: 'AlignedContig'
    index: int


@dataclass(frozen=True)
class InitialStrip:
    contig: 'AlignedContig'
    q_st: int
    q_ei: int


@dataclass(frozen=True)
class StitchCut:
    left: 'AlignedContig'
    right: 'AlignedContig'
    left_overlap: 'AlignedContig'
    right_overlap: 'AlignedContig'
    left_remainder: 'AlignedContig'
    right_remainder: 'AlignedContig'


@dataclass(frozen=True)
class Overlap:
    left: 'AlignedContig'
    right: 'AlignedContig'
    left_overlap: 'AlignedContig'
    right_overlap: 'AlignedContig'
    left_remainder: 'AlignedContig'
    right_remainder: 'AlignedContig'
    left_take: 'AlignedContig'
    right_take: 'AlignedContig'
    concordance: List[Fraction]
    average: Fraction
    cut_point: int
    cut_point_scaled: Fraction


@dataclass(frozen=True)
class NoOverlap:
    contig: 'AlignedContig'


@dataclass(frozen=True)
class Stitch:
    left: 'AlignedContig'
    right: 'AlignedContig'
    result: 'AlignedContig'


@dataclass(frozen=True)
class Drop:
    contig: 'AlignedContig'
    covering: List['AlignedContig']


@dataclass(frozen=True)
class IgnoreGap:
    contig: 'AlignedContig'
    gap: 'CigarHit'


@dataclass(frozen=True)
class SplitGap:
    contig: 'AlignedContig'
    gap: 'CigarHit'
    left: 'AlignedContig'
    right: 'AlignedContig'


@dataclass(frozen=True)
class Intro:
    contig: 'GenotypedContig'


@dataclass(frozen=True)
class FinalCombine:
    contigs: List['AlignedContig']
    result: 'AlignedContig'


AlignmentEvent = Union[NoRef, InitialHit, ZeroHits, StrandConflict, ReverseComplement, HitNumber, ConnectedHit]
ModifyEvent = Union[LStrip, RStrip]
EventType = Union[Cut, ModifyEvent, Munge, Combine, AlignmentEvent, InitialStrip, StitchCut, Overlap, NoOverlap, Stitch, Drop, IgnoreGap, SplitGap, Intro, FinalCombine]
