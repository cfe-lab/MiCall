from typing import Union, List
from dataclasses import dataclass
from fractions import Fraction


@dataclass
class Cut:
    original: 'Contig'
    left: 'Contig'
    right: 'Contig'


@dataclass
class LStrip:
    original: 'AlignedContig'
    result: 'AlignedContig'


@dataclass
class RStrip:
    original: 'AlignedContig'
    result: 'AlignedContig'


@dataclass
class Munge:
    left: 'AlignedContig'
    right: 'AlignedContig'
    result: 'AlignedContig'


@dataclass
class Combine:
    contigs: List['AlignedContig']
    result: 'AlignedContig'


@dataclass
class NoRef:
    contig: 'GenotypedContig'


@dataclass
class ZeroHits:
    contig: 'GenotypedContig'


@dataclass
class StrandConflict:
    contig: 'GenotypedContig'


@dataclass
class HitNumber:
    contig: 'GenotypedContig'
    initial: object
    connected: object


@dataclass
class ReverseComplement:
    contig: 'GenotypedContig'
    result: 'GenotypedContig'


@dataclass
class Hit:
    contig: 'GenotypedContig'
    part: 'AlignedContig'
    index: int


@dataclass
class StitchCut:
    left: 'AlignedContig'
    right: 'AlignedContig'
    left_overlap: 'AlignedContig'
    right_overlap: 'AlignedContig'
    left_remainder: 'AlignedContig'
    right_remainder: 'AlignedContig'


@dataclass
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


@dataclass
class NoOverlap:
    contig: 'AlignedContig'


@dataclass
class Stitch:
    left: 'AlignedContig'
    right: 'AlignedContig'
    result: 'AlignedContig'


@dataclass
class Drop:
    contig: 'AlignedContig'
    covering: List['AlignedContig']


@dataclass
class IgnoreGap:
    contig: 'AlignedContig'
    gap: 'CigarHit'


@dataclass
class SplitGap:
    contig: 'AlignedContig'
    gap: 'CigarHit'
    left: 'AlignedContig'
    right: 'AlignedContig'


@dataclass
class Intro:
    contig: 'GenotypedContig'


@dataclass
class FinalCombine:
    contigs: List['AlignedContig']
    result: 'AlignedContig'


AlignmentEvent = Union[NoRef, ZeroHits, StrandConflict, HitNumber, ReverseComplement, Hit]
ModifyEvent = Union[LStrip, RStrip]
EventType = Union[Cut, ModifyEvent, Munge, Combine, AlignmentEvent, StitchCut, Overlap, NoOverlap, Stitch, Drop, IgnoreGap, SplitGap, Intro, FinalCombine]
