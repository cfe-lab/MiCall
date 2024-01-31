from typing import Iterable, Optional, Tuple, List, Dict, Union, Literal, TypeVar, Callable, Set
from collections import deque, defaultdict
from dataclasses import dataclass, replace
from math import ceil, floor
from mappy import Aligner
from functools import cached_property, reduce
from itertools import accumulate, takewhile, tee, islice, chain
from gotoh import align_it
from queue import LifoQueue
from Bio import Seq
import logging
from contextvars import ContextVar, Context
from contextlib import contextmanager
from fractions import Fraction

from micall.utils.cigar_tools import Cigar, connect_cigar_hits, CigarHit
from micall.utils.consensus_aligner import CigarActions
import micall.utils.contig_stitcher_events as events

T = TypeVar("T")
logger = logging.getLogger(__name__)

class StitcherContext:
    def __init__(self) -> None:
        self.name_generator_state: int = 0
        self.nameset: Set[str] = set()
        self.events: List[events.EventType] = []

    def generate_new_name(self) -> str:
        while True:
            self.name_generator_state += 1
            name = f"c{self.name_generator_state}"
            if name not in self.nameset:
                self.nameset.add(name)
                return name

    def emit(self, event: events.EventType) -> None:
        self.events.append(event)


    @staticmethod
    @contextmanager
    def fresh():
        ctx = StitcherContext()
        token = context.set(ctx)
        try:
            yield ctx
        finally:
            context.reset(token)


context: ContextVar[StitcherContext] = ContextVar("StitcherContext")


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

    def cut_query(self, cut_point: float) -> Tuple['GenotypedContig', 'GenotypedContig']:
        """ Cuts query sequence in two parts with cut_point between them. """

        cut_point = max(0, cut_point)
        left = replace(self, name=context.get().generate_new_name(), seq=self.seq[:ceil(cut_point)])
        right = replace(self, name=context.get().generate_new_name(), seq=self.seq[ceil(cut_point):])
        return (left, right)


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


    def cut_reference(self, cut_point: float) -> Tuple['AlignedContig', 'AlignedContig']:
        """ Cuts this alignment in two parts with cut_point between them. """

        alignment_left, alignment_right = self.alignment.cut_reference(cut_point)
        left = replace(self, name=context.get().generate_new_name(), alignment=alignment_left)
        right = replace(self, name=context.get().generate_new_name(), alignment=alignment_right)

        logger.debug("Created contigs %r at %s and %r at %s by cutting %r at %s at cut point = %s.",
                     left.name, left.alignment, right.name, right.alignment,
                     self.name, self.alignment, round(cut_point, 1))
        context.get().emit(events.Cut(self, left, right))

        return (left, right)


    def lstrip(self) -> 'AlignedContig':
        """
        Trims the query sequence of the contig from its beginning up to the start of the
        alignment. The CIGAR alignment is also updated to reflect the trimming.
        """

        alignment = self.alignment.lstrip_reference().lstrip_query()
        q_remainder, query = self.cut_query(alignment.q_st - 0.5)
        alignment = alignment.translate(0, -1 * alignment.q_st)
        result = AlignedContig.make(query, alignment, self.strand)
        logger.debug("Doing lstrip of %r at %s (len %s) resulted in %r at %s (len %s).",
                     self.name, self.alignment, len(self.seq),
                     result.name, result.alignment, len(result.seq))
        context.get().emit(events.LStrip(self, result))
        return result


    def rstrip(self) -> 'AlignedContig':
        """
        Trims the query sequence of the contig from its end based on the end of the
        alignment. The CIGAR alignment is also updated to reflect the trimming.
        """

        alignment = self.alignment.rstrip_reference().rstrip_query()
        query, q_remainder = self.cut_query(alignment.q_ei + 0.5)
        result = AlignedContig.make(query, alignment, self.strand)
        logger.debug("Doing rstrip of %r at %s (len %s) resulted in %r at %s (len %s).",
                     self.name, self.alignment, len(self.seq),
                     result.name, result.alignment, len(result.seq))
        context.get().emit(events.RStrip(self, result))
        return result


    def overlaps(self, other) -> bool:
        def intervals_overlap(x, y):
            return x[0] <= y[1] and x[1] >= y[0]

        if self.group_ref != other.group_ref:
            return False

        return intervals_overlap((self.alignment.r_st, self.alignment.r_ei),
                                 (other.alignment.r_st, other.alignment.r_ei))


    def munge(self, other: 'AlignedContig') -> 'AlignedContig':
        """
        Combines two adjacent contigs into a single contig by joining their
        query sequences and alignments.
        """

        match_fraction = min(self.match_fraction, other.match_fraction)
        ref_name = max([self, other], key=lambda x: x.alignment.ref_length).ref_name
        query = GenotypedContig(seq=self.seq + other.seq,
                                name=context.get().generate_new_name(),
                                ref_name=ref_name,
                                group_ref=self.group_ref,
                                ref_seq=self.ref_seq,
                                match_fraction=match_fraction)

        self_alignment = self.alignment
        other_alignment = \
            other.alignment.translate(
                query_delta=(-1 * other.alignment.q_st + self.alignment.q_ei + 1),
                reference_delta=0)
        alignment = self_alignment.connect(other_alignment)

        assert self.strand == other.strand
        ret = AlignedContig.make(query=query, alignment=alignment, strand=self.strand)
        logger.debug("Munged contigs %r at %s with %r at %s resulting in %r at %s.",
                     self.name, self.alignment, other.name, other.alignment, ret.name, ret.alignment)
        context.get().emit(events.Munge(self, other, ret))
        return ret


def sliding_window(sequence: Iterable[T]) -> Iterable[Tuple[Optional[T], T, Optional[T]]]:
    """
    Generate a three-element sliding window of a sequence.

    Each element generated contains a tuple with the previous item (None if the first item),
    the current item, and the next item (None if the last item) in the sequence.
    """

    a, b, c = tee(sequence, 3)
    prevs = chain([None], a)
    nexts = chain(islice(c, 1, None), [None])
    return zip(prevs, b, nexts)


def combine_contigs(parts: List[AlignedContig]) -> AlignedContig:
    """
    Combine a list of contigs into a single AlignedContig by trimming and merging overlapping parts.

    Left-trimming and right-trimming occur at any shared overlapping points
    between adjacent parts. AlignedContig.munge() is used to combine contiguous parts without overlap.
    """

    stripped_parts = []
    for prev_part, part, next_part in sliding_window(parts):
        if prev_part is not None:
            part = part.lstrip()
        if next_part is not None:
            part = part.rstrip()
        stripped_parts.append(part)

    ret = reduce(AlignedContig.munge, stripped_parts)
    logger.debug("Created a frankenstein %r at %s (len %s) from %s.",
                 ret.name, ret.alignment, len(ret.seq),
                 [f"{x.name!r} at {x.alignment} (len {len(x.seq)})" for x in stripped_parts])
    context.get().emit(events.Combine(stripped_parts, ret))
    return ret


def align_to_reference(contig: GenotypedContig) -> Iterable[GenotypedContig]:
    """
    Align a single Contig to its reference sequence, producing potentially multiple aligned contigs.

    If the reference sequence (ref_seq) is unavailable, the contig is returned unaltered.
    Otherwise, alignments are performed and contigs corresponding to each alignment are yielded.
    """

    if contig.ref_seq is None:
        logger.debug("Contig %r not aligned - no reference.", contig.name)
        context.get().emit(events.NoRef(contig))
        yield contig
        return

    aligner = Aligner(seq=contig.ref_seq, preset='map-ont')
    alignments = list(aligner.map(contig.seq))
    hits_array: List[Tuple[CigarHit, Literal["forward", "reverse"]]] = \
        [(CigarHit(Cigar(x.cigar),
                   min(x.r_st, x.r_en - 1), max(x.r_st, x.r_en - 1),
                   min(x.q_st, x.q_en - 1), max(x.q_st, x.q_en - 1)),
          "forward" if x.strand == 1 else "reverse") for x in alignments]

    connected = connect_cigar_hits([hit for hit, strand in hits_array]) if hits_array else []
    if not connected:
        logger.debug("Contig %r not aligned - backend's choice.", contig.name)
        context.get().emit(events.ZeroHits(contig))
        yield contig
        return

    if len(set(map(lambda p: p[1], hits_array))) > 1:
        logger.debug("Discarding contig %r because it aligned both in forward and reverse sense.", contig.name)
        context.get().emit(events.StrandConflict(contig))
        yield contig
        return

    logger.debug("Contig %r produced %s aligner hits. After connecting them, the number became %s.",
                 contig.name, len(hits_array), len(connected))
    context.get().emit(events.HitNumber(contig, hits_array, connected))

    strand = hits_array[0][1]
    if strand == "reverse":
        rc = str(Seq(contig.seq).reverse_complement())
        new_contig = replace(contig, seq=rc)
        logger.debug("Reverse complemented contig %r.", contig.name)
        context.get().emit(events.ReverseComplement(contig, new_contig))
        contig = new_contig

    for i, single_hit in enumerate(connected):
        query = replace(contig, name=context.get().generate_new_name())
        part = AlignedContig.make(query, single_hit, strand)

        logger.debug("Part %r of contig %r aligned as %r at %s%s.", i, contig.name,
                     part.name,part.alignment, " (rev)" if strand == "reverse" else "")
        context.get().emit(events.Hit(contig, part, i))
        yield part


def strip_conflicting_mappings(contigs: Iterable[GenotypedContig]) -> Iterable[GenotypedContig]:
    contigs = list(contigs)
    names = {contig.name: contig for contig in contigs}

    def get_indexes(name: str) -> Tuple[int, int]:
        contig = names[name]
        if isinstance(contig, AlignedContig):
            return (contig.alignment.q_st, contig.alignment.r_st)
        else:
            return (-1, -1)

    reference_sorted = list(sorted(names.keys(), key=lambda name: get_indexes(name)[1]))
    query_sorted = list(sorted(names.keys(), key=lambda name: get_indexes(name)[0]))

    def is_out_of_order(name: str) -> bool:
        return reference_sorted.index(name) != query_sorted.index(name)

    sorted_by_query = sorted(contigs, key=lambda contig: get_indexes(contig.name))
    for prev_contig, contig, next_contig in sliding_window(sorted_by_query):
        if isinstance(contig, AlignedContig):
            original = contig
            start = prev_contig.alignment.q_ei + 1 if isinstance(prev_contig, AlignedContig) else 0
            end = next_contig.alignment.q_st - 1 if isinstance(next_contig, AlignedContig) else len(contig.seq) - 1

            if prev_contig is not None or is_out_of_order(original.name):
                contig = contig.lstrip()
                context.get().emit(events.InitialStrip(original, start, original.alignment.q_st - 1))
            if next_contig is not None or is_out_of_order(original.name):
                contig = contig.rstrip()
                context.get().emit(events.InitialStrip(original, original.alignment.q_ei + 1, end))

        yield contig


def align_all_to_reference(contigs: Iterable[GenotypedContig]) -> Iterable[GenotypedContig]:
    """
    Align multiple contigs to their respective reference sequences.

    Applies align_to_reference to each contig in the given collection,
    flattening the result into a single list.
    """

    groups = map(align_to_reference, contigs)
    groups = map(strip_conflicting_mappings, groups)
    for group in groups:
        yield from group


def align_queries(seq1: str, seq2: str) -> Tuple[str, str]:
    """ Globally align two query sequences against each other and return the resulting aligned sequences in MSA format. """

    gap_open_penalty = 15
    gap_extend_penalty = 3
    use_terminal_gap_penalty = 1
    aseq1, aseq2, score = \
        align_it(
            seq1, seq2,
            gap_open_penalty,
            gap_extend_penalty,
            use_terminal_gap_penalty)

    return aseq1, aseq2


def find_all_overlapping_contigs(self, aligned_contigs):
    """"
    Yield all contigs from a collection that overlap with a given contig.
    Contigs are considered overlapping if they have overlapping intervals on the same reference genome.
    """

    for other in aligned_contigs:
        if self.overlaps(other):
            yield other


def find_overlapping_contig(self, aligned_contigs):
    """
    Find the single contig in a collection that overlaps the most with a given contig.
    It returns the contig with the maximum overlapped reference length with the given contig (self).
    """

    every = find_all_overlapping_contigs(self, aligned_contigs)
    return max(every, key=lambda other: other.alignment.ref_length if other else 0, default=None)


def calculate_concordance(left: str, right: str) -> List[Fraction]:
    """
    Calculate concordance for two given sequences using a sliding average.

    The function compares the two strings character by character, simultaneously from
    both left to right and right to left, calculating a score that represents a moving
    average of matches at each position. If characters match at a given position,
    a score of 1 is added; otherwise, a score of 0 is added. The score is then
    averaged with the previous scores using a weighted sliding average where the
    current score has a weight of 1/3 and the accumulated score has a weight of 2/3.
    This sliding average score is halved and then processed again, but in reverse direction.

    :param left: string representing first sequence
    :param right: string representing second sequence
    :return: list representing concordance ratio for each position
    """

    if len(left) != len(right):
        raise ValueError("Can only calculate concordance for same sized sequences")

    result: List[Fraction] = [Fraction(0)] * len(left)

    def slide(start, end):
        scores_sum = Fraction(0)
        inputs = list(zip(left, right))
        increment = 1 if start <= end else -1

        for i in range(start, end, increment):
            (a, b) = inputs[i]
            current = Fraction(1) if a == b else Fraction(0)
            scores_sum = (scores_sum * 2 / 3 + current * 1 / 3)
            result[i] += scores_sum / 2

    # Slide forward, then in reverse, adding the scores at each position.
    slide(0, len(left))
    slide(len(left) - 1, -1)

    return result


def disambiguate_concordance(concordance: List[float]) -> Iterable[Tuple[float, int]]:
    for i, x in enumerate(concordance):
        global_rank = i if i < len(concordance) / 2 else len(concordance) - i - 1
        yield (x, global_rank)


def concordance_to_cut_points(left_overlap, right_overlap, aligned_left, aligned_right, concordance):
    """ Determine optimal cut points for stitching based on sequence concordance in the overlap region. """

    concordance_d = list(disambiguate_concordance(concordance))
    sorted_concordance_indexes = sorted(range(len(concordance)), key=lambda i: concordance_d[i])
    remove_dashes = lambda s: ''.join(c for c in s if c != '-')

    for max_concordance_index in reversed(sorted_concordance_indexes):
        aligned_left_q_index = len(remove_dashes(aligned_left[:max_concordance_index]))
        aligned_right_q_index = right_overlap.alignment.query_length - len(remove_dashes(aligned_right[max_concordance_index:])) + 1
        aligned_left_r_index = left_overlap.alignment.coordinate_mapping.query_to_ref.left_max(aligned_left_q_index)
        if aligned_left_r_index is None:
            aligned_left_r_index = left_overlap.alignment.r_st - 1
        aligned_right_r_index = right_overlap.alignment.coordinate_mapping.query_to_ref.right_min(aligned_right_q_index)
        if aligned_right_r_index is None:
            aligned_right_r_index = right_overlap.alignment.r_ei + 1
        if aligned_right_r_index > aligned_left_r_index:
            return (aligned_left_r_index + 0.5, aligned_right_r_index - 0.5, max_concordance_index)

    return (left_overlap.alignment.r_st - 1 + 0.5, right_overlap.alignment.r_ei + 1 - 0.5, 0)


def stitch_2_contigs(left, right):
    """
    Stitch two contigs together into a single coherent contig.

    The function handles the overlap by cutting both contigs into segments, aligning the
    overlapping segments, and then choosing the optimal stitching points based on sequence
    concordance. Non-overlapping segments are retained as is.
    """

    # Cut in 4 parts.
    left_remainder, left_overlap = left.cut_reference(right.alignment.r_st - 0.5)
    right_overlap, right_remainder = right.cut_reference(left.alignment.r_ei + 0.5)
    left_overlap = left_overlap.rstrip().lstrip()
    right_overlap = right_overlap.lstrip().rstrip()
    left_remainder = left_remainder.rstrip()
    right_remainder = right_remainder.lstrip()

    logger.debug("Stitching %r at %s (len %s) with %r at %s (len %s)."
                 " The left_overlap %r is at %s (len %s)"
                 " and the right_overlap %r is at %s (len %s).",
                 left.name, left.alignment, len(left.seq),
                 right.name, right.alignment, len(right.seq),
                 left_overlap.name, left_overlap.alignment, len(left_overlap.seq),
                 right_overlap.name, right_overlap.alignment, len(right_overlap.seq))
    context.get().emit(events.StitchCut(left, right, left_overlap, right_overlap, left_remainder, right_remainder))

    # Align overlapping parts, then recombine based on concordance.
    aligned_left, aligned_right = align_queries(left_overlap.seq, right_overlap.seq)
    concordance = calculate_concordance(aligned_left, aligned_right)
    aligned_left_cutpoint, aligned_right_cutpoint, max_concordance_index = \
        concordance_to_cut_points(left_overlap, right_overlap, aligned_left, aligned_right, concordance)
    left_overlap_take, left_overlap_drop = left_overlap.cut_reference(aligned_left_cutpoint)
    right_overlap_drop, right_overlap_take = right_overlap.cut_reference(aligned_right_cutpoint)

    # Log it.
    average_concordance = Fraction(sum(concordance) / (len(concordance) or 1))
    concordance_str = ', '.join(map(lambda x: str(int(round(x * 100)) / 100), concordance))
    cut_point_location_scaled = max_concordance_index / (((len(concordance) or 1) - 1) or 1)
    logger.debug("Created overlap contigs %r at %s and %r at %s based on parts of %r and %r, with avg. concordance %s%%, cut point at %s%%, and full concordance [%s].",
                 left_overlap_take.name, left_overlap.alignment, right_overlap_take.name, right_overlap_take.alignment,
                 left.name, right.name, round(average_concordance * 100),
                 round(cut_point_location_scaled * 100), concordance_str)
    context.get().emit(events.Overlap(left, right, left_overlap, right_overlap,
                                      left_remainder, right_remainder, left_overlap_take,
                                      right_overlap_take, concordance, average_concordance,
                                      max_concordance_index, cut_point_location_scaled))

    return combine_contigs([left_remainder, left_overlap_take, right_overlap_take, right_remainder])


def combine_overlaps(contigs: List[AlignedContig]) -> Iterable[AlignedContig]:
    """"
    Repeatedly combine all overlapping aligned contigs into an iterable collection of contiguous AlignedContigs.
    It proceeds by iterating through sorted contigs and stitching any overlapping ones until none are left.
    """

    # Going left-to-right through aligned contigs.
    contigs = list(sorted(contigs, key=lambda x: x.alignment.r_st))
    while contigs:
        current = contigs.pop(0)

        # Find overlap. If there isn't one - we are done with the current contig.
        overlapping_contig = find_overlapping_contig(current, contigs)
        if not overlapping_contig:
            logger.debug("Nothing overlaps with %r.", current.name)
            context.get().emit(events.NoOverlap(current))
            yield current
            continue

        # Replace two contigs by their stitched version, then loop with it.
        new_contig = stitch_2_contigs(current, overlapping_contig)
        contigs.remove(overlapping_contig)
        contigs.insert(0, new_contig)

        logger.debug("Stitching %r with %r results in %r at %s (len %s).",
                     current.name, overlapping_contig.name,
                     new_contig.name, new_contig.alignment, len(new_contig.seq))
        context.get().emit(events.Stitch(current, overlapping_contig, new_contig))


def merge_intervals(intervals: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    """
    Merge overlapping and adjacent intervals.
    Note that intervals are inclusive.

    :param intervals: A list of intervals [start, end] where 'start' and 'end' are integers.
    :return: A list of merged intervals.
    """

    if not intervals:
        return []

    # Sort intervals by their starting values
    sorted_intervals = sorted(intervals, key=lambda x: x[0])

    merged_intervals = [sorted_intervals[0]]
    for current in sorted_intervals[1:]:
        current_start, current_end = current
        last_start, last_end = merged_intervals[-1]
        if current_start <= last_end + 1:
            # Extend the last interval if there is an overlap or if they are adjacent
            merged_intervals[-1] = (min(last_start, current_start), max(last_end, current_end))
        else:
            # Add the current interval if there is no overlap
            merged_intervals.append(current)

    return merged_intervals


def find_covered_contig(contigs: List[AlignedContig]) -> Tuple[Optional[AlignedContig], List[AlignedContig]]:
    """
    Find and return the first contig that is completely covered by other contigs.

    :param contigs: List of all aligned contigs to be considered.
    :return: An AlignedContig if there is one completely covered by others, None otherwise.
    """

    def calculate_cumulative_coverage(contigs) -> List[Tuple[int, int]]:
        intervals = [(contig.alignment.r_st, contig.alignment.r_ei) for contig in contigs]
        merged_intervals = merge_intervals(intervals)
        return merged_intervals

    for current in contigs:
        current_interval = (current.alignment.r_st, current.alignment.r_ei)

        # Create a map of cumulative coverage for contigs
        overlaping_contigs = [x for x in contigs if x != current and x.overlaps(current)]
        cumulative_coverage = calculate_cumulative_coverage(overlaping_contigs)

        # Check if the current contig is covered by the cumulative coverage intervals
        if any((cover_interval[0] <= current_interval[0] and cover_interval[1] >= current_interval[1])
               for cover_interval in cumulative_coverage):
            return current, overlaping_contigs

    return None, []


def drop_completely_covered(contigs: List[AlignedContig]) -> List[AlignedContig]:
    """ Filter out all contigs that are contained within other contigs. """

    contigs = contigs[:]
    while contigs:
        covered, covering = find_covered_contig(contigs)
        if covered:
            contigs.remove(covered)
            logger.debug("Droped contig %r as it is completely covered by these contigs: %s.",
                         covered.name, ", ".join(repr(x.name) for x in covering))
            context.get().emit(events.Drop(covered, covering))
        else:
            break

    return contigs


def split_contigs_with_gaps(contigs: List[AlignedContig]) -> List[AlignedContig]:
    """
    Split contigs at large gaps if those gaps are covered by other contigs in the list.

    A gap within a contig is considered large based on a pre-defined threshold. If another contig aligns
    within that gap's range, the contig is split into two around the midpoint of the gap.
    """

    def covered_by(gap, other):
        # Check if any 1 reference coordinate in gap is mapped in other.
        gap_coords = gap.coordinate_mapping.ref_to_query.domain
        cover_coords = set(other.alignment.coordinate_mapping.ref_to_query.keys())
        return not gap_coords.isdisjoint(cover_coords)

    def covered(contig, gap):
        return any(covered_by(gap, other) for other in contigs if other != contig)

    def significant(gap):
        return gap.ref_length > 5

    def try_split(contig):
        for gap in contig.alignment.deletions():
            if not significant(gap):
                # Really we do not want to split on every little deletion
                # because that would mean that we would need to stitch
                # overlaps around them.
                # And we are likely to lose quality with every stitching operation.
                # By skipping we assert that this gap is aligner's fault.
                logger.debug("Ignored insignificant gap of %r, %s.", contig.name, gap)
                context.get().emit(events.IgnoreGap(contig, gap))
                continue

            if covered(contig, gap):
                midpoint = gap.r_st + (gap.r_ei - gap.r_st) / 2 + contig.alignment.epsilon
                left_part, right_part = contig.cut_reference(midpoint)
                left_part = left_part.rstrip()
                right_part = right_part.lstrip()

                contigs.remove(contig)
                contigs.append(left_part)
                contigs.append(right_part)
                process_queue.put(right_part)

                logger.debug("Split contig %r at %s around its gap at [%s, %s]->[%s, %s]. "
                             "Left part: %r at %s, "
                             "right part: %r at %s.",
                             contig.name, contig.alignment,
                             gap.q_st, gap.q_ei, gap.r_st, gap.r_ei,
                             left_part.name, left_part.alignment,
                             right_part.name, right_part.alignment)
                context.get().emit(events.SplitGap(contig, gap, left_part, right_part))
                return

    process_queue: LifoQueue = LifoQueue()
    for contig in contigs: process_queue.put(contig)

    while not process_queue.empty():
        contig = process_queue.get()
        try_split(contig)

    return contigs


def stitch_contigs(contigs: Iterable[GenotypedContig]) -> Iterable[GenotypedContig]:
    contigs = list(contigs)
    for contig in contigs:
        logger.debug("Introduced contig %r (seq = %s) of ref %r, group_ref %r (seq = %s), and length %s.",
                     contig.name, contig.seq, contig.ref_name,
                     contig.group_ref, contig.ref_seq, len(contig.seq))
        context.get().emit(events.Intro(contig))
        context.get().nameset.add(contig.name)

    maybe_aligned = list(align_all_to_reference(contigs))

    # Contigs that did not align do not need any more processing
    yield from (x for x in maybe_aligned if not isinstance(x, AlignedContig))
    aligned = [x for x in maybe_aligned if isinstance(x, AlignedContig)]

    aligned = split_contigs_with_gaps(aligned)
    aligned = drop_completely_covered(aligned)
    yield from combine_overlaps(aligned)


GroupRef = str

def stitch_consensus(contigs: Iterable[GenotypedContig]) -> Iterable[GenotypedContig]:
    contigs = list(stitch_contigs(contigs))
    consensus_parts: Dict[GroupRef, List[AlignedContig]] = defaultdict(list)

    for contig in contigs:
        if isinstance(contig, AlignedContig) and contig.strand == "forward":
            consensus_parts[contig.group_ref].append(contig)
        else:
            yield contig

    def combine(group_ref):
        contigs = sorted(consensus_parts[group_ref], key=lambda x: x.alignment.r_st)
        result = combine_contigs(contigs)
        logger.debug("Combining these contigs for final output for %r: %s.",
                     group_ref, [f"{x.name!r} at {x.alignment} (len {len(x.seq)})" for x in contigs])
        context.get().emit(events.FinalCombine(contigs, result))
        return result

    yield from map(combine, consensus_parts)


def main(args):
    import argparse
    from micall.core.denovo import write_contig_refs # TODO(vitalik): move denovo stuff here.

    parser = argparse.ArgumentParser()
    parser.add_argument('contigs', type=argparse.FileType('r'))
    parser.add_argument('stitched_contigs', type=argparse.FileType('w'))
    parser.add_argument('--plot')
    verbosity_group = parser.add_mutually_exclusive_group()
    verbosity_group.add_argument('--verbose', action='store_true', help='Increase output verbosity')
    verbosity_group.add_argument('--no-verbose', action='store_true', help='Normal output verbosity', default=True)
    verbosity_group.add_argument('--debug', action='store_true', help='Maximum output verbosity')
    verbosity_group.add_argument('--quiet', action='store_true', help='Minimize output verbosity')
    args = parser.parse_args(args)

    if args.quiet:
        logger.setLevel(logging.ERROR)
    elif args.verbose:
        logger.setLevel(logging.INFO)
    elif args.debug:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.WARN)

    logging.basicConfig(level=logger.level)
    with StitcherContext.fresh():
        write_contig_refs(args.contigs.name, args.stitched_contigs, stitcher_plot_path=args.plot)
        args.contigs.close()
        args.stitched_contigs.close()


if __name__ == '__main__':
    import sys
    main(sys.argv[1:])
