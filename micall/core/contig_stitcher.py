from typing import Iterable, Optional, Tuple, List, Dict
from collections import deque, defaultdict
from dataclasses import dataclass
from math import ceil
from mappy import Aligner
from functools import cached_property, reduce
from itertools import accumulate, takewhile
from gotoh import align_it
from queue import LifoQueue
from math import floor
import logging

from micall.utils.cigar_tools import Cigar, connect_cigar_hits, CigarHit
from micall.utils.consensus_aligner import CigarActions
from micall.utils.structured_logger import register_structured_logger


logger = logging.getLogger(__name__)
register_structured_logger(logger)

@dataclass
class Contig:
    name: str
    seq: str


@dataclass
class GenotypedContig(Contig):
    ref_name: str
    group_ref: str
    ref_seq: Optional[str]          # The sequence of self.group_ref. None in cases where the reference organism is unknown.
    match_fraction: float           # Approximated overall concordance between `seq` and `ref_seq`. It is calculated by BLAST as qcovhsp/100, where qcovhsp means Query Coverage Per HSP.

    def cut_query(self, cut_point: float) -> Tuple['GenotypedContig', 'GenotypedContig']:
        """
        Cuts this alignment in two parts with cut_point between them.
        Reference sequence is kept untouched.
        """

        cut_point = max(0, cut_point)
        match_fraction = self.match_fraction
        left = GenotypedContig(name=f'left({self.name})',
                               seq=self.seq[:ceil(cut_point)],
                               ref_seq=self.ref_seq,
                               ref_name=self.ref_name,
                               group_ref=self.group_ref,
                               match_fraction=match_fraction)
        right = GenotypedContig(name=f'right({self.name})',
                                seq=self.seq[ceil(cut_point):],
                                ref_seq=self.ref_seq,
                                ref_name=self.ref_name,
                                group_ref=self.group_ref,
                                match_fraction=match_fraction)

        return (left, right)


@dataclass
class AlignedContig(GenotypedContig):
    query: GenotypedContig
    alignment: CigarHit

    def __init__(self, query: GenotypedContig, alignment: CigarHit):
        self.query = query
        self.alignment = alignment
        super().__init__(
            seq = query.seq,
            name = query.name,
            ref_name = query.ref_name,
            group_ref = query.group_ref,
            ref_seq = query.ref_seq,
            match_fraction = query.match_fraction)


    def cut_reference(self, cut_point: float) -> Tuple['AlignedContig', 'AlignedContig']:
        """ Cuts this alignment in two parts with cut_point between them. """

        alignment_left, alignment_right = self.alignment.cut_reference(cut_point)
        query_left, query_right = self.query.cut_query(alignment_left.q_ei + 0.5)
        alignment_right = alignment_right.translate(0, -1 * alignment_right.q_st)

        return (AlignedContig(query_left, alignment_left),
                AlignedContig(query_right, alignment_right))


    def lstrip_query(self) -> 'AlignedContig':
        alignment = self.alignment.lstrip_query()
        q_remainder, query = self.query.cut_query(alignment.q_st - 0.5)
        alignment = alignment.translate(0, -1 * alignment.q_st)
        return AlignedContig(query, alignment)


    def rstrip_query(self) -> 'AlignedContig':
        alignment = self.alignment.rstrip_query()
        query, q_remainder = self.query.cut_query(alignment.q_ei + 0.5)
        return AlignedContig(query, alignment)


    def overlaps(self, other) -> bool:
        def intervals_overlap(x, y):
            return x[0] <= y[1] and x[1] >= y[0]

        if self.group_ref != other.group_ref:
            return False

        return intervals_overlap((self.alignment.r_st, self.alignment.r_ei),
                                 (other.alignment.r_st, other.alignment.r_ei))


class SyntheticContig(AlignedContig):
    """
    Contig that is not really aligned, but its boundaries are known.
    It is created as a result of overlaps between the real contigs.
    """
    # TODO(vitalik): maybe it is worth to realign overlaps to get rid of this special-case class.

    def __init__(self, query: GenotypedContig, r_st: int, r_ei: int):
        alignment = CigarHit.from_default_alignment(r_st=r_st, r_ei=r_ei,
                                                    q_st=0, q_ei=len(query.seq)-1)
        super().__init__(query, alignment)


    def cut_reference(self, cut_point: float):
        raise NotImplementedError("SyntheticContigs cannot be cut because they are not properly aligned")


    def lstrip_query(self):
        return self


    def rstrip_query(self):
        return self


class FrankensteinContig(AlignedContig):
    """
    Assembled of parts that were not even aligned together,
    and of some parts that were not aligned at all.
    Yet its self.seq string looks like a real contig.
    """

    def __init__(self, parts: List[AlignedContig]):
        if len(parts) == 0:
            raise ValueError("Empty Frankenstei do not exist")

        # Flatten any possible Frankenstein parts
        self.parts: List[AlignedContig] = \
            [subpart for part in parts for subpart in
             (part.parts if isinstance(part, FrankensteinContig) else [part])]

        aligned = reduce(FrankensteinContig.munge, self.parts)

        super().__init__(aligned.query, aligned.alignment)


    def cut_reference(self, cut_point: float) -> Tuple['FrankensteinContig', 'FrankensteinContig']:
        # Search for the part that needs to be cut
        left_parts = list(takewhile(lambda part: cut_point >= part.alignment.r_ei + 1, self.parts))
        target_part = self.parts[len(left_parts)]
        right_parts = self.parts[len(left_parts) + 1:]

        # Cut the target part and add its pieces to left and right.
        target_part_left, target_part_right = target_part.cut_reference(cut_point)
        left = FrankensteinContig(left_parts + [target_part_left])
        right = FrankensteinContig([target_part_right] + right_parts)

        return (left, right)


    def lstrip_query(self):
        return FrankensteinContig([self.parts[0].lstrip_query()] + self.parts[1:])


    def rstrip_query(self):
        return FrankensteinContig(self.parts[:-1] + [self.parts[-1].rstrip_query()])


    @staticmethod
    def munge(left: AlignedContig, right: AlignedContig) -> AlignedContig:
        query_seq = left.rstrip_query().seq + right.lstrip_query().seq
        match_fraction = min(left.match_fraction, right.match_fraction)
        ref_name = max([left, right], key=lambda x: x.alignment.ref_length).ref_name
        query = GenotypedContig(seq=query_seq,
                                name=f'{left.name}+{right.name}',
                                ref_name=ref_name,
                                group_ref=left.group_ref,
                                ref_seq=left.ref_seq,
                                match_fraction=match_fraction)

        left_alignment = left.alignment
        right_alignment = \
            right.alignment.translate(
                query_delta=(-1 * right.alignment.q_st + left.alignment.q_ei + 1),
                reference_delta=0)
        alignment = left_alignment.connect(right_alignment)

        return AlignedContig(query, alignment)


def align_to_reference(contig) -> Iterable[GenotypedContig]:
    if contig.ref_seq is None:
        logger.info("Contig %r not aligned - no reference.", contig.name,
                    extra={"action": "alignment", "type": "noref", "contig": contig})
        yield contig
        return

    aligner = Aligner(seq=contig.ref_seq, preset='map-ont')
    alignments = list(aligner.map(contig.seq))
    if not alignments:
        logger.info("Contig %r not aligned - backend choice.", contig.name,
                    extra={"action": "alignment", "type": "zerohits", "contig": contig})
        yield contig
        return

    hits_array = [CigarHit(Cigar.coerce(x.cigar), x.r_st, x.r_en - 1, x.q_st, x.q_en - 1) for x in alignments]
    connected = connect_cigar_hits(hits_array)

    logger.info("Contig %r aligned in %s parts.", contig.name, len(connected),
                extra={"action": "alignment", "type": "hitnumber",
                       "contig": contig, "n": len(connected)})

    def logpart(i, part):
        logger.info("Part %r of contig %r aligned as [%s, %s]->[%s, %s].",
                    i, contig.name, part.q_st, part.q_ei, part.r_st, part.r_ei,
                    extra={"action": "alignment", "type": "hit",
                           "contig": contig, "part": part, "i": i})
        logger.debug("Part %r of contig %r aligned as %s.", i, contig.name, part)

    if len(connected) == 1:
        logpart(0, connected[0])
        yield AlignedContig(query=contig, alignment=connected[0])
        return

    for i, single_hit in enumerate(connected):
        logpart(i, single_hit)
        query = GenotypedContig(name=f'part({i}, {contig.name})',
                                seq=contig.seq,
                                ref_name=contig.ref_name,
                                group_ref=contig.group_ref,
                                ref_seq=contig.ref_seq,
                                match_fraction=contig.match_fraction)
        yield AlignedContig(query=query, alignment=single_hit)


def align_all_to_reference(contigs):
    return [contig for parts in map(align_to_reference, contigs) for contig in parts]


def align_queries(seq1: str, seq2: str) -> Tuple[str, str]:
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
    for other in aligned_contigs:
        if self.overlaps(other):
            yield other


def find_overlapping_contig(self, aligned_contigs):
    every = find_all_overlapping_contigs(self, aligned_contigs)
    return max(every, key=lambda other: other.alignment.ref_length if other else 0,
               default=None)


def calculate_concordance(left: str, right: str) -> List[float]:
    """
    Calculate concordance for two given sequences using a sliding window method.

    The function compares the two strings from both left to right and then right to left,
    calculating for each position the ratio of matching characters in a window around the
    current position.

    It's required that the input strings are of the same length.

    :param left: string representing first sequence
    :param right: string representing second sequence
    :return: list representing concordance ratio for each position
    """

    if len(left) != len(right):
        raise ValueError("Can only calculate concordance for same sized sequences")

    result: List[float] = [0] * len(left)

    def slide(left, right):
        window_size = 30
        scores = deque([0] * window_size, maxlen=window_size)
        scores_sum = 0

        for i, (a, b) in enumerate(zip(left, right)):
            current = a == b
            scores_sum -= scores.popleft()
            scores_sum += current
            scores.append(current)
            result[i] += (scores_sum / window_size) / 2

    # Slide forward, then in reverse, adding the scores at each position.
    slide(left, right)
    slide(reversed(left), reversed(right))

    return result


def stitch_2_contigs(left, right):
    # Cut in 4 parts.
    left_remainder, left_overlap = left.cut_reference(right.alignment.r_st - 0.5)
    right_overlap, right_remainder = right.cut_reference(left.alignment.r_ei + 0.5)
    left_overlap = left_overlap.rstrip_query()
    right_overlap = right_overlap.lstrip_query()

    # Align overlapping parts, then recombine based on concordance.
    aligned_left, aligned_right = align_queries(left_overlap.seq, right_overlap.seq)
    concordance = calculate_concordance(aligned_left, aligned_right)
    max_concordance_index = max(range(len(concordance)), key=lambda i: concordance[i])
    aligned_left_part = aligned_left[:max_concordance_index]
    aligned_right_part = aligned_right[max_concordance_index:]
    overlap_seq = ''.join(c for c in aligned_left_part + aligned_right_part if c != '-')

    average_concordance = sum(concordance) / (len(concordance) or 1)
    logger.debug("Average concordance between overlapping parts of %r and %r is %s (full is %s).",
                 left.name, right.name, average_concordance, concordance,
                 extra={"action": "concordance", "left": left, "right": right,
                        "value": concordance, "avg": average_concordance})

    # Return something that can be fed back into the loop.
    match_fraction = min(left.match_fraction, right.match_fraction)
    ref_name = max([left, right], key=lambda x: x.alignment.ref_length).ref_name
    overlap_query = GenotypedContig(name=f'overlap({left.name},{right.name})',
                                    ref_name=ref_name,
                                    seq=overlap_seq, group_ref=left.group_ref,
                                    ref_seq=left.ref_seq, match_fraction=match_fraction)
    overlap_contig = SyntheticContig(overlap_query,
                                     r_st=left_overlap.alignment.r_st,
                                     r_ei=right_overlap.alignment.r_ei)

    return FrankensteinContig([left_remainder, overlap_contig, right_remainder])


def combine_overlaps(contigs: List[AlignedContig]) -> Iterable[AlignedContig]:
    # Going left-to-right through aligned contigs.
    contigs = list(sorted(contigs, key=lambda x: x.alignment.r_st))
    while contigs:
        current = contigs.pop(0)

        # Find overlap. If there isn't one - we are done with the current contig.
        overlapping_contig = find_overlapping_contig(current, contigs)
        if not overlapping_contig:
            logger.info("Nothing overlaps with %r.",
                        current.name, extra={"action": "nooverlap", "contig": current})
            yield current
            continue

        # Replace two contigs by their stitched version, then loop with it.
        new_contig = stitch_2_contigs(current, overlapping_contig)
        contigs.remove(overlapping_contig)
        contigs.insert(0, new_contig)

        logger.info("Stitching %r with %r results in %r at [%s,%s]->[%s,%s].",
                    current.name, overlapping_contig.name,
                    new_contig.name, new_contig.alignment.q_st, new_contig.alignment.q_ei,
                    new_contig.alignment.r_st, new_contig.alignment.r_ei,
                    extra={"action": "stitch", "result": new_contig,
                           "left": current, "right": overlapping_contig})
        logger.debug("Stitching %r with %r results in %r at %s.",
                     current.name, overlapping_contig.name,
                     new_contig.name, new_contig.alignment)


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


def find_covered_contig(contigs: List[AlignedContig]) -> Optional[AlignedContig]:
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
        other_contigs = [x for x in contigs if x != current and x.group_ref == current.group_ref]
        cumulative_coverage = calculate_cumulative_coverage(other_contigs)

        # Check if the current contig is covered by the cumulative coverage intervals
        if any((cover_interval[0] <= current_interval[0] and cover_interval[1] >= current_interval[1])
               for cover_interval in cumulative_coverage):
            return current

    return None


def drop_completely_covered(contigs: List[AlignedContig]) -> List[AlignedContig]:
    """ Filter out all contigs that are contained within other contigs. """

    contigs = contigs[:]
    while contigs:
        covered = find_covered_contig(contigs)
        if covered:
            contigs.remove(covered)
            logger.info("Droped contig %r as it is completely covered by other contigs.",
                        covered.name, extra={"action": "drop", "contig": covered})
        else:
            break

    return contigs


def split_contigs_with_gaps(contigs: List[AlignedContig]) -> List[AlignedContig]:
    def covered_by(gap, other):
        # Check if any 1 reference coordinate in gap is mapped in other.
        gap_coords = gap.coordinate_mapping.ref_to_query.domain
        cover_coords = set(other.alignment.coordinate_mapping.ref_to_query.keys())
        return not gap_coords.isdisjoint(cover_coords)

    def covered(contig, gap):
        return any(covered_by(gap, other) for other in contigs
                   if other != contig)

    def significant(gap):
        return gap.ref_length > 5

    def try_split(contig):
        for gap in contig.alignment.gaps():
            if not significant(gap):
                # Really we do not want to split on every little deletion
                # because that would mean that we would need to stitch
                # overlaps around them.
                # And we are likely to lose quality with every stitching operation.
                # By skipping we assert that this gap is aligner's fault.
                logger.debug("Ignored insignificant gap of %r, %s.", contig.name, gap,
                             extra={"action": "ignoregap", "contig": contig, "gap": gap})
                continue

            if covered(contig, gap):
                midpoint = gap.r_st + (gap.r_ei - gap.r_st) // 2 + contig.alignment.epsilon
                left_part, right_part = contig.cut_reference(midpoint)
                left_part = left_part.rstrip_query()
                right_part = right_part.lstrip_query()

                contigs.remove(contig)
                contigs.append(left_part)
                contigs.append(right_part)
                process_queue.put(right_part)

                logger.info("Split contig %r around its gap at [%s, %s]->[%s, %s]. "
                            "Left part: %r at [%s, %s]->[%s, %s], "
                            "right part: %r at [%s, %s]->[%s, %s].",
                            contig.name, gap.q_st, gap.q_ei, gap.r_st, gap.r_ei,
                            left_part.name, left_part.alignment.q_st, left_part.alignment.q_ei,
                            left_part.alignment.r_st, left_part.alignment.r_ei,
                            right_part.name, right_part.alignment.q_st, right_part.alignment.q_ei,
                            right_part.alignment.r_st, right_part.alignment.r_ei,
                            extra={"action": "splitgap", "contig": contig,
                                   "gap": gap, "left": left_part, "right": right_part})
                return

    process_queue: LifoQueue = LifoQueue()
    for contig in contigs: process_queue.put(contig)

    while not process_queue.empty():
        contig = process_queue.get()
        try_split(contig)

    return contigs


def stitch_contigs(contigs: Iterable[GenotypedContig]) -> Iterable[AlignedContig]:
    contigs = list(contigs)
    for contig in contigs:
        logger.info("Introduced contig %r of ref %r, group_ref %r, and length %s.",
                    contig.name, contig.ref_name, contig.group_ref, len(contig.seq),
                    extra={"action": "intro", "contig": contig})

    maybe_aligned = align_all_to_reference(contigs)

    # Contigs that did not align do not need any more processing
    yield from (x for x in maybe_aligned if not isinstance(x, AlignedContig))
    aligned: List[AlignedContig] = \
        [x for x in maybe_aligned if isinstance(x, AlignedContig)]

    aligned = split_contigs_with_gaps(aligned)
    aligned = drop_completely_covered(aligned)
    yield from combine_overlaps(aligned)


GroupRef = str

def stitch_consensus(contigs: Iterable[GenotypedContig]) -> Iterable[GenotypedContig]:
    contigs = list(stitch_contigs(contigs))
    consensus_parts: Dict[GroupRef, List[AlignedContig]] = defaultdict(list)

    for contig in contigs:
        if isinstance(contig, AlignedContig):
            consensus_parts[contig.group_ref].append(contig)
        else:
            yield contig

    def combine(contigs):
        contigs = sorted(contigs, key=lambda x: x.alignment.r_st)
        return FrankensteinContig(contigs)

    yield from map(combine, consensus_parts.values())
