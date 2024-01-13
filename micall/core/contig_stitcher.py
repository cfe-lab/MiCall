from typing import Iterable, Optional, Tuple, List, Dict, Union, Literal
from collections import deque, defaultdict
from dataclasses import dataclass
from math import ceil, floor
from mappy import Aligner
from functools import cached_property, reduce
from itertools import accumulate, takewhile
from gotoh import align_it
from queue import LifoQueue
import logging

from micall.utils.cigar_tools import Cigar, connect_cigar_hits, CigarHit
from micall.utils.consensus_aligner import CigarActions


logger = logging.getLogger(__name__)


name_generator_state = 0
def generate_new_name():
    global name_generator_state
    name_generator_state += 1
    return f"c{name_generator_state}"


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
        left_name = generate_new_name()
        left = GenotypedContig(name=left_name,
                               seq=self.seq[:ceil(cut_point)],
                               ref_seq=self.ref_seq,
                               ref_name=self.ref_name,
                               group_ref=self.group_ref,
                               match_fraction=match_fraction)
        right_name = generate_new_name()
        right = GenotypedContig(name=right_name,
                                seq=self.seq[ceil(cut_point):],
                                ref_seq=self.ref_seq,
                                ref_name=self.ref_name,
                                group_ref=self.group_ref,
                                match_fraction=match_fraction)

        return (left, right)


    def rename(self, new_name: str) -> 'GenotypedContig':
        return GenotypedContig(
            name=new_name,
            seq=self.seq,
            ref_name=self.ref_name,
            group_ref=self.group_ref,
            ref_seq=self.ref_seq,
            match_fraction=self.match_fraction)


@dataclass
class AlignedContig(GenotypedContig):
    query: GenotypedContig
    alignment: CigarHit
    reverse: bool

    def __init__(self,
                 query: GenotypedContig,
                 alignment: CigarHit,
                 reverse: bool):
        self.query = query
        self.alignment = alignment
        self.reverse = reverse
        super().__init__(
            seq=query.seq,
            name=query.name,
            ref_name=query.ref_name,
            group_ref=query.group_ref,
            ref_seq=query.ref_seq,
            match_fraction=query.match_fraction)


    def modify(self, query: GenotypedContig, alignment: CigarHit) -> 'AlignedContig':
        if query.seq == self.query.seq and alignment == self.alignment:
            return self
        return AlignedContig(
            reverse=self.reverse,
            query=query,
            alignment=alignment)


    def cut_reference(self, cut_point: float) -> Tuple['AlignedContig', 'AlignedContig']:
        """ Cuts this alignment in two parts with cut_point between them. """

        alignment_left, alignment_right = self.alignment.cut_reference(cut_point)
        left_query = self.query.rename(generate_new_name())
        right_query = self.query.rename(generate_new_name())
        left = self.modify(left_query, alignment_left)
        right = self.modify(right_query, alignment_right)

        logger.debug("Created contigs %r at %s and %r at %s by cutting %r.",
                     left.name, left.alignment, right.name, right.alignment, self.name,
                     extra={"action": "cut", "original": self,
                            "left": left, "right": right})

        return (left, right)


    def lstrip_query(self) -> 'AlignedContig':
        alignment = self.alignment.lstrip_query()
        q_remainder, query = self.query.cut_query(alignment.q_st - 0.5)
        alignment = alignment.translate(0, -1 * alignment.q_st)
        result = self.modify(query, alignment)
        logger.debug("Contig %r morfed into contig %r, so %s became %s",
                     self.name, result.name, self.alignment, result.alignment,
                     extra={"action": "modify", "original": self, "result": result})
        return result


    def rstrip_query(self) -> 'AlignedContig':
        alignment = self.alignment.rstrip_query()
        query, q_remainder = self.query.cut_query(alignment.q_ei + 0.5)
        result = self.modify(query, alignment)
        logger.debug("Contig %r morfed into contig %r, so %s became %s",
                     self.name, result.name, self.alignment, result.alignment,
                     extra={"action": "modify", "original": self, "result": result})
        return result


    def overlaps(self, other) -> bool:
        def intervals_overlap(x, y):
            return x[0] <= y[1] and x[1] >= y[0]

        if self.group_ref != other.group_ref:
            return False

        return intervals_overlap((self.alignment.r_st, self.alignment.r_ei),
                                 (other.alignment.r_st, other.alignment.r_ei))


    def munge(self, other: 'AlignedContig') -> 'AlignedContig':
        query_seq = self.rstrip_query().seq + other.lstrip_query().seq
        match_fraction = min(self.match_fraction, other.match_fraction)
        ref_name = max([self, other], key=lambda x: x.alignment.ref_length).ref_name
        query = GenotypedContig(seq=query_seq,
                                name=generate_new_name(),
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

        assert self.reverse == other.reverse
        ret = AlignedContig(reverse=self.reverse, query=query, alignment=alignment)
        logger.debug("Munged contigs %r at %s with %r at %s resulting in %r at %s.",
                     self.name, self.alignment, other.name, other.alignment,
                     ret.name, ret.alignment, extra={"action": "munge", "left": self,
                                                     "right": other, "result": ret})
        return ret


def combine_contigs(parts: List[AlignedContig]) -> AlignedContig:
    ret = reduce(AlignedContig.munge, parts)
    logger.debug("Created a frankenstein %r at %s (len %s) from %s.",
                 ret.name, ret.alignment, len(ret.seq),
                 [f"{x.name!r} at {x.alignment} (len {len(x.seq)})" for x in parts],
                 extra={"action": "frankenstein", "contigs": parts, "result": ret})
    return ret


def align_to_reference(contig) -> Iterable[GenotypedContig]:
    if contig.ref_seq is None:
        logger.info("Contig %r not aligned - no reference.", contig.name,
                    extra={"action": "alignment", "type": "noref", "contig": contig})
        yield contig
        return

    aligner = Aligner(seq=contig.ref_seq, preset='map-ont')
    alignments = list(aligner.map(contig.seq))
    hits_array = [(CigarHit(Cigar(x.cigar), x.r_st, x.r_en - 1, x.q_st, x.q_en - 1), x.strand == -1)
                  for x in alignments]
    reversed_alignments = [alignment for alignment, is_rev in hits_array if is_rev]
    alignments = [alignment for alignment, is_rev in hits_array if not is_rev]

    logger.info("Contig %r produced %s reverse-complement alignments.",
                contig.name, len(reversed_alignments),
                extra={"action": "alignment", "type": "reversenumber",
                       "contig": contig, "n": len(reversed_alignments)})

    connected = connect_cigar_hits(alignments) if alignments else []

    logger.info("Contig %r produced %s forward alignments.", contig.name, len(connected),
                extra={"action": "alignment", "type": "hitnumber",
                       "contig": contig, "n": len(connected)})

    def logpart(i, part, is_rev):
        logger.info("Part %r of contig %r aligned as %r at [%s, %s]->[%s, %s]%s.",
                    i, contig.name, part.name, part.alignment.q_st,
                    part.alignment.q_ei, part.alignment.r_st, part.alignment.r_ei,
                    " (rev)" if is_rev else "",
                    extra={"action": "alignment", "type": "hit",
                           "contig": contig, "part": part, "i": i})
        logger.debug("Part %r of contig %r aligned as %r at %s%s.", i, contig.name,
                     part.name, part.alignment, " (rev)" if is_rev else "")

    def make_aligned(query, alignment, is_rev):
        return AlignedContig(
            query=query,
            alignment=alignment,
            reverse=is_rev)

    to_return = connected + reversed_alignments
    if len(to_return) == 0:
        logger.info("Contig %r not aligned - backend choice.", contig.name,
                    extra={"action": "alignment", "type": "zerohits", "contig": contig})
        yield contig
        return

    if len(to_return) == 1:
        is_rev = to_return[0] in reversed_alignments
        part = make_aligned(contig, to_return[0], is_rev)
        logpart(0, part, is_rev)
        yield part
        return

    for i, single_hit in enumerate(to_return):
        query = GenotypedContig(name=generate_new_name(),
                                seq=contig.seq,
                                ref_name=contig.ref_name,
                                group_ref=contig.group_ref,
                                ref_seq=contig.ref_seq,
                                match_fraction=contig.match_fraction)
        is_rev = single_hit in reversed_alignments
        part = make_aligned(query, single_hit, is_rev)
        logpart(i, part, is_rev)
        yield part


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
    current position. So position holds a moving avarage score.

    It's required that the input strings are of the same length.

    :param left: string representing first sequence
    :param right: string representing second sequence
    :return: list representing concordance ratio for each position
    """

    if len(left) != len(right):
        raise ValueError("Can only calculate concordance for same sized sequences")

    result: List[float] = [0] * len(left)

    def slide(start, end):
        window_size = 30
        scores = deque([0] * window_size, maxlen=window_size)
        scores_sum = 0
        inputs = list(zip(left, right))
        increment = 1 if start <= end else -1

        for i in range(start, end, increment):
            (a, b) = inputs[i]
            current = a == b
            scores_sum -= scores.popleft()
            scores_sum += current
            scores.append(current)
            result[i] += (scores_sum / window_size) / 2

    # Slide forward, then in reverse, adding the scores at each position.
    slide(0, len(left))
    slide(len(left) - 1, -1)

    return result


def stitch_2_contigs(left, right):
    # Cut in 4 parts.
    left_remainder, left_overlap = left.cut_reference(right.alignment.r_st - 0.5)
    right_overlap, right_remainder = right.cut_reference(left.alignment.r_ei + 0.5)
    left_overlap = left_overlap.rstrip_query().lstrip_query()
    right_overlap = right_overlap.lstrip_query().rstrip_query()
    left_remainder = left_remainder.rstrip_query()
    right_remainder = right_remainder.lstrip_query()

    logger.debug("Stitching %r at %s (len %s) with %r at %s (len %s)."
                 " The left_overlap %r is at %s (len %s)"
                 " and the right_overlap %r is at %s (len %s).",
                 left.name, left.alignment, len(left.seq),
                 right.name, right.alignment, len(right.seq),
                 left_overlap.name, left_overlap.alignment, len(left_overlap.seq),
                 right_overlap.name, right_overlap.alignment, len(right_overlap.seq),
                 extra={"action": "stitchcut", "left": left, "right": right,
                        "left_overlap": left_overlap, "right_overlap": right_overlap,
                        "left_remainder": left_remainder, "right_remainder": right_remainder})

    # Align overlapping parts, then recombine based on concordance.
    aligned_left, aligned_right = align_queries(left_overlap.seq, right_overlap.seq)
    concordance = calculate_concordance(aligned_left, aligned_right)
    valuator = lambda i: (concordance[i], i if i < len(concordance) / 2 else len(concordance) - i - 1)
    max_concordance_index = max(range(len(concordance)), key=valuator)

    # Return something that can be fed back into the loop.
    without_dashes = lambda s: ''.join(c for c in s if c != '-')
    aligned_left_q_index = len(without_dashes(aligned_left[:max_concordance_index]))
    aligned_right_q_index = right_overlap.alignment.query_length - len(without_dashes(aligned_right[max_concordance_index:])) + 1
    aligned_left_r_index = left_overlap.alignment.coordinate_mapping.query_to_ref.left_max(aligned_left_q_index)
    if aligned_left_r_index is None:
        aligned_left_r_index = left_overlap.alignment.r_st - 1
    aligned_right_r_index = right_overlap.alignment.coordinate_mapping.query_to_ref.right_min(aligned_right_q_index)
    if aligned_right_r_index is None:
        aligned_right_r_index = right_overlap.alignment.r_ei + 1
    left_overlap_take, left_overlap_drop = left_overlap.cut_reference(aligned_left_r_index + 0.5)
    right_overlap_drop, right_overlap_take = right_overlap.cut_reference(aligned_right_r_index - 0.5)

    # Log it.
    average_concordance = sum(concordance) / (len(concordance) or 1)
    concordance_str = ', '.join(map(lambda x: str(round(x, 2)), concordance)),
    cut_point_location_scaled = max_concordance_index / (((len(concordance) or 1) - 1) or 1)
    logger.debug("Created overlap contigs %r at %s and %r at %s based on parts of %r and %r, with avg. concordance %s%%, cut point at %s%%, and full concordance [%s].",
                 left_overlap_take.name, left_overlap.alignment, right_overlap_take.name, right_overlap_take.alignment,
                 left.name, right.name, round(average_concordance * 100),
                 round(cut_point_location_scaled * 100), concordance_str,
                 extra={"action": "overlap", "left": left, "right": right,
                        "left_remainder": left_remainder, "right_remainder": right_remainder,
                        "left_overlap": left_overlap, "right_original": right_overlap,
                        "left_take": left_overlap_take, "right_take": right_overlap_take,
                        "concordance": concordance, "avg": average_concordance,
                        "cut_point": max_concordance_index,
                        "cut_point_scaled": cut_point_location_scaled})

    return combine_contigs([left_remainder, left_overlap_take, right_overlap_take, right_remainder])


def combine_overlaps(contigs: List[AlignedContig]) -> Iterable[AlignedContig]:
    # Going left-to-right through aligned contigs.
    contigs = list(sorted(contigs, key=lambda x: x.alignment.r_st))
    while contigs:
        current = contigs.pop(0)

        # Find overlap. If there isn't one - we are done with the current contig.
        overlapping_contig = find_overlapping_contig(current, contigs)
        if not overlapping_contig:
            logger.info("Nothing overlaps with %r.", current.name,
                        extra={"action": "nooverlap", "contig": current})
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
        logger.debug("Stitching %r with %r results in %r at %s (len %s).",
                     current.name, overlapping_contig.name,
                     new_contig.name, new_contig.alignment, len(new_contig.seq))


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
        logger.debug("Introduced contig %r (seq = %s) of ref %r, group_ref %r (seq = %s), and length %s.",
                     contig.name, contig.seq, contig.ref_name,
                     contig.group_ref, contig.ref_seq, len(contig.seq))

    aligned = align_all_to_reference(contigs)

    # Contigs that did not align do not need any more processing
    yield from (x for x in aligned if not isinstance(x, AlignedContig))
    aligned = [x for x in aligned if isinstance(x, AlignedContig)]

    # Contigs aligned in reverse do not need any more processing
    yield from (x for x in aligned if x.reverse)
    aligned = [x for x in aligned if not x.reverse]

    aligned = split_contigs_with_gaps(aligned)
    aligned = drop_completely_covered(aligned)
    yield from combine_overlaps(aligned)


GroupRef = str

def stitch_consensus(contigs: Iterable[GenotypedContig]) -> Iterable[GenotypedContig]:
    contigs = list(stitch_contigs(contigs))
    consensus_parts: Dict[GroupRef, List[AlignedContig]] = defaultdict(list)

    for contig in contigs:
        if isinstance(contig, AlignedContig) and not contig.reverse:
            consensus_parts[contig.group_ref].append(contig)
        else:
            yield contig

    def combine(group_ref):
        contigs = sorted(consensus_parts[group_ref], key=lambda x: x.alignment.r_st)
        ret = combine_contigs(contigs)
        logger.info("Combined these contigs for final output for %r: %s,"
                     " resulting in %r at [%s, %s]->[%s, %s].", group_ref,
                     [repr(x.name) for x in contigs],
                     ret.name, ret.alignment.q_st, ret.alignment.q_ei,
                     ret.alignment.r_st, ret.alignment.r_ei,
                     extra={"action": "finalcombine", "contigs": contigs, "result": ret})
        return ret

    yield from map(combine, consensus_parts)


def main(args):
    import argparse
    from micall.core.denovo import write_contig_refs # TODO(vitalik): move denovo stuff here.

    parser = argparse.ArgumentParser()
    parser.add_argument('contigs', type=argparse.FileType('r'))
    parser.add_argument('stitched_contigs', type=argparse.FileType('w'))
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

    write_contig_refs(args.contigs.name, args.stitched_contigs)
    args.contigs.close()
    args.stitched_contigs.close()


if __name__ == '__main__':
    import sys
    main(sys.argv[1:])
