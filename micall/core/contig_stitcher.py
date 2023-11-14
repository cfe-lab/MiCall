from typing import Iterable, Optional, Tuple, List
from collections import deque
from dataclasses import dataclass
from mappy import Aligner
from functools import cached_property, reduce
from itertools import accumulate, takewhile
from gotoh import align_it
from queue import LifoQueue
from math import floor

from micall.utils.cigar_tools import Cigar, connect_cigar_hits, CigarHit
from micall.utils.consensus_aligner import CigarActions


@dataclass
class Contig:
    name: str
    seq: str


@dataclass
class GenotypedContig(Contig):
    ref_name: str
    ref_seq: str
    matched_fraction: Optional[float] # Approximated overall concordance between `seq` and `ref_seq`.

    def align_to_reference(self):
        aligner = Aligner(seq=self.ref_seq, preset='map-ont')
        alignments = list(aligner.map(self.seq))
        if not alignments:
            return self

        hits_array = [CigarHit(x.cigar, x.r_st, x.r_en - 1, x.q_st, x.q_en - 1) for x in alignments]
        single_cigar_hit = connect_cigar_hits(hits_array)
        return AlignedContig(query=self, alignment=single_cigar_hit)


class AlignedContig(GenotypedContig):

    def __init__(self, query: GenotypedContig, alignment: CigarHit):
        self.alignment = alignment
        self.query = query

        ref_msa, query_msa = self.alignment.to_msa(self.query.ref_seq, self.query.seq)
        seq = ''.join((c for c in query_msa if c != '-'))

        super().__init__(
            seq = seq,
            name = query.name,
            ref_name = query.ref_name,
            ref_seq = query.ref_seq,
            matched_fraction = query.matched_fraction)


    def cut_reference(self, cut_point: float) -> Tuple['AlignedContig', 'AlignedContig']:
        """ Cuts this alignment in two parts with cut_point between them. """

        alignment_left, alignment_right = self.alignment.cut_reference(cut_point)
        return (AlignedContig(self.query, alignment_left),
                AlignedContig(self.query, alignment_right))


    def lstrip_query(self) -> 'AlignedContig':
        alignment = self.alignment.lstrip_query()
        return AlignedContig(self.query, alignment)


    def rstrip_query(self) -> 'AlignedContig':
        alignment = self.alignment.rstrip_query()
        return AlignedContig(self.query, alignment)


    def overlaps(self, other) -> bool:
        def intervals_overlap(x, y):
            return x[0] <= y[1] and x[1] >= y[0]

        if self.ref_name != other.ref_name:
            return False

        return intervals_overlap((self.alignment.r_st, self.alignment.r_ei),
                                 (other.alignment.r_st, other.alignment.r_ei))


    def contains(self, other) -> bool:
        def interval_contains(x, y):
            return x[0] <= y[0] and x[1] >= y[1]

        if self.ref_name != other.ref_name:
            return False

        return interval_contains((self.alignment.r_st, self.alignment.r_ei),
                                 (other.alignment.r_st, other.alignment.r_ei))


    def gaps(self) -> Iterable[CigarHit]:
        return self.alignment.gaps()


class SyntheticContig(AlignedContig):
    def __init__(self, query: GenotypedContig, r_st: int, r_ei: int):
        alignment = CigarHit.from_default_alignment(r_st=r_st, r_ei=r_ei,
                                                    q_st=0, q_ei=len(query.seq)-1)
        super().__init__(query, alignment)


    def cut_reference(self, cut_point: float):
        raise NotImplementedError("SyntheticContigs cannot be cut because they are not properly aligned")


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
        self.parts = [subpart for part in parts for subpart in
                      (part.parts if isinstance(part, FrankensteinContig) else [part])]

        aligned = reduce(FrankensteinContig.munge, self.parts)

        super().__init__(aligned.query, aligned.alignment)


    def cut_reference(self, cut_point: float) -> 'FrankensteinContig':
        # The cut_reference version of super() works here.
        # But it loses information about parts,
        # and does not check if the cut is legal
        # i.e. whether it slices a SyntheticContig.

        # Search for the part that needs to be cut:
        left_parts = list(takewhile(lambda part: cut_point >= part.alignment.r_ei + 1, self.parts))
        target_part = self.parts[len(left_parts)]
        right_parts = self.parts[len(left_parts) + 1:]

        target_part_left, target_part_right = target_part.cut_reference(cut_point)
        left = FrankensteinContig(left_parts + [target_part_left])
        right = FrankensteinContig([target_part_right] + right_parts)

        return (left, right)


    @staticmethod
    def munge(left: AlignedContig, right: AlignedContig) -> AlignedContig:
        left_query_seq = left.query.seq[0:left.alignment.q_ei + 1]
        right_query_seq = right.query.seq[right.alignment.q_st:]
        query_seq = left_query_seq + right_query_seq

        left_alignment = left.alignment
        right_alignment = \
            right.alignment.translate(
                query_delta=(-1 * right.alignment.q_st + len(left_query_seq)),
                reference_delta=0)
        alignment = left_alignment + right_alignment

        query = GenotypedContig(seq=query_seq,
                                name=f'{left.name}+{right.name}',
                                ref_name=left.ref_name,
                                ref_seq=left.ref_seq,
                                matched_fraction=None)
        return AlignedContig(query, alignment)


def align_equal(seq1: str, seq2: str) -> Tuple[str, str]:
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
    return max(every, key=lambda other: other.alignment.r_ei - other.alignment.r_st if other else 0,
               default=None)


def calculate_concordance(left: str, right: str) -> List[float]:
    """
    Calculate concordance for two given sequences using a sliding window method.

    The function compares the two strings from both left to right and then right to left,
    calculating for each position the ratio of matching characters in a window around the
    current position (10 characters to the left and right).

    It's required that the input strings are of the same length.

    :param left: string representing first sequence
    :param right: string representing second sequence
    :return: list representing concordance ratio for each position
    """

    result = [0] * len(left)

    assert len(left) == len(right), "Can only calculate concordance for same sized sequences"

    def slide(left, right):
        window_size = 10
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

    # Align overlapping parts, then recombine based on concordance.
    aligned_left, aligned_right = align_equal(left_overlap.seq, right_overlap.seq)
    concordance = calculate_concordance(aligned_left, aligned_right)
    max_concordance_index = max(range(len(concordance)),
                                key=lambda i: concordance[i])
    aligned_left_part = aligned_left[:max_concordance_index]
    aligned_right_part = aligned_right[max_concordance_index:]
    overlap_seq = ''.join(c for c in aligned_left_part + aligned_right_part if c != '-')

    # Return something that can be fed back into the loop.
    overlap_query = GenotypedContig(name=f'overlap({left.name},{right.name})',
                                    seq=overlap_seq, ref_name=left.ref_name,
                                    ref_seq=left.ref_seq, matched_fraction=None)
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
            yield current
            continue

        # Replace two contigs by their stitched version, then loop with it.
        new_contig = stitch_2_contigs(current, overlapping_contig)
        contigs.remove(overlapping_contig)
        contigs.insert(0, new_contig)


def merge_intervals(intervals: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    """
    Merge overlapping and adjacent intervals.
    Note that intervals are inclusive.

    :param intervals: A list of intervals [start, end] where 'start' and 'end' are integers.
    :eturn: A list of merged intervals.
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
        other_contigs = [x for x in contigs if x != current and x.ref_name == current.ref_name]
        cumulative_coverage = calculate_cumulative_coverage(other_contigs)

        # Check if the current contig is covered by the cumulative coverage intervals
        if any((cover_interval[0] <= current_interval[0] and cover_interval[1] >= current_interval[1])
               for cover_interval in cumulative_coverage):
            return current


def drop_completely_covered(contigs: List[AlignedContig]) -> List[AlignedContig]:
    """ Filter out all contigs that are contained within other contigs. """

    while contigs:
        covered = find_covered_contig(contigs)
        if covered:
            contigs.remove(covered)
        else:
            break

    return contigs


def split_contigs_with_gaps(contigs: List[AlignedContig]) -> Iterable[AlignedContig]:
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
        for gap in contig.gaps():
            if not significant(gap):
                # Really we do not want to split on every little deletion
                # because that would mean that we would need to stitch
                # overlaps around them.
                # And we are likely to lose quality with every stitching operation.
                # By skipping we assert that this gap is aligner's fault.
                continue

            if covered(contig, gap):
                midpoint = gap.r_st + (gap.r_ei - gap.r_st) / 2 + contig.alignment.epsilon
                left_part, right_part = contig.cut_reference(midpoint)
                left_part = left_part.rstrip_query()
                right_part = right_part.lstrip_query()

                contigs.remove(contig)
                contigs.append(left_part)
                contigs.append(right_part)
                process_queue.put(right_part)
                return

    process_queue = LifoQueue()
    for contig in contigs: process_queue.put(contig)

    while not process_queue.empty():
        contig = process_queue.get()
        try_split(contig)

    return contigs


def stitch_contigs(contigs: Iterable[GenotypedContig]):
    maybe_aligned = list(map(GenotypedContig.align_to_reference, contigs))

    # Contigs that did not align do not need any more processing
    yield from (x for x in maybe_aligned if not isinstance(x, AlignedContig))
    aligned = [x for x in maybe_aligned if isinstance(x, AlignedContig)]

    aligned = split_contigs_with_gaps(aligned)
    aligned = drop_completely_covered(aligned)

    yield from combine_overlaps(aligned)
