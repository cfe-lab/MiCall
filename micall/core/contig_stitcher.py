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


def drop_completely_covered(contigs: List[AlignedContig]) -> List[AlignedContig]:
    """ Filter out all contigs that are contained within other contigs. """

    # TODO: filter out if covered by multiple contigs
    # TODO: split contigs that have big gaps in them first, otherwise they will cover too much.

    def find_most_covered(contigs) -> Optional[AlignedContig]:
        for current in contigs:
            if any(x for x in contigs if x != current and x.contains(current)):
                return current

    while contigs:
        most_covered = find_most_covered(contigs)
        if most_covered:
            contigs.remove(most_covered)
        else:
            break

    return contigs


def split_contigs_with_gaps(contigs: List[AlignedContig]) -> Iterable[AlignedContig]:
    def covered_by(gap, contig):
        # TODO(vitalik): implement the more precise check
        possible_reference_coordinates = set(range(gap.r_st, gap.r_ei + 1))
        return possible_reference_coordinates \
                .issubset(contig.alignment.coordinate_mapping.reference_coordinates())

    def covered(contig, gap):
        return any(covered_by(gap, other) for other in contigs
                   if other != contig)

    def gap_boundaries(gap):
        midpoint = gap.r_st + (gap.r_ei - gap.r_st) / 2
        left_slice, right_slice = contig.cut_reference(floor(midpoint) + 0.5)
        left_closest_query = left_slice.alignment.coordinate_mapping.ref_to_closest_query(midpoint)
        right_closest_query = right_slice.alignment.coordinate_mapping.ref_to_closest_query(midpoint)
        left_closest_ref = left_slice.alignment.coordinate_mapping.query_to_ref(left_closest_query)
        right_closest_ref = right_slice.alignment.coordinate_mapping.query_to_ref(right_closest_query)
        return (left_closest_ref, right_closest_ref)

    def try_split(contig):
        for gap in contig.gaps():
            if covered(contig, gap):
                left_closest_ref, right_closest_ref = gap_boundaries(gap)
                left_part, left_gap = contig.cut_reference(left_closest_ref + contig.alignment.epsilon)
                right_gap, right_part = contig.cut_reference(right_closest_ref - contig.alignment.epsilon)

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
