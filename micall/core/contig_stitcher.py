from typing import Iterable, Optional, Tuple, List
from collections import deque
from dataclasses import dataclass
from mappy import Aligner
from functools import cached_property
from gotoh import align_it

from micall.utils.cigar_tools import connect_cigar_hits, CigarHit


@dataclass
class Contig:
    name: str
    seq: str


@dataclass
class GenotypedContig(Contig):
    ref_name: str
    ref_seq: str
    matched_fraction: Optional[float] # Approximated overall concordance between `seq` and `ref_seq`.

    @property
    def contig(self):
        return self

@dataclass
class AlignedContig:
    contig: GenotypedContig
    alignment: CigarHit

    def cut_reference(self, cut_point: float) -> Tuple['AlignedContig', 'AlignedContig']:
        """ Cuts this alignment in two parts with cut_point between them. """

        alignment_left, alignment_right = self.alignment.cut_reference(cut_point)
        return (AlignedContig(self.contig, alignment_left),
                AlignedContig(self.contig, alignment_right))


    @cached_property
    def msa(self):
        return self.alignment.to_msa(self.contig.ref_seq, self.contig.seq)


    @cached_property
    def seq(self):
        seq_left, ref_seq_left = self.msa
        return ''.join((c for c in ref_seq_left if c != '-'))


class FrankensteinContig(AlignedContig):
    """ Assembled of parts that were not even aligned together,
    and of some parts that were not aligned at all.
    Yet its .seq string looks like a real contig. """

    def __init__(self, parts: List[GenotypedContig]):
        self.parts = [subpart for part in parts for subpart in
                      (part.parts if isinstance(part, FrankensteinContig) else [part])]

        name = '+'.join(map(lambda acontig: acontig.contig.name, self.parts))
        ref = self.parts[0].contig
        contig = GenotypedContig(name=name, seq=self.seq,
                                 ref_name=ref.ref_name,
                                 ref_seq=ref.ref_seq,
                                 matched_fraction=ref.matched_fraction)

        alignment = connect_cigar_hits([part.alignment for part in self.parts
                                        if isinstance(part, AlignedContig)])

        super().__init__(contig, alignment)


    @cached_property
    def seq(self):
        return ''.join(map(lambda part: part.seq, self.parts))



def align_to_reference(contig: GenotypedContig):
    aligner = Aligner(seq=contig.ref_seq, preset='map-ont')
    alignments = list(aligner.map(contig.seq))
    if not alignments:
        return contig

    hits_array = [CigarHit(x.cigar, x.r_st, x.r_en - 1, x.q_st, x.q_en - 1) for x in alignments]
    single_cigar_hit = connect_cigar_hits(hits_array)
    return AlignedContig(contig=contig, alignment=single_cigar_hit)


def align_equal(seq1, seq2) -> Tuple[str, str]:
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


def interval_contains(x, y):
    """ Check if interval (x0, x1) contains interval (y0, y1). """
    return x[0] <= y[0] and x[1] >= y[1]


def find_all_overlapping_contigs(self, aligned_contigs):
    for other in aligned_contigs:
        if self.contig.ref_name != other.contig.ref_name:
            continue

        if self.alignment.overlaps(other.alignment):
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
            result[i] += scores_sum / window_size

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
    overlap_contig = GenotypedContig(name=f'overlap({left.contig.name},{right.contig.name})',
                                     seq=overlap_seq, ref_name=left.contig.ref_name,
                                     ref_seq=left.contig.ref_seq, matched_fraction=None)
    return FrankensteinContig([left_remainder, overlap_contig, right_remainder])


def stitch_contigs(contigs: Iterable[GenotypedContig]):
    aligned = list(map(align_to_reference, contigs))

    # Contigs that did not align do not need any more processing
    stitched = yield from (x for x in aligned if not isinstance(x, AlignedContig))
    aligned = [x for x in aligned if isinstance(x, AlignedContig)]

    # Going left-to-right through aligned contigs.
    aligned = list(sorted(aligned, key=lambda x: x.alignment.r_st))
    while aligned:
        current = aligned.pop(0)

        # Filter out all contigs that are contained within the current one.
        # TODO: actually filter out if covered by multiple contigs
        # TODO: split contigs that have big gaps in them first, otherwise they will cover too much.
        aligned = [x for x in aligned if not \
                   interval_contains((current.alignment.r_st, current.alignment.r_ei),
                                     (x.alignment.r_st, x.alignment.r_ei))]

        # Find overlap. If there isn't one - we are done with the current contig.
        overlapping_contig = find_overlapping_contig(current, aligned)
        if not overlapping_contig:
            yield current
            continue

        # Replace two contigs by their stitched version, then loop with it.
        new_contig = stitch_2_contigs(current, overlapping_contig)
        aligned.remove(overlapping_contig)
        aligned.insert(0, new_contig)
