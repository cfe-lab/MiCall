from typing import Tuple, List, Sequence, Optional, Iterable, Iterator
from dataclasses import dataclass
from operator import attrgetter
from itertools import groupby

from aligntools import CigarActions, Cigar, CigarHit, connect_cigar_hits
from gotoh import align_it
from mappy import Aligner
import mappy


#
# Alignments with deletions larger than MAX_GAP_SIZE
# will be split around those deletions into multiple
# separate alignments.
#
MAX_GAP_SIZE = 600  # TODO: make this smaller?


@dataclass(frozen=True)
class Alignment:
    """
    Our representation of mappy's Alignment object.
    """

    ctg: str
    ctg_len: int
    r_st: int
    r_en: int
    strand: int
    q_st: int
    q_en: int
    mapq: int
    cigar: Sequence[Tuple[int, CigarActions]]
    cigar_str: str

    @staticmethod
    def coerce(obj: object) -> 'Alignment':
        if isinstance(obj, Alignment):
            return obj
        elif isinstance(obj, mappy.Alignment):
            cigar: List[Tuple[int, CigarActions]] = []
            for (size, action) in obj.cigar:
                cigar.append((size, CigarActions(action)))

            return Alignment(ctg=obj.ctg,
                             ctg_len=obj.ctg_len,
                             r_st=obj.r_st, r_en=obj.r_en,
                             strand=obj.strand,
                             q_st=obj.q_st, q_en=obj.q_en,
                             mapq=obj.mapq,
                             cigar=cigar,
                             cigar_str=obj.cigar_str,
                             )
        else:
            raise TypeError(f"Cannot coerce from {obj!r}.")

    def to_cigar_hit(self) -> CigarHit:
        return CigarHit(Cigar(self.cigar),
                        r_st=self.r_st, r_ei=self.r_en - 1,
                        q_st=self.q_st, q_ei=self.q_en - 1)

    @staticmethod
    def from_cigar_hit(hit: CigarHit, ctg: str, ctg_len: int, strand: int, mapq: int) -> 'Alignment':
        return Alignment(ctg=ctg,
                         ctg_len=ctg_len,
                         r_st=hit.r_st, r_en=hit.r_ei + 1,
                         strand=strand,
                         q_st=hit.q_st, q_en=hit.q_ei + 1,
                         mapq=mapq,
                         cigar=list(hit.cigar._data),
                         cigar_str=str(hit.cigar),
                         )


def align_gotoh(coordinate_seq: str, consensus: str) -> Optional[Alignment]:
    gap_open_penalty = 15
    gap_extend_penalty = 3
    use_terminal_gap_penalty = 1
    assert '&' not in consensus, "Consensus contains forbidden character '&'"
    consensus = ''.join('&' if x == '-' else x for x in consensus)
    aligned_coordinate, aligned_consensus, score = align_it(
        coordinate_seq,
        consensus,
        gap_open_penalty,
        gap_extend_penalty,
        use_terminal_gap_penalty)

    if min(len(coordinate_seq), len(consensus)) < score:
        cigar = Cigar.from_msa(aligned_coordinate, aligned_consensus)
        hit = CigarHit(cigar,
                       q_st=0, q_ei=len(consensus)-1,
                       r_st=0, r_ei=len(coordinate_seq)-1)
        hit = hit.lstrip_query().lstrip_reference().rstrip_query().rstrip_reference()
        return Alignment.from_cigar_hit(
            hit,
            ctg='N/A',
            ctg_len=len(coordinate_seq),
            strand=1,
            mapq=0)
    else:
        return None


def alignment_quality(alignment: Alignment) -> Tuple[int, ...]:
    cigar = Cigar(alignment.cigar)
    mlen = sum(1 for action in cigar.iterate_operations()
               if action == CigarActions.MATCH)
    return (alignment.mapq * cigar.query_length, mlen, cigar.query_length)


def connect_alignments(alignments: Iterable[Alignment]) -> Iterator[Alignment]:
    stranded = groupby(alignments, key=lambda x: (x.strand, x.ctg, x.ctg_len))
    for (strand, ctg, ctg_len), group_iter in stranded:
        group = list(group_iter)
        hits = list(map(Alignment.to_cigar_hit, group))
        connected_hits = connect_cigar_hits(hits)
        mapq = min(x.mapq for x in group)
        for hit in connected_hits:
            yield Alignment.from_cigar_hit(hit,
                                           ctg=ctg, ctg_len=ctg_len,
                                           strand=strand, mapq=mapq)


def collect_big_gaps_cut_points(alignment: Alignment) -> Iterator[float]:
    hit = alignment.to_cigar_hit()
    for deletion in hit.deletions():
        if deletion.ref_length > MAX_GAP_SIZE:
            midpoint = deletion.r_st + deletion.ref_length / 2
            yield int(midpoint) + hit.epsilon


def cut_hit_into_multiple_parts(hit: CigarHit, cut_points: Iterable[float]) -> Iterator[CigarHit]:
    for cut_point in cut_points:
        left, right = hit.cut_reference(cut_point)
        left = left.rstrip_reference()
        right = right.lstrip_reference()
        yield left
        hit = right
    yield hit


def split_around_big_gaps(alignments: Iterable[Alignment]) -> Iterator[Alignment]:
    for alignment in alignments:
        cut_points = list(collect_big_gaps_cut_points(alignment))
        if cut_points:
            hit = alignment.to_cigar_hit()
            for part in cut_hit_into_multiple_parts(hit, cut_points):
                yield Alignment.from_cigar_hit(part,
                                               ctg=alignment.ctg,
                                               ctg_len=alignment.ctg_len,
                                               strand=alignment.strand,
                                               mapq=alignment.mapq)
        else:
            yield alignment


def align_consensus(coordinate_seq: str, consensus: str) -> Tuple[List[Alignment], str]:
    aligner = Aligner(seq=coordinate_seq, bw=500, bw_long=500, preset='map-ont')
    mappy_alignments: List[mappy.Alignment] = list(aligner.map(consensus))
    if mappy_alignments or 10_000 < len(consensus):
        algorithm = 'minimap2'
        alignments = [Alignment.coerce(alignment)
                      for alignment in mappy_alignments
                      if alignment.is_primary]

        # Following code will connect non-overlapping alignments
        # that mappy outputs sometimes.
        # It will also drop overlapping (in query coords) alignments.
        # We are sorting the alignments before connect in order
        # to drop the lowest quality contigs in case they overlap with
        # higher quality alignments.
        alignments.sort(key=alignment_quality)
        alignments = list(connect_alignments(reversed(alignments)))
    else:
        algorithm = 'gotoh'
        gotoh_alignment = align_gotoh(coordinate_seq, consensus)
        if gotoh_alignment:
            alignments = [gotoh_alignment]
        else:
            alignments = []

    alignments = list(split_around_big_gaps(alignments))
    alignments.sort(key=attrgetter('q_st'))
    return (alignments, algorithm)
