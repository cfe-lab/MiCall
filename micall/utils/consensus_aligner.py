from typing import Dict, List, Optional, Set, Iterator, Iterable, Tuple
from dataclasses import dataclass, replace
from itertools import count, groupby
from operator import attrgetter
import csv
import os
import logging
from aligntools import CigarActions, Cigar, CigarHit, connect_nonoverlapping_cigar_hits

from gotoh import align_it, align_it_aa
from mappy import Aligner
import mappy

from micall.core.project_config import ProjectConfig
from micall.utils.report_amino import SeedAmino, ReportAmino, ReportNucleotide, SeedNucleotide
from micall.utils.translation import translate
from micall.utils.alignment import Alignment

logger = logging.getLogger(__name__)

# A section between reading frame shifts must be at least this long to be split.
MINIMUM_READING_FRAME_SHIFT = 30
# Minimum aminos in a section between reading frame shifts to be split.
MINIMUM_AMINO_ALIGNMENT = 10
# Most codons in an insertion or deletion that is still aligned in amino acids.
MAXIMUM_AMINO_GAP = 10


#
# Alignments with deletions larger than MAX_GAP_SIZE
# will be split around those deletions into multiple
# separate alignments.
#
MAX_GAP_SIZE = 600  # TODO: make this smaller?



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
        cigar = cigar.relax()  # turn '=' and 'X' into 'M'.
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


def connect_alignments(alignments: Iterable[Alignment]) -> Iterator[Alignment]:
    stranded = groupby(alignments, key=lambda x: (x.strand, x.ctg, x.ctg_len))
    for (strand, ctg, ctg_len), group_iter in stranded:
        group = list(group_iter)
        hits = list(map(Alignment.to_cigar_hit, group))
        connected_hits = connect_nonoverlapping_cigar_hits(hits)
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


def align_aminos(reference: str,
                 query: str,
                 gap_open: int = 40,
                 gap_extend: int = 10,
                 use_terminal_gap_penalty: int = 1):
    """ Align two amino acid sequences.

    """
    aligned_ref, aligned_query, score = align_it_aa(
        reference,
        query,
        gap_open,
        gap_extend,
        use_terminal_gap_penalty)
    return aligned_ref, aligned_query, score


def map_amino_sequences(from_seq: str, to_seq: str):
    from_aligned, to_aligned, _score = align_aminos(from_seq, to_seq)
    seq_map = {}
    from_index = to_index = 0
    for from_aa, to_aa in zip(from_aligned, to_aligned):
        if (to_index < len(to_seq) and
                from_index < len(from_seq) and
                to_aa == to_seq[to_index]):
            seq_map[from_index] = to_index
            to_index += 1
        if (from_index < len(from_seq) and
                from_aa == from_seq[from_index]):
            from_index += 1
    return seq_map


class ConsensusAligner:
    def __init__(self,
                 projects: ProjectConfig,
                 alignments_file=None,
                 unmerged_alignments_file=None,
                 intermediate_alignments_file=None,
                 overall_alignments_file=None,
                 seed_concordance_file=None,
                 contig_name=None):
        self.projects = projects
        self.coordinate_name = self.consensus = self.amino_consensus = ''
        self.algorithm = ''
        self.consensus_offset = 0
        self.alignments: List[Alignment] = []
        self.reading_frames: Dict[int, List[SeedAmino]] = {}
        self.seed_nucs: List[SeedNucleotide] = []
        self.amino_alignments: List[AminoAlignment] = []
        self.contig_name = contig_name

        # consensus nucleotide positions that were inserts
        self.inserts: Set[int] = set()

        if alignments_file is not None:
            self.alignments_writer = self._create_alignments_writer(alignments_file)
        else:
            self.alignments_writer = None

        if unmerged_alignments_file is not None:
            self.unmerged_alignments_writer = self._create_alignments_writer(unmerged_alignments_file)
        else:
            self.unmerged_alignments_writer = None

        if intermediate_alignments_file is not None:
            self.intermediate_alignments_writer = self._create_alignments_writer(intermediate_alignments_file)
        else:
            self.intermediate_alignments_writer = None

        if overall_alignments_file is not None:
            columns = ["coordinate_name",
                       "contig",
                       "query_start",
                       "query_end",
                       "consensus_offset",
                       "ref_start",
                       "ref_end",
                       "cigar_str"]
            self.overall_alignments_writer = \
                self._create_alignments_writer(overall_alignments_file, different_columns=columns)
        else:
            self.overall_alignments_writer = None

        if seed_concordance_file is not None:
            columns = ["seed_name",
                       "contig",
                       "region",
                       "pct_concordance",
                       "pct_covered"]
            self.seed_concordance_writer = self._create_alignments_writer(seed_concordance_file,
                                                                          different_columns=columns)
        else:
            self.seed_concordance_writer = None

    @staticmethod
    def _create_alignments_writer(alignments_file, different_columns=None):
        columns = different_columns or ["coordinate_name",
                                        "contig",
                                        "action",
                                        "query_start",
                                        "query_end",
                                        "ref_start",
                                        "ref_end",
                                        "reading_frame",
                                        "ref_amino_start",
                                        "aligned_query",
                                        "aligned_ref"]
        writer = csv.DictWriter(alignments_file, columns, lineterminator=os.linesep)
        pos = alignments_file.tell()
        if pos == 0:
            writer.writeheader()
        return writer

    def start_contig(self,
                     coordinate_name: Optional[str] = None,
                     consensus: Optional[str] = None,
                     reading_frames: Optional[Dict[int, List[SeedAmino]]] = None):
        self.clear()

        if consensus:
            self.consensus = consensus
        elif reading_frames:
            self.reading_frames = {
                frame_offset: list(seed_aminos)
                for frame_offset, seed_aminos in reading_frames.items()}
            seed_aminos = self.reading_frames[0]
            self.seed_nucs = [seed_nuc
                              for seed_amino in seed_aminos
                              for seed_nuc in seed_amino.nucleotides]
            aligned_consensus = ''.join(nuc.get_consensus('MAX') or '?'
                                        for nuc in self.seed_nucs)
            consensus = aligned_consensus.lstrip('?')
            self.consensus_offset = len(aligned_consensus) - len(consensus)
            self.consensus = consensus.rstrip('?')
        if not coordinate_name:
            return
        self.coordinate_name = coordinate_name
        try:
            coordinate_seq = self.projects.getGenotypeReference(coordinate_name)
        except KeyError:
            coordinate_seq = self.projects.getReference(coordinate_name)

        self.alignments, self.algorithm = align_consensus(coordinate_seq, self.consensus)

        if self.overall_alignments_writer is not None:
            for alignment in self.alignments:
                row = {"coordinate_name": self.coordinate_name,
                       "contig": self.contig_name,
                       "query_start": alignment.q_st,
                       "query_end": alignment.q_en,
                       "consensus_offset": self.consensus_offset,
                       "ref_start": alignment.r_st,
                       "ref_end": alignment.r_en,
                       "cigar_str": alignment.cigar_str}
                self.overall_alignments_writer.writerow(row)

    def find_amino_alignments(self,
                              start_pos: int,
                              end_pos: int,
                              repeat_pos: Optional[int],
                              skip_pos: Optional[int],
                              amino_ref: Optional[str]):
        translations = {
            reading_frame: translate(
                '-'*(reading_frame + self.consensus_offset) +
                self.consensus +
                '-'*((3-len(self.consensus) -
                      reading_frame -
                      self.consensus_offset) % 3))
            for reading_frame in range(3)}
        self.amino_alignments.clear()
        query_progress = 0
        for alignment in self.alignments:
            amino_sections = []
            query_end = alignment.q_st + self.consensus_offset
            ref_end = alignment.r_st
            for cigar_index, (size, action) in enumerate(alignment.cigar):
                ref_reading_frame = (ref_end - (start_pos-1)) % 3
                conseq_codon_start = query_end - ref_reading_frame
                reading_frame = -conseq_codon_start % 3
                amino_alignment = AminoAlignment(query_end,
                                                 query_end,
                                                 ref_end,
                                                 ref_end,
                                                 action,
                                                 reading_frame)
                if action == CigarActions.MATCH:
                    amino_alignment.query_end += size
                    amino_alignment.ref_end += size
                elif action == CigarActions.INSERT:
                    amino_alignment.query_end += size
                else:
                    assert action == CigarActions.DELETE
                    amino_alignment.ref_end += size
                query_end = amino_alignment.query_end
                ref_end = amino_alignment.ref_end
                if amino_alignment.has_overlap(start_pos, end_pos):
                    start_shift = max(0,
                                      start_pos-1 - amino_alignment.ref_start,
                                      query_progress -
                                      amino_alignment.query_start -
                                      self.consensus_offset)
                    end_shift = max(0, amino_alignment.ref_end - end_pos)
                    amino_alignment.ref_start += start_shift
                    amino_alignment.ref_end -= end_shift
                    if action != CigarActions.DELETE:
                        amino_alignment.query_start += start_shift
                        amino_alignment.query_end -= end_shift
                    alignment_size = amino_alignment.ref_end - amino_alignment.ref_start
                    if action == CigarActions.MATCH and alignment_size <= 0:
                        pass
                    else:
                        amino_sections.append(amino_alignment)
                query_progress = query_end
            if repeat_pos is not None:
                for i, amino_alignment in enumerate(amino_sections):
                    if amino_alignment.action != CigarActions.MATCH:
                        continue
                    ref_start = amino_alignment.ref_start
                    ref_end = amino_alignment.ref_end
                    # make sure that there are at least 3 nucleotides in each split alignment
                    if ref_start <= repeat_pos-3 and repeat_pos+3 <= ref_end:
                        offset = repeat_pos - ref_start
                        query_start = amino_alignment.query_start
                        amino_alignment2 = replace(
                            amino_alignment,
                            ref_start=repeat_pos,
                            ref_end=ref_end + 1,
                            query_start=query_start + offset-1,
                            reading_frame=(amino_alignment.reading_frame+1) % 3)
                        amino_alignment.ref_end = repeat_pos
                        amino_alignment.query_end = query_start + offset
                        amino_sections.insert(i+1, amino_alignment2)
                        break
            if skip_pos is not None:
                for i, amino_alignment in enumerate(amino_sections):
                    if amino_alignment.action != CigarActions.MATCH:
                        continue
                    ref_start = amino_alignment.ref_start
                    ref_end = amino_alignment.ref_end
                    # make sure that there are at least 3 nucleotides in each split alignment
                    if ref_start <= skip_pos-4 and skip_pos+3 <= ref_end:
                        offset = skip_pos - ref_start
                        query_start = amino_alignment.query_start
                        query_end = amino_alignment.query_end
                        amino_alignment2 = replace(
                            amino_alignment,
                            ref_start=skip_pos,
                            query_start=query_start + offset - 1,
                            query_end=query_end+1,
                            reading_frame=(amino_alignment.reading_frame+1) % 3)
                        amino_alignment.ref_end = skip_pos - 1
                        amino_alignment.query_end = query_start + offset
                        amino_sections.insert(i+1, amino_alignment2)
                        break
            for amino_alignment in amino_sections:
                if amino_alignment.action != CigarActions.MATCH:
                    continue
                amino_alignment.find_reading_frame(amino_ref,
                                                   start_pos,
                                                   translations)

            if self.unmerged_alignments_writer is not None:
                self.write_alignments_file(amino_sections, self.unmerged_alignments_writer)

            amino_sections = self.combine_alignments(amino_sections, amino_ref, start_pos, translations)

            if self.intermediate_alignments_writer is not None:
                self.write_alignments_file(amino_sections, self.intermediate_alignments_writer)

            # do a second pass to combine alignments
            amino_sections = self.combine_alignments(amino_sections, amino_ref, start_pos, translations)

            self.amino_alignments.extend(amino_sections)

        if self.alignments_writer is not None:
            self.write_alignments_file(self.amino_alignments, self.alignments_writer)

    def clear(self):
        self.coordinate_name = self.consensus = self.amino_consensus = ''
        self.algorithm = ''
        self.alignments.clear()
        self.seed_nucs.clear()
        self.inserts.clear()

    def write_alignments_file(self, amino_alignments, alignments_writer):
        for alignment in amino_alignments:
            row = {"action": CigarActions(alignment.action).name,
                   "contig": self.contig_name,
                   "query_start": alignment.query_start,
                   "query_end": alignment.query_end,
                   "ref_start": alignment.ref_start,
                   "ref_end": alignment.ref_end,
                   "aligned_query": alignment.aligned_query,
                   "aligned_ref": alignment.aligned_ref,
                   "reading_frame": alignment.reading_frame,
                   "ref_amino_start": alignment.ref_amino_start,
                   "coordinate_name": self.coordinate_name}
            alignments_writer.writerow(row)

    @staticmethod
    def combine_alignments(amino_sections, amino_ref, start_pos, translations):
        for i in range(len(amino_sections) - 2, 0, -1):
            amino_alignment = amino_sections[i]
            if amino_alignment.action == CigarActions.MATCH:
                continue
            size = max(amino_alignment.ref_end - amino_alignment.ref_start,
                       amino_alignment.query_end - amino_alignment.query_start)
            prev_alignment = amino_sections[i - 1]
            next_alignment = amino_sections[i + 1]
            can_align = (MINIMUM_AMINO_ALIGNMENT <= min(
                prev_alignment.amino_size,
                next_alignment.amino_size))
            has_frame_shift = (prev_alignment.reading_frame !=
                               next_alignment.reading_frame)
            has_big_gap = MAXIMUM_AMINO_GAP < size // 3
            if can_align and (has_frame_shift or has_big_gap):
                # Both neighbours are big enough, so we have a choice.
                # Either there's a frame shift or a big gap, so keep
                # dividing around this indel.
                continue
            # Merge the two sections on either side of this indel.
            amino_sections.pop(i + 1)
            amino_sections.pop(i)
            new_alignment = AminoAlignment(prev_alignment.query_start,
                                           next_alignment.query_end,
                                           prev_alignment.ref_start,
                                           next_alignment.ref_end,
                                           CigarActions.MATCH,
                                           prev_alignment.reading_frame)
            amino_sections[i - 1] = new_alignment
            new_alignment.find_reading_frame(amino_ref,
                                             start_pos,
                                             translations)
        return amino_sections

    def report_region(
            self,
            start_pos: int,
            end_pos: int,
            report_nucleotides: List[ReportNucleotide],
            report_aminos: Optional[List[ReportAmino]] = None,
            repeat_position: Optional[int] = None,
            skip_position: Optional[int] = None,
            amino_ref: Optional[str] = None):
        """ Add read counts to report counts for a section of the reference.

        :param start_pos: 1-based position of first nucleotide to report in
            reference coordinates.
        :param end_pos: 1-based position of last nucleotide to report in
            reference coordinates.
        :param report_nucleotides: list of nucleotide objects to collect counts.
            Will be extended to required size, if needed.
        :param report_aminos: list of amino objects to collect counts.
            Will be extended to required size, if needed.
        :param repeat_position: if not None, repeat the nucleotide at this
            1-based position in reference coordinates, but only in report_aminos.
        :param skip_position: if not None, skip the nucleotide at this
            1-based position in reference coordinates, but only in report_aminos.
        :param amino_ref: amino acid sequence to align with, or None if only
            nucleotides are reported.
        """
        self.inserts = set()
        if not self.alignments:
            # Take a wild guess at the consensus to report the failure.
            self.amino_consensus = ''.join(
                seed_amino.get_consensus() or '?'
                for seed_amino in self.reading_frames[0])
            return

        self.amino_consensus = ''
        nuc_count = end_pos - start_pos + 1
        report_nucleotides.extend(ReportNucleotide(i+1)
                                  for i in range(nuc_count))
        if report_aminos is None:
            self.build_nucleotide_report(start_pos,
                                         end_pos,
                                         report_nucleotides)
        elif amino_ref is not None:
            report_aminos.extend(ReportAmino(SeedAmino(None), i + 1)
                                 for i in range(len(amino_ref)))
            self.build_amino_report(start_pos,
                                    end_pos,
                                    report_nucleotides,
                                    report_aminos,
                                    repeat_position,
                                    skip_position,
                                    amino_ref)

    def get_seed_nuc(self, consensus_nuc_index: int):
        seed_amino = self.reading_frames[0][consensus_nuc_index // 3]
        return seed_amino.nucleotides[consensus_nuc_index % 3]

    def get_deletion_coverage(self, consensus_nuc_index):
        try:
            prev_nuc = self.get_seed_nuc(consensus_nuc_index)
            next_nuc = self.get_seed_nuc(consensus_nuc_index + 1)
            coverage = min(prev_nuc.get_coverage(), next_nuc.get_coverage())
        except IndexError:
            coverage = 0
            logger.warning(f"Could not get deletion coverage for consensus index {consensus_nuc_index}")
        return coverage

    def build_amino_report(self,
                           start_pos: int,
                           end_pos: int,
                           report_nucleotides: List[ReportNucleotide],
                           report_aminos: Optional[List[ReportAmino]] = None,
                           repeat_position: Optional[int] = None,
                           skip_position: Optional[int] = None,
                           amino_ref: Optional[str] = None):
        """ Add read counts to report counts for a section of the reference.

        Used for regions that translate to amino acids.

        :param start_pos: 1-based position of first nucleotide to report in
            reference coordinates.
        :param end_pos: 1-based position of last nucleotide to report in
            reference coordinates.
        :param report_nucleotides: list of nucleotide objects to collect counts.
            Will be extended to required size, if needed.
        :param report_aminos: list of amino objects to collect counts.
            Will be extended to required size, if needed.
        :param repeat_position: if not None, repeat the nucleotide at this
            1-based position in reference coordinates, but only in report_aminos.
        :param skip_position: if not None, skip the nucleotide at this
            1-based position in reference coordinates, but only in report_aminos.
        :param amino_ref: amino acid sequence to align with, or None if only
            nucleotides are reported.
        """
        self.find_amino_alignments(start_pos,
                                   end_pos,
                                   repeat_position,
                                   skip_position,
                                   amino_ref)
        has_skipped_nucleotide = False
        if skip_position is not None:
            reading_frame1 = reading_frame2 = None
            for alignment in self.amino_alignments:
                if alignment.ref_end == skip_position-1:
                    reading_frame1 = alignment.reading_frame
                elif alignment.ref_start == skip_position:
                    reading_frame2 = alignment.reading_frame
            if reading_frame1 is not None and reading_frame2 is not None:
                if reading_frame1 != reading_frame2:
                    has_skipped_nucleotide = True
        for amino_alignment in self.amino_alignments:
            if amino_alignment.action == CigarActions.DELETE:
                self.count_deletion(amino_alignment,
                                    amino_ref,
                                    report_aminos,
                                    report_nucleotides,
                                    start_pos,
                                    repeat_position,
                                    skip_position)
            elif amino_alignment.action == CigarActions.INSERT:
                self.count_insertion(amino_alignment, start_pos, end_pos)
            else:
                assert amino_alignment.action == CigarActions.MATCH
                self.count_match(amino_alignment,
                                 amino_ref,
                                 report_aminos,
                                 report_nucleotides,
                                 start_pos,
                                 repeat_position,
                                 skip_position,
                                 has_skipped_nucleotide)

    @staticmethod
    def update_report_amino(coord_index: int,
                            report_aminos: List[ReportAmino],
                            report_nucleotides: List[ReportNucleotide],
                            seed_amino: SeedAmino,
                            start_pos: int,
                            repeat_position: Optional[int] = None,
                            skip_position: Optional[int] = None,
                            skipped_nuc: Optional[SeedAmino] =None):
        report_amino = report_aminos[coord_index]
        report_amino.seed_amino.add(seed_amino)
        ref_nuc_pos = coord_index * 3 + start_pos
        for codon_nuc_index, seed_nuc in enumerate(
                seed_amino.nucleotides):
            if len(report_nucleotides) <= coord_index * 3 + codon_nuc_index:
                if repeat_position is None:
                    continue
                elif repeat_position is not None and repeat_position >= ref_nuc_pos:
                    continue

            report_nuc_index = coord_index * 3 + codon_nuc_index
            if repeat_position is not None and start_pos < repeat_position:
                if repeat_position < ref_nuc_pos:
                    report_nuc_index -= 1
                if repeat_position == ref_nuc_pos - 1 and codon_nuc_index == 0:
                    continue
            if skip_position is not None:
                if skip_position == start_pos + report_nuc_index and skipped_nuc is not None:
                    report_nuc = report_nucleotides[report_nuc_index]
                    report_nuc.seed_nucleotide.add(skipped_nuc)
                if skip_position-1 < start_pos + report_nuc_index:
                    report_nuc_index += 1
            report_nuc = report_nucleotides[report_nuc_index]
            report_nuc.seed_nucleotide.add(seed_nuc)

    def count_deletion(self,
                       amino_alignment,
                       amino_ref,
                       report_aminos,
                       report_nucleotides,
                       start_pos,
                       repeat_position,
                       skip_position):
        coverage = self.get_deletion_coverage(
            amino_alignment.query_start)
        seed_amino = SeedAmino(None)
        seed_amino.count_aminos('---', coverage)

        # Start index within region
        del_start_nuc_index = amino_alignment.ref_start - start_pos + 1
        del_end_nuc_index = amino_alignment.ref_end - start_pos
        # Calculate amino index. amino_alignment is already trimmed to
        # be between start_pos and end_pos, but HIV's vpr has an extra
        # base that can roll over to an extra amino acid. Trim to length
        # of amino reference.
        del_start_amino_index = min(del_start_nuc_index // 3,
                                    len(amino_ref) - 1)
        del_end_amino_index = min(del_end_nuc_index // 3,
                                  len(amino_ref) - 1)
        for coord_index in count(del_start_amino_index):
            if coord_index == del_start_amino_index:
                start_codon_nuc = del_start_nuc_index % 3
                trimmed_seed_amino = SeedAmino(None)
                trimmed_seed_amino.add(seed_amino, start_codon_nuc)
                trimmed_seed_amino.consensus_nuc_index = (
                        amino_alignment.query_start -
                        start_codon_nuc)
            elif coord_index == del_end_amino_index:
                end_codon_nuc = del_end_nuc_index % 3
                trimmed_seed_amino = SeedAmino(None)
                trimmed_seed_amino.add(seed_amino, 0, end_codon_nuc)
            else:
                trimmed_seed_amino = seed_amino
            self.update_report_amino(coord_index,
                                     report_aminos,
                                     report_nucleotides,
                                     trimmed_seed_amino,
                                     start_pos,
                                     repeat_position=repeat_position,
                                     skip_position=skip_position)
            if coord_index == del_end_amino_index:
                break

    def count_insertion(self, amino_alignment, start_pos, end_pos):
        ref_nuc_index = amino_alignment.ref_start
        consensus_nuc_index = amino_alignment.query_start
        section_size = (amino_alignment.query_end - amino_alignment.query_start) // 3
        for _ in range(section_size):
            if start_pos - 1 <= ref_nuc_index < end_pos:
                self.inserts.add(consensus_nuc_index)
            consensus_nuc_index += 3

    def count_match(self,
                    amino_alignment,
                    amino_ref,
                    report_aminos,
                    report_nucleotides,
                    start_pos,
                    repeat_position,
                    skip_position,
                    has_skipped_nucleotide):
        region_seed_aminos = self.reading_frames[amino_alignment.reading_frame]
        coord2conseq = amino_alignment.map_amino_sequences()
        if not coord2conseq:
            # coord2conseq is empty (alignment was too small / at region boundary) - there is nothing to do
            return

        coordinate_inserts = {seed_amino.consensus_nuc_index
                              for seed_amino in region_seed_aminos}
        prev_conseq_index = None
        prev_consensus_nuc_index = None
        max_conseq_index = max(coord2conseq.values())
        for coord_index in range(len(amino_ref)):
            conseq_index = coord2conseq.get(coord_index)
            if conseq_index is None:
                seed_amino = SeedAmino(None)
                if (prev_conseq_index is not None and
                        prev_conseq_index < max_conseq_index):
                    prev_seed_amino = region_seed_aminos[prev_conseq_index]
                    prev_count = sum(prev_seed_amino.counts.values())
                    prev_count += prev_seed_amino.deletions
                    next_seed_amino = region_seed_aminos[prev_conseq_index + 1]
                    next_count = sum(next_seed_amino.counts.values())
                    next_count += next_seed_amino.deletions
                    min_count = min(prev_count, next_count)
                    seed_amino.deletions = min_count
                    for nuc in seed_amino.nucleotides:
                        nuc.count_nucleotides('-', min_count)
            else:
                seed_amino = region_seed_aminos[conseq_index]
                consensus_nuc_index = seed_amino.consensus_nuc_index
                if consensus_nuc_index < amino_alignment.query_start:
                    start_codon_nuc = max(0,
                                          amino_alignment.query_start -
                                          consensus_nuc_index)
                    end_codon_nuc = 2
                    seed_amino2 = SeedAmino(None)
                    seed_amino2.add(seed_amino,
                                    start_codon_nuc,
                                    end_codon_nuc)
                    seed_amino = seed_amino2
                if amino_alignment.query_end < consensus_nuc_index + 3:
                    start_codon_nuc = 0
                    end_codon_nuc = (amino_alignment.query_end -
                                     consensus_nuc_index - 1)
                    seed_amino2 = SeedAmino(None)
                    seed_amino2.add(seed_amino,
                                    start_codon_nuc,
                                    end_codon_nuc)
                    seed_amino = seed_amino2
                if prev_conseq_index is None:
                    coordinate_inserts = {i
                                          for i in coordinate_inserts
                                          if i >= consensus_nuc_index}
                prev_conseq_index = conseq_index

            if seed_amino.consensus_nuc_index is not None:
                coordinate_inserts.remove(seed_amino.consensus_nuc_index)
                prev_consensus_nuc_index = seed_amino.consensus_nuc_index
            if skip_position is not None and coord_index == (skip_position - start_pos) // 3 and \
                    amino_alignment.ref_start < skip_position - 1 <= amino_alignment.ref_end:
                conseq_pos = coord2conseq.get((skip_position - 1 - start_pos) // 3)
                if conseq_pos is not None:
                    if has_skipped_nucleotide:
                        skipped_nuc = \
                            self.reading_frames[amino_alignment.reading_frame][conseq_pos + 1].nucleotides[0]
                    else:
                        skipped_nuc = SeedNucleotide()
                        coverage = self.get_deletion_coverage(conseq_pos)
                        skipped_nuc.count_nucleotides('-', coverage)
                else:
                    skipped_nuc = SeedNucleotide()
            else:
                skipped_nuc = None
            self.update_report_amino(coord_index,
                                     report_aminos,
                                     report_nucleotides,
                                     seed_amino,
                                     start_pos,
                                     repeat_position=repeat_position,
                                     skip_position=skip_position,
                                     skipped_nuc=skipped_nuc, )
        if prev_consensus_nuc_index is None:
            coordinate_inserts.clear()
        else:
            coordinate_inserts = {i
                                  for i in coordinate_inserts
                                  if i <= prev_consensus_nuc_index}
        self.inserts.update(coordinate_inserts)

    def build_nucleotide_report(self,
                                start_pos: int,
                                end_pos: int,
                                report_nucleotides: List[ReportNucleotide]):
        """ Add read counts to report counts for a section of the reference.

        Used for regions that don't translate to amino acids.

        :param start_pos: 1-based position of first nucleotide to report in
            reference coordinates.
        :param end_pos: 1-based position of last nucleotide to report in
            reference coordinates.
        :param report_nucleotides: list of nucleotide objects to collect counts.
            Will be extended to required size, if needed.
        """
        for alignment in self.alignments:
            ref_nuc_index = alignment.r_st
            consensus_nuc_index = alignment.q_st + self.consensus_offset
            source_amino_index = 0
            seed_aminos = self.reading_frames[0]
            for section_size, section_action in alignment.cigar:
                if section_action == CigarActions.INSERT:
                    for _ in range(section_size):
                        if start_pos - 1 <= ref_nuc_index < end_pos:
                            self.inserts.add(consensus_nuc_index)
                        consensus_nuc_index += 1
                    continue
                if section_action == CigarActions.DELETE:
                    coverage = self.get_deletion_coverage(consensus_nuc_index)
                    seed_nuc = SeedNucleotide()
                    seed_nuc.count_nucleotides('-', coverage)
                    for _ in range(section_size):
                        if start_pos - 1 <= ref_nuc_index < end_pos:
                            target_nuc_index = ref_nuc_index - start_pos + 1
                            report_nuc = report_nucleotides[target_nuc_index]
                            report_nuc.seed_nucleotide.add(seed_nuc)
                        ref_nuc_index += 1
                    continue
                assert section_action == CigarActions.MATCH, section_action
                for _ in range(section_size):
                    if start_pos - 1 <= ref_nuc_index < end_pos:
                        target_nuc_index = ref_nuc_index - start_pos + 1
                        while True:
                            seed_amino = seed_aminos[source_amino_index]
                            codon_nuc_index = (consensus_nuc_index -
                                               seed_amino.consensus_nuc_index)
                            if codon_nuc_index < 3:
                                break
                            source_amino_index += 1
                        seed_nuc = seed_amino.nucleotides[codon_nuc_index]
                        report_nuc = report_nucleotides[target_nuc_index]
                        report_nuc.seed_nucleotide.add(seed_nuc)
                    ref_nuc_index += 1
                    consensus_nuc_index += 1

    def seed_concordance(self, seed_name, projects, seed_coordinates, excluded_regions, included_regions=None):
        if self.seed_concordance_writer is None:
            return
        seed_ref = self.projects.getReference(seed_name)
        seed_aligner = mappy.Aligner(seq=seed_ref, bw=500, bw_long=500, preset='map-ont')
        seed_alignments: List[mappy.Alignment] = list(seed_aligner.map(self.consensus))

        regions = projects.getCoordinateReferences(seed_name)
        for region in regions:
            if not self.projects.isAmino(region) or region in excluded_regions:
                continue
            if included_regions and region not in included_regions:
                continue
            coordinates = seed_coordinates[seed_name][region]
            try:
                start_pos = coordinates['start']
                end_pos = coordinates['end']
            except KeyError:
                continue
            self.region_seed_concordance(region, seed_name, seed_alignments, seed_ref, start_pos, end_pos)

    def coord_concordance(self, half_window_size: int = 10) -> List[float]:
        coord_alignments = self.alignments
        try:
            coord_ref = self.projects.getGenotypeReference(self.coordinate_name)
        except KeyError:
            coord_ref = self.projects.getReference(self.coordinate_name)
        query_matches = [0] * len(self.consensus)
        concordance_list: List[float] = [0] * len(self.consensus)

        for alignment in coord_alignments:
            ref_progress = alignment.r_st
            query_progress = alignment.q_st
            for cigar_index, (size, action) in enumerate(alignment.cigar):
                if action == CigarActions.INSERT:
                    query_progress += size
                elif action == CigarActions.DELETE:
                    ref_progress += size
                else:
                    assert action == CigarActions.MATCH
                    for pos in range(0, size):
                        ref_pos = ref_progress + pos
                        query_pos = query_progress + pos
                        if self.consensus[query_pos] == coord_ref[ref_pos]:
                            query_matches[query_pos] = 1
                        else:
                            query_matches[query_pos] = 0
                    ref_progress += size
                    query_progress += size

        for pos in range(0, len(self.consensus)):
            start = max(0, pos - half_window_size)
            end = min(len(self.consensus), pos + half_window_size)
            concordance = sum(query_matches[start:end])
            normalized_concordance = concordance / (end-start)
            concordance_list[pos] = normalized_concordance

        return concordance_list

    def region_seed_concordance(self, region, seed_name, seed_alignments, seed_ref, start_pos, end_pos):
        length_aligned = end_pos - start_pos
        nuc_agreements = [0] * length_aligned
        nuc_covered = [0] * length_aligned
        for alignment in seed_alignments:
            ref_progress = alignment.r_st
            query_progress = alignment.q_st + self.consensus_offset
            for cigar_index, (size, action) in enumerate(alignment.cigar):
                if action == CigarActions.INSERT:
                    query_progress += size
                elif action == CigarActions.DELETE:
                    ref_progress += size
                else:
                    assert action == CigarActions.MATCH
                    amino_alignment = AminoAlignment(query_progress,
                                                     query_progress + size,
                                                     ref_progress,
                                                     ref_progress + size,
                                                     action,
                                                     0)
                    query_progress += size
                    ref_progress += size
                    if amino_alignment.has_overlap(start_pos, end_pos):
                        start_shift = max(0, start_pos - amino_alignment.ref_start)
                        end_shift = max(0, amino_alignment.ref_end - end_pos)
                        amino_alignment.ref_start += start_shift
                        amino_alignment.ref_end -= end_shift
                        amino_alignment.query_start += start_shift
                        amino_alignment.query_end -= end_shift
                        match_size = amino_alignment.ref_end - amino_alignment.ref_start
                        for pos in range(0, match_size):
                            query_nuc = self.consensus[amino_alignment.query_start - self.consensus_offset + pos]
                            seed_nuc = seed_ref[amino_alignment.ref_start + pos]
                            offset_pos = pos + amino_alignment.ref_start - start_pos
                            if query_nuc == seed_nuc:
                                nuc_agreements[offset_pos] = 1
                            nuc_covered[offset_pos] = 1
        covered = 100 * sum(nuc_covered) / length_aligned
        try:
            concordance_covered = 100 * sum(nuc_agreements) / sum(nuc_covered)
        except ZeroDivisionError:
            concordance_covered = 0.0
        region_row = dict(seed_name=seed_name,
                          region=region,
                          contig=self.contig_name,
                          pct_concordance=concordance_covered,
                          pct_covered=covered)
        self.seed_concordance_writer.writerow(region_row)


@dataclass
class AminoAlignment:
    """ Part of a full alignment, with aligned aminos.
     Positions are 0-based with non-inclusive ends."""
    query_start: int
    query_end: int
    ref_start: int
    ref_end: int
    action: CigarActions
    reading_frame: int
    query: Optional[str] = None  # Amino sequence
    ref: Optional[str] = None  # Amino sequence
    aligned_query: Optional[str] = None
    aligned_ref: Optional[str] = None
    ref_amino_start: Optional[int] = None

    def has_overlap(self, start_pos: int, end_pos: int) -> bool:
        before_end = self.ref_start < end_pos
        after_start = start_pos <= self.ref_end
        return before_end == after_start

    def find_reading_frame(self, amino_ref, start_pos, translations):
        ref_amino_start = (self.ref_start - start_pos + 1) // 3
        ref_amino_end = (self.ref_end - start_pos + 3) // 3
        self.ref_amino_start = ref_amino_start

        ref_section = amino_ref[ref_amino_start:ref_amino_end]
        # Exclamation marks avoid extra gaps near the ends.
        assert '!' not in ref_section
        ref_to_align = '!' + ref_section + '!'
        max_score = (ref_amino_end - ref_amino_start) // 2

        for reading_frame in range(3):
            query_amino_start = (self.query_start + reading_frame) // 3
            query_amino_end = (self.query_end + reading_frame + 2) // 3
            translation = translations[reading_frame]
            query_section = translation[query_amino_start:query_amino_end]
            aref, aquery, score = align_aminos(ref_to_align, query_section)
            aref = aref.replace('!', '-')
            if max_score < score:
                max_score = score
                self.reading_frame = reading_frame
                # Record alignments after trimming exclamation marks.
                self.aligned_query = aquery
                self.aligned_ref = aref

            if self.reading_frame == reading_frame:
                # Either had a good match, or we're falling back to best guess.
                self.query = query_section
                self.ref = ref_section
        if self.aligned_ref is None:
            # Falling back to best guess.
            self.aligned_ref = self.ref
            self.aligned_query = self.query

    @property
    def size(self):
        return max(self.ref_end - self.ref_start,
                   self.query_end - self.query_start)

    @property
    def amino_size(self):
        return (self.size + 2) // 3

    def map_amino_sequences(self) -> Dict[int, int]:
        """ Map reference amino indexes to query amino indexes. """

        assert self.aligned_ref is not None, "For this operation, aligned_ref must not be None"
        assert self.aligned_query is not None, "For this operation, aligned_query must not be None"
        assert self.query is not None, "For this operation, query must not be None"
        assert self.ref is not None, "For this operation, ref must not be None"
        assert self.ref_amino_start is not None, "For this operation, ref_amino_start must not be None"

        seq_map = {}
        query_offset = (self.query_start + self.reading_frame) // 3
        ref_index = query_index = 0
        for from_aa, to_aa in zip(self.aligned_ref, self.aligned_query):
            if (query_index < len(self.query) and
                    ref_index < len(self.ref) and
                    to_aa == self.query[query_index]):
                seq_map[ref_index+self.ref_amino_start] = query_index+query_offset
                query_index += 1
            if (ref_index < len(self.ref) and
                    from_aa == self.ref[ref_index]):
                ref_index += 1
        return seq_map
