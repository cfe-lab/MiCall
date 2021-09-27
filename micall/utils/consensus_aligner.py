import typing
from dataclasses import dataclass, replace
from enum import IntEnum
from itertools import count
from operator import attrgetter

from gotoh import align_it, align_it_aa
from mappy import Alignment, Aligner

from micall.core.project_config import ProjectConfig
from micall.utils.report_amino import SeedAmino, ReportAmino, ReportNucleotide, SeedNucleotide

# A section between reading frame shifts must be at least this long to be split.
from micall.utils.translation import translate

MINIMUM_READING_FRAME_SHIFT = 30
# Minimum aminos in a section between reading frame shifts to be split.
MINIMUM_AMINO_ALIGNMENT = 10
# Most codons in an insertion or deletion that is still aligned in amino acids.
MAXIMUM_AMINO_GAP = 10

CigarActions = IntEnum(
    'CigarActions',
    'MATCH INSERT DELETE SKIPPED SOFT_CLIPPED HARD_CLIPPED',
    start=0)


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


class AlignmentWrapper(Alignment):
    init_fields = (
        'ctg ctg_len r_st r_en strand q_st q_en mapq cigar is_primary mlen '
        'blen NM trans_strand read_num cs MD').split()

    @classmethod
    def wrap(cls, source: Alignment, **overrides):
        """ Wrap an Alignment object to make it easier to compare and display.

        Mostly used when testing.
        """
        args = [getattr(source, field_name)
                for field_name in cls.init_fields]
        for name, value in overrides.items():
            i = cls.init_fields.index(name)
            args[i] = value
        return cls(*args)

    # noinspection PyPep8Naming
    def __new__(cls,
                ctg='',
                ctg_len=0,
                r_st=0,
                r_en=0,
                strand=1,
                q_st=0,
                q_en=0,
                mapq=0,
                cigar: typing.Iterable[typing.List[int]] = tuple(),
                is_primary=True,
                mlen=0,
                blen=0,
                NM=0,
                trans_strand=0,
                read_num=1,
                cs='',
                MD=''):
        """ Create an instance.

        :param ctg: name of the reference sequence the query is mapped to
        :param ctg_len: total length of the reference sequence
        :param r_st and r_en: start and end positions on the reference
        :param strand: +1 if on the forward strand; -1 if on the reverse strand
        :param q_st and q_en: start and end positions on the query
        :param mapq: mapping quality
        :param cigar: CIGAR returned as an array of shape (n_cigar,2). The two
            numbers give the length and the operator of each CIGAR operation.
        :param is_primary: if the alignment is primary (typically the best and
            the first to generate)
        :param mlen: length of the matching bases in the alignment, excluding
            ambiguous base matches.
        :param blen: length of the alignment, including both alignment matches
            and gaps but excluding ambiguous bases.
        :param NM: number of mismatches, gaps and ambiguous positions in the
            alignment
        :param trans_strand: transcript strand. +1 if on the forward strand; -1
            if on the reverse strand; 0 if unknown
        :param read_num: read number that the alignment corresponds to; 1 for
            the first read and 2 for the second read
        :param cs: the cs tag.
        :param MD: the MD tag as in the SAM format. It is an empty string unless
            the MD argument is applied when calling mappy.Aligner.map().
        """
        cigar = list(cigar)
        if not mlen:
            mlen = min(q_en-q_st, r_en-r_st)
        if not blen:
            blen = max(q_en-q_st, r_en-r_st)
        if not cigar:
            cigar = [[max(q_en-q_st, r_en-r_st), CigarActions.MATCH]]
        return super().__new__(cls,
                               ctg,
                               ctg_len,
                               r_st,
                               r_en,
                               strand,
                               q_st,
                               q_en,
                               mapq,
                               cigar,
                               is_primary,
                               mlen,
                               blen,
                               NM,
                               trans_strand,
                               read_num-1,
                               cs,
                               MD)

    def __eq__(self, other: Alignment):
        for field_name in self.init_fields:
            self_value = getattr(self, field_name)
            other_value = getattr(other, field_name)
            if self_value != other_value:
                return False
        return True

    def __repr__(self):
        return (f'AlignmentWrapper({self.ctg!r}, {self.ctg_len}, '
                f'{self.r_st}, {self.r_en}, {self.strand}, '
                f'{self.q_st}, {self.q_en})')


def clip_seed_aminos(seed_aminos: typing.List[SeedAmino],
                     region_start_consensus_amino_index: int,
                     region_end_consensus_amino_index: int,
                     start_codon_nuc_index: int,
                     end_codon_nuc_index: int):
    """ Extract a section of seed aminos for counting.

    Handles boundary codons by clipping as needed.
    :param seed_aminos: source seed aminos to copy
    :param region_start_consensus_amino_index: first seed amino to copy
    :param region_end_consensus_amino_index: one past the last one to copy
    :param start_codon_nuc_index: first codon index within the first seed amino
        that should be copied
    :param end_codon_nuc_index: last codon index to copy within the last codon
    """
    assert 0 <= region_start_consensus_amino_index
    assert 0 <= region_end_consensus_amino_index
    assert 0 <= start_codon_nuc_index
    assert 0 <= end_codon_nuc_index
    result = seed_aminos[region_start_consensus_amino_index:
                         region_end_consensus_amino_index]
    if start_codon_nuc_index != 0:
        old_start_amino = result[0]
        start_amino = SeedAmino(old_start_amino.consensus_nuc_index)
        start_amino.add(old_start_amino, start_nuc=start_codon_nuc_index)
        result[0] = start_amino
    if end_codon_nuc_index != 2:
        old_end_amino = result[-1]
        end_amino = SeedAmino(old_end_amino.consensus_nuc_index)
        end_amino.add(old_end_amino, end_nuc=end_codon_nuc_index)
        result[-1] = end_amino
    return result


class ConsensusAligner:
    def __init__(self, projects: ProjectConfig):
        self.projects = projects
        self.coordinate_name = self.consensus = self.amino_consensus = ''
        self.algorithm = ''
        self.consensus_offset = 0
        self.alignments: typing.List[Alignment] = []
        self.reading_frames: typing.List[typing.List[SeedAmino]] = []
        self.seed_nucs: typing.List[SeedNucleotide] = []
        self.amino_alignments: typing.List[AminoAlignment] = []

        # consensus nucleotide positions that were inserts
        self.inserts: typing.Set[int] = set()

    def start_contig(self,
                     coordinate_name: str = None,
                     consensus: str = None,
                     reading_frames: typing.Dict[
                         int,
                         typing.List[SeedAmino]] = None):
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
        coordinate_seq = self.projects.getReference(coordinate_name)
        aligner = Aligner(seq=coordinate_seq, preset='map-ont')
        self.alignments = list(aligner.map(self.consensus))
        if self.alignments or 10_000 < len(self.consensus):
            self.algorithm = 'minimap2'
        else:
            self.algorithm = 'gotoh'
            self.align_gotoh(coordinate_seq, self.consensus)
        self.alignments = [alignment
                           for alignment in self.alignments
                           if alignment.is_primary]
        self.alignments.sort(key=attrgetter('q_st'))

    def align_gotoh(self, coordinate_seq, consensus):
        gap_open_penalty = 15
        gap_extend_penalty = 3
        use_terminal_gap_penalty = 1
        aligned_coordinate, aligned_consensus, score = align_it(
            coordinate_seq,
            consensus,
            gap_open_penalty,
            gap_extend_penalty,
            use_terminal_gap_penalty)
        if min(len(coordinate_seq), len(consensus)) < score:
            ref_start = len(aligned_consensus) - len(aligned_consensus.lstrip('-'))
            aligned_consensus: str = aligned_consensus[ref_start:]
            aligned_coordinate: str = aligned_coordinate[ref_start:]
            aligned_consensus = aligned_consensus.rstrip('-')
            ref_index = ref_start
            consensus_index = 0
            cigar = []
            for ref_nuc, nuc in zip(aligned_coordinate, aligned_consensus):
                expected_nuc = consensus[consensus_index]
                ref_index += 1
                consensus_index += 1
                expected_action = CigarActions.MATCH
                if nuc == '-' and nuc != expected_nuc:
                    expected_action = CigarActions.DELETE
                    consensus_index -= 1
                if ref_nuc == '-':
                    expected_action = CigarActions.INSERT
                    ref_index -= 1
                if cigar and cigar[-1][1] == expected_action:
                    cigar[-1][0] += 1
                else:
                    cigar.append([1, expected_action])
            self.alignments.append(AlignmentWrapper(
                'N/A',
                len(coordinate_seq),
                ref_start,
                ref_index,
                q_st=0,
                q_en=consensus_index,
                cigar=cigar))

    def find_amino_alignments(self,
                              start_pos: int,
                              end_pos: int,
                              repeat_pos: typing.Optional[int],
                              skip_pos: typing.Optional[int],
                              amino_ref: str):
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
                    amino_alignment.ref_start += start_shift
                    amino_alignment.query_start += start_shift
                    end_shift = max(0, amino_alignment.ref_end - end_pos)
                    amino_alignment.ref_end -= end_shift
                    amino_alignment.query_end -= end_shift
                    amino_sections.append(amino_alignment)
                query_progress = query_end
            if repeat_pos is not None or skip_pos is not None:
                for i, amino_alignment in enumerate(amino_sections):
                    if amino_alignment.action != CigarActions.MATCH:
                        continue
                    ref_start = amino_alignment.ref_start
                    ref_end = amino_alignment.ref_end
                    notable_pos = repeat_pos if repeat_pos is not None else skip_pos
                    if ref_start < notable_pos-1 < ref_end:
                        offset = notable_pos - ref_start
                        query_start = amino_alignment.query_start
                        amino_alignment2 = replace(
                            amino_alignment,
                            ref_start=notable_pos,
                            ref_end=ref_end+1,
                            query_start=query_start + offset-1,
                            reading_frame=(amino_alignment.reading_frame+1) % 3)
                        if repeat_pos is not None:
                            amino_alignment2.ref_end = ref_end+1
                        amino_alignment.ref_end = notable_pos
                        amino_alignment.query_end = query_start + offset
                        amino_sections.insert(i+1, amino_alignment2)
                        break
            for amino_alignment in amino_sections:
                if amino_alignment.action != CigarActions.MATCH:
                    continue
                amino_alignment.find_reading_frame(amino_ref,
                                                   start_pos,
                                                   translations)
            for i in range(len(amino_sections)-2, 0, -1):
                amino_alignment = amino_sections[i]
                if amino_alignment.action == CigarActions.MATCH:
                    continue
                size = max(amino_alignment.ref_end-amino_alignment.ref_start,
                           amino_alignment.query_end-amino_alignment.query_start)
                prev_alignment = amino_sections[i-1]
                next_alignment = amino_sections[i+1]
                can_align = (MINIMUM_AMINO_ALIGNMENT <= min(
                    prev_alignment.amino_size,
                    next_alignment.amino_size))
                if can_align and (size % 3 != 0 or MAXIMUM_AMINO_GAP < size//3):
                    # Keep dividing around this indel.
                    continue
                # Merge the two sections on either side of this indel.
                amino_sections.pop(i+1)
                amino_sections.pop(i)
                new_alignment = AminoAlignment(prev_alignment.query_start,
                                               next_alignment.query_end,
                                               prev_alignment.ref_start,
                                               next_alignment.ref_end,
                                               CigarActions.MATCH,
                                               prev_alignment.reading_frame)
                amino_sections[i-1] = new_alignment
                new_alignment.find_reading_frame(amino_ref,
                                                 start_pos,
                                                 translations)
            self.amino_alignments.extend(amino_sections)

    def clear(self):
        self.coordinate_name = self.consensus = self.amino_consensus = ''
        self.algorithm = ''
        self.alignments.clear()
        self.seed_nucs.clear()
        self.inserts.clear()

    def report_region(
            self,
            start_pos: int,
            end_pos: int,
            report_nucleotides: typing.List[ReportNucleotide],
            report_aminos: typing.List[ReportAmino] = None,
            repeat_position: int = None,
            skip_position: int = None,
            amino_ref: str = None):
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
        else:
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
        prev_nuc = self.get_seed_nuc(consensus_nuc_index)
        next_nuc = self.get_seed_nuc(consensus_nuc_index + 1)
        coverage = min(prev_nuc.get_coverage(), next_nuc.get_coverage())
        return coverage

    def build_amino_report(self,
                           start_pos: int,
                           end_pos: int,
                           report_nucleotides: typing.List[ReportNucleotide],
                           report_aminos: typing.List[ReportAmino] = None,
                           repeat_position: int = None,
                           skip_position: int = None,
                           amino_ref: str = None):
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
        if skip_position is not None:
            has_skipped_nucleotide = False
            reading_frame1 = reading_frame2 = None
            for alignment in self.amino_alignments:
                if alignment.ref_end == skip_position:
                    reading_frame1 = alignment.reading_frame
                elif alignment.ref_start == skip_position:
                    reading_frame2 = alignment.reading_frame
            if reading_frame1 is not None and reading_frame2 is not None:
                if reading_frame1 != reading_frame2:
                    has_skipped_nucleotide = True
        for amino_alignment in self.amino_alignments:
            if amino_alignment.action == CigarActions.DELETE:
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
                                            len(amino_ref)-1)
                del_end_amino_index = min(del_end_nuc_index // 3,
                                          len(amino_ref)-1)
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
                                             repeat_position)
                    if coord_index == del_end_amino_index:
                        break

            if amino_alignment.action != CigarActions.MATCH:
                continue
            region_seed_aminos = self.reading_frames[amino_alignment.reading_frame]
            coord2conseq = amino_alignment.map_amino_sequences()

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
                        next_seed_amino = region_seed_aminos[prev_conseq_index+1]
                        next_count = sum(next_seed_amino.counts.values())
                        next_count += next_seed_amino.deletions
                        min_count = min(prev_count, next_count)
                        seed_amino.deletions = min_count
                        for nuc in seed_amino.nucleotides:
                            nuc.count_nucleotides('-', min_count)
                else:
                    seed_amino = region_seed_aminos[conseq_index]
                    if (seed_amino.consensus_nuc_index <
                            amino_alignment.query_start):
                        start_codon_nuc = max(0,
                                              amino_alignment.query_start -
                                              seed_amino.consensus_nuc_index)
                        end_codon_nuc = 2
                        seed_amino2 = SeedAmino(None)
                        seed_amino2.add(seed_amino,
                                        start_codon_nuc,
                                        end_codon_nuc)
                        seed_amino = seed_amino2
                    if (amino_alignment.query_end <
                            seed_amino.consensus_nuc_index+3):
                        start_codon_nuc = 0
                        end_codon_nuc = (amino_alignment.query_end -
                                         seed_amino.consensus_nuc_index - 1)
                        seed_amino2 = SeedAmino(None)
                        seed_amino2.add(seed_amino,
                                        start_codon_nuc,
                                        end_codon_nuc)
                        seed_amino = seed_amino2
                    if prev_conseq_index is None:
                        coordinate_inserts = {i
                                              for i in coordinate_inserts
                                              if i >= seed_amino.consensus_nuc_index}
                    prev_conseq_index = conseq_index

                if seed_amino.consensus_nuc_index is not None:
                    coordinate_inserts.remove(seed_amino.consensus_nuc_index)
                    prev_consensus_nuc_index = seed_amino.consensus_nuc_index
                if skip_position is not None and coord_index == (skip_position-start_pos)//3 and amino_alignment.ref_end == skip_position and has_skipped_nucleotide:
                    skipped_nuc = self.reading_frames[amino_alignment.reading_frame][coord2conseq[(skip_position-start_pos)//3]+1].nucleotides[0]
                else:
                    skipped_nuc = None
                self.update_report_amino(coord_index,
                                         report_aminos,
                                         report_nucleotides,
                                         seed_amino,
                                         start_pos,
                                         repeat_position=repeat_position,
                                         skip_position=skip_position,
                                         skipped_nuc=skipped_nuc,)
            if prev_consensus_nuc_index is None:
                coordinate_inserts.clear()
            else:
                coordinate_inserts = {i
                                      for i in coordinate_inserts
                                      if i <= prev_consensus_nuc_index}
            self.inserts.update(coordinate_inserts)

    @staticmethod
    def update_report_amino(coord_index: int,
                            report_aminos: typing.List[ReportAmino],
                            report_nucleotides: typing.List[ReportNucleotide],
                            seed_amino: SeedAmino,
                            start_pos: int,
                            repeat_position: int = None,
                            skip_position: int = None,
                            skipped_nuc = None,):
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
                    report_nuc = report_nucleotides[report_nuc_index+1]
                    report_nuc.seed_nucleotide.add(skipped_nuc)
                if skip_position < start_pos + report_nuc_index:
                    report_nuc_index += 1
            report_nuc = report_nucleotides[report_nuc_index]
            report_nuc.seed_nucleotide.add(seed_nuc)

    def build_nucleotide_report(self,
                                start_pos: int,
                                end_pos: int,
                                report_nucleotides: typing.List[ReportNucleotide]):
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
                    consensus_nuc_index += section_size
                    # TODO: Record insertion positions.
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


@dataclass
class AminoAlignment:
    """ Part of a full alignment, with aligned aminos. """
    query_start: int
    query_end: int
    ref_start: int
    ref_end: int
    action: CigarActions
    reading_frame: int
    query: str = None  # Amino sequence
    ref: str = None  # Amino sequence
    aligned_query: str = None
    aligned_ref: str = None
    ref_amino_start: int = None

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

    def map_amino_sequences(self) -> typing.Dict[int, int]:
        """ Map reference amino indexes to query amino indexes. """
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
