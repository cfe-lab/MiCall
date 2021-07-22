import typing
from dataclasses import dataclass
from enum import IntEnum
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

    def split_alignments(self):
        """ Split alignments into sections that can be translated to aminos. """

        # Remove overlaps between alignments.
        max_query_pos = -1
        for alignment_num, alignment in enumerate(self.alignments):
            query_start = alignment.q_st
            if query_start < max_query_pos:
                offset = max_query_pos - alignment.q_st
                new_cigar = [section[:] for section in alignment.cigar]
                new_cigar[0][0] -= offset
                new_alignment = AlignmentWrapper.wrap(alignment,
                                                      q_st=max_query_pos,
                                                      r_st=alignment.r_st+offset,
                                                      cigar=new_cigar,
                                                      # Unneeded fields => -1.
                                                      mlen=-1,
                                                      blen=-1,
                                                      NM=-1)
                self.alignments[alignment_num] = new_alignment
            max_query_pos = alignment.q_en

        # Split alignments around frame shifts and big deletions.
        new_alignments = []
        for alignment in self.alignments:
            breakpoints = [AlignmentBreakpoint(alignment.q_st,
                                               alignment.q_st,
                                               alignment.r_st,
                                               alignment.r_st,
                                               0,
                                               0)]
            query_end = alignment.q_st
            ref_end = alignment.r_st
            old_reading_frame = (alignment.r_st - alignment.q_st) % 3
            for cigar_index, (size, action) in enumerate(alignment.cigar):
                if action == CigarActions.MATCH:
                    query_end += size
                    ref_end += size
                elif action == CigarActions.INSERT:
                    query_end += size
                    new_reading_frame = (ref_end - query_end) % 3
                    if old_reading_frame != new_reading_frame:
                        breakpoints.append(AlignmentBreakpoint(query_end-size,
                                                               query_end,
                                                               ref_end,
                                                               ref_end,
                                                               cigar_index,
                                                               cigar_index+1))
                        old_reading_frame = new_reading_frame
                else:
                    assert action == CigarActions.DELETE
                    new_reading_frame = (ref_end + size - query_end)
                    if new_reading_frame != old_reading_frame or 6 < size:
                        breakpoints.append(AlignmentBreakpoint(query_end,
                                                               query_end,
                                                               ref_end,
                                                               ref_end+size,
                                                               cigar_index,
                                                               cigar_index+1))
                        old_reading_frame = new_reading_frame
                    ref_end += size
            breakpoints.append(AlignmentBreakpoint(alignment.q_en,
                                                   alignment.q_en,
                                                   alignment.r_en,
                                                   alignment.r_en,
                                                   len(alignment.cigar),
                                                   len(alignment.cigar)))
            filtered_breakpoints = [breakpoints[0]]
            prev_breakpoints = iter(breakpoints)
            mid_breakpoints = iter(breakpoints)
            next_breakpoints = iter(breakpoints)
            next(mid_breakpoints)
            next(next_breakpoints)
            next(next_breakpoints)
            for prev_breakpoint, mid_breakpoint, next_breakpoint in zip(
                    prev_breakpoints, mid_breakpoints, next_breakpoints):
                prev_size = (mid_breakpoint.from_query_pos -
                             prev_breakpoint.to_query_pos)
                next_size = (next_breakpoint.from_query_pos -
                             mid_breakpoint.to_query_pos)
                if MINIMUM_READING_FRAME_SHIFT <= min(prev_size, next_size):
                    filtered_breakpoints.append(mid_breakpoint)
            filtered_breakpoints.append(breakpoints[-1])
            if len(filtered_breakpoints) == 2:
                new_alignments.append(alignment)
            else:
                prev_breakpoints = iter(filtered_breakpoints)
                next_breakpoints = iter(filtered_breakpoints)
                next(next_breakpoints)
                for prev_breakpoint, next_breakpoint in zip(
                        prev_breakpoints, next_breakpoints):
                    new_alignment = AlignmentWrapper.wrap(
                        alignment,
                        q_st=prev_breakpoint.to_query_pos,
                        q_en=next_breakpoint.from_query_pos,
                        r_st=prev_breakpoint.to_ref_pos,
                        r_en=next_breakpoint.from_ref_pos,
                        cigar=alignment.cigar[prev_breakpoint.to_cigar:
                                              next_breakpoint.from_cigar],
                        # Unneeded fields => -1.
                        mlen=-1,
                        blen=-1,
                        NM=-1)
                    new_alignments.append(new_alignment)
        self.alignments = new_alignments

    def find_amino_alignments(self,
                              start_pos: int,
                              end_pos: int,
                              amino_ref: str):
        translations = {
            reading_frame: translate(
                '-'*reading_frame +
                self.consensus +
                '-'*((3-len(self.consensus)-reading_frame) % 3))
            for reading_frame in range(3)}
        self.amino_alignments.clear()
        for alignment in self.alignments:
            amino_sections = []
            query_end = alignment.q_st
            ref_end = alignment.r_st
            for cigar_index, (size, action) in enumerate(alignment.cigar):
                reading_frame = (ref_end - query_end) % 3
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
                    start_shift = max(0, start_pos-1 - amino_alignment.ref_start)
                    amino_alignment.ref_start += start_shift
                    amino_alignment.query_start += start_shift
                    end_shift = max(0, amino_alignment.ref_end - end_pos)
                    amino_alignment.ref_end -= end_shift
                    amino_alignment.query_end -= end_shift
                    amino_sections.append(amino_alignment)
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
                if can_align and (size % 3 != 0 or size//3 <= MAXIMUM_AMINO_GAP):
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
                                    amino_ref)

    def build_amino_report(self,
                           start_pos: int,
                           end_pos: int,
                           report_nucleotides: typing.List[ReportNucleotide],
                           report_aminos: typing.List[ReportAmino] = None,
                           repeat_position: int = None,
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
        :param amino_ref: amino acid sequence to align with, or None if only
            nucleotides are reported.
        """
        self.find_amino_alignments(start_pos, end_pos, amino_ref)
        for amino_alignment in self.amino_alignments:
            if amino_alignment.action != CigarActions.MATCH:
                continue
            region_seed_aminos = self.reading_frames[amino_alignment.reading_frame]
            coord2conseq = amino_alignment.map_amino_sequences()

            coordinate_inserts = {seed_amino.consensus_nuc_index
                                  for seed_amino in region_seed_aminos}
            prev_conseq_index = None
            prev_consensus_nuc_index = None
            for coord_index in range(len(amino_ref)):
                conseq_index = coord2conseq.get(coord_index)
                if conseq_index is None:
                    seed_amino = SeedAmino(None)
                    if (prev_conseq_index is not None and
                            prev_conseq_index+1 < len(region_seed_aminos)):
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
                    if prev_conseq_index is None:
                        coordinate_inserts = {i
                                              for i in coordinate_inserts
                                              if i >= seed_amino.consensus_nuc_index}
                    prev_conseq_index = conseq_index

                report_amino = report_aminos[coord_index]
                report_amino.seed_amino.add(seed_amino)
                if seed_amino.consensus_nuc_index is not None:
                    coordinate_inserts.remove(seed_amino.consensus_nuc_index)
                    prev_consensus_nuc_index = seed_amino.consensus_nuc_index
                ref_nuc_pos = coord_index*3 + start_pos
                for codon_nuc_index, seed_nuc in enumerate(
                        seed_amino.nucleotides):
                    if len(report_nucleotides) <= coord_index*3 + codon_nuc_index:
                        continue

                    report_nuc_index = coord_index*3 + codon_nuc_index
                    if repeat_position is not None and start_pos < repeat_position:
                        if repeat_position < ref_nuc_pos:
                            report_nuc_index -= 1
                        if repeat_position == ref_nuc_pos-1 and codon_nuc_index == 0:
                            continue
                    report_nuc = report_nucleotides[report_nuc_index]
                    report_nuc.seed_nucleotide.add(seed_nuc)
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
                    # TODO: Increase deletion counts.
                    ref_nuc_index += section_size
                    continue
                assert section_action == CigarActions.MATCH, section_action
                section_nuc_index = 0
                while section_nuc_index < section_size:
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
                    section_nuc_index += 1


@dataclass
class AlignmentBreakpoint:
    """ Describe the pieces to remove from an alignment. """
    from_query_pos: int
    to_query_pos: int
    from_ref_pos: int
    to_ref_pos: int
    from_cigar: int
    to_cigar: int


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
        before_end = self.ref_start <= end_pos
        after_start = start_pos <= self.ref_end
        return before_end == after_start

    def find_reading_frame(self, amino_ref, start_pos, translations):
        ref_amino_start = (self.ref_start - start_pos + 1) // 3
        ref_amino_end = (self.ref_end - start_pos + 3) // 3
        self.ref_amino_start = ref_amino_start

        ref_section = amino_ref[ref_amino_start:ref_amino_end]
        if ref_amino_end - ref_amino_start < MINIMUM_AMINO_ALIGNMENT:
            # Too small to usefully align.
            possible_frames = [self.reading_frame]
            max_score = 0
        else:
            possible_frames = list(range(3))
            # Exclamation marks avoid extra gaps near the ends.
            assert '!' not in ref_section
            ref_section = '!' + ref_section + '!'
            max_score = (ref_amino_end - ref_amino_start) // 2

        for reading_frame in range(3):
            query_amino_start = (self.query_start + reading_frame) // 3
            query_amino_end = (self.query_end + reading_frame + 2) // 3
            translation = translations[reading_frame]
            query_section = translation[query_amino_start:query_amino_end]
            if len(possible_frames) == 1:
                aref = ref_section
                aquery = query_section
                score = 1
            else:
                aref, aquery, score = align_aminos(ref_section, query_section)
                aref = aref.replace('!', '-')
                ref_section = ref_section[1:-1]
            if max_score < score:
                max_score = score
                self.reading_frame = reading_frame

                # Record alignments after trimming exclamation marks.
                self.query = query_section
                self.ref = ref_section
                self.aligned_query = aquery
                self.aligned_ref = aref

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
