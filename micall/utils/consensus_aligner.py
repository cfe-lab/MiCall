import typing
from enum import IntEnum
from operator import attrgetter

from gotoh import align_it, align_it_aa
from mappy import Alignment, Aligner

from micall.core.project_config import ProjectConfig
from micall.utils.report_amino import SeedAmino, ReportAmino, ReportNucleotide, SeedNucleotide

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
            cigar = [[blen, CigarActions.MATCH]]
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


class ConsensusAligner:
    def __init__(self, projects: ProjectConfig):
        self.projects = projects
        self.coordinate_name = self.consensus = self.algorithm = ''
        self.consensus_offset = 0
        self.alignments: typing.List[Alignment] = []
        self.reading_frames: typing.List[typing.List[SeedAmino]] = []
        self.seed_nucs: typing.List[SeedNucleotide] = []

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
        self.alignments = list(aligner.map(consensus))
        if self.alignments or 10_000 < len(consensus):
            self.algorithm = 'minimap2'
        else:
            self.algorithm = 'gotoh'
            self.align_gotoh(coordinate_seq, consensus)
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
        if len(consensus) < score:
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

    def clear(self):
        self.coordinate_name = self.consensus = self.algorithm = ''
        self.alignments.clear()
        self.seed_nucs.clear()

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
        if not self.alignments:
            return

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
        if repeat_position is None:
            repeat_index = None
        else:
            repeat_index = repeat_position - start_pos

        for alignment in self.alignments:
            ref_nuc_index = alignment.r_st
            consensus_nuc_index = alignment.q_st + self.consensus_offset
            frame_offset = -1
            region_start_consensus_amino_index = -1
            region_end_consensus_amino_index = -1
            source_amino_index = 0
            seed_aminos = self.reading_frames[0]
            region_seed_aminos: typing.List[SeedAmino] = []
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
                if frame_offset < 0:
                    # Skip over insertions at the start of the consensus before
                    # choosing the right reading frame.
                    frame_offset = (ref_nuc_index - start_pos + 1 -
                                    consensus_nuc_index) % 3
                    seed_aminos = self.reading_frames[frame_offset]
                section_nuc_index = 0
                while section_nuc_index < section_size:
                    if start_pos - 1 <= ref_nuc_index < end_pos:
                        while True:
                            seed_amino = seed_aminos[source_amino_index]
                            codon_nuc_index = (consensus_nuc_index -
                                               seed_amino.consensus_nuc_index)
                            if codon_nuc_index < 3:
                                break
                            source_amino_index += 1
                        if region_start_consensus_amino_index < 0:
                            region_start_consensus_amino_index = source_amino_index
                        region_end_consensus_amino_index = source_amino_index + 1
                        if ref_nuc_index + 1 == repeat_position:
                            region_seed_aminos = seed_aminos[
                                                 region_start_consensus_amino_index:
                                                 region_end_consensus_amino_index]
                            region_start_consensus_amino_index = -1
                            frame_offset = (frame_offset + 1) % 3
                            seed_aminos = self.reading_frames[frame_offset]
                    ref_nuc_index += 1
                    consensus_nuc_index += 1
                    section_nuc_index += 1
            region_seed_aminos.extend(
                seed_aminos[region_start_consensus_amino_index:
                            region_end_consensus_amino_index])
            amino_consensus = ''.join(seed_amino.get_consensus() or '?'
                                      for seed_amino in region_seed_aminos)
            aref, aconseq, score = align_aminos(amino_ref, amino_consensus)

            consensus_amino_index = ref_amino_index = 0
            consensus_nuc_index = -1
            repeat_count = 0
            for consensus_amino, ref_amino in zip(aconseq, aref):
                if consensus_amino == '-':
                    seed_amino = None
                else:
                    seed_amino = region_seed_aminos[consensus_amino_index]
                    consensus_amino_index += 1
                if ref_amino == '-':
                    report_amino = report_codon = None
                else:
                    report_amino = report_aminos[ref_amino_index]
                    codon_start = ref_amino_index*3
                    if repeat_index is not None and repeat_index < codon_start:
                        codon_start -= 1
                    report_codon = report_nucleotides[codon_start:codon_start+3]
                    ref_amino_index += 1
                if seed_amino is not None and report_amino is not None:
                    report_amino.seed_amino.add(seed_amino)
                    for seed_nuc, report_nuc in zip(seed_amino.nucleotides,
                                                    report_codon):
                        if seed_nuc.consensus_index <= consensus_nuc_index:
                            repeat_count += 1
                            continue
                        consensus_nuc_index = seed_nuc.consensus_index
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
