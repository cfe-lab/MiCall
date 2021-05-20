import typing
from enum import IntEnum
from operator import attrgetter

from gotoh import align_it
from mappy import Alignment, Aligner

from micall.core import aln2counts
from micall.core.project_config import ProjectConfig

CigarActions = IntEnum(
    'CigarActions',
    'MATCH INSERT DELETE SKIPPED SOFT_CLIPPED HARD_CLIPPED',
    start=0)


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
        self.seed_nucs: typing.List[aln2counts.SeedNucleotide] = []

    def align(self,
              coordinate_name: str = None,
              consensus: str = None,
              seed_aminos: typing.Iterable['aln2counts.SeedAmino'] = None):
        self.clear()
        if consensus:
            self.consensus = consensus
        elif seed_aminos:
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
        aligner = Aligner(seq=coordinate_seq, preset='map-ont', best_n=5)
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
                ref_index += 1
                consensus_index += 1
                expected_action = CigarActions.MATCH
                if nuc == '-':
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
