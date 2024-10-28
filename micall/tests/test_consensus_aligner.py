import math
import typing
from typing import Iterable, Tuple
import random
from io import StringIO
from pytest import approx

from micall.core.aln2counts import SeedAmino, ReportAmino
from micall.utils.consensus_aligner import ConsensusAligner, Alignment, AminoAlignment
from aligntools import CigarActions, Cigar
from micall.core.project_config import ProjectConfig

# noinspection PyUnresolvedReferences
from micall.tests.test_remap import load_projects
from micall.tests.utils import fixed_random_seed
from micall.utils.report_amino import ReportNucleotide


def mutate_sequence(rate, seq):
    def mutate(x):
        if random.random() >= rate:
            return x

        while True:
            y = random.choice(['A', 'C', 'G', 'T'])
            if y != x: return y

    with fixed_random_seed(42):
        return ''.join(map(mutate, seq))


def assert_alignments(aligner: ConsensusAligner,
                      *expected_alignments: Alignment):
    __tracebackhide__ = True
    wrapped_alignments = tuple(Alignment.coerce(alignment)
                               for alignment in aligner.alignments)
    if repr(wrapped_alignments) != repr(expected_alignments):
        assert wrapped_alignments == expected_alignments
    for i, (wrapped_alignment, expected_alignment) in enumerate(
            zip(wrapped_alignments, expected_alignments)):
        for field_name in dir(expected_alignment):
            if callable(getattr(expected_alignment, field_name)) or field_name.startswith('_'):
                continue
            wrapped = (i, field_name, getattr(wrapped_alignment, field_name))
            expected = (i, field_name, getattr(expected_alignment, field_name))
            assert wrapped == expected
    assert wrapped_alignments == expected_alignments


def assert_consensus_nuc_indexes(
        report_aminos: typing.List[ReportAmino],
        ref_positions_in_consensus: typing.List[typing.Optional[int]],
        start_pos: int,
        end_pos: int):
    __hide_traceback_frame__ = True
    amino_count = (end_pos - start_pos + 1) // 3
    expected_consensus_nuc_indexes: typing.List[typing.Optional[int]] = (
            [None] * amino_count)
    ref_positions_map = {pos: i
                         for i, pos in enumerate(ref_positions_in_consensus)}
    for pos in range(start_pos, end_pos+1, 3):
        for nuc_pos in range(pos, pos+3):
            consensus_nuc_index = ref_positions_map.get(nuc_pos)
            if consensus_nuc_index is not None:
                consensus_nuc_index -= nuc_pos - pos
                break
        else:
            continue
        amino_num = (pos - start_pos) // 3
        expected_consensus_nuc_indexes[amino_num] = consensus_nuc_index
    consensus_nuc_indexes = [report_amino.seed_amino.consensus_nuc_index
                             for report_amino in report_aminos]
    assert consensus_nuc_indexes == expected_consensus_nuc_indexes


def create_reading_frames(consensus: str) -> typing.Dict[int,
                                                         typing.List[SeedAmino]]:
    reading_frames = {}
    for frame_offset in range(3):
        shifted_consensus = ' '*frame_offset + consensus
        amino_count = math.ceil(len(shifted_consensus)/3)
        seed_aminos = []
        for i in range(amino_count):
            codon = shifted_consensus[i*3:(i+1)*3]
            seed_amino = SeedAmino(i*3-frame_offset)
            seed_amino.count_aminos(codon, 1)
            seed_aminos.append(seed_amino)
        reading_frames[frame_offset] = seed_aminos
    return reading_frames


def make_alignment(
        ctg='',
        ctg_len=0,
        r_st=0,
        r_en=0,
        strand=1,
        q_st=0,
        q_en=0,
        mapq=0,
        cigar: Iterable[Tuple[int, CigarActions]] = tuple(),
        cigar_str=None) -> Alignment:

    cigar = list(cigar)
    if not cigar:
        cigar = [(max(q_en-q_st, r_en-r_st), CigarActions.MATCH)]
    if cigar_str is None:
        cigar_str = str(Cigar(cigar))

    return Alignment(ctg=ctg,
                     ctg_len=ctg_len,
                     r_st=r_st,
                     r_en=r_en,
                     strand=strand,
                     q_st=q_st,
                     q_en=q_en,
                     mapq=mapq,
                     cigar=cigar,
                     cigar_str=cigar_str)


def test_create_reading_frames():
    reading_frames = create_reading_frames('AAACCCTTTGGG')

    assert len(reading_frames) == 3
    assert len(reading_frames[0]) == 4
    assert len(reading_frames[1]) == 5
    assert reading_frames[0][0].consensus_nuc_index == 0
    assert reading_frames[0][1].consensus_nuc_index == 3
    assert reading_frames[0][0].codon_counts == {'AAA': 1}
    assert reading_frames[1][0].consensus_nuc_index == -1
    assert reading_frames[1][0].codon_counts == {' AA': 1}


def test_alignment_repr():
    alignment = make_alignment('R1', 0, 1001, 1100, 1, 1, 100)

    assert repr(alignment) == "Alignment(ctg='R1', ctg_len=0, r_st=1001, r_en=1100, strand=1, q_st=1, q_en=100, mapq=0, cigar=[(99, CigarActions.MATCH)], cigar_str='99M')"


def test_start_contig(projects):
    seed_name = 'SARS-CoV-2-seed'
    seed_seq = projects.getReference(seed_name)
    consensus = seed_seq[1000:2000]
    expected_alignment = make_alignment(ctg='N/A',
                                   ctg_len=len(seed_seq),
                                   r_st=1000,
                                   r_en=2000,
                                   q_st=0,
                                   q_en=1000,
                                   mapq=60)
    aligner = ConsensusAligner(projects)

    aligner.start_contig(seed_name, consensus)

    assert aligner.consensus == consensus
    assert_alignments(aligner, expected_alignment)
    assert aligner.algorithm == 'minimap2'


def test_start_contig_multiple_sections(projects):
    seed_name = 'SARS-CoV-2-seed'
    seed_seq = projects.getReference(seed_name)
    consensus = seed_seq[6000:6500] + seed_seq[3000:3500] + seed_seq[1000:2000]
    expected_alignments = [make_alignment(ctg='N/A',
                                     ctg_len=len(seed_seq),
                                     r_st=6000,
                                     r_en=6500,
                                     q_st=0,
                                     q_en=500,
                                     mapq=60),
                           make_alignment(ctg='N/A',
                                     ctg_len=len(seed_seq),
                                     r_st=3000,
                                     r_en=3500,
                                     q_st=500,
                                     q_en=1000,
                                     mapq=60),
                           make_alignment(ctg='N/A',
                                     ctg_len=len(seed_seq),
                                     r_st=1000,
                                     r_en=2000,
                                     q_st=1000,
                                     q_en=2000,
                                     mapq=60)]
    aligner = ConsensusAligner(projects)

    aligner.start_contig(seed_name, consensus)

    assert_alignments(aligner, *expected_alignments)


def test_start_contig_overlapping_sections(projects):
    """ Similar sections can fool minimap2 into reporting a section twice.

    In this example, positions 1-60 of the read map to pos 4441-4500 of the
    reference. Positions 55-120 of the read map to pos 2995-3060 of the ref.
    Since positions 55-60 are in both alignments, remove them from the second
    one.
    """
    seed_name = 'HIV1-B-FR-K03455-seed'
    seed_seq = projects.getReference(seed_name)
    consensus = seed_seq[4440:4500] + seed_seq[3000:3060]
    reading_frames = create_reading_frames(consensus)
    int_ref = projects.getReference('INT')
    rt_ref = projects.getReference('RT')
    aligner = ConsensusAligner(projects)

    aligner.start_contig(seed_name, reading_frames=reading_frames)

    rt_aminos: typing.List[ReportAmino] = []
    rt_nucleotides: typing.List[ReportNucleotide] = []
    aligner.report_region(2550,
                          4229,
                          rt_nucleotides,
                          rt_aminos,
                          amino_ref=rt_ref)
    int_aminos: typing.List[ReportAmino] = []
    int_nucleotides: typing.List[ReportNucleotide] = []
    aligner.report_region(4230,
                          5096,
                          int_nucleotides,
                          int_aminos,
                          amino_ref=int_ref)

    assert_consensus_nuc_indexes(int_aminos,
                                 list(range(4441, 4501)) +
                                 list(range(3001, 3061)),
                                 4230,
                                 5096)
    assert_consensus_nuc_indexes(rt_aminos,
                                 list(range(4441, 4501)) +
                                 list(range(3001, 3061)),
                                 2550,
                                 4229)


# noinspection DuplicatedCode
def test_start_contig_big_deletion(projects):
    """ Split alignments around big deletions. """
    seed_name = 'HIV1-B-FR-K03455-seed'
    seed_seq = projects.getReference(seed_name)
    consensus = seed_seq[789:1282] + seed_seq[1863:2278]
    reading_frames = create_reading_frames(consensus)
    amino_ref = projects.getReference('HIV1B-gag')
    aligner = ConsensusAligner(projects)

    aligner.start_contig(seed_name, reading_frames=reading_frames)

    report_aminos: typing.List[ReportAmino] = []
    report_nucleotides: typing.List[ReportNucleotide] = []
    aligner.report_region(790,
                          2292,
                          report_nucleotides,
                          report_aminos,
                          amino_ref=amino_ref)

    assert_consensus_nuc_indexes(report_aminos,
                                 list(range(790, 1283)) +
                                 list(range(1864, 2279)),
                                 790,
                                 2292)
    assert report_aminos[163].seed_amino.get_consensus() == 'F'
    assert report_aminos[164].seed_amino.get_consensus() == '?'
    assert report_aminos[163].seed_amino.counts['F'] == 1
    assert report_aminos[164].seed_amino.deletions == 1


# noinspection DuplicatedCode
def test_start_contig_insert_and_big_deletion(projects):
    """ Split alignments around a big deletion after an insert. """
    seed_name = 'HIV1-B-FR-K03455-seed'
    seed_seq = projects.getReference(seed_name)
    consensus = (seed_seq[789:900] +
                 'ACTGAC' +
                 seed_seq[900:1282] +
                 seed_seq[1863:2278])
    reading_frames = create_reading_frames(consensus)
    amino_ref = projects.getReference('HIV1B-gag')
    aligner = ConsensusAligner(projects)

    aligner.start_contig(seed_name, reading_frames=reading_frames)

    report_aminos: typing.List[ReportAmino] = []
    report_nucleotides: typing.List[ReportNucleotide] = []
    aligner.report_region(790,
                          2292,
                          report_nucleotides,
                          report_aminos,
                          amino_ref=amino_ref)

    assert_consensus_nuc_indexes(report_aminos,
                                 list(range(790, 901)) +
                                 [None] * 6 +
                                 list(range(901, 1283)) +
                                 list(range(1864, 2279)),
                                 790,
                                 2292)


# noinspection DuplicatedCode
def test_start_contig_frame_change_insert(projects):
    """ Split alignments when the reading frame changes. """
    seed_name = 'HIV1-B-FR-K03455-seed'
    seed_seq = projects.getReference(seed_name)
    consensus = (seed_seq[800:900] +
                 'C' +
                 seed_seq[900:1000])
    reading_frames = create_reading_frames(consensus)
    amino_ref = projects.getReference('HIV1B-gag')
    aligner = ConsensusAligner(projects)

    aligner.start_contig(seed_name, reading_frames=reading_frames)
    report_aminos: typing.List[ReportAmino] = []
    report_nucleotides: typing.List[ReportNucleotide] = []
    aligner.report_region(790,
                          2292,
                          report_nucleotides,
                          report_aminos,
                          amino_ref=amino_ref)

    assert_consensus_nuc_indexes(report_aminos,
                                 list(range(801, 901)) +
                                 [None] +
                                 list(range(901, 1001)),
                                 790,
                                 2292)


# noinspection DuplicatedCode
def test_start_contig_frame_change_delete(projects):
    seed_name = 'HIV1-B-FR-K03455-seed'
    seed_seq = projects.getReference(seed_name)
    consensus = (seed_seq[800:899] + seed_seq[900:1000])
    reading_frames = create_reading_frames(consensus)
    amino_ref = projects.getReference('HIV1B-gag')
    aligner = ConsensusAligner(projects)

    aligner.start_contig(seed_name, reading_frames=reading_frames)
    report_aminos: typing.List[ReportAmino] = []
    report_nucleotides: typing.List[ReportNucleotide] = []
    aligner.report_region(790,
                          2292,
                          report_nucleotides,
                          report_aminos,
                          amino_ref=amino_ref)

    assert_consensus_nuc_indexes(report_aminos,
                                 list(range(801, 900)) +
                                 list(range(901, 1001)),
                                 790,
                                 2292)


# noinspection DuplicatedCode
def test_start_contig_frame_change_delete_across_vpr_boundary(projects):
    """ Because vpr has 1 base inserted in HXB2, amino ref is shorter. """
    seed_name = 'HIV1-B-FR-K03455-seed'
    seed_seq = projects.getReference(seed_name)
    consensus = (seed_seq[5500:5771] + seed_seq[5772:5849] + seed_seq[5850:6200])
    reading_frames = create_reading_frames(consensus)
    amino_ref = projects.getReference('HIV1B-vpr')
    aligner = ConsensusAligner(projects)

    aligner.start_contig(seed_name, reading_frames=reading_frames)
    report_aminos: typing.List[ReportAmino] = []
    report_nucleotides: typing.List[ReportNucleotide] = []
    aligner.report_region(5559,
                          5850,
                          report_nucleotides,
                          report_aminos,
                          amino_ref=amino_ref)

    # The expected positions look cleaner, because of the single base insertion
    # in HXB2 at 5772.
    assert_consensus_nuc_indexes(report_aminos,
                                 list(range(5501, 5847)),
                                 5559,
                                 5849)


# noinspection DuplicatedCode
def test_start_contig_close_frame_changes(projects):
    """ The reading frame changes and then changes back within 30 bases. """
    seed_name = 'HIV1-B-FR-K03455-seed'
    seed_seq = projects.getReference(seed_name)
    consensus = (seed_seq[800:899] + seed_seq[900:920] + seed_seq[922:1000])
    reading_frames = create_reading_frames(consensus)
    amino_ref = projects.getReference('HIV1B-gag')
    aligner = ConsensusAligner(projects)

    aligner.start_contig(seed_name, reading_frames=reading_frames)
    report_aminos: typing.List[ReportAmino] = []
    report_nucleotides: typing.List[ReportNucleotide] = []
    aligner.report_region(790,
                          2292,
                          report_nucleotides,
                          report_aminos,
                          amino_ref=amino_ref)

    assert_consensus_nuc_indexes(report_aminos,
                                 list(range(801, 904)) +
                                 list(range(907, 1001)),
                                 790,
                                 2292)


def test_start_contig_short_consensus(projects):
    """ Consensus is too short for minimap2, fall back to Gotoh. """
    seed_name = 'SARS-CoV-2-seed'
    seed_seq = projects.getReference(seed_name)
    start = 1560
    end = 1617
    consensus = seed_seq[start:end]
    expected_alignment = make_alignment(ctg='N/A',
                                        ctg_len=len(seed_seq),
                                        r_st=start,
                                        r_en=end,
                                        q_st=0,
                                        q_en=end-start)
    aligner = ConsensusAligner(projects)

    aligner.start_contig(seed_name, consensus)

    assert_alignments(aligner, expected_alignment)
    assert aligner.algorithm == 'gotoh'


def test_start_contig_deletion_minimap2(projects):
    seed_name = 'SARS-CoV-2-seed'
    seed_seq = projects.getReference(seed_name)
    consensus = seed_seq[2000:2030] + seed_seq[2031:2060]
    expected_alignment = make_alignment(ctg='N/A',
                                        ctg_len=len(seed_seq),
                                        r_st=2000,
                                        r_en=2060,
                                        q_st=0,
                                        q_en=59,
                                        mapq=9,
                                        cigar=[(30, CigarActions.MATCH),
                                               (1, CigarActions.DELETE),
                                               (29, CigarActions.MATCH)])
    aligner = ConsensusAligner(projects)

    aligner.start_contig(seed_name, consensus)

    assert_alignments(aligner, expected_alignment)
    assert aligner.algorithm == 'minimap2'


def test_start_contig_big_deletion_minimap2(projects):
    seed_name = 'HCV-1a'
    seed_seq = projects.getReference(seed_name)
    seed_seq = mutate_sequence(seq=seed_seq, rate=0.04)
    consensus = seed_seq[290:983] + seed_seq[3000:9269]

    expected_alignment = [make_alignment(ctg='N/A',
                                         ctg_len=len(seed_seq),
                                         r_st=290,
                                         r_en=983,
                                         q_st=0,
                                         q_en=693,
                                         mapq=60,
                                         cigar=[(693, CigarActions.MATCH)]),
                          make_alignment(ctg='N/A',
                                         ctg_len=len(seed_seq),
                                         r_st=3000,
                                         r_en=9269,
                                         q_st=693,
                                         q_en=6962,
                                         mapq=60,
                                         cigar=[(6269, CigarActions.MATCH)])]

    aligner = ConsensusAligner(projects)

    aligner.start_contig(seed_name, consensus)

    assert_alignments(aligner, *expected_alignment)
    assert aligner.algorithm == 'minimap2'


def test_start_contig_deletion_gotoh(projects):
    seed_name = 'SARS-CoV-2-seed'
    seed_seq = projects.getReference(seed_name)
    consensus = seed_seq[2000:2030] + seed_seq[2031:2050]
    expected_alignment = make_alignment(ctg='N/A',
                                        ctg_len=len(seed_seq),
                                        r_st=2000,
                                        r_en=2050,
                                        q_st=0,
                                        q_en=49,
                                        mapq=0,
                                        cigar=[(30, CigarActions.MATCH),
                                               (1, CigarActions.DELETE),
                                               (19, CigarActions.MATCH)])
    aligner = ConsensusAligner(projects)

    aligner.start_contig(seed_name, consensus)

    assert_alignments(aligner, expected_alignment)
    assert aligner.algorithm == 'gotoh'


def test_start_contig_matched_deletion_gotoh(projects):
    """ The consensus contains a deletion, but that aligns as a match. """
    seed_name = 'SARS-CoV-2-seed'
    seed_seq = projects.getReference(seed_name)
    consensus = seed_seq[2000:2030] + '-' + seed_seq[2031:2050]
    expected_alignment = make_alignment(ctg='N/A',
                                        ctg_len=len(seed_seq),
                                        r_st=2000,
                                        r_en=2050,
                                        q_st=0,
                                        q_en=50,
                                        mapq=0,
                                        cigar=[(50, CigarActions.MATCH)])
    aligner = ConsensusAligner(projects)

    aligner.start_contig(seed_name, consensus)

    assert aligner.alignments[0].cigar_str == '50M'
    assert_alignments(aligner, expected_alignment)
    assert aligner.algorithm == 'gotoh'


def test_start_contig_insertion_minimap2(projects):
    seed_name = 'SARS-CoV-2-seed'
    seed_seq = projects.getReference(seed_name)
    consensus = seed_seq[2000:2030] + 'ACT' + seed_seq[2030:2060]
    expected_alignment = make_alignment(ctg='N/A',
                                        ctg_len=len(seed_seq),
                                        r_st=2000,
                                        r_en=2060,
                                        q_st=0,
                                        q_en=63,
                                        mapq=8,
                                        cigar=[(30, CigarActions.MATCH),
                                               (3, CigarActions.INSERT),
                                               (30, CigarActions.MATCH)])
    aligner = ConsensusAligner(projects)

    aligner.start_contig(seed_name, consensus)

    assert_alignments(aligner, expected_alignment)
    assert aligner.algorithm == 'minimap2'


def test_start_contig_insertion_gotoh(projects):
    seed_name = 'SARS-CoV-2-seed'
    seed_seq = projects.getReference(seed_name)
    consensus = seed_seq[2000:2030] + 'T' + seed_seq[2030:2050]
    expected_alignment = make_alignment(ctg='N/A',
                                        ctg_len=len(seed_seq),
                                        r_st=2000,
                                        r_en=2050,
                                        q_st=0,
                                        q_en=51,
                                        mapq=0,
                                        cigar=[(30, CigarActions.MATCH),
                                               (1, CigarActions.INSERT),
                                               (20, CigarActions.MATCH)])
    aligner = ConsensusAligner(projects)

    aligner.start_contig(seed_name, consensus)

    assert_alignments(aligner, expected_alignment)
    assert aligner.algorithm == 'gotoh'


def test_start_contig_with_only_primary_matches(projects):
    """ HXB2 has very similar sections in 5' LTR and 3' LTR. """
    seed_name = 'HIV1-B-FR-K03455-seed'
    seed_seq = projects.getReference(seed_name)
    consensus = seed_seq[:500]
    expected_alignment = make_alignment(ctg='N/A',
                                        ctg_len=len(seed_seq),
                                        r_st=0,
                                        r_en=500,
                                        q_st=0,
                                        q_en=500,
                                        mapq=60)
    aligner = ConsensusAligner(projects)

    aligner.start_contig(seed_name, consensus)

    assert_alignments(aligner, expected_alignment)
    assert aligner.algorithm == 'minimap2'


# noinspection DuplicatedCode
def test_start_contig_reading_frames(projects):
    expected_consensus = 'AAACCCGGG'
    reading_frames = create_reading_frames(expected_consensus)
    seed_name = 'HCV-6t'
    seed_seq = projects.getReference(seed_name)
    expected_alignment = make_alignment(ctg='N/A',
                                        ctg_len=len(seed_seq),
                                        r_st=4798,
                                        r_en=4807,
                                        q_st=0,
                                        q_en=9,
                                        mapq=0)
    aligner = ConsensusAligner(projects)

    aligner.start_contig(seed_name, reading_frames=reading_frames)

    assert aligner.consensus == expected_consensus
    assert_alignments(aligner, expected_alignment)
    assert len(aligner.seed_nucs) == 9


def test_clear(projects):
    seed_name = 'SARS-CoV-2-seed'
    seed_seq = projects.getReference(seed_name)
    consensus = seed_seq[1000:2000]
    aligner = ConsensusAligner(projects)
    aligner.start_contig(seed_name, consensus)

    aligner.clear()

    assert aligner.consensus == ''
    assert aligner.algorithm == ''
    assert len(aligner.alignments) == 0


# noinspection DuplicatedCode
def test_start_contig_twice(projects):
    consensus1 = 'AAACCCGGG'
    reading_frames = create_reading_frames(consensus1)
    seed_name = 'HCV-6t'
    seed_seq = projects.getReference(seed_name)
    consensus2 = seed_seq[:100]
    aligner = ConsensusAligner(projects)
    aligner.start_contig(seed_name, reading_frames=reading_frames)

    aligner.start_contig(seed_name, consensus2)

    assert aligner.consensus == consensus2
    assert not aligner.seed_nucs


def test_start_contig_without_seed_name(projects):
    expected_consensus = 'ATACGCGTG'
    reading_frames = create_reading_frames(expected_consensus)
    seed_name = None
    aligner = ConsensusAligner(projects)
    aligner.start_contig(seed_name, reading_frames=reading_frames)

    assert aligner.consensus == expected_consensus


def test_start_contig_consensus_offset(projects):
    reading_frames = create_reading_frames('   CCCGGG')
    seed_name = 'HCV-6t'
    expected_consensus = 'CCCGGG'
    aligner = ConsensusAligner(projects)

    aligner.start_contig(seed_name, reading_frames=reading_frames)

    assert aligner.consensus == expected_consensus
    assert aligner.consensus_offset == 3


# noinspection DuplicatedCode
def test_report_region(projects):
    seed_name = 'SARS-CoV-2-seed'
    seed_seq = projects.getReference(seed_name)
    # ORF7b runs from 1-based positions 27756 to 27887.
    consensus = seed_seq[27000:27999]
    reading_frames = create_reading_frames(consensus)
    amino_ref = projects.getReference('SARS-CoV-2-ORF7b')
    aligner = ConsensusAligner(projects)
    aligner.start_contig(seed_name, reading_frames=reading_frames)

    report_aminos: typing.List[ReportAmino] = []
    report_nucleotides: typing.List[ReportNucleotide] = []
    aligner.report_region(27756,
                          27887,
                          report_nucleotides,
                          report_aminos,
                          amino_ref=amino_ref)

    assert len(report_aminos) == 44
    assert len(report_nucleotides) == 132  # 27887-27756+1
    assert report_aminos[0].seed_amino.consensus_nuc_index == 755  # currently, 758
    assert report_nucleotides[0].seed_nucleotide.consensus_index == 755


# noinspection DuplicatedCode
def test_report_region_nucleotides_only(projects):
    seed_name = 'SARS-CoV-2-seed'
    seed_seq = projects.getReference(seed_name)
    # TRS-B-8 runs from 1-based positions 28260 to 28273: ACGAACAAACTAAA
    consensus = seed_seq[28000:28999]
    reading_frames = create_reading_frames(consensus)
    aligner = ConsensusAligner(projects)
    aligner.start_contig(seed_name, reading_frames=reading_frames)

    report_nucleotides: typing.List[ReportNucleotide] = []
    aligner.report_region(28260, 28273, report_nucleotides)

    assert len(report_nucleotides) == 14  # 28273-28260+1
    assert report_nucleotides[0].seed_nucleotide.consensus_index == 259


# noinspection DuplicatedCode
def test_report_region_no_overlap(projects):
    seed_name = 'SARS-CoV-2-seed'
    seed_seq = projects.getReference(seed_name)
    consensus = seed_seq[27000:27999]
    reading_frames = create_reading_frames(consensus)
    amino_ref = 'F' * 100
    aligner = ConsensusAligner(projects)
    aligner.start_contig(seed_name, reading_frames=reading_frames)

    report_aminos: typing.List[ReportAmino] = []
    report_nucleotides: typing.List[ReportNucleotide] = []
    aligner.report_region(10_001,
                          10_300,
                          report_nucleotides,
                          report_aminos,
                          amino_ref=amino_ref)

    assert len(report_aminos) == 100
    assert len(report_nucleotides) == 300
    assert report_aminos[0].seed_amino.consensus_nuc_index is None
    assert report_nucleotides[0].seed_nucleotide.consensus_index is None


# noinspection DuplicatedCode
def test_report_region_after_start(projects):
    seed_name = 'SARS-CoV-2-seed'
    seed_seq = projects.getReference(seed_name)
    consensus = seed_seq[28000:28999]
    reading_frames = create_reading_frames(consensus)
    amino_ref = 'PVSYSLLF*M'
    aligner = ConsensusAligner(projects)
    aligner.start_contig(seed_name, reading_frames=reading_frames)

    report_aminos: typing.List[ReportAmino] = []
    report_nucleotides: typing.List[ReportNucleotide] = []
    aligner.report_region(27998,
                          28027,
                          report_nucleotides,
                          report_aminos,
                          amino_ref=amino_ref)

    assert len(report_aminos) == 10
    assert report_aminos[0].seed_amino.consensus_nuc_index is None
    assert report_aminos[1].seed_amino.consensus_nuc_index == 0
    assert report_aminos[-1].seed_amino.consensus_nuc_index == 24


# noinspection DuplicatedCode
def test_report_region_before_end(projects):
    seed_name = 'SARS-CoV-2-seed'
    seed_seq = projects.getReference(seed_name)
    consensus = seed_seq[28000:28999]
    reading_frames = create_reading_frames(consensus)
    amino_ref = 'QQQGQ'
    aligner = ConsensusAligner(projects)
    aligner.start_contig(seed_name, reading_frames=reading_frames)

    report_aminos: typing.List[ReportAmino] = []
    report_nucleotides: typing.List[ReportNucleotide] = []
    aligner.report_region(28991,
                          29005,
                          report_nucleotides,
                          report_aminos,
                          amino_ref=amino_ref)

    assert len(report_aminos) == 5
    assert report_aminos[0].seed_amino.consensus_nuc_index == 990
    assert report_aminos[-1].seed_amino.consensus_nuc_index is None


# noinspection DuplicatedCode
def test_report_region_with_repeated_nucleotide(projects):
    seed_name = 'SARS-CoV-2-seed'
    seed_seq = projects.getReference(seed_name)
    # nsp12 runs from 1-based positions 13442 to 16236, with repeat at 13468.
    consensus = seed_seq[13000:13500]
    reading_frames = create_reading_frames(consensus)
    amino_ref = projects.getReference('SARS-CoV-2-nsp12')
    aligner = ConsensusAligner(projects)
    aligner.start_contig(seed_name, reading_frames=reading_frames)

    report_aminos: typing.List[ReportAmino] = []
    report_nucleotides: typing.List[ReportNucleotide] = []
    aligner.report_region(13442,
                          16236,
                          report_nucleotides,
                          report_aminos,
                          repeat_position=13468,
                          amino_ref=amino_ref)

    amino_consensus = ''.join(report_amino.seed_amino.get_consensus()
                              for report_amino in report_aminos[:21])
    assert amino_consensus == 'SADAQSFLNRVCGVSAARLT-'
    nuc_consensus = ''.join(report_nucleotide.seed_nucleotide.get_consensus('MAX')
                            for report_nucleotide in report_nucleotides)
    assert nuc_consensus == seed_seq[13441:13500]


# noinspection DuplicatedCode
def test_alignments_file(projects):
    alignments_file = StringIO()
    unmerged_alignments_file = StringIO()
    intermediate_alignments_file = StringIO()
    _ = ConsensusAligner(projects,
                         alignments_file=alignments_file,
                         unmerged_alignments_file=unmerged_alignments_file,
                         intermediate_alignments_file=intermediate_alignments_file)
    expected_text = """\
coordinate_name,contig,action,query_start,query_end,ref_start,ref_end,reading_frame,\
ref_amino_start,aligned_query,aligned_ref
"""

    assert alignments_file.getvalue() == expected_text
    assert unmerged_alignments_file.getvalue() == expected_text
    assert intermediate_alignments_file.getvalue() == expected_text


# noinspection DuplicatedCode
def test_overall_alignments_file(projects):
    overall_alignments_file = StringIO()
    _ = ConsensusAligner(projects, overall_alignments_file=overall_alignments_file)
    expected_text = """\
coordinate_name,contig,query_start,query_end,consensus_offset,ref_start,ref_end,cigar_str
"""

    assert overall_alignments_file.getvalue() == expected_text


# noinspection DuplicatedCode
def test_write_alignment(projects):
    alignments_file = StringIO()
    aligner = ConsensusAligner(projects, alignments_file=alignments_file, contig_name='test_contig')
    aligner.coordinate_name = "test"
    amino_alignments = [AminoAlignment(0, 10, 20, 30, 0, 1, aligned_query="ATA", aligned_ref="ACA", ref_amino_start=10)]

    expected_text = """\
coordinate_name,contig,action,query_start,query_end,ref_start,ref_end,reading_frame,\
ref_amino_start,aligned_query,aligned_ref
test,test_contig,MATCH,0,10,20,30,1,10,ATA,ACA
"""

    aligner.write_alignments_file(amino_alignments, aligner.alignments_writer)

    assert alignments_file.getvalue() == expected_text


# noinspection DuplicatedCode
def test_write_multiple_alignments(projects):
    alignments_file = StringIO()
    aligner = ConsensusAligner(projects, alignments_file=alignments_file, contig_name='test_contig')
    aligner.coordinate_name = "test"
    amino_alignments = [AminoAlignment(0, 10, 20, 30, 0, 1, aligned_query="ATA", aligned_ref="ACA", ref_amino_start=10)]
    aligner.write_alignments_file(amino_alignments, aligner.alignments_writer)
    aligner2 = ConsensusAligner(projects, alignments_file=alignments_file, contig_name='test_contig_2')
    aligner2.coordinate_name = "2ndtest"
    amino_alignments2 = [AminoAlignment(0, 20, 40, 60, 0, 2, aligned_query="ABC", aligned_ref="DEF", ref_amino_start=4)]
    aligner2.write_alignments_file(amino_alignments2, aligner2.alignments_writer)

    expected_text = """\
coordinate_name,contig,action,query_start,query_end,ref_start,ref_end,reading_frame,\
ref_amino_start,aligned_query,aligned_ref
test,test_contig,MATCH,0,10,20,30,1,10,ATA,ACA
2ndtest,test_contig_2,MATCH,0,20,40,60,2,4,ABC,DEF
"""

    assert alignments_file.getvalue() == expected_text


def test_count_coord_concordance():
    """Check the concordance for a 27 base consensus with a 27 base match."""
    projects = ProjectConfig()
    projects.load(StringIO("""\
    {"genotype_references": {"test-region": {"is_nucleotide": true,"reference": ["AGATTTCGATGATTCAGAAGATAAGCA"]}}}
    """))
    aligner = ConsensusAligner(projects)
    aligner.consensus = "AGATTTCGATGATTCAGAAGATAAGCA"
    aligner.coordinate_name = 'test-region'
    aligner.alignments = [make_alignment(r_st=0, r_en=27, q_st=0, q_en=27, cigar=[(27, CigarActions.MATCH)])]

    expected_concordance_list = [1.0]*len(aligner.consensus)

    concordance_list = aligner.coord_concordance()

    assert concordance_list == approx(expected_concordance_list)


# noinspection DuplicatedCode
def test_count_coord_concordance_mismatch():
    """Check the concordance for a 27 base consensus with a 27 base match, but two nucleotides mismatch."""
    projects = ProjectConfig()
    projects.load(StringIO("""\
    {"genotype_references": {"test-region": {"is_nucleotide": true,"reference": ["AGATTTCGATGATTCAGAAGATAAGCA"]}}}
    """))
    aligner = ConsensusAligner(projects)
    aligner.consensus = "AGATTTCGATGATTCAGAAGATTTGCA"
    # changed nucs:                            ^^
    aligner.coordinate_name = 'test-region'
    aligner.alignments = [make_alignment(r_st=0, r_en=27, q_st=0, q_en=27, cigar=[(27, CigarActions.MATCH)])]

    # At the end of the consensus, the size of the averaging window for the concordance decreases from 20 to 11.
    # The concordance therefore decreases from 18/20 to 9/11
    expected_concordance_list = [1.0]*13 + [19/20, 18/20, 18/20, 18/20, 18/20] +\
                                [17/19, 16/18, 15/17, 14/16, 13/15, 12/14, 11/13, 10/12, 9/11]

    concordance_list = aligner.coord_concordance()

    assert concordance_list == approx(expected_concordance_list)


# noinspection DuplicatedCode
def test_count_coord_concordance_short_match():
    """Check the concordance for a 27 base consensus with a 15 base match."""
    projects = ProjectConfig()
    projects.load(StringIO("""\
    {"genotype_references": {"test-region": {"is_nucleotide": true,"reference": ["AGATTTCGATGATTCAGAAGATTTGCATTT"]}}}
    """))
    aligner = ConsensusAligner(projects)
    aligner.consensus = "AGATTTCGATGATTCTCTTCTAAACGT"
    # last match position:             ^
    aligner.coordinate_name = 'test-region'
    aligner.alignments = [make_alignment(r_st=0, r_en=15, q_st=0, q_en=15, cigar=[(15, CigarActions.MATCH)])]
    # We start out with 100% match for the first 6 positions
    expected_concordance_list = [1.0] * 6
    # After that, the averaging window (whose size is still increasing) starts to slide past the match:
    # 15 matches / 16 positions to 15 matches / 19 positions
    expected_concordance_list += [15/16, 15/17, 15/18, 15/19]
    # We then reach the full size of the averaging window (20), and the concordance decreases from 15/20
    # in 1 base intervals down to 8 bases of the match left in the window (8/20).
    expected_concordance_list += [15/20, 14/20, 13/20, 12/20, 11/20, 10/20, 9/20, 8/20]
    # After that, the averaging window starts to become smaller: 7/19 to 1/13,
    # until no more bases of the match are left in the window.
    expected_concordance_list += [7/19, 6/18, 5/17, 4/16, 3/15, 2/14, 1/13, 0, 0]

    concordance_list = aligner.coord_concordance()

    assert concordance_list == approx(expected_concordance_list)


# noinspection DuplicatedCode
def test_count_coord_concordance_two_matches():
    """Check the concordance for a 30 base consensus with two matches and a gap."""
    projects = ProjectConfig()
    projects.load(StringIO("""\
    {"genotype_references": {"test-region": {"is_nucleotide": true,"reference": ["AGATTTCGATGATTCAGAAGATTTGCATTT"]}}}
    """))
    aligner = ConsensusAligner(projects)
    aligner.consensus = "AGATTTCGATGATTCAGAAGATTTGCATTT"
    aligner.coordinate_name = 'test-region'
    aligner.alignments = [make_alignment(r_st=0, r_en=12, q_st=0, q_en=12, cigar=[(12, CigarActions.MATCH)]),
                          make_alignment(r_st=15, r_en=30, q_st=15, q_en=30, cigar=[(15, CigarActions.MATCH)])]

    expected_concordance_list = [1.0] * 3 + [12/13, 12/14, 12/15, 13/16, 14/17, 15/18, 16/19] + [17/20]*11 + \
                                [16/19, 15/18, 15/17, 15/16] + [1.0]*5

    concordance_list = aligner.coord_concordance()

    assert concordance_list == approx(expected_concordance_list)


# noinspection DuplicatedCode
def test_count_coord_concordance_with_insertion():
    """Check the concordance for a 30 base consensus with a 27 base match and an insertion of three."""
    projects = ProjectConfig()
    projects.load(StringIO("""\
    {"genotype_references": {"test-region": {"is_nucleotide": true,"reference": ["AGATTTCGATGATTCAGAAGATTTGCA"]}}}
    """))
    aligner = ConsensusAligner(projects)
    aligner.consensus = "AGATTTCGACCCTGATTCAGAAGATTTGCA"
    # insertion:                  ^^^
    aligner.coordinate_name = 'test-region'
    aligner.alignments = [make_alignment(r_st=0, r_en=27, q_st=0, q_en=30, cigar=[(9, CigarActions.MATCH),
                                                                                  (3, CigarActions.INSERT),
                                                                                  (18, CigarActions.MATCH)])]
    # the window size increases from 10 to 20, while the averaging window slides over the insertion
    expected_concordance_list = [9/10, 9/11, 9/12, 10/13, 11/14, 12/15, 13/16, 14/17, 15/18, 16/19]
    # for 10 positions in the middle, the insertion is included in the full window size fo 20
    expected_concordance_list += [17/20]*10
    # the averaging window then slides over the insertion and starts becoming smaller
    expected_concordance_list += [18/20, 18/19] + [1.0]*8

    concordance_list = aligner.coord_concordance()

    assert concordance_list == approx(expected_concordance_list)


# noinspection DuplicatedCode
def test_count_coord_concordance_with_deletion():
    """Check the concordance for a 24 base consensus with a 27 base match and a deletion of three."""
    projects = ProjectConfig()
    projects.load(StringIO("""\
    {"genotype_references": {"test-region": {"is_nucleotide": true,"reference": ["AGATTTCGATGATTCAGAAGATTTGCA"]}}}
    """))
    aligner = ConsensusAligner(projects)
    aligner.consensus = "AGATTTCGATTCAGAAGATTTGCA"
    # deletion behind this pos:  ^
    aligner.coordinate_name = 'test-region'
    aligner.alignments = [make_alignment(r_st=0, r_en=27, q_st=0, q_en=30, cigar=[(9, CigarActions.MATCH),
                                                                                  (3, CigarActions.DELETE),
                                                                                  (15, CigarActions.MATCH)])]
    # the deletion does not decrease the concordance
    expected_concordance_list = [1.0]*len(aligner.consensus)

    concordance_list = aligner.coord_concordance()

    assert concordance_list == approx(expected_concordance_list)


# noinspection DuplicatedCode
def test_count_seed_region_concordance(projects):
    seed_concordance_file = StringIO()
    aligner = ConsensusAligner(projects,
                               seed_concordance_file=seed_concordance_file,
                               contig_name='test-contig')
    aligner.consensus = "AGATTTCGATGATTCAGAAGATTTGCA"
    seed_name = 'test-seed'
    seed_ref = "AGATTTCGATGATTCAGAAGATTTGCA"
    region = 'test-region'
    seed_alignments = [make_alignment(r_st=0, r_en=27, q_st=0, q_en=27, cigar=[(27, CigarActions.MATCH)])]

    expected_file = """\
seed_name,contig,region,pct_concordance,pct_covered
test-seed,test-contig,test-region,100.0,100.0
"""

    aligner.region_seed_concordance(region, seed_name, seed_alignments, seed_ref, 0, 27)

    assert seed_concordance_file.getvalue() == expected_file


# noinspection DuplicatedCode
def test_count_seed_region_concordance_mismatch(projects):
    seed_concordance_file = StringIO()
    aligner = ConsensusAligner(projects,
                               seed_concordance_file=seed_concordance_file,
                               contig_name='test-contig')
    aligner.consensus = "AGATTTCGATGATTCAGACCCCCCGCATGA"
    # mismatch:                            ^^^^^^
    seed_name = 'test-seed'
    seed_ref = "AGATTTCGATGATTCAGAAGATTTGCATGA"
    region = 'test-region'
    seed_alignments = [make_alignment(r_st=0, r_en=30, q_st=0, q_en=30, cigar=[(30, CigarActions.MATCH)])]

    expected_file = """\
seed_name,contig,region,pct_concordance,pct_covered
test-seed,test-contig,test-region,80.0,100.0
"""

    aligner.region_seed_concordance(region, seed_name, seed_alignments, seed_ref, 0, 30)

    assert seed_concordance_file.getvalue() == expected_file


# noinspection DuplicatedCode
def test_count_seed_region_concordance_seed_not_aligned(projects):
    seed_concordance_file = StringIO()
    aligner = ConsensusAligner(projects,
                               seed_concordance_file=seed_concordance_file,
                               contig_name='test-contig')
    aligner.consensus = "AGATTTCGATGATTCAGAAGATTTGCATGA"
    seed_name = 'test-seed'
    seed_ref = "AGATTTCGATGATTCAGAAGATTTGCATGA"
    region = 'test-region'
    seed_alignments = [make_alignment(r_st=0, r_en=15, q_st=0, q_en=15, cigar=[(15, CigarActions.MATCH)])]

    expected_file = """\
seed_name,contig,region,pct_concordance,pct_covered
test-seed,test-contig,test-region,100.0,50.0
"""

    aligner.region_seed_concordance(region, seed_name, seed_alignments, seed_ref, 0, 30)

    assert seed_concordance_file.getvalue() == expected_file


# noinspection DuplicatedCode
def test_count_seed_region_concordance_larger_match(projects):
    seed_concordance_file = StringIO()
    aligner = ConsensusAligner(projects,
                               seed_concordance_file=seed_concordance_file,
                               contig_name='test-contig')
    aligner.consensus = "AGATTTCGATGATTCAGAAGATTTGCATGA"
    seed_name = 'test-seed'
    seed_ref = "AGATTTCGATGATTCAGAAGATTTGCATGA"
    region = 'test-region'
    seed_alignments = [make_alignment(r_st=0, r_en=30, q_st=0, q_en=30, cigar=[(30, CigarActions.MATCH)])]

    expected_file = """\
seed_name,contig,region,pct_concordance,pct_covered
test-seed,test-contig,test-region,100.0,100.0
"""

    aligner.region_seed_concordance(region, seed_name, seed_alignments, seed_ref, 0, 15)

    assert seed_concordance_file.getvalue() == expected_file


# noinspection DuplicatedCode
def test_count_seed_region_concordance_insertion(projects):
    seed_concordance_file = StringIO()
    aligner = ConsensusAligner(projects,
                               seed_concordance_file=seed_concordance_file,
                               contig_name='test-contig')
    aligner.consensus = "AGATTTCGACCCTGATTCAGAAGATTTGCA"
    # insert here:                ^^^
    seed_name = 'test-seed'
    seed_ref = "AGATTTCGATGATTCAGAAGATTTGCA"
    region = 'test-region'
    seed_alignments = [make_alignment(r_st=0, r_en=27, q_st=0, q_en=30, cigar=[(9, CigarActions.MATCH),
                                                                               (3, CigarActions.INSERT),
                                                                               (18, CigarActions.MATCH)])]

    expected_file = """\
seed_name,contig,region,pct_concordance,pct_covered
test-seed,test-contig,test-region,100.0,100.0
"""

    aligner.region_seed_concordance(region, seed_name, seed_alignments, seed_ref, 0, 27)

    assert seed_concordance_file.getvalue() == expected_file


# noinspection DuplicatedCode
def test_count_seed_region_concordance_deletion(projects):
    seed_concordance_file = StringIO()
    aligner = ConsensusAligner(projects,
                               seed_concordance_file=seed_concordance_file,
                               contig_name='test-contig')
    aligner.consensus = "AGATTTCGATTCAGAAGATTTGCATGA"
    # deletion after this nuc:   ^
    seed_name = 'test-seed'
    seed_ref = "AGATTTCGATGATTCAGAAGATTTGCATGA"
    region = 'test-region'
    seed_alignments = [make_alignment(r_st=0, r_en=30, q_st=0, q_en=27, cigar=[(9, CigarActions.MATCH),
                                                                               (3, CigarActions.DELETE),
                                                                               (18, CigarActions.MATCH)])]

    expected_file = """\
seed_name,contig,region,pct_concordance,pct_covered
test-seed,test-contig,test-region,100.0,90.0
"""

    aligner.region_seed_concordance(region, seed_name, seed_alignments, seed_ref, 0, 30)

    assert seed_concordance_file.getvalue() == expected_file
