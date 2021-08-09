import math
import typing

from micall.core.aln2counts import SeedAmino, ReportAmino
from micall.utils.consensus_aligner import ConsensusAligner, AlignmentWrapper, CigarActions

# noinspection PyUnresolvedReferences
from micall.tests.test_remap import load_projects
from micall.utils.report_amino import ReportNucleotide


def assert_alignments(aligner: ConsensusAligner,
                      *expected_alignments: AlignmentWrapper):
    __tracebackhide__ = True
    wrapped_alignments = tuple(AlignmentWrapper.wrap(alignment)
                               for alignment in aligner.alignments)
    if repr(wrapped_alignments) != repr(expected_alignments):
        assert wrapped_alignments == expected_alignments
    for i, (wrapped_alignment, expected_alignment) in enumerate(
            zip(wrapped_alignments, expected_alignments)):
        for field_name in AlignmentWrapper.init_fields:
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
    alignment = AlignmentWrapper('R1', 0, 1001, 1100, 1, 1, 100)

    assert repr(alignment) == "AlignmentWrapper('R1', 0, 1001, 1100, 1, 1, 100)"


def test_wrap_overrides():
    alignment1 = AlignmentWrapper(r_st=100, r_en=200)
    alignment2 = AlignmentWrapper.wrap(alignment1, r_en=300, blen=200, cigar=[])
    expected_alignment = AlignmentWrapper(r_st=100, r_en=300, cigar=[])

    assert alignment2 == expected_alignment


def test_start_contig(projects):
    seed_name = 'SARS-CoV-2-seed'
    seed_seq = projects.getReference(seed_name)
    consensus = seed_seq[1000:2000]
    expected_alignment = AlignmentWrapper(ctg='N/A',
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
    expected_alignments = [AlignmentWrapper(ctg='N/A',
                                            ctg_len=len(seed_seq),
                                            r_st=6000,
                                            r_en=6500,
                                            q_st=0,
                                            q_en=500,
                                            mapq=60),
                           AlignmentWrapper(ctg='N/A',
                                            ctg_len=len(seed_seq),
                                            r_st=3000,
                                            r_en=3500,
                                            q_st=500,
                                            q_en=1000,
                                            mapq=60),
                           AlignmentWrapper(ctg='N/A',
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
                                 list(range(801, 900)) +
                                 list(range(903, 1001)),
                                 790,
                                 2292)


def test_start_contig_short_consensus(projects):
    """ Consensus is too short for minimap2, fall back to Gotoh. """
    seed_name = 'SARS-CoV-2-seed'
    seed_seq = projects.getReference(seed_name)
    start = 1560
    end = 1617
    consensus = seed_seq[start:end]
    expected_alignment = AlignmentWrapper(ctg='N/A',
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
    expected_alignment = AlignmentWrapper(ctg='N/A',
                                          ctg_len=len(seed_seq),
                                          r_st=2000,
                                          r_en=2060,
                                          q_st=0,
                                          q_en=59,
                                          mapq=9,
                                          cigar=[[30, CigarActions.MATCH],
                                                 [1, CigarActions.DELETE],
                                                 [29, CigarActions.MATCH]],
                                          NM=1)
    aligner = ConsensusAligner(projects)

    aligner.start_contig(seed_name, consensus)

    assert_alignments(aligner, expected_alignment)
    assert aligner.algorithm == 'minimap2'


def test_start_contig_deletion_gotoh(projects):
    seed_name = 'SARS-CoV-2-seed'
    seed_seq = projects.getReference(seed_name)
    consensus = seed_seq[2000:2030] + seed_seq[2031:2050]
    expected_alignment = AlignmentWrapper(ctg='N/A',
                                          ctg_len=len(seed_seq),
                                          r_st=2000,
                                          r_en=2050,
                                          q_st=0,
                                          q_en=49,
                                          mapq=0,
                                          cigar=[[30, CigarActions.MATCH],
                                                 [1, CigarActions.DELETE],
                                                 [19, CigarActions.MATCH]],
                                          NM=0)
    aligner = ConsensusAligner(projects)

    aligner.start_contig(seed_name, consensus)

    assert_alignments(aligner, expected_alignment)
    assert aligner.algorithm == 'gotoh'


def test_start_contig_matched_deletion_gotoh(projects):
    """ The consensus contains a deletion, but that aligns as a match. """
    seed_name = 'SARS-CoV-2-seed'
    seed_seq = projects.getReference(seed_name)
    consensus = seed_seq[2000:2030] + '-' + seed_seq[2031:2050]
    expected_alignment = AlignmentWrapper(ctg='N/A',
                                          ctg_len=len(seed_seq),
                                          r_st=2000,
                                          r_en=2050,
                                          q_st=0,
                                          q_en=50,
                                          mapq=0,
                                          cigar=[[50, CigarActions.MATCH]],
                                          NM=0)
    aligner = ConsensusAligner(projects)

    aligner.start_contig(seed_name, consensus)

    assert aligner.alignments[0].cigar_str == '50M'
    assert_alignments(aligner, expected_alignment)
    assert aligner.algorithm == 'gotoh'


def test_start_contig_insertion_minimap2(projects):
    seed_name = 'SARS-CoV-2-seed'
    seed_seq = projects.getReference(seed_name)
    consensus = seed_seq[2000:2030] + 'ACT' + seed_seq[2030:2060]
    expected_alignment = AlignmentWrapper(ctg='N/A',
                                          ctg_len=len(seed_seq),
                                          r_st=2000,
                                          r_en=2060,
                                          q_st=0,
                                          q_en=63,
                                          mapq=8,
                                          cigar=[[30, CigarActions.MATCH],
                                                 [3, CigarActions.INSERT],
                                                 [30, CigarActions.MATCH]],
                                          NM=3)
    aligner = ConsensusAligner(projects)

    aligner.start_contig(seed_name, consensus)

    assert_alignments(aligner, expected_alignment)
    assert aligner.algorithm == 'minimap2'


def test_start_contig_insertion_gotoh(projects):
    seed_name = 'SARS-CoV-2-seed'
    seed_seq = projects.getReference(seed_name)
    consensus = seed_seq[2000:2030] + 'T' + seed_seq[2030:2050]
    expected_alignment = AlignmentWrapper(ctg='N/A',
                                          ctg_len=len(seed_seq),
                                          r_st=2000,
                                          r_en=2050,
                                          q_st=0,
                                          q_en=51,
                                          mapq=0,
                                          cigar=[[30, CigarActions.MATCH],
                                                 [1, CigarActions.INSERT],
                                                 [20, CigarActions.MATCH]],
                                          NM=0)
    aligner = ConsensusAligner(projects)

    aligner.start_contig(seed_name, consensus)

    assert_alignments(aligner, expected_alignment)
    assert aligner.algorithm == 'gotoh'


def test_start_contig_with_only_primary_matches(projects):
    """ HXB2 has very similar sections in 5' LTR and 3' LTR. """
    seed_name = 'HIV1-B-FR-K03455-seed'
    seed_seq = projects.getReference(seed_name)
    consensus = seed_seq[:500]
    expected_alignment = AlignmentWrapper(ctg='N/A',
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
    expected_alignment = AlignmentWrapper(ctg='N/A',
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
