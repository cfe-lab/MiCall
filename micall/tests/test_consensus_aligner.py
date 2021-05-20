from micall.core.aln2counts import SeedAmino
from micall.utils.consensus_aligner import ConsensusAligner, AlignmentWrapper, CigarActions

# noinspection PyUnresolvedReferences
from micall.tests.test_remap import load_projects


def assert_alignments(aligner: ConsensusAligner,
                      *expected_alignments: AlignmentWrapper):
    __tracebackhide__ = True
    wrapped_alignments = tuple(AlignmentWrapper.wrap(alignment)
                               for alignment in aligner.alignments)
    assert wrapped_alignments == expected_alignments


def test_alignment_repr():
    alignment = AlignmentWrapper('R1', 0, 1001, 1100, 1, 1, 100)

    assert repr(alignment) == "AlignmentWrapper('R1', 0, 1001, 1100, 1, 1, 100)"


def test_wrap_overrides():
    alignment1 = AlignmentWrapper(r_st=100, r_en=200)
    alignment2 = AlignmentWrapper.wrap(alignment1, r_en=300, blen=200, cigar=[])
    expected_alignment = AlignmentWrapper(r_st=100, r_en=300, cigar=[])

    assert alignment2 == expected_alignment


def test_align(projects):
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

    aligner.align(seed_name, consensus)

    assert aligner.consensus == consensus
    assert_alignments(aligner, expected_alignment)
    assert aligner.algorithm == 'minimap2'


def test_align_multiple_sections(projects):
    seed_name = 'SARS-CoV-2-seed'
    seed_seq = projects.getReference(seed_name)
    consensus = seed_seq[3000:3500] + seed_seq[1000:2000] + seed_seq[5000:5500]
    expected_alignments = [AlignmentWrapper(ctg='N/A',
                                            ctg_len=len(seed_seq),
                                            r_st=3000,
                                            r_en=3500,
                                            q_st=0,
                                            q_en=500,
                                            mapq=60),
                           AlignmentWrapper(ctg='N/A',
                                            ctg_len=len(seed_seq),
                                            r_st=1000,
                                            r_en=2000,
                                            q_st=500,
                                            q_en=1500,
                                            mapq=60),
                           AlignmentWrapper(ctg='N/A',
                                            ctg_len=len(seed_seq),
                                            r_st=5000,
                                            r_en=5500,
                                            q_st=1500,
                                            q_en=2000,
                                            mapq=60)]
    aligner = ConsensusAligner(projects)

    aligner.align(seed_name, consensus)

    assert_alignments(aligner, *expected_alignments)


def test_align_short_consensus(projects):
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

    aligner.align(seed_name, consensus)

    assert_alignments(aligner, expected_alignment)
    assert aligner.algorithm == 'gotoh'


def test_align_deletion_minimap2(projects):
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

    aligner.align(seed_name, consensus)

    assert_alignments(aligner, expected_alignment)
    assert aligner.algorithm == 'minimap2'


def test_align_deletion_gotoh(projects):
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

    aligner.align(seed_name, consensus)

    assert_alignments(aligner, expected_alignment)
    assert aligner.algorithm == 'gotoh'


def test_align_insertion_minimap2(projects):
    seed_name = 'SARS-CoV-2-seed'
    seed_seq = projects.getReference(seed_name)
    consensus = seed_seq[2000:2030] + 'T' + seed_seq[2030:2060]
    expected_alignment = AlignmentWrapper(ctg='N/A',
                                          ctg_len=len(seed_seq),
                                          r_st=2000,
                                          r_en=2060,
                                          q_st=0,
                                          q_en=61,
                                          mapq=9,
                                          cigar=[[30, CigarActions.MATCH],
                                                 [1, CigarActions.INSERT],
                                                 [30, CigarActions.MATCH]],
                                          NM=1)
    aligner = ConsensusAligner(projects)

    aligner.align(seed_name, consensus)

    assert_alignments(aligner, expected_alignment)
    assert aligner.algorithm == 'minimap2'


def test_align_insertion_gotoh(projects):
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

    aligner.align(seed_name, consensus)

    assert_alignments(aligner, expected_alignment)
    assert aligner.algorithm == 'gotoh'


def test_align_with_only_primary_matches(projects):
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

    aligner.align(seed_name, consensus)

    assert_alignments(aligner, expected_alignment)
    assert aligner.algorithm == 'minimap2'


# noinspection DuplicatedCode
def test_align_seed_aminos(projects):
    seed_aminos = [SeedAmino(i) for i in range(3)]
    seed_aminos[0].count_aminos('AAA', 1)
    seed_aminos[1].count_aminos('CCC', 1)
    seed_aminos[2].count_aminos('GGG', 1)
    seed_name = 'HCV-6t'
    seed_seq = projects.getReference(seed_name)
    expected_consensus = 'AAACCCGGG'
    expected_alignment = AlignmentWrapper(ctg='N/A',
                                          ctg_len=len(seed_seq),
                                          r_st=4798,
                                          r_en=4807,
                                          q_st=0,
                                          q_en=9,
                                          mapq=0)
    aligner = ConsensusAligner(projects)

    aligner.align(seed_name, seed_aminos=seed_aminos)

    assert aligner.consensus == expected_consensus
    assert_alignments(aligner, expected_alignment)
    assert len(aligner.seed_nucs) == 9


def test_clear(projects):
    seed_name = 'SARS-CoV-2-seed'
    seed_seq = projects.getReference(seed_name)
    consensus = seed_seq[1000:2000]
    aligner = ConsensusAligner(projects)
    aligner.align(seed_name, consensus)

    aligner.clear()

    assert aligner.consensus == ''
    assert aligner.algorithm == ''
    assert len(aligner.alignments) == 0


# noinspection DuplicatedCode
def test_align_twice(projects):
    seed_aminos = [SeedAmino(i) for i in range(3)]
    seed_aminos[0].count_aminos('AAA', 1)
    seed_aminos[1].count_aminos('CCC', 1)
    seed_aminos[2].count_aminos('GGG', 1)
    seed_name = 'HCV-6t'
    seed_seq = projects.getReference(seed_name)
    consensus2 = seed_seq[:100]
    aligner = ConsensusAligner(projects)
    aligner.align(seed_name, seed_aminos=seed_aminos)

    aligner.align(seed_name, consensus2)

    assert aligner.consensus == consensus2
    assert not aligner.seed_nucs


def test_align_without_seed_name(projects):
    seed_aminos = [SeedAmino(i) for i in range(3)]
    seed_aminos[0].count_aminos('ATA', 1)
    seed_aminos[1].count_aminos('CGC', 1)
    seed_aminos[2].count_aminos('GTG', 1)
    expected_consensus = 'ATACGCGTG'
    seed_name = None
    aligner = ConsensusAligner(projects)
    aligner.align(seed_name, seed_aminos=seed_aminos)

    assert aligner.consensus == expected_consensus
    assert len(aligner.seed_nucs) == 9


def test_align_consensus_offset(projects):
    seed_aminos = [SeedAmino(i) for i in range(3)]
    seed_aminos[1].count_aminos('CCC', 1)
    seed_aminos[2].count_aminos('GGG', 1)
    seed_name = 'HCV-6t'
    expected_consensus = 'CCCGGG'
    aligner = ConsensusAligner(projects)

    aligner.align(seed_name, seed_aminos=seed_aminos)

    assert aligner.consensus == expected_consensus
    assert aligner.consensus_offset == 3
    assert len(aligner.seed_nucs) == 9
