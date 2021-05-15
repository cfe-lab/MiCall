# noinspection PyUnresolvedReferences
from random import randrange, choice

from micall.utils.consensus_aligner import ConsensusAligner, AlignmentWrapper, CigarActions


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
    assert len(aligner.alignments) == 1
    alignment, = aligner.alignments
    alignment = AlignmentWrapper.wrap(alignment)
    assert alignment == expected_alignment
    assert aligner.algorithm == 'minimap2'


def test_align_multiple_sections(projects):
    seed_name = 'SARS-CoV-2-seed'
    seed_seq = projects.getReference(seed_name)
    consensus = seed_seq[3000:4000] + seed_seq[1000:2000]
    expected_alignments = [AlignmentWrapper(ctg='N/A',
                                            ctg_len=len(seed_seq),
                                            r_st=3000,
                                            r_en=4000,
                                            q_st=0,
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

    aligner.align(seed_name, consensus)

    assert aligner.consensus == consensus
    alignments = [AlignmentWrapper.wrap(alignment)
                  for alignment in aligner.alignments]
    assert alignments == expected_alignments


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

    assert len(aligner.alignments) == 1
    alignment, = aligner.alignments
    alignment = AlignmentWrapper.wrap(alignment, mapq=0)
    assert alignment == expected_alignment
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

    assert len(aligner.alignments) == 1
    alignment, = aligner.alignments
    assert alignment == expected_alignment
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

    assert len(aligner.alignments) == 1
    alignment, = aligner.alignments
    assert alignment == expected_alignment
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

    assert len(aligner.alignments) == 1
    alignment, = aligner.alignments
    alignment = AlignmentWrapper.wrap(alignment)
    assert alignment == expected_alignment
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

    assert len(aligner.alignments) == 1
    alignment, = aligner.alignments
    assert alignment == expected_alignment
    assert aligner.algorithm == 'gotoh'
