from micall.core.aln2counts import SeedNucleotide, MAX_CUTOFF, FIRST_CUTOFF


def test_single_read():
    """ Read a single nucleotide, and report on counts.
    Columns are:       A,C,G,T
    """
    nuc_seq = 'C'
    expected_counts = '0,8,0,0'
    nuc = SeedNucleotide()
    
    nuc.count_nucleotides(nuc_seq, 8)
    counts = nuc.get_report()

    assert expected_counts == counts


def test_consensus_no_mixes():
    nuc = SeedNucleotide()
    nuc.count_nucleotides('C', 1)
    consensus_max = nuc.get_consensus(MAX_CUTOFF)
    consensus_mix = nuc.get_consensus(0.1)

    expected_consensus = 'C'
    assert expected_consensus == consensus_max
    assert expected_consensus == consensus_mix


def test_consensus_mixed():
    nuc = SeedNucleotide()
    nuc.count_nucleotides('C', 2)
    nuc.count_nucleotides('T', 1)
    consensus_max = nuc.get_consensus(MAX_CUTOFF)
    consensus_mix = nuc.get_consensus(0.1)

    expected_consensus_max = 'C'
    expected_consensus_mix = 'Y'
    assert expected_consensus_max == consensus_max
    assert expected_consensus_mix == consensus_mix


def test_consensus_mixed_three():
    nuc = SeedNucleotide()
    nuc.count_nucleotides('C', 2)
    nuc.count_nucleotides('T', 1)
    nuc.count_nucleotides('G', 1)
    consensus_max = nuc.get_consensus(MAX_CUTOFF)
    consensus_mix = nuc.get_consensus(0.1)

    expected_consensus_max = 'C'
    expected_consensus_mix = 'B'  # B is a mix of T, G, and C
    assert expected_consensus_max == consensus_max
    assert expected_consensus_mix == consensus_mix


def test_consensus_mixed_all():
    nuc = SeedNucleotide()
    nuc.count_nucleotides('C', 2)
    nuc.count_nucleotides('T', 1)
    nuc.count_nucleotides('G', 1)
    nuc.count_nucleotides('A', 1)
    consensus_max = nuc.get_consensus(MAX_CUTOFF)
    consensus_mix = nuc.get_consensus(0.1)

    expected_consensus_max = 'C'
    expected_consensus_mix = 'N'  # All four are reported as N
    assert expected_consensus_max == consensus_max
    assert expected_consensus_mix == consensus_mix


def test_consensus_mixed_max():
    nuc = SeedNucleotide()
    nuc.count_nucleotides('C', 2)
    nuc.count_nucleotides('T', 2)
    nuc.count_nucleotides('G', 1)
    consensus_first = nuc.get_consensus(FIRST_CUTOFF)
    consensus_max = nuc.get_consensus(MAX_CUTOFF)
    consensus_mix = nuc.get_consensus(0.1)

    expected_consensus_first = 'C'  # C and T tie for max, alphabetical first=C
    expected_consensus_max = 'Y'  # C and T tie for max, mix is Y
    expected_consensus_mix = 'B'  # C, T, and G mix is B
    assert expected_consensus_first == consensus_first
    assert expected_consensus_max == consensus_max
    assert expected_consensus_mix == consensus_mix


def test_consensus_cutoff():
    nuc = SeedNucleotide()
    nuc.count_nucleotides('C', 2)
    nuc.count_nucleotides('T', 1)
    consensus_mix = nuc.get_consensus(0.5)

    expected_consensus = 'C'  # T was below the cutoff
    assert expected_consensus == consensus_mix


def test_consensus_cutoff_at_boundary():
    nuc = SeedNucleotide()
    nuc.count_nucleotides('C', 9000)
    nuc.count_nucleotides('T', 1000)
    consensus_mix = nuc.get_consensus(0.1)

    expected_consensus = 'Y'  # T was at the cutoff
    assert expected_consensus == consensus_mix


def test_consensus_cutoff_below_boundary():
    nuc = SeedNucleotide()
    nuc.count_nucleotides('C', 9001)
    nuc.count_nucleotides('T', 999)
    consensus_mix = nuc.get_consensus(0.1)

    expected_consensus = 'C'  # T was below the cutoff
    assert expected_consensus == consensus_mix


def test_consensus_mixed_with_poor_quality():
    nuc = SeedNucleotide()
    nuc.count_nucleotides('N', 99)
    nuc.count_nucleotides('T', 1)
    consensus_max = nuc.get_consensus(MAX_CUTOFF)
    consensus_mix_one_pct = nuc.get_consensus(0.01)
    consensus_mix_ten_pct = nuc.get_consensus(0.10)

    expected_consensus_max = 'T'  # N always overruled
    expected_consensus_mix_one_pct = 'T'
    expected_consensus_mix_ten_pct = 'T'
    assert expected_consensus_max == consensus_max
    assert expected_consensus_mix_one_pct == consensus_mix_one_pct
    assert expected_consensus_mix_ten_pct == consensus_mix_ten_pct


def test_consensus_mixed_with_gap():
    nuc = SeedNucleotide()
    nuc.count_nucleotides('-', 99)
    nuc.count_nucleotides('T', 1)
    consensus_max = nuc.get_consensus(MAX_CUTOFF)
    consensus_mix_one_pct = nuc.get_consensus(0.01)
    consensus_mix_ten_pct = nuc.get_consensus(0.10)

    expected_consensus_max = '-'  # most common
    expected_consensus_mix_one_pct = 't'  # mix of both
    expected_consensus_mix_ten_pct = '-'  # only deletions
    assert expected_consensus_max == consensus_max
    assert expected_consensus_mix_one_pct == consensus_mix_one_pct
    assert expected_consensus_mix_ten_pct == consensus_mix_ten_pct


def test_consensus_mixed_with_gap_and_poor_quality():
    nuc = SeedNucleotide()
    nuc.count_nucleotides('N', 3)
    nuc.count_nucleotides('-', 2)
    nuc.count_nucleotides('T', 1)
    consensus_max = nuc.get_consensus(MAX_CUTOFF)
    consensus_mix = nuc.get_consensus(0.1)

    expected_consensus_max = '-'
    expected_consensus_mix = 't'
    assert expected_consensus_max == consensus_max
    assert expected_consensus_mix == consensus_mix


def test_consensus_poor_quality_only():
    nuc = SeedNucleotide()
    nuc.count_nucleotides('N', 1)
    consensus_max = nuc.get_consensus(MAX_CUTOFF)
    consensus_mix = nuc.get_consensus(0.1)

    expected_consensus_max = 'N'
    expected_consensus_mix = 'N'
    assert expected_consensus_max == consensus_max
    assert expected_consensus_mix == consensus_mix


def test_consensus_mixed_gap_and_poor_quality_only():
    nuc = SeedNucleotide()
    nuc.count_nucleotides('N', 3)
    nuc.count_nucleotides('-', 2)
    consensus_max = nuc.get_consensus(MAX_CUTOFF)
    consensus_mix = nuc.get_consensus(0.1)

    expected_consensus_max = '-'
    expected_consensus_mix = '-'
    assert expected_consensus_max == consensus_max
    assert expected_consensus_mix == consensus_mix


def test_consensus_all_below_cutoff():
    nuc = SeedNucleotide()
    nuc.count_nucleotides('C', 101)
    nuc.count_nucleotides('T', 100)
    nuc.count_nucleotides('G', 99)
    consensus_max = nuc.get_consensus(MAX_CUTOFF)
    consensus_mix = nuc.get_consensus(0.5)

    expected_consensus_max = 'C'
    expected_consensus_mix = 'N'
    assert expected_consensus_max == consensus_max
    assert expected_consensus_mix == consensus_mix


def test_consensus_between_reads():
    """Lower-case n represents the gap between forward and reverse reads.

    Should not be counted in consensus totals"""
    nuc = SeedNucleotide()
    nuc.count_nucleotides('C', 9)
    nuc.count_nucleotides('T', 1)
    nuc.count_nucleotides('n', 2)
    consensus_mix = nuc.get_consensus(0.1)

    expected_consensus = 'Y'
    assert expected_consensus == consensus_mix


def test_consensus_missing_positions():
    """ Positions that are never read are ignored in the consensus. """

    # No counts added
    nuc = SeedNucleotide()

    consensus_max = nuc.get_consensus(MAX_CUTOFF)
    consensus_mix = nuc.get_consensus(0.1)

    expected_consensus = ''
    assert expected_consensus == consensus_max
    assert expected_consensus == consensus_mix


def test_overlap():
    nuc = SeedNucleotide()
    nuc.count_nucleotides('T', 4)
    other = SeedNucleotide()
    other.count_nucleotides('C', 5)
    expected_counts = {'T': 4}
    expected_v3_overlap = 5

    nuc.count_overlap(other)

    assert expected_counts == nuc.counts
    assert expected_v3_overlap == nuc.v3_overlap


def test_add():
    nuc = SeedNucleotide()
    nuc.count_nucleotides('T', 4)
    nuc.count_nucleotides('C', 1)
    nuc.clip_count = 10
    nuc.insertion_count = 7

    other = SeedNucleotide()
    other.count_nucleotides('C', 5)
    other.clip_count = 9
    other.insertion_count = 8
    expected_counts = {'T': 4, 'C': 6}
    expected_clip_count = 19
    expected_insertion_count = 15

    nuc.add(other)

    assert expected_counts == nuc.counts
    assert expected_clip_count == nuc.clip_count
    assert expected_insertion_count == nuc.insertion_count
