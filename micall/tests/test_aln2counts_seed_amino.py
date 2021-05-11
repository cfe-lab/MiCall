from micall.core.aln2counts import SeedAmino


def test_single_read():
    """ Read a single codon, and report on counts.
    Columns are:       A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*
    """
    nuc_seq = 'AAA'  # -> K
    expected_counts = '0,0,0,0,0,0,0,0,8,0,0,0,0,0,0,0,0,0,0,0,0'
    amino = SeedAmino(None)

    amino.count_aminos(nuc_seq, 8)
    counts = amino.get_report()

    assert expected_counts == counts


def test_different_codon():
    """ Read two different codons, and report on counts.
    Columns are:       A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*
    """
    nuc_seq1 = 'AAA'  # -> K
    nuc_seq2 = 'GGG'  # -> G
    expected_counts = '0,0,0,0,0,5,0,0,8,0,0,0,0,0,0,0,0,0,0,0,0'
    amino = SeedAmino(None)

    amino.count_aminos(nuc_seq1, 8)
    amino.count_aminos(nuc_seq2, 5)
    counts = amino.get_report()

    assert expected_counts == counts


def test_same_amino_acid():
    """ Read same codon twice, and report on counts.
    Columns are:       A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*
    """
    nuc_seq1 = 'AAA'  # -> K
    nuc_seq2 = 'AAG'  # -> K
    expected_counts = '0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0'
    amino = SeedAmino(None)

    amino.count_aminos(nuc_seq1, 4)
    amino.count_aminos(nuc_seq2, 5)
    counts = amino.get_report()

    assert expected_counts == counts


def test_nucleotides():
    nuc_seq1 = 'AAA'  # -> K
    nuc_seq2 = 'AAG'  # -> K
    expected_nuc_counts = '4,0,5,0'
    amino = SeedAmino(None)

    amino.count_aminos(nuc_seq1, 4)
    amino.count_aminos(nuc_seq2, 5)
    counts = amino.nucleotides[2].get_report()

    assert expected_nuc_counts == counts

  
def test_consensus():
    nuc_seq1 = 'AAA'  # -> K
    nuc_seq2 = 'GGG'  # -> G
    expected_consensus = 'G'
    amino = SeedAmino(None)

    amino.count_aminos(nuc_seq1, 4)
    amino.count_aminos(nuc_seq2, 5)
    consensus = amino.get_consensus()

    assert expected_consensus == consensus


def test_consensus_mixture():
    nuc_seq1 = 'AAA'  # -> K
    nuc_seq2 = 'GGG'  # -> G
    nuc_seq3 = 'TTT'  # -> F
    allowed_consensus_values = ('G', 'K')
    amino = SeedAmino(None)

    amino.count_aminos(nuc_seq1, 4)
    amino.count_aminos(nuc_seq2, 4)
    amino.count_aminos(nuc_seq3, 3)
    consensus = amino.get_consensus()

    assert consensus in allowed_consensus_values


def test_consensus_partial():
    nuc_seq1 = 'AAA'  # -> K
    nuc_seq2 = '-GG'  # -> ?
    expected_consensus = 'K'  # Complete reads always override partials.
    amino = SeedAmino(None)

    amino.count_aminos(nuc_seq1, 4)
    amino.count_aminos(nuc_seq2, 10)
    consensus = amino.get_consensus()

    assert consensus == expected_consensus


def test_consensus_with_no_reads():
    amino = SeedAmino(None)

    consensus = amino.get_consensus()

    assert consensus == '-'


def test_consensus_with_only_partials():
    nuc_seq = '-GG'  # -> ?
    expected_consensus = '?'  # No complete reads, only partials.
    amino = SeedAmino(None)

    amino.count_aminos(nuc_seq, 4)
    consensus = amino.get_consensus()

    assert consensus == expected_consensus


def test_missing_data():
    """ Lower-case n represents a gap between the forward and reverse reads. """

    nuc_seq = 'CTn'
    expected_consensus = '?'
    amino = SeedAmino(None)

    amino.count_aminos(nuc_seq, 1)
    consensus = amino.get_consensus()

    assert expected_consensus == consensus


def test_ambiguous_data():
    """If a read is ambiguous, don't count it toward consensus."""

    nuc_seq1 = 'Cnn'  # -> ?
    nuc_seq2 = 'AAA'  # -> K
    expected_consensus = 'K'
    amino = SeedAmino(None)

    amino.count_aminos(nuc_seq1, 9)
    amino.count_aminos(nuc_seq2, 1)
    consensus = amino.get_consensus()

    assert expected_consensus == consensus


def test_overlap():
    amino = SeedAmino(None)
    amino.count_aminos('GGG', 4)
    other = SeedAmino(consensus_nuc_index=7)
    other.count_aminos('TAG', 5)
    expected_counts = {'G': 4}
    expected_v3_overlap = 5

    amino.count_overlap(other)

    assert expected_counts == amino.counts
    assert expected_v3_overlap == amino.v3_overlap
    assert expected_v3_overlap == amino.nucleotides[0].v3_overlap


def test_overlap_partial_codon():
    amino = SeedAmino(None)
    amino.count_aminos('GGG', 4)
    other = SeedAmino(consensus_nuc_index=7)
    other.count_aminos('TA', 5)
    expected_counts = {'G': 4}
    expected_v3_overlap = 5

    amino.count_overlap(other)

    assert expected_counts == amino.counts
    assert expected_v3_overlap == amino.v3_overlap
    assert expected_v3_overlap == amino.nucleotides[0].v3_overlap


def test_add():
    amino = SeedAmino(None)
    amino.count_aminos('GGG', 4)  # => G
    other = SeedAmino(consensus_nuc_index=7)
    other.count_aminos('AAA', 5)  # => K
    expected_counts = {'G': 4, 'K': 5}
    expected_nucleotide_counts = {'G': 4, 'A': 5}

    amino.add(other)

    assert amino.counts == expected_counts
    assert amino.nucleotides[0].counts == expected_nucleotide_counts
    assert amino.consensus_nuc_index is None


def test_add_first_time():
    amino = SeedAmino(None)
    other = SeedAmino(consensus_nuc_index=7)
    other.count_aminos('AAA', 5)  # => K
    expected_counts = {'K': 5}

    amino.add(other)

    assert amino.counts == expected_counts
    assert amino.consensus_nuc_index == 7


def test_amino_repeat_nuc0():
    """ Repeat first nucleotide in the codon.
    Columns are:       A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*
    """
    nuc_seq1 = 'GAG'  # -> E
    nuc_seq2 = 'GAC'  # -> D
    # repeat = 'GGA'  # -> G
    expected_counts = '0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0'
    amino = SeedAmino(None)

    amino.count_aminos(nuc_seq1, 4)
    amino.count_aminos(nuc_seq2, 5)
    counts = amino.apply_repeat(0).get_report()

    assert counts == expected_counts


def test_amino_repeat_nuc1():
    """ Repeat second nucleotide in the codon.
    Columns are:       A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*
    """
    nuc_seq1 = 'GAG'  # -> E
    nuc_seq2 = 'GAC'  # -> D
    # repeat = 'GAA'  # -> E
    expected_counts = '0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0'
    amino = SeedAmino(None)

    amino.count_aminos(nuc_seq1, 4)
    amino.count_aminos(nuc_seq2, 5)
    counts = amino.apply_repeat(1).get_report()

    assert counts == expected_counts


def test_amino_repeat_nuc2():
    """ Repeat last nucleotide in the codon.
    Columns are:       A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*
    """
    nuc_seq1 = 'GAG'  # -> E
    nuc_seq2 = 'GAC'  # -> D
    expected_counts = '0,0,5,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0'
    amino = SeedAmino(None)

    amino.count_aminos(nuc_seq1, 4)
    amino.count_aminos(nuc_seq2, 5)
    counts = amino.apply_repeat(2).get_report()

    assert counts == expected_counts
