from csv import DictReader
from StringIO import StringIO

from micall.utils.find_chimera import find_consensus, map_sequence


def test_map_sequence_subsequence():
    source_seq = 'ACTGATGC'
    dest_seq = 'TGAT'
    expected_positions = [3, 4, 5, 6]

    positions = map_sequence(source_seq, dest_seq)

    assert expected_positions == positions


def test_map_sequence_identical():
    source_seq = 'ACTG'
    dest_seq = 'ACTG'
    expected_positions = [1, 2, 3, 4]

    positions = map_sequence(source_seq, dest_seq)

    assert expected_positions == positions


def test_map_sequence_substitution():
    source_seq = 'ACTGATGC'
    dest_seq = 'TCAT'
    expected_positions = [3, 4, 5, 6]

    positions = map_sequence(source_seq, dest_seq)

    assert expected_positions == positions


def test_map_sequence_deletion():
    source_seq = 'ACTGATGC'
    dest_seq = 'ACTATGC'
    expected_positions = [1, 2, 3, 5, 6, 7, 8]

    positions = map_sequence(source_seq, dest_seq)

    assert expected_positions == positions


def test_map_sequence_insertion():
    source_seq = 'ACTGATGC'
    dest_seq = 'ACTGCATGC'
    expected_positions = [1, 2, 3, 4, None, 5, 6, 7, 8]

    positions = map_sequence(source_seq, dest_seq)

    assert expected_positions == positions


def test_map_sequence_low_quality():
    source_seq = 'AGAGCGAACCGATTC'
    dest_seq = 'NNNNNNNNNNNATTC'
    expected_positions = 11 * [None] + [12, 13, 14, 15]

    positions = map_sequence(source_seq, dest_seq)

    assert expected_positions == positions


def test_find_consensus():
    nuc_csv = StringIO("""\
seed,region,q-cutoff,query.nuc.pos,refseq.nuc.pos,A,C,G,T
R1-seed,R1,15,1,1,9,0,0,0
R1-seed,R1,15,2,2,0,9,0,0
R1-seed,R1,15,3,3,0,0,9,0
R1-seed,R1,15,4,4,0,0,0,9
""")
    expected_consensus = 'ACGT'
    rows = list(DictReader(nuc_csv))
    consensus = find_consensus(rows, 'R1-seed', 'R1')

    assert consensus == expected_consensus


def test_find_consensus_with_mixture():
    nuc_csv = StringIO("""\
seed,region,q-cutoff,query.nuc.pos,refseq.nuc.pos,A,C,G,T
R1-seed,R1,15,1,1,6,0,3,0
R1-seed,R1,15,2,2,0,9,0,0
R1-seed,R1,15,3,3,0,0,9,0
R1-seed,R1,15,4,4,0,0,0,9
""")
    expected_consensus = 'ACGT'
    rows = list(DictReader(nuc_csv))
    consensus = find_consensus(rows, 'R1-seed', 'R1')

    assert consensus == expected_consensus


def test_find_consensus_with_gap():
    nuc_csv = StringIO("""\
seed,region,q-cutoff,query.nuc.pos,refseq.nuc.pos,A,C,G,T
R1-seed,R1,15,1,1,9,0,0,0
R1-seed,R1,15,2,2,0,9,0,0
R1-seed,R1,15,3,3,0,0,0,0
R1-seed,R1,15,4,4,0,0,0,9
""")
    expected_consensus = 'ACNT'
    rows = list(DictReader(nuc_csv))
    consensus = find_consensus(rows, 'R1-seed', 'R1')

    assert consensus == expected_consensus


def test_find_consensus_with_offset():
    nuc_csv = StringIO("""\
seed,region,q-cutoff,query.nuc.pos,refseq.nuc.pos,A,C,G,T
R1-seed,R1,15,1,11,9,0,0,0
R1-seed,R1,15,2,12,0,9,0,0
R1-seed,R1,15,3,13,0,0,9,0
R1-seed,R1,15,4,14,0,0,0,9
""")
    expected_consensus = 'NNNNNNNNNNACGT'
    rows = list(DictReader(nuc_csv))
    consensus = find_consensus(rows, 'R1-seed', 'R1')

    assert consensus == expected_consensus


def test_find_consensus_one_region():
    nuc_csv = StringIO("""\
seed,region,q-cutoff,query.nuc.pos,refseq.nuc.pos,A,C,G,T
R1-seed,R1,15,1,1,9,0,0,0
R1-seed,R1,15,2,2,0,9,0,0
R1-seed,R1,15,3,3,0,0,9,0
R1-seed,R1,15,4,4,0,0,0,9
R1-seed,R2,15,100,1,0,0,0,9
""")
    expected_consensus = 'ACGT'
    rows = list(DictReader(nuc_csv))
    consensus = find_consensus(rows, 'R1-seed', 'R1')

    assert consensus == expected_consensus


def test_find_consensus_one_seed():
    nuc_csv = StringIO("""\
seed,region,q-cutoff,query.nuc.pos,refseq.nuc.pos,A,C,G,T
R1-seed,R1,15,1,1,9,0,0,0
R1-seed,R1,15,2,2,0,9,0,0
R1-seed,R1,15,3,3,0,0,9,0
R1-seed,R1,15,4,4,0,0,0,9
R2-seed,R1,15,1,1,0,0,0,9
""")
    expected_consensus = 'ACGT'
    rows = list(DictReader(nuc_csv))
    consensus = find_consensus(rows, 'R1-seed', 'R1')

    assert consensus == expected_consensus
