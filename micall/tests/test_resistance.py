import typing
from csv import DictReader, DictWriter
from io import StringIO
from unittest import TestCase

import pytest

from micall.resistance.resistance import read_aminos, write_resistance, \
    select_reported_regions, AminoList, filter_aminos, load_asi, get_genotype, \
    combine_aminos, write_consensus, create_consensus_writer, write_failures, \
    write_nuc_mutations, HIVDB_VERSION


@pytest.fixture(scope="module")
def asi_algorithms():
    return load_asi()


def assert_aminos_lists_match(actual_lists: typing.List[AminoList],
                              expected_lists: typing.List[AminoList]):
    if expected_lists and actual_lists and (len(expected_lists) == len(actual_lists)):
        for i, (expected, actual) in enumerate(zip(expected_lists, actual_lists)):
            assert actual.aminos == expected.aminos, i
    assert actual_lists == expected_lists


class SelectReportedRegionsTest(TestCase):
    def test_all_regions(self):
        choices = ['PR', 'RT', 'IN']
        original_regions = {'PR', 'RT', 'IN'}
        expected_regions = original_regions.copy()

        selected_regions = select_reported_regions(choices, original_regions)

        self.assertEqual(expected_regions, selected_regions)

    def test_some_regions(self):
        choices = ['PR', 'IN']
        original_regions = {'PR', 'RT', 'IN'}
        expected_regions = {'PR', 'IN'}

        selected_regions = select_reported_regions(choices, original_regions)

        self.assertEqual(expected_regions, selected_regions)

    def test_combined_regions(self):
        choices = ['PR_RT']
        original_regions = {'PR', 'RT', 'IN'}
        expected_regions = {'PR', 'RT'}

        selected_regions = select_reported_regions(choices, original_regions)

        self.assertEqual(expected_regions, selected_regions)


def format_rows(column_names, rows):
    display = StringIO()
    writer = DictWriter(display, column_names)
    writer.writeheader()
    writer.writerows(rows)
    return display.getvalue()


def check_combination(
        amino_csv,
        midi_amino_csv,
        expected_csv,
        expected_failures=None):
    """ Combine two amino count CSV files, and check result.

    :param amino_csv: open file with main amino counts
    :param midi_amino_csv: open file with MIDI amino counts, or the same
        file as amino_csv, in which case it will be ignored.
    :param expected_csv: text expected to find in combined CSV
    :param expected_failures: None for no failures, or {(region, is_midi): message}
    """
    if expected_failures is None:
        expected_failures = {}
    failures = {}
    combined_rows = list(combine_aminos(amino_csv,
                                        midi_amino_csv,
                                        failures))

    assert_csv_rows(combined_rows, expected_csv)
    assert failures == expected_failures


def assert_csv_rows(rows: typing.Iterable[dict], expected_csv: str):
    expected_reader = DictReader(StringIO(expected_csv))
    expected_rows = list(expected_reader)
    assert (format_rows(expected_reader.fieldnames, rows) ==
            format_rows(expected_reader.fieldnames, expected_rows))


def test_combine_aminos_simple():
    amino_csv = StringIO("""\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,g2p_overlap,coverage
HCV-1a,HCV1A-H77-NS3,15,1,1,0,0,0,0,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS3,15,7,181,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5a,15,1,1,0,0,0,0,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5a,15,7,101,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,1,1,0,0,0,0,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,4,2,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,7,228,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
""")
    amino_csv.name = 'main_amino.csv'
    expected_csv = amino_csv.getvalue()  # No change.

    check_combination(amino_csv, amino_csv, expected_csv)


def test_combine_aminos_low_average():
    amino_csv = StringIO("""\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,g2p_overlap,coverage
HCV-1a,HCV1A-H77-NS3,15,1,1,0,0,0,0,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS3,15,7,181,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5a,15,1,1,0,0,0,0,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5a,15,7,101,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,1,1,0,0,0,0,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,4,2,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
""")
    amino_csv.name = 'main_amino.csv'
    expected_csv = """\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,g2p_overlap,coverage
HCV-1a,HCV1A-H77-NS3,15,1,1,0,0,0,0,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS3,15,7,181,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5a,15,1,1,0,0,0,0,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5a,15,7,101,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
"""
    expected_failures = {('HCV-1a', 'HCV1A-H77-NS5b', False): 'low average coverage'}

    check_combination(amino_csv, amino_csv, expected_csv, expected_failures)


def test_combine_aminos_ns5b_window():
    amino_csv = StringIO("""\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,g2p_overlap,coverage
HCV-1a,HCV1A-H77-NS3,15,1,1,0,0,0,0,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS3,15,7,181,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5a,15,1,1,0,0,0,0,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5a,15,7,101,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,1,1,0,0,0,0,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,4,2,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,7,229,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
""")
    amino_csv.name = 'main_amino.csv'
    expected_csv = """\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,g2p_overlap,coverage
HCV-1a,HCV1A-H77-NS3,15,1,1,0,0,0,0,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS3,15,7,181,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5a,15,1,1,0,0,0,0,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5a,15,7,101,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
"""
    expected_failures = {('HCV-1a', 'HCV1A-H77-NS5b', False): 'low average coverage'}

    check_combination(amino_csv, amino_csv, expected_csv, expected_failures)


def test_combine_aminos_ns3_and_ns5a_windows():
    amino_csv = StringIO("""\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,g2p_overlap,coverage
HCV-1a,HCV1A-H77-NS3,15,1,1,0,0,0,0,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS3,15,7,182,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5a,15,1,1,0,0,0,0,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5a,15,7,102,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,1,1,0,0,0,0,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,4,2,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,7,228,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
""")
    amino_csv.name = 'main_amino.csv'
    expected_csv = """\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,g2p_overlap,coverage
HCV-1a,HCV1A-H77-NS5b,15,1,1,0,0,0,0,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,4,2,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,7,228,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
"""
    expected_failures = {('HCV-1a', 'HCV1A-H77-NS3', False): 'low average coverage',
                         ('HCV-1a', 'HCV1A-H77-NS5a', False): 'low average coverage'}

    check_combination(amino_csv, amino_csv, expected_csv, expected_failures)


def test_combine_aminos_ns3_and_ns5_missing_tails():
    amino_csv = StringIO("""\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,g2p_overlap,coverage
HCV-1a,HCV1A-H77-NS3,15,1,1,0,0,0,0,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS3,15,7,180,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5a,15,1,1,0,0,0,0,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5a,15,7,100,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,1,1,0,0,0,0,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,4,2,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,7,228,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
""")
    amino_csv.name = 'main_amino.csv'
    expected_csv = """\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,g2p_overlap,coverage
HCV-1a,HCV1A-H77-NS5b,15,1,1,0,0,0,0,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,4,2,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,7,228,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
"""
    expected_failures = {
        ('HCV-1a', 'HCV1A-H77-NS3', False): 'not enough high-coverage amino acids',
        ('HCV-1a', 'HCV1A-H77-NS5a', False): 'not enough high-coverage amino acids'}

    check_combination(amino_csv, amino_csv, expected_csv, expected_failures)


def test_combine_aminos_ns3_and_ns5a_missing_heads():
    amino_csv = StringIO("""\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,g2p_overlap,coverage
HCV-1a,HCV1A-H77-NS3,15,1,2,0,0,0,0,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS3,15,7,181,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5a,15,1,2,0,0,0,0,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5a,15,7,101,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,1,1,0,0,0,0,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,4,2,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,7,228,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
""")
    amino_csv.name = 'main_amino.csv'
    expected_csv = """\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,g2p_overlap,coverage
HCV-1a,HCV1A-H77-NS5b,15,1,1,0,0,0,0,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,4,2,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,7,228,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
"""
    expected_failures = {
        ('HCV-1a', 'HCV1A-H77-NS3', False): 'not enough high-coverage amino acids',
        ('HCV-1a', 'HCV1A-H77-NS5a', False): 'not enough high-coverage amino acids'}

    check_combination(amino_csv, amino_csv, expected_csv, expected_failures)


def test_combine_aminos_ns5b():
    amino_csv = StringIO("""\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,g2p_overlap,coverage
HCV-1a,HCV1A-H77-NS5b,15,1,1,0,0,0,0,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,4,2,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,7,228,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
""")
    amino_csv.name = 'main_amino.csv'
    midi_amino_csv = StringIO("""\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,g2p_overlap,coverage
HCV-1a,HCV1A-H77-NS5b,15,1,558,0,0,0,0,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,4,559,0,0,0,0,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,7,560,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,10,561,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
""")
    midi_amino_csv.name = 'midi_amino.csv'
    expected_csv = """\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,g2p_overlap,coverage
HCV-1a,HCV1A-H77-NS5b,15,1,1,0,0,0,0,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,4,2,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,7,228,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,1,558,0,0,0,0,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,4,559,0,0,0,0,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,7,560,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,10,561,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
"""

    check_combination(amino_csv, midi_amino_csv, expected_csv)


def test_combine_aminos_ns5b_multiple_seeds():
    amino_csv = StringIO("""\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,g2p_overlap,coverage
HCV-1a,HCV1A-H77-NS5b,15,1,1,0,0,0,0,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,4,2,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,7,228,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
""")
    amino_csv.name = 'main_amino.csv'
    midi_amino_csv = StringIO("""\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,g2p_overlap,coverage
HCV-1a,HCV1A-H77-NS5b,15,1,558,0,0,0,0,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,1,559,0,0,0,0,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,4,560,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,7,561,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1b,HCV1B-Con1-NS5b,15,1,558,0,0,0,0,0,0,0,0,20000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,20000
HCV-1b,HCV1B-Con1-NS5b,15,1,559,0,0,0,0,0,0,0,0,20000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,20000
HCV-1b,HCV1B-Con1-NS5b,15,4,560,0,0,0,0,20000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,20000
HCV-1b,HCV1B-Con1-NS5b,15,7,561,20000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,20000
""")
    midi_amino_csv.name = 'midi_amino.csv'
    expected_csv = """\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,g2p_overlap,coverage
HCV-1a,HCV1A-H77-NS5b,15,1,1,0,0,0,0,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,4,2,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,7,228,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,1,558,0,0,0,0,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,1,559,0,0,0,0,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,4,560,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,7,561,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1b,HCV1B-Con1-NS5b,15,1,558,0,0,0,0,0,0,0,0,20000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,20000
HCV-1b,HCV1B-Con1-NS5b,15,1,559,0,0,0,0,0,0,0,0,20000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,20000
HCV-1b,HCV1B-Con1-NS5b,15,4,560,0,0,0,0,20000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,20000
HCV-1b,HCV1B-Con1-NS5b,15,7,561,20000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,20000
"""
    expected_failures = {('HCV-1b', 'HCV1B-Con1-NS5b', False): 'nothing mapped'}

    check_combination(amino_csv, midi_amino_csv, expected_csv, expected_failures)


def test_combine_aminos_ns5b_multiple_subtypes():
    """ Main sample's subtype overrides MIDI subtype when they disagree. """
    amino_csv = StringIO("""\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,g2p_overlap,coverage
HCV-4c,HCV4-ED43-NS5b,15,1,1,0,0,0,0,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-4c,HCV4-ED43-NS5b,15,4,2,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-4c,HCV4-ED43-NS5b,15,7,228,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
""")
    amino_csv.name = 'main_amino.csv'
    midi_amino_csv = StringIO("""\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,g2p_overlap,coverage
HCV-4m,HCV4-ED43-NS5b,15,1,558,0,0,0,0,0,0,0,0,20000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,20000
HCV-4m,HCV4-ED43-NS5b,15,1,559,0,0,0,0,0,0,0,0,20000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,20000
HCV-4m,HCV4-ED43-NS5b,15,4,560,0,0,0,0,20000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,20000
HCV-4m,HCV4-ED43-NS5b,15,7,561,20000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,20000
""")
    midi_amino_csv.name = 'midi_amino.csv'
    expected_csv = """\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,g2p_overlap,coverage
HCV-4c,HCV4-ED43-NS5b,15,1,1,0,0,0,0,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-4c,HCV4-ED43-NS5b,15,4,2,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-4c,HCV4-ED43-NS5b,15,7,228,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-4c,HCV4-ED43-NS5b,15,1,558,0,0,0,0,0,0,0,0,20000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,20000
HCV-4c,HCV4-ED43-NS5b,15,1,559,0,0,0,0,0,0,0,0,20000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,20000
HCV-4c,HCV4-ED43-NS5b,15,4,560,0,0,0,0,20000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,20000
HCV-4c,HCV4-ED43-NS5b,15,7,561,20000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,20000
"""

    check_combination(amino_csv, midi_amino_csv, expected_csv)


def test_combine_aminos_low_coverage_in_midi():
    amino_csv = StringIO("""\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,g2p_overlap,coverage
HCV-1a,HCV1A-H77-NS5b,15,1,1,0,0,0,0,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,4,2,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,7,228,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
""")
    amino_csv.name = 'main_amino.csv'
    midi_amino_csv = StringIO("""\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,g2p_overlap,coverage
HCV-1a,HCV1A-H77-NS5b,15,4,230,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,4,559,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,4,560,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,7,561,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
""")
    midi_amino_csv.name = 'midi_amino.csv'
    expected_csv = """\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,g2p_overlap,coverage
HCV-1a,HCV1A-H77-NS5b,15,1,1,0,0,0,0,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,4,2,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,7,228,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
"""
    expected_failures = {('HCV-1a', 'HCV1A-H77-NS5b', True): 'low average coverage'}

    check_combination(amino_csv, midi_amino_csv, expected_csv, expected_failures)


def test_combine_aminos_no_coverage_in_midi():
    amino_csv = StringIO("""\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,g2p_overlap,coverage
HCV-1a,HCV1A-H77-NS5b,15,1,1,0,0,0,0,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,4,2,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,7,228,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
""")
    amino_csv.name = 'main_amino.csv'
    midi_amino_csv = StringIO("""\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,g2p_overlap,coverage
""")
    midi_amino_csv.name = 'midi_amino.csv'
    expected_csv = """\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,g2p_overlap,coverage
HCV-1a,HCV1A-H77-NS5b,15,1,1,0,0,0,0,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,4,2,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,7,228,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
"""
    expected_failures = {('HCV-1a', 'HCV1A-H77-NS5b', True): 'nothing mapped'}

    check_combination(amino_csv, midi_amino_csv, expected_csv, expected_failures)


def test_combine_aminos_no_coverage_in_ns5b():
    amino_csv = StringIO("""\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,g2p_overlap,coverage
""")
    amino_csv.name = 'main_amino.csv'
    midi_amino_csv = StringIO("""\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,g2p_overlap,coverage
HCV-1a,HCV1A-H77-NS5b,15,1,558,0,0,0,0,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,1,559,0,0,0,0,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,4,560,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,7,561,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
""")
    midi_amino_csv.name = 'midi_amino.csv'
    expected_csv = """\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,g2p_overlap,coverage
HCV-1a,HCV1A-H77-NS5b,15,1,558,0,0,0,0,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,1,559,0,0,0,0,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,4,560,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,7,561,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
"""
    expected_failures = {('HCV-1a', 'HCV1A-H77-NS5b', False): 'nothing mapped'}

    check_combination(amino_csv, midi_amino_csv, expected_csv, expected_failures)


def test_combine_aminos_midi_other_regions():
    amino_csv = StringIO("""\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,g2p_overlap,coverage
HCV-1a,HCV1A-H77-NS5b,15,1,1,0,0,0,0,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,4,2,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,7,228,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
""")
    amino_csv.name = 'main_amino.csv'
    midi_amino_csv = StringIO("""\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,g2p_overlap,coverage
HCV-1a,HCV1A-H77-NS3,15,1,1,0,0,0,0,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS3,15,7,181,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5a,15,1,1,0,0,0,0,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5a,15,7,101,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
""")
    midi_amino_csv.name = 'midi_amino.csv'
    expected_csv = """\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,g2p_overlap,coverage
HCV-1a,HCV1A-H77-NS5b,15,1,1,0,0,0,0,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,4,2,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,7,228,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
"""
    expected_failures = {('HCV-1a', 'HCV1A-H77-NS5b', True): 'nothing mapped'}

    check_combination(amino_csv, midi_amino_csv, expected_csv, expected_failures)


def test_combine_aminos_combine_ignores_early_midi():
    amino_csv = StringIO("""\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,g2p_overlap,coverage
HCV-1a,HCV1A-H77-NS5b,15,1,1,0,0,0,0,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,4,2,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,7,228,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
""")
    amino_csv.name = 'main_amino.csv'
    midi_amino_csv = StringIO("""\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,g2p_overlap,coverage
HCV-1a,HCV1A-H77-NS5b,15,1,226,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,1,558,0,0,0,0,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,1,559,0,0,0,0,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,4,560,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,7,561,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
""")
    midi_amino_csv.name = 'midi_amino.csv'
    expected_csv = """\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,g2p_overlap,coverage
HCV-1a,HCV1A-H77-NS5b,15,1,1,0,0,0,0,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,4,2,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,7,228,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,1,558,0,0,0,0,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,1,559,0,0,0,0,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,4,560,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,7,561,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
"""

    check_combination(amino_csv, midi_amino_csv, expected_csv)


def test_combine_aminos_combine_takes_higher_coverage():
    amino_csv = StringIO("""\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,g2p_overlap,coverage
HCV-1a,HCV1A-H77-NS5b,15,1,1,0,0,0,0,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,4,2,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,7,228,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,10,229,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
""")
    amino_csv.name = 'main_amino.csv'
    midi_amino_csv = StringIO("""\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,g2p_overlap,coverage
HCV-1a,HCV1A-H77-NS5b,15,1,228,0,0,10001,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10001
HCV-1a,HCV1A-H77-NS5b,15,4,229,0,0,9999,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9999
HCV-1a,HCV1A-H77-NS5b,15,7,558,0,0,0,0,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,10,559,0,0,0,0,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,13,560,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,16,561,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
""")
    midi_amino_csv.name = 'midi_amino.csv'
    expected_csv = """\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,g2p_overlap,coverage
HCV-1a,HCV1A-H77-NS5b,15,1,1,0,0,0,0,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,4,2,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,1,228,0,0,10001,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10001
HCV-1a,HCV1A-H77-NS5b,15,10,229,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,7,558,0,0,0,0,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,10,559,0,0,0,0,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,13,560,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,16,561,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
"""

    check_combination(amino_csv, midi_amino_csv, expected_csv)


def test_combine_aminos_ignores_main_tail():
    amino_csv = StringIO("""\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,g2p_overlap,coverage
HCV-1a,HCV1A-H77-NS5b,15,1,1,0,0,0,0,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,4,2,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,7,228,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,10,337,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
""")
    amino_csv.name = 'main_amino.csv'
    midi_amino_csv = StringIO("""\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,g2p_overlap,coverage
HCV-1a,HCV1A-H77-NS5b,15,7,558,0,0,0,0,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,10,559,0,0,0,0,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,13,560,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,16,561,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
""")
    midi_amino_csv.name = 'midi_amino.csv'
    expected_csv = """\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,g2p_overlap,coverage
HCV-1a,HCV1A-H77-NS5b,15,1,1,0,0,0,0,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,4,2,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,7,228,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,7,558,0,0,0,0,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,10,559,0,0,0,0,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,13,560,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,16,561,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
"""

    check_combination(amino_csv, midi_amino_csv, expected_csv)


def test_combine_aminos_ignores_main_tail_with_midi():
    amino_csv = StringIO("""\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,g2p_overlap,coverage
HCV-1a,HCV1A-H77-NS5b,15,1,1,0,0,0,0,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,4,2,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,7,228,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,10,337,10001,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10001
""")
    amino_csv.name = 'main_amino.csv'
    midi_amino_csv = StringIO("""\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,g2p_overlap,coverage
HCV-1a,HCV1A-H77-NS5b,15,4,337,0,0,0,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,7,558,0,0,0,0,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,10,559,0,0,0,0,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,13,560,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,16,561,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
""")
    midi_amino_csv.name = 'midi_amino.csv'
    expected_csv = """\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,g2p_overlap,coverage
HCV-1a,HCV1A-H77-NS5b,15,1,1,0,0,0,0,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,4,2,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,7,228,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,4,337,0,0,0,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,7,558,0,0,0,0,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,10,559,0,0,0,0,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,13,560,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,16,561,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
"""

    check_combination(amino_csv, midi_amino_csv, expected_csv)


def test_combine_aminos_ignores_main_tail_without_midi():
    amino_csv = StringIO("""\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,g2p_overlap,coverage
HCV-1a,HCV1A-H77-NS5b,15,1,1,0,0,0,0,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,4,2,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,7,228,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,10,337,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
""")
    amino_csv.name = 'main_amino.csv'
    midi_amino_csv = amino_csv
    expected_csv = """\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,g2p_overlap,coverage
HCV-1a,HCV1A-H77-NS5b,15,1,1,0,0,0,0,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,4,2,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,7,228,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
"""

    check_combination(amino_csv, midi_amino_csv, expected_csv)


def test_combine_aminos_midi_only():
    amino_csv = StringIO("""\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,g2p_overlap,coverage
HCV-1a,HCV1A-H77-NS5b,15,1,1,0,0,0,0,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
""")
    amino_csv.name = 'main_amino.csv'
    midi_amino_csv = StringIO("""\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,g2p_overlap,coverage
HCV-1a,HCV1A-H77-NS5b,15,1,558,0,0,0,0,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,1,559,0,0,0,0,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,4,560,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,7,561,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
""")
    midi_amino_csv.name = 'midi_amino.csv'
    expected_csv = """\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,g2p_overlap,coverage
HCV-1a,HCV1A-H77-NS5b,15,1,558,0,0,0,0,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,1,559,0,0,0,0,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,4,560,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,7,561,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
"""
    expected_failures = {('HCV-1a', 'HCV1A-H77-NS5b', False): 'low average coverage'}

    check_combination(amino_csv, midi_amino_csv, expected_csv, expected_failures)


def test_combine_aminos_midi_only_from_random_primer():
    """ Random primer samples can cover the whole range. """
    amino_csv = StringIO("""\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,g2p_overlap,coverage
HCV-1a,HCV1A-H77-NS5b,15,1,558,0,0,0,0,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,4,559,0,0,0,0,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,7,560,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,10,561,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
""")
    amino_csv.name = 'amino.csv'
    expected_csv = """\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,g2p_overlap,coverage
HCV-1a,HCV1A-H77-NS5b,15,1,558,0,0,0,0,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,4,559,0,0,0,0,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,7,560,0,0,0,0,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,10,561,10000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10000
"""
    expected_failures = {}

    check_combination(amino_csv, amino_csv, expected_csv, expected_failures)


def test_write_nuc_mutations_not_sars():
    """ Only SARS-CoV-2 samples report nucleotide mutations. """
    nuc_csv = StringIO("""\
seed,region,q-cutoff,query.nuc.pos,refseq.nuc.pos,A,C,G,T,N,del,ins,clip,v3_overlap,coverage
HCV-1a,HCV1A-H77-NS5b,15,2101,919,0,0,10000,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,2102,920,0,0,10000,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,2103,921,0,0,0,10000,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,2104,922,0,10000,0,0,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,2105,923,0,0,0,10000,0,0,0,0,0,10000
HCV-1a,HCV1A-H77-NS5b,15,2106,924,0,10000,0,0,0,0,0,0,0,10000
""")
    expected_nuc_mutations_csv = """\
seed,region,wt,refseq_nuc_pos,var,prevalence,ref_genome_pos
"""
    nuc_mutations_csv = StringIO()

    write_nuc_mutations(nuc_csv, nuc_mutations_csv)

    assert nuc_mutations_csv.getvalue() == expected_nuc_mutations_csv


def test_write_nuc_mutation_complete():
    nuc_csv = StringIO("""\
seed,region,q-cutoff,query.nuc.pos,refseq.nuc.pos,genome.pos,\
A,C,G,T,N,del,ins,clip,v3_overlap,coverage
SARS-CoV-2-seed,SARS-CoV-2-ORF3a,15,25393,1,25393,100,0,0,0,0,0,0,0,0,100
SARS-CoV-2-seed,SARS-CoV-2-ORF3a,15,25394,2,25394,0,0,0,100,0,0,0,0,0,100
SARS-CoV-2-seed,SARS-CoV-2-ORF3a,15,25395,3,25395,0,100,0,0,0,0,0,0,0,100
""")
    expected_nuc_mutations_csv = """\
seed,region,wt,refseq_nuc_pos,var,prevalence,ref_genome_pos
SARS-CoV-2-seed,SARS-CoV-2-ORF3a,G,3,C,1.0,25395
"""
    nuc_mutations_csv = StringIO()

    write_nuc_mutations(nuc_csv, nuc_mutations_csv)

    assert nuc_mutations_csv.getvalue() == expected_nuc_mutations_csv


def test_write_nuc_mutations_prevalence():
    nuc_csv = StringIO("""\
seed,region,q-cutoff,query.nuc.pos,refseq.nuc.pos,genome.pos,\
A,C,G,T,N,del,ins,clip,v3_overlap,coverage
SARS-CoV-2-seed,SARS-CoV-2-ORF3a,15,25393,1,25393,100,0,0,0,0,0,0,0,0,100
SARS-CoV-2-seed,SARS-CoV-2-ORF3a,15,25394,2,25394,0,0,5,95,0,0,0,0,0,100
SARS-CoV-2-seed,SARS-CoV-2-ORF3a,15,25395,3,25395,0,4,96,0,0,0,0,0,0,100
""")
    expected_nuc_mutations_csv = """\
seed,region,wt,refseq_nuc_pos,var,prevalence,ref_genome_pos
SARS-CoV-2-seed,SARS-CoV-2-ORF3a,T,2,G,0.05,25394
"""
    nuc_mutations_csv = StringIO()

    write_nuc_mutations(nuc_csv, nuc_mutations_csv)

    assert nuc_mutations_csv.getvalue() == expected_nuc_mutations_csv


def test_write_nuc_mutations_multiple():
    nuc_csv = StringIO("""\
seed,region,q-cutoff,query.nuc.pos,refseq.nuc.pos,genome.pos,\
A,C,G,T,N,del,ins,clip,v3_overlap,coverage
SARS-CoV-2-seed,SARS-CoV-2-ORF3a,15,25393,1,25393,100,0,0,0,0,0,0,0,0,100
SARS-CoV-2-seed,SARS-CoV-2-ORF3a,15,25394,2,25394,0,5,6,89,0,0,0,0,0,100
SARS-CoV-2-seed,SARS-CoV-2-ORF3a,15,25395,3,25395,0,4,96,0,0,0,0,0,0,100
""")
    expected_nuc_mutations_csv = """\
seed,region,wt,refseq_nuc_pos,var,prevalence,ref_genome_pos
SARS-CoV-2-seed,SARS-CoV-2-ORF3a,T,2,C,0.05,25394
SARS-CoV-2-seed,SARS-CoV-2-ORF3a,T,2,G,0.06,25394
"""
    nuc_mutations_csv = StringIO()

    write_nuc_mutations(nuc_csv, nuc_mutations_csv)

    assert nuc_mutations_csv.getvalue() == expected_nuc_mutations_csv


def test_write_nuc_mutations_no_coverage():
    nuc_csv = StringIO("""\
seed,region,q-cutoff,query.nuc.pos,refseq.nuc.pos,genome.pos,\
A,C,G,T,N,del,ins,clip,v3_overlap,coverage
SARS-CoV-2-seed,SARS-CoV-2-ORF3a,15,25393,1,25393,100,0,0,0,0,0,0,0,0,100
SARS-CoV-2-seed,SARS-CoV-2-ORF3a,15,25394,2,25394,0,0,0,0,0,0,0,0,0,0
SARS-CoV-2-seed,SARS-CoV-2-ORF3a,15,25395,3,25395,0,0,100,0,0,0,0,0,0,100
""")
    expected_nuc_mutations_csv = """\
seed,region,wt,refseq_nuc_pos,var,prevalence,ref_genome_pos
"""
    nuc_mutations_csv = StringIO()

    write_nuc_mutations(nuc_csv, nuc_mutations_csv)

    assert nuc_mutations_csv.getvalue() == expected_nuc_mutations_csv


def test_read_aminos_simple():
    amino_csv = DictReader(StringIO("""\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,g2p_overlap,coverage
R1-seed,R1,15,1,1,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
R1-seed,R1,15,4,2,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
R1-seed,R1,15,7,3,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
"""))
    min_fraction = 0.2
    expected_aminos = [AminoList('R1',
                                 [{'K': 1.0}, {'F': 1.0}, {'A': 1.0}],
                                 None,
                                 'R1-seed')]

    aminos = list(read_aminos(amino_csv, min_fraction))

    assert expected_aminos == aminos


def test_read_aminos_hcv():
    amino_csv = DictReader(StringIO("""\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,g2p_overlap,coverage
HCV-1a,HCV1A-H77-NS5b,15,1,1,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
HCV-1a,HCV1A-H77-NS5b,15,4,2,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
HCV-1a,HCV1A-H77-NS5b,15,7,3,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
"""))
    min_fraction = 0.2
    expected_aminos = [AminoList('HCV1A-H77-NS5b',
                                 [{'K': 1.0}, {'F': 1.0}, {'A': 1.0}],
                                 '1A',
                                 'HCV-1a')]

    aminos = list(read_aminos(amino_csv, min_fraction))

    assert expected_aminos == aminos


def test_read_aminos_mixtures():
    amino_csv = DictReader(StringIO("""\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,g2p_overlap,coverage
R1-seed,R1,15,1,1,1,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10
R1-seed,R1,15,4,2,2,0,0,0,8,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10
R1-seed,R1,15,,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
"""))
    min_fraction = 0.2
    expected_aminos = [AminoList('R1',
                                 [{'K': 0.9}, {'F': 0.8, 'A': 0.2}, {}],
                                 None,
                                 'R1-seed')]

    aminos = list(read_aminos(amino_csv, min_fraction))

    assert expected_aminos == aminos


def test_read_aminos_no_coverage():
    amino_csv = DictReader(StringIO("""\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,g2p_overlap,coverage
R1-seed,R1,15,1,1,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
R1-seed,R1,15,4,2,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
R1-seed,R1,15,,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
"""))
    min_fraction = 0.2
    expected_aminos = [AminoList('R1',
                                 [{'K': 1.0}, {'F': 1.0}, {}],
                                 None,
                                 'R1-seed')]

    aminos = list(read_aminos(amino_csv, min_fraction))

    assert expected_aminos == aminos


def test_read_aminos_insertion_before_coverage():
    amino_csv = DictReader(StringIO("""\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,g2p_overlap,coverage
R1-seed,R1,15,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9,0,0,0
R1-seed,R1,15,4,2,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
R1-seed,R1,15,7,3,0,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
"""))
    min_fraction = 0.2
    expected_aminos = [AminoList('R1',
                                 [{}, {'F': 1.0}, {'L': 1.0}],
                                 None,
                                 'R1-seed')]

    aminos = list(read_aminos(amino_csv, min_fraction))

    assert expected_aminos == aminos


def test_read_aminos_multiple_regions():
    amino_csv = DictReader(StringIO("""\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,g2p_overlap,coverage
R1-seed,R1,15,1,1,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
R1-seed,R1,15,4,2,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
R1-seed,R1,15,,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
R2-seed,R2,15,1,1,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
R2-seed,R2,15,4,2,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
R2-seed,R2,15,7,3,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
"""))
    min_fraction = 0.2
    expected_aminos = [AminoList('R1',
                                 [{'K': 1.0}, {'F': 1.0}, {}],
                                 None,
                                 'R1-seed'),
                       AminoList('R2',
                                 [{'K': 1.0}, {'F': 1.0}, {'C': 1.0}],
                                 None,
                                 'R2-seed')]

    aminos = list(read_aminos(amino_csv, min_fraction))

    assert expected_aminos == aminos


def test_read_aminos_reported_regions():
    amino_csv = DictReader(StringIO("""\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,g2p_overlap,coverage
R1-seed,R1,15,1,1,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
R1-seed,R1,15,4,2,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
R1-seed,R1,15,,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
R2-seed,R2,15,1,1,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
R2-seed,R2,15,4,2,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
R2-seed,R2,15,7,3,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
R3-seed,R3,15,1,1,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
R3-seed,R3,15,4,2,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
R3-seed,R3,15,7,3,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
"""))
    min_fraction = 0.2
    reported_regions = {'R1', 'R2'}
    expected_aminos = [AminoList('R1',
                                 [{'K': 1.0}, {'F': 1.0}, {}],
                                 None,
                                 'R1-seed'),
                       AminoList('R2',
                                 [{'K': 1.0}, {'F': 1.0}, {'C': 1.0}],
                                 None,
                                 'R2-seed')]

    aminos = list(read_aminos(amino_csv,
                              min_fraction,
                              reported_regions=reported_regions))

    assert expected_aminos == aminos


def test_read_aminos_low_coverage():
    amino_csv = DictReader(StringIO("""\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,g2p_overlap,coverage
R1-seed,R1,15,1,1,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
R1-seed,R1,15,4,2,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
R1-seed,R1,15,7,3,0,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4
R2-seed,R2,15,1,1,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
R2-seed,R2,15,4,2,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
R2-seed,R2,15,7,3,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
R3-seed,R3,15,1,1,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
R3-seed,R3,15,4,2,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
R3-seed,R3,15,7,3,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
"""))
    min_fraction = 0.2
    min_coverage = 5
    reported_regions = {'R1', 'R2'}
    expected_aminos = [AminoList('R1',
                                 [{'K': 1.0}, {'F': 1.0}, {}],
                                 None,
                                 'R1-seed'),
                       AminoList('R2',
                                 [{'K': 1.0}, {'F': 1.0}, {'C': 1.0}],
                                 None,
                                 'R2-seed')]

    aminos = list(read_aminos(amino_csv,
                              min_fraction,
                              reported_regions=reported_regions,
                              min_coverage=min_coverage))

    assert expected_aminos == aminos


def test_read_aminos_good_average_coverage():
    amino_csv = DictReader(StringIO("""\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,g2p_overlap,coverage
HCV-1a,HCV-Foo-NS5a,15,1,100,0,0,0,0,0,0,0,0,505,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,505
HCV-1a,HCV-Foo-NS5a,15,4,101,0,0,0,0,404,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,404
"""))
    min_fraction = 0.2
    min_coverage = 9
    expected_aminos = [AminoList('HCV-Foo-NS5a',
                                 [{}]*99 + [{'K': 1.0}, {'F': 1.0}],
                                 '1A',
                                 'HCV-1a')]

    aminos = list(read_aminos(amino_csv,
                              min_fraction,
                              min_coverage=min_coverage))

    assert expected_aminos == aminos


def test_read_aminos_resistant_even_with_missing_midi():
    amino_csv = DictReader(StringIO("""\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,g2p_overlap,coverage
HCV-1a,HCV1A-H77-NS5b,15,1,395,0,0,0,0,0,5000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,5000
HCV-1a,HCV1A-H77-NS5b,15,4,396,0,0,0,0,5000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,5000
"""))
    min_fraction = 0.2
    min_coverage = 9
    expected_aminos = [AminoList('HCV1A-H77-NS5b',
                                 [{}]*394 + [{'G': 1.0}, {'F': 1.0}],
                                 '1A',
                                 'HCV-1a')]

    aminos = list(read_aminos(amino_csv,
                              min_fraction,
                              min_coverage=min_coverage))

    assert expected_aminos == aminos


def test_read_aminos_missing_region():
    amino_csv = DictReader(StringIO("""\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,g2p_overlap,coverage
R2-seed,R2,15,1,1,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
R2-seed,R2,15,4,2,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
R2-seed,R2,15,7,3,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
"""))
    min_fraction = 0.2
    reported_regions = {'R1', 'R2'}
    expected_aminos = [AminoList('R2',
                                 [{'K': 1.0}, {'F': 1.0}, {'C': 1.0}],
                                 None,
                                 'R2-seed')]

    aminos = list(read_aminos(amino_csv,
                              min_fraction,
                              reported_regions=reported_regions))

    assert expected_aminos == aminos


def test_read_aminos_deletions():
    amino_csv = DictReader(StringIO("""\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,g2p_overlap,coverage
R1-seed,R1,15,1,1,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,10
R1-seed,R1,15,4,2,0,0,0,0,8,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,10
R1-seed,R1,15,,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
"""))
    min_fraction = 0.2
    expected_aminos = [AminoList('R1',
                                 [{'K': 0.9}, {'F': 0.8, 'd': 0.2}, {}],
                                 None,
                                 'R1-seed')]

    aminos = list(read_aminos(amino_csv, min_fraction))

    assert expected_aminos == aminos


def test_read_aminos_insertions():
    amino_csv = DictReader(StringIO("""\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,g2p_overlap,coverage
R1-seed,R1,15,1,1,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,9
R1-seed,R1,15,4,2,0,0,0,0,8,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,8
R1-seed,R1,15,,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
"""))
    min_fraction = 0.2
    expected_aminos = [AminoList('R1',
                                 [{'K': 1.0}, {'F': 1.0, 'i': 0.25}, {}],
                                 None,
                                 'R1-seed')]

    aminos = list(read_aminos(amino_csv, min_fraction))

    assert expected_aminos == aminos


def test_read_aminos_stop_codons():
    amino_csv = DictReader(StringIO("""\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,g2p_overlap,coverage
R1-seed,R1,15,1,1,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,10
R1-seed,R1,15,4,2,0,0,0,0,8,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,10
R1-seed,R1,15,,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
"""))
    min_fraction = 0.2
    expected_aminos = [AminoList('R1',
                                 [{'K': 0.9}, {'F': 0.8}, {}],
                                 None,
                                 'R1-seed')]

    aminos = list(read_aminos(amino_csv, min_fraction))

    assert expected_aminos == aminos


def test_read_aminos_missing_position():
    amino_csv = DictReader(StringIO("""\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,g2p_overlap,coverage
R1-seed,R1,15,4,2,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
R1-seed,R1,15,7,3,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
"""))
    min_fraction = 0.2
    expected_aminos = [AminoList('R1',
                                 [{}, {'F': 1.0}, {'A': 1.0}],
                                 None,
                                 'R1-seed')]

    aminos = list(read_aminos(amino_csv, min_fraction))

    assert expected_aminos == aminos


def test_read_aminos_pos_55_resistant(asi_algorithms):
    amino_lines = [
        'seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,'
        'A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,'
        'partial,del,ins,clip,g2p_overlap,coverage'] + [
        f'HCV-1a,HCV1A-H77-NS3,15,{i},{i},20,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,20'
        for i in range(1, 632)]
    amino_lines[55] = \
        'HCV-1a,HCV1A-H77-NS3,15,55,55,8,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,8'
    amino_csv = DictReader(amino_lines)
    min_fraction = 0.2
    min_coverage = 9
    expected_aminos = [AminoList('HCV1A-H77-NS3',
                                 [{'A': 1.0}]*54 + [{}] + [{'A': 1.0}]*576,
                                 '1A',
                                 'HCV-1a',
                                 False)]

    aminos = list(read_aminos(amino_csv,
                              min_fraction,
                              min_coverage=min_coverage,
                              algorithms=asi_algorithms))

    assert expected_aminos == aminos


def test_read_aminos_pos_55_not_resistant(asi_algorithms):
    amino_lines = [
        'seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,'
        'A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,'
        'partial,del,ins,clip,g2p_overlap,coverage'] + [
        f'HCV-1b,HCV1B-Con1-NS3,15,{i},{i},20,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,20'
        for i in range(1, 632)]
    amino_lines[55] = \
        'HCV-1b,HCV1B-Con1-NS3,15,55,55,8,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,8'
    amino_csv = DictReader(amino_lines)
    min_fraction = 0.2
    min_coverage = 9
    expected_aminos = [AminoList('HCV1B-Con1-NS3',
                                 [{'A': 1.0}]*54 + [{}] + [{'A': 1.0}]*576,
                                 '1B',
                                 'HCV-1b',
                                 True)]

    aminos = list(read_aminos(amino_csv,
                              min_fraction,
                              min_coverage=min_coverage,
                              algorithms=asi_algorithms))

    assert expected_aminos == aminos


def test_read_aminos_no_midi(asi_algorithms):
    amino_lines = [
        'seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,'
        'A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,'
        'partial,del,ins,clip,g2p_overlap,coverage'] + [
        f'HCV-1b,HCV1B-Con1-NS5b,15,{i},{i},8,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,8'
        for i in range(1, 201)] + [
        f'HCV-1b,HCV1B-Con1-NS5b,15,{i},{i},20,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,20'
        for i in range(201, 337)]
    amino_csv = DictReader(amino_lines)
    min_fraction = 0.2
    min_coverage = 9
    no_counts: typing.Dict[str, float] = {}
    expected_aminos = [AminoList('HCV1B-Con1-NS5b',
                                 [no_counts] * 200 +
                                 [{'A': 1.0}]*136 + [no_counts] * 256,
                                 '1B',
                                 'HCV-1b',
                                 False)]

    aminos = list(read_aminos(amino_csv,
                              min_fraction,
                              min_coverage=min_coverage,
                              algorithms=asi_algorithms))

    assert aminos == expected_aminos


def test_read_aminos_ns5b_missing_142(asi_algorithms):
    """ NS5b still reports if it has coverage on all the MIDI positions. """
    amino_lines = [
        'seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,'
        'A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,'
        'partial,del,ins,clip,g2p_overlap,coverage'] + [
        f'HCV-1b,HCV1B-Con1-NS5b,15,{i},{i},20,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,20'
        for i in range(1, 593)]
    amino_lines[142] = \
        'HCV-1b,HCV1B-Con1-NS5b,15,142,142,8,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,8'
    amino_csv = DictReader(amino_lines)
    min_fraction = 0.2
    min_coverage = 9
    expected_aminos = [AminoList('HCV1B-Con1-NS5b',
                                 [{'A': 1.0}]*141 + [{}] + [{'A': 1.0}]*450,
                                 '1B',
                                 'HCV-1b',
                                 True)]

    aminos = list(read_aminos(amino_csv,
                              min_fraction,
                              min_coverage=min_coverage,
                              algorithms=asi_algorithms))

    assert_aminos_lists_match(expected_aminos, aminos)


def test_read_aminos_ns5b_missing_355(asi_algorithms):
    """ NS5b still reports if it has coverage on all the whole genome positions. """
    amino_lines = [
        'seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,'
        'A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,'
        'partial,del,ins,clip,g2p_overlap,coverage'] + [
        f'HCV-1b,HCV1B-Con1-NS5b,15,{i},{i},20,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,20'
        for i in range(1, 593)]
    amino_lines[355] = \
        'HCV-1b,HCV1B-Con1-NS5b,15,355,355,8,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,8'
    amino_csv = DictReader(amino_lines)
    min_fraction = 0.2
    min_coverage = 9
    expected_aminos = [AminoList('HCV1B-Con1-NS5b',
                                 [{'A': 1.0}]*354 + [{}] + [{'A': 1.0}]*237,
                                 '1B',
                                 'HCV-1b',
                                 True)]

    aminos = list(read_aminos(amino_csv,
                              min_fraction,
                              min_coverage=min_coverage,
                              algorithms=asi_algorithms))

    assert_aminos_lists_match(expected_aminos, aminos)


def test_read_aminos_ns5b_missing_142_and_355(asi_algorithms):
    """ NS5b does not report if it is missing MIDI and whole genome positions. """
    amino_lines = [
        'seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,'
        'A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,'
        'partial,del,ins,clip,g2p_overlap,coverage'] + [
        f'HCV-1b,HCV1B-Con1-NS5b,15,{i},{i},20,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,20'
        for i in range(1, 593)]
    amino_lines[142] = \
        'HCV-1b,HCV1B-Con1-NS5b,15,142,142,8,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,8'
    amino_lines[355] = \
        'HCV-1b,HCV1B-Con1-NS5b,15,355,355,8,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,8'
    amino_csv = DictReader(amino_lines)
    min_fraction = 0.2
    min_coverage = 9
    expected_mixtures = [{'A': 1.0}] * 592
    expected_mixtures[141] = {}
    expected_mixtures[354] = {}
    expected_aminos = [AminoList('HCV1B-Con1-NS5b',
                                 expected_mixtures,
                                 '1B',
                                 'HCV-1b',
                                 False)]

    aminos = list(read_aminos(amino_csv,
                              min_fraction,
                              min_coverage=min_coverage,
                              algorithms=asi_algorithms))

    assert_aminos_lists_match(expected_aminos, aminos)


def test_read_aminos_ns5b_short(asi_algorithms):
    """ Low coverage at the end is no longer reported in amino.csv. """
    amino_lines = [
        'seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,'
        'A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,'
        'partial,del,ins,clip,g2p_overlap,coverage'] + [
        f'HCV-1b,HCV1B-Con1-NS5b,15,{i},{i},20,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,20'
        for i in range(1, 401)]
    amino_csv = DictReader(amino_lines)
    min_fraction = 0.2
    min_coverage = 9
    no_counts: typing.Dict[str, float] = {}
    expected_mixtures = [{'A': 1.0}] * 400 + [no_counts] * 192
    expected_aminos = [AminoList('HCV1B-Con1-NS5b',
                                 expected_mixtures,
                                 '1B',
                                 'HCV-1b',
                                 True)]

    aminos = list(read_aminos(amino_csv,
                              min_fraction,
                              min_coverage=min_coverage,
                              algorithms=asi_algorithms))

    assert_aminos_lists_match(aminos, expected_aminos)


def test_filter_aminos_all(asi_algorithms):
    all_aminos = [AminoList('PR',
                            [{'K': 1.0}, {'F': 1.0}, {'A': 1.0}],
                            None,
                            'HIV1B-seed',
                            True),
                  AminoList('RT',
                            [{'C': 1.0}, {'A': 1.0}, {'N': 1.0}],
                            None,
                            'HIV1B-seed',
                            True),
                  AminoList('INT',
                            [{'E': 1.0}, {'A': 1.0}, {'T': 1.0}],
                            None,
                            'HIV1B-seed',
                            True)]
    expected_aminos = [AminoList('CA', [{}]*249, None, 'HIV1B-seed'),
                       AminoList('INT',
                                 [{'E': 1.0}, {'A': 1.0}, {'T': 1.0}],
                                 None,
                                 'HIV1B-seed',
                                 True),
                       AminoList('PR',
                                 [{'K': 1.0}, {'F': 1.0}, {'A': 1.0}],
                                 None,
                                 'HIV1B-seed',
                                 True),
                       AminoList('RT',
                                 [{'C': 1.0}, {'A': 1.0}, {'N': 1.0}],
                                 None,
                                 'HIV1B-seed',
                                 True)]

    aminos = filter_aminos(all_aminos, asi_algorithms)

    assert expected_aminos == aminos


def test_filter_aminos_missing(asi_algorithms):
    all_aminos = [AminoList('RT',
                            [{'K': 1.0}, {'F': 1.0}, {'A': 1.0}],
                            None,
                            'HIV1B-seed',
                            True)]
    expected_aminos = [AminoList('CA', [{}]*249, None, 'HIV1B-seed'),
                       AminoList('INT', [{}]*289, None, 'HIV1B-seed'),
                       AminoList('PR', [{}]*99, None, 'HIV1B-seed'),
                       AminoList('RT',
                                 [{'K': 1.0}, {'F': 1.0}, {'A': 1.0}],
                                 None,
                                 'HIV1B-seed',
                                 True)]

    aminos = filter_aminos(all_aminos, asi_algorithms)

    assert expected_aminos == aminos


def test_filter_aminos_missing_with_genotypes(asi_algorithms):
    all_aminos = [AminoList('HCV1A-H77-NS5b',
                            [{'K': 1.0}, {'F': 1.0}, {'A': 1.0}],
                            '1A',
                            'HCV-1a',
                            True)]
    expected_aminos = [AminoList('HCV1A-H77-NS3', [{}]*631, '1A', 'HCV-1a'),
                       AminoList('HCV1A-H77-NS5a', [{}]*448, '1A', 'HCV-1a'),
                       AminoList('HCV1A-H77-NS5b',
                                 [{'K': 1.0}, {'F': 1.0}, {'A': 1.0}],
                                 '1A',
                                 'HCV-1a',
                                 True)]

    aminos = filter_aminos(all_aminos, asi_algorithms)

    assert expected_aminos == aminos


def test_filter_aminos_missing_with_some_genotypes(asi_algorithms):
    all_aminos = [AminoList('RT',
                            [{'K': 1.0}, {'F': 1.0}, {'A': 1.0}],
                            None,
                            'HIV1B-seed',
                            True),
                  AminoList('HCV1A-H77-NS5b',
                            [{'K': 1.0}, {'F': 1.0}, {'A': 1.0}],
                            '1A',
                            'HCV-1a',
                            True)]
    expected_aminos = [AminoList('CA', [{}]*249, None, 'HIV1B-seed'),
                       AminoList('HCV1A-H77-NS3', [{}] * 631, '1A', 'HCV-1a'),
                       AminoList('HCV1A-H77-NS5a', [{}] * 448, '1A', 'HCV-1a'),
                       AminoList('HCV1A-H77-NS5b',
                                 [{'K': 1.0}, {'F': 1.0}, {'A': 1.0}],
                                 '1A',
                                 'HCV-1a',
                                 True),
                       AminoList('INT', [{}]*289, None, 'HIV1B-seed'),
                       AminoList('PR', [{}]*99, None, 'HIV1B-seed'),
                       AminoList('RT',
                                 [{'K': 1.0}, {'F': 1.0}, {'A': 1.0}],
                                 None,
                                 'HIV1B-seed',
                                 True)]

    aminos = filter_aminos(all_aminos, asi_algorithms)

    assert expected_aminos == aminos


def test_filter_aminos_exclude(asi_algorithms):
    all_aminos = [AminoList('PR',
                            [{'K': 1.0}, {'F': 1.0}, {'A': 1.0}],
                            None,
                            'HIV1B-seed',
                            True),
                  AminoList('RT',
                            [{'K': 1.0}, {'F': 1.0}, {'A': 1.0}],
                            None,
                            'HIV1B-seed',
                            True),
                  AminoList('INT',
                            [{'K': 1.0}, {'F': 1.0}, {'A': 1.0}],
                            None,
                            'HIV1B-seed',
                            True),
                  AminoList('HCV1A-H77-NS3',
                            [{'K': 1.0}, {'F': 1.0}, {'A': 1.0}],
                            '1A',
                            'HCV-1a',
                            False)]
    expected_aminos = [AminoList('CA', [{}]*249, None, 'HIV1B-seed'),
                       AminoList('INT',
                                 [{'K': 1.0}, {'F': 1.0}, {'A': 1.0}],
                                 None,
                                 'HIV1B-seed',
                                 True),
                       AminoList('PR',
                                 [{'K': 1.0}, {'F': 1.0}, {'A': 1.0}],
                                 None,
                                 'HIV1B-seed',
                                 True),
                       AminoList('RT',
                                 [{'K': 1.0}, {'F': 1.0}, {'A': 1.0}],
                                 None,
                                 'HIV1B-seed',
                                 True)]

    aminos = filter_aminos(all_aminos, asi_algorithms)

    assert expected_aminos == aminos


def test_filter_aminos_no_good(asi_algorithms):
    all_aminos = [AminoList('PR', [{}, {}, {}], None, 'HIV1B-seed'),
                  AminoList('RT', [{}, {}, {}], None, 'HIV1B-seed'),
                  AminoList('HCV1A-H77-NS3', [{}, {}, {}], '1A', 'HCV-1a')]
    expected_aminos = []

    aminos = filter_aminos(all_aminos, asi_algorithms)

    assert expected_aminos == aminos


def test_write_failures_empty():
    failures = {}
    expected_fail_csv = """\
seed,region,coord_region,genotype,reason
"""
    check_failure_report(failures, expected_fail_csv)


def test_write_failures_one():
    failures = {('R1-seed', 'R1', False): 'low average coverage'}
    expected_fail_csv = """\
seed,region,coord_region,genotype,reason
R1-seed,R1,R1,,low average coverage
"""
    check_failure_report(failures, expected_fail_csv)


def test_write_failures_two():
    failures = {('R2-seed', 'R2', False): 'not enough high-coverage amino acids',
                ('R1-seed', 'R1', False): 'low average coverage'}
    expected_fail_csv = """\
seed,region,coord_region,genotype,reason
R1-seed,R1,R1,,low average coverage
R2-seed,R2,R2,,not enough high-coverage amino acids
"""
    check_failure_report(failures, expected_fail_csv)


def test_write_failures_midi():
    failures = {('R1-seed', 'R1', False): 'not enough high-coverage amino acids',
                ('R1-seed', 'R1', True): 'low average coverage'}
    expected_fail_csv = """\
seed,region,coord_region,genotype,reason
R1-seed,R1,R1,,not enough high-coverage amino acids
R1-seed,R1,R1,,MIDI: low average coverage
"""
    check_failure_report(failures, expected_fail_csv)


def check_failure_report(failures, expected_fail_csv):
    filtered_aminos = [AminoList('R1',
                                 [{'F': 1.0}, {'A': 1.0}],
                                 None,
                                 'R1-seed'),
                       AminoList('R2',
                                 [{'V': 1.0}, {'L': 1.0}],
                                 None,
                                 'R2-seed')]
    fail_csv = StringIO()

    write_failures(failures, filtered_aminos, fail_csv)

    assert expected_fail_csv == fail_csv.getvalue()


def check_failure_report_nothing_mapped(failures, expected_fail_csv):
    filtered_aminos = [AminoList('R1',
                                 [{}, {}],
                                 None,
                                 'R1-seed'),
                       AminoList('R2',
                                 [{'V': 1.0}, {'L': 1.0}],
                                 None,
                                 'R2-seed')]
    fail_csv = StringIO()
    write_failures(failures, filtered_aminos, fail_csv)
    assert expected_fail_csv == fail_csv.getvalue()


def test_write_failures_nothing_mapped_in_region():
    failures = {}
    expected_fail_csv = """\
seed,region,coord_region,genotype,reason
R1-seed,R1,R1,,nothing mapped
"""
    check_failure_report_nothing_mapped(failures, expected_fail_csv)


def test_write_failures_nothing_mapped_conflict():
    failures = {('R1-seed', 'R1', False): 'low average coverage'}
    expected_fail_csv = """\
seed,region,coord_region,genotype,reason
R1-seed,R1,R1,,low average coverage
"""
    check_failure_report_nothing_mapped(failures, expected_fail_csv)


def test_write_failures_extra_region():
    failures = {('R3-seed', 'R3', False): 'low average coverage'}
    expected_fail_csv = """\
seed,region,coord_region,genotype,reason
"""
    check_failure_report(failures, expected_fail_csv)


def test_write_failures_hcv():
    failures = {('HCV-1a', 'HCV1A-H77-NS3', False): 'low average coverage',
                ('HCV-1b', 'HCV1B-Con1-NS5b-NS3', False): 'low average coverage'}
    expected_fail_csv = """\
seed,region,coord_region,genotype,reason
HCV-1a,NS3,HCV1A-H77-NS3,1A,low average coverage
"""
    filtered_aminos = [AminoList('HCV1A-H77-NS3',
                                 [{'F': 1.0}, {'A': 1.0}],
                                 '1A',
                                 'HCV-1a'),
                       AminoList('HCV1A-H77-NS5a',
                                 [{'V': 1.0}, {'L': 1.0}],
                                 '1A',
                                 'HCV-1a')]
    fail_csv = StringIO()

    write_failures(failures, filtered_aminos, fail_csv)

    assert expected_fail_csv == fail_csv.getvalue()


def test_write_failures_no_mapped_regions():
    failures = {('HCV-1a', 'HCV1A-H77-NS3', False): 'low average coverage',
                ('HCV-1b', 'HCV1B-Con1-NS3', False): 'low average coverage'}
    expected_fail_csv = """\
seed,region,coord_region,genotype,reason
HCV-1a,NS3,HCV1A-H77-NS3,1A,low average coverage
HCV-1b,NS3,HCV1B-Con1-NS3,1B,low average coverage
"""
    filtered_aminos = []
    fail_csv = StringIO()

    write_failures(failures, filtered_aminos, fail_csv)

    assert expected_fail_csv == fail_csv.getvalue()


class WriteResistanceTest(TestCase):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.algorithms = load_asi()

    def test_simple(self):
        rt_aminos = [{c: 1.0} for c in self.algorithms[None].stds['RT']]
        rt_aminos[40] = {'L': 1.0}  # pos 41
        aminos = [AminoList('RT', rt_aminos, None, 'HIV1B-seed')]
        resistance_csv = StringIO()
        mutations_csv = StringIO()
        resistance_consensus_csv = StringIO()
        expected_resistance = f"""\
region,drug_class,drug,drug_name,level,level_name,score,genotype,seed,coord_region,version
RT,NRTI,3TC,lamivudine,1,Susceptible,0.0,,HIV1B-seed,RT,{HIVDB_VERSION}
RT,NRTI,ABC,abacavir,1,Susceptible,5.0,,HIV1B-seed,RT,{HIVDB_VERSION}
RT,NRTI,AZT,zidovudine,3,Low-Level Resistance,15.0,,HIV1B-seed,RT,{HIVDB_VERSION}
RT,NRTI,D4T,stavudine,3,Low-Level Resistance,15.0,,HIV1B-seed,RT,{HIVDB_VERSION}
RT,NRTI,DDI,didanosine,2,Susceptible,10.0,,HIV1B-seed,RT,{HIVDB_VERSION}
RT,NRTI,FTC,emtricitabine,1,Susceptible,0.0,,HIV1B-seed,RT,{HIVDB_VERSION}
RT,NRTI,TDF,tenofovir,1,Susceptible,5.0,,HIV1B-seed,RT,{HIVDB_VERSION}
RT,NNRTI,DOR,doravirine,1,Susceptible,0.0,,HIV1B-seed,RT,{HIVDB_VERSION}
RT,NNRTI,DPV,dapirivine,1,Susceptible,0.0,,HIV1B-seed,RT,{HIVDB_VERSION}
RT,NNRTI,EFV,efavirenz,1,Susceptible,0.0,,HIV1B-seed,RT,{HIVDB_VERSION}
RT,NNRTI,ETR,etravirine,1,Susceptible,0.0,,HIV1B-seed,RT,{HIVDB_VERSION}
RT,NNRTI,NVP,nevirapine,1,Susceptible,0.0,,HIV1B-seed,RT,{HIVDB_VERSION}
RT,NNRTI,RPV,rilpivirine,1,Susceptible,0.0,,HIV1B-seed,RT,{HIVDB_VERSION}
"""
        expected_mutations = f"""\
drug_class,mutation,prevalence,genotype,region,seed,coord_region,version
NRTI,M41L,1.0,,RT,HIV1B-seed,RT,{HIVDB_VERSION}
"""
        expected_resistance_consensus = f"""\
seed,region,coord_region,version,offset,sequence
HIV1B-seed,RT,RT,{HIVDB_VERSION},0,\
PISPIETVPVKLKPGMDGPKVKQWPLTEEKIKALVEICTELEKEGKISKIGPENPYNTPVFAIKK\
KDSTKWRKLVDFRELNKRTQDFWEVQLGIPHPAGLKKKKSVTVLDVGDAYFSVPLDKDFRKYTAF\
TIPSINNETPGIRYQYNVLPQGWKGSPAIFQSSMTKILEPFRKQNPDIVIYQYMDDLYVGSDLEI\
GQHRTKIEELRQHLLRWGFTTPDKKHQKEPPFLWMGYELHPDKWTVQPIVLPEKDSWTVNDIQKL\
VGKLNWASQIYAGIKVKQLCKLLRGTKALTEVIPLTEEAELELAENREILKEPVHGVYYDPSKDL\
IAEIQKQGQGQWTYQIYQEPFKNLKTGKYARMRGAHTNDVKQLTEAVQKIATESIVIWGKTPKFK\
LPIQKETWEAWWTEYWQATWIPEWEFVNTPPLVKLWYQLEKEPIVGAETFYVDGAANRETKLGKA\
GYVTDRGRQKVVSLTDTTNQKTELQAIHLALQDSGLEVNIVTDSQYALGIIQAQPDKSESELVSQ\
IIEQLIKKEKVYLAWVPAHKGIGGNEQVDKLVSAGIRKVL
"""

        write_resistance(aminos,
                         resistance_csv,
                         mutations_csv,
                         self.algorithms,
                         resistance_consensus_csv)

        assert expected_resistance == resistance_csv.getvalue()
        assert expected_mutations == mutations_csv.getvalue()
        assert expected_resistance_consensus == resistance_consensus_csv.getvalue()

    def test_low_coverage(self):
        rt_aminos = [{c: 1.0} for c in self.algorithms[None].stds['RT']]
        rt_aminos[40] = {'L': 1.0}  # pos 41
        aminos = [AminoList('IN', [{}]*289, None, 'HIV1B-seed'),
                  AminoList('PR', [{}]*99, None, 'HIV1B-seed'),
                  AminoList('RT', rt_aminos, None, 'HIV1B-seed')]
        resistance_csv = StringIO()
        mutations_csv = StringIO()
        expected_resistance = f"""\
region,drug_class,drug,drug_name,level,level_name,score,genotype,seed,coord_region,version
IN,INSTI,BIC,bictegravir,0,Sequence does not meet quality-control standards,0.0,,HIV1B-seed,IN,{HIVDB_VERSION}
IN,INSTI,CAB,cabotegravir,0,Sequence does not meet quality-control standards,0.0,,HIV1B-seed,IN,{HIVDB_VERSION}
IN,INSTI,DTG,dolutegravir,0,Sequence does not meet quality-control standards,0.0,,HIV1B-seed,IN,{HIVDB_VERSION}
IN,INSTI,EVG,elvitegravir,0,Sequence does not meet quality-control standards,0.0,,HIV1B-seed,IN,{HIVDB_VERSION}
IN,INSTI,RAL,raltegravir,0,Sequence does not meet quality-control standards,0.0,,HIV1B-seed,IN,{HIVDB_VERSION}
PR,PI,ATV/r,atazanavir/r,0,Sequence does not meet quality-control standards,0.0,,HIV1B-seed,PR,{HIVDB_VERSION}
PR,PI,DRV/r,darunavir/r,0,Sequence does not meet quality-control standards,0.0,,HIV1B-seed,PR,{HIVDB_VERSION}
PR,PI,FPV/r,fosamprenavir/r,0,Sequence does not meet quality-control standards,0.0,,HIV1B-seed,PR,{HIVDB_VERSION}
PR,PI,IDV/r,indinavir/r,0,Sequence does not meet quality-control standards,0.0,,HIV1B-seed,PR,{HIVDB_VERSION}
PR,PI,LPV/r,lopinavir/r,0,Sequence does not meet quality-control standards,0.0,,HIV1B-seed,PR,{HIVDB_VERSION}
PR,PI,NFV,nelfinavir,0,Sequence does not meet quality-control standards,0.0,,HIV1B-seed,PR,{HIVDB_VERSION}
PR,PI,SQV/r,saquinavir/r,0,Sequence does not meet quality-control standards,0.0,,HIV1B-seed,PR,{HIVDB_VERSION}
PR,PI,TPV/r,tipranavir/r,0,Sequence does not meet quality-control standards,0.0,,HIV1B-seed,PR,{HIVDB_VERSION}
RT,NRTI,3TC,lamivudine,1,Susceptible,0.0,,HIV1B-seed,RT,{HIVDB_VERSION}
RT,NRTI,ABC,abacavir,1,Susceptible,5.0,,HIV1B-seed,RT,{HIVDB_VERSION}
RT,NRTI,AZT,zidovudine,3,Low-Level Resistance,15.0,,HIV1B-seed,RT,{HIVDB_VERSION}
RT,NRTI,D4T,stavudine,3,Low-Level Resistance,15.0,,HIV1B-seed,RT,{HIVDB_VERSION}
RT,NRTI,DDI,didanosine,2,Susceptible,10.0,,HIV1B-seed,RT,{HIVDB_VERSION}
RT,NRTI,FTC,emtricitabine,1,Susceptible,0.0,,HIV1B-seed,RT,{HIVDB_VERSION}
RT,NRTI,TDF,tenofovir,1,Susceptible,5.0,,HIV1B-seed,RT,{HIVDB_VERSION}
RT,NNRTI,DOR,doravirine,1,Susceptible,0.0,,HIV1B-seed,RT,{HIVDB_VERSION}
RT,NNRTI,DPV,dapirivine,1,Susceptible,0.0,,HIV1B-seed,RT,{HIVDB_VERSION}
RT,NNRTI,EFV,efavirenz,1,Susceptible,0.0,,HIV1B-seed,RT,{HIVDB_VERSION}
RT,NNRTI,ETR,etravirine,1,Susceptible,0.0,,HIV1B-seed,RT,{HIVDB_VERSION}
RT,NNRTI,NVP,nevirapine,1,Susceptible,0.0,,HIV1B-seed,RT,{HIVDB_VERSION}
RT,NNRTI,RPV,rilpivirine,1,Susceptible,0.0,,HIV1B-seed,RT,{HIVDB_VERSION}
"""
        expected_mutations = f"""\
drug_class,mutation,prevalence,genotype,region,seed,coord_region,version
NRTI,M41L,1.0,,RT,HIV1B-seed,RT,{HIVDB_VERSION}
"""

        write_resistance(aminos, resistance_csv, mutations_csv, self.algorithms)

        assert expected_resistance == resistance_csv.getvalue()
        assert expected_mutations == mutations_csv.getvalue()

    def test_mixture(self):
        rt_aminos = [{c: 1.0} for c in self.algorithms[None].stds['RT']]
        rt_aminos[40] = {'L': 0.3, 'd': 0.7}  # pos 41
        aminos = [AminoList('RT', rt_aminos, None, 'HIV1B-seed')]
        resistance_csv = StringIO()
        mutations_csv = StringIO()
        expected_resistance = f"""\
region,drug_class,drug,drug_name,level,level_name,score,genotype,seed,coord_region,version
RT,NRTI,3TC,lamivudine,1,Susceptible,0.0,,HIV1B-seed,RT,{HIVDB_VERSION}
RT,NRTI,ABC,abacavir,1,Susceptible,5.0,,HIV1B-seed,RT,{HIVDB_VERSION}
RT,NRTI,AZT,zidovudine,3,Low-Level Resistance,15.0,,HIV1B-seed,RT,{HIVDB_VERSION}
RT,NRTI,D4T,stavudine,3,Low-Level Resistance,15.0,,HIV1B-seed,RT,{HIVDB_VERSION}
RT,NRTI,DDI,didanosine,2,Susceptible,10.0,,HIV1B-seed,RT,{HIVDB_VERSION}
RT,NRTI,FTC,emtricitabine,1,Susceptible,0.0,,HIV1B-seed,RT,{HIVDB_VERSION}
RT,NRTI,TDF,tenofovir,1,Susceptible,5.0,,HIV1B-seed,RT,{HIVDB_VERSION}
RT,NNRTI,DOR,doravirine,1,Susceptible,0.0,,HIV1B-seed,RT,{HIVDB_VERSION}
RT,NNRTI,DPV,dapirivine,1,Susceptible,0.0,,HIV1B-seed,RT,{HIVDB_VERSION}
RT,NNRTI,EFV,efavirenz,1,Susceptible,0.0,,HIV1B-seed,RT,{HIVDB_VERSION}
RT,NNRTI,ETR,etravirine,1,Susceptible,0.0,,HIV1B-seed,RT,{HIVDB_VERSION}
RT,NNRTI,NVP,nevirapine,1,Susceptible,0.0,,HIV1B-seed,RT,{HIVDB_VERSION}
RT,NNRTI,RPV,rilpivirine,1,Susceptible,0.0,,HIV1B-seed,RT,{HIVDB_VERSION}
"""
        expected_mutations = f"""\
drug_class,mutation,prevalence,genotype,region,seed,coord_region,version
NRTI,M41L,0.3,,RT,HIV1B-seed,RT,{HIVDB_VERSION}
"""

        write_resistance(aminos, resistance_csv, mutations_csv, self.algorithms)

        assert expected_resistance == resistance_csv.getvalue()
        assert expected_mutations == mutations_csv.getvalue()

    def test_hcv(self):
        aminos = [AminoList('HCV6-EUHK2-NS5b',
                            [{'T': 1.0}] * 592,
                            '6',
                            'HCV-6a')]
        resistance_csv = StringIO()
        mutations_csv = StringIO()
        expected_resistance = """\
region,drug_class,drug,drug_name,level,level_name,score,genotype,seed,coord_region,version
NS5b,NS5b,DSV,Dasabuvir,2,Not Indicated,0.0,6,HCV-6a,HCV6-EUHK2-NS5b,1.7
NS5b,NS5b,SOF-EPC,Sofosbuvir in Epclusa,5,Resistance Likely,8.0,6,HCV-6a,HCV6-EUHK2-NS5b,1.7
NS5b,NS5b,SOF-HAR,Sofosbuvir in Harvoni,2,Not Indicated,0.0,6,HCV-6a,HCV6-EUHK2-NS5b,1.7
"""
        expected_mutations = """\
drug_class,mutation,prevalence,genotype,region,seed,coord_region,version
NS5b,L159T,1.0,6,NS5b,HCV-6a,HCV6-EUHK2-NS5b,1.7
NS5b,S282T,1.0,6,NS5b,HCV-6a,HCV6-EUHK2-NS5b,1.7
NS5b,C316T,1.0,6,NS5b,HCV-6a,HCV6-EUHK2-NS5b,1.7
NS5b,L320T,1.0,6,NS5b,HCV-6a,HCV6-EUHK2-NS5b,1.7
NS5b,V321T,1.0,6,NS5b,HCV-6a,HCV6-EUHK2-NS5b,1.7
"""

        write_resistance(aminos, resistance_csv, mutations_csv, self.algorithms)

        self.assertEqual(expected_resistance, resistance_csv.getvalue())
        self.assertEqual(expected_mutations, mutations_csv.getvalue())

    def test_hcv_missing_midi(self):
        ns5b_aminos = [{c: 1.0}
                       for c in self.algorithms[None].stds['HCV1B-Con1-NS5b']]
        ns5b_aminos[315] = {'Y': 1.0}  # C316Y mutation => Resistance Likely
        ns5b_aminos = ns5b_aminos[:336]  # Missing MIDI
        aminos = [AminoList('HCV1B-Con1-NS5b', ns5b_aminos, '1B', 'HCV-1b')]
        resistance_csv = StringIO()
        mutations_csv = StringIO()
        expected_resistance = """\
region,drug_class,drug,drug_name,level,level_name,score,genotype,seed,coord_region,version
NS5b,NS5b,DSV,Dasabuvir,5,Resistance Likely,8.0,1B,HCV-1b,HCV1B-Con1-NS5b,1.7
NS5b,NS5b,SOF-EPC,Sofosbuvir in Epclusa,0,Sequence does not meet quality-control standards,\
0.0,1B,HCV-1b,HCV1B-Con1-NS5b,1.7
NS5b,NS5b,SOF-HAR,Sofosbuvir in Harvoni,0,Sequence does not meet quality-control standards,\
0.0,1B,HCV-1b,HCV1B-Con1-NS5b,1.7
"""
        expected_mutations = """\
drug_class,mutation,prevalence,genotype,region,seed,coord_region,version
NS5b,C316Y,1.0,1B,NS5b,HCV-1b,HCV1B-Con1-NS5b,1.7
"""

        write_resistance(aminos, resistance_csv, mutations_csv, self.algorithms)

        self.assertEqual(expected_resistance, resistance_csv.getvalue())
        self.assertEqual(expected_mutations, mutations_csv.getvalue())

    def test_hcv_ns3_coverage_gap(self):
        self.maxDiff = None
        ns3_aminos = [{c: 1.0}
                      for c in self.algorithms[None].stds['HCV1A-H77-NS3']]
        # Avoid wild-type resistance to Simeprevir
        ns3_aminos[247] = {'I': 1.0}
        ns3_aminos[146] = {'S': 1.0}

        ns3_aminos[53] = {}  # Coverage gap at position 54
        aminos = [AminoList('HCV1A-H77-NS3', ns3_aminos, '1A', 'HCV-1a')]
        resistance_csv = StringIO()
        mutations_csv = StringIO()
        expected_resistance = """\
region,drug_class,drug,drug_name,level,level_name,score,genotype,seed,coord_region,version
NS3,NS3,ASV,Asunaprevir,0,Sequence does not meet quality-control standards,0.0,1A,HCV-1a,HCV1A-H77-NS3,1.7
NS3,NS3,GLE,Glecaprevir,1,Likely Susceptible,0.0,1A,HCV-1a,HCV1A-H77-NS3,1.7
NS3,NS3,GZR,Grazoprevir,1,Likely Susceptible,0.0,1A,HCV-1a,HCV1A-H77-NS3,1.7
NS3,NS3,PTV,Paritaprevir,1,Likely Susceptible,0.0,1A,HCV-1a,HCV1A-H77-NS3,1.7
NS3,NS3,SPV,Simeprevir,0,Sequence does not meet quality-control standards,0.0,1A,HCV-1a,HCV1A-H77-NS3,1.7
NS3,NS3,VOX,Voxilaprevir,1,Likely Susceptible,0.0,1A,HCV-1a,HCV1A-H77-NS3,1.7
"""
        expected_mutations = """\
drug_class,mutation,prevalence,genotype,region,seed,coord_region,version
"""

        write_resistance(aminos, resistance_csv, mutations_csv, self.algorithms)

        self.assertEqual(expected_mutations, mutations_csv.getvalue())
        self.assertEqual(expected_resistance, resistance_csv.getvalue())

    def test_hcv_low_coverage(self):
        aminos = []
        self.check_low_coverage_reports(aminos)

    def check_low_coverage_reports(self, aminos):
        resistance_csv = StringIO()
        mutations_csv = StringIO()
        expected_resistance = """\
region,drug_class,drug,drug_name,level,level_name,score,genotype,seed,coord_region,version
"""
        expected_mutations = """\
drug_class,mutation,prevalence,genotype,region,seed,coord_region,version
"""

        write_resistance(aminos, resistance_csv, mutations_csv, self.algorithms)

        self.assertEqual(expected_resistance, resistance_csv.getvalue())
        self.assertEqual(expected_mutations, mutations_csv.getvalue())

    def test_hcv_mostly_low_coverage(self):
        aminos = [AminoList('HCV1B-Con1-NS5b',
                            [{'T': 1.0}] + [{}] * 319 + [{'I': 1.0}] + [{}] * 271,
                            '1B',
                            'HCV-1b')]
        self.check_low_coverage_reports(aminos)


class GenotypeTest(TestCase):
    def test_hiv(self):
        self.assertEqual(None, get_genotype('HIV1-B-FR-KF716496-seed'))

    def test_hcv2(self):
        self.assertEqual('2', get_genotype('HCV-2c'))

    def test_hcv6(self):
        self.assertEqual('6', get_genotype('HCV-6o'))

    def test_hcv6e(self):
        self.assertEqual('6E', get_genotype('HCV-6e'))

    def test_hcv1a(self):
        self.assertEqual('1A', get_genotype('HCV-1a'))

    def test_hcv1b(self):
        self.assertEqual('1B', get_genotype('HCV-1b'))

    def test_hcv1c(self):
        self.assertEqual('1A', get_genotype('HCV-1c'))


class WriteConsensusTest(TestCase):
    def setUp(self):
        self.version = '1.0'

    def test_simple(self):
        resistance_consensus_csv = StringIO()
        writer = create_consensus_writer(resistance_consensus_csv)
        amino_list = AminoList('RT',
                               [{'A': 1.0}, {'R': 1.0}, {'M': 1.0}, {'S': 1.0}],
                               None,
                               'HIV1B-seed')
        expected_resistance_consensus_csv = """\
seed,region,coord_region,version,offset,sequence
HIV1B-seed,RT,RT,1.0,0,ARMS
"""
        write_consensus(writer, amino_list, self.version)

        self.assertEqual(expected_resistance_consensus_csv,
                         resistance_consensus_csv.getvalue())

    def test_gap(self):
        resistance_consensus_csv = StringIO()
        writer = create_consensus_writer(resistance_consensus_csv)
        amino_list = AminoList('RT',
                               [{'A': 1.0}, {}, {'M': 1.0}, {'S': 1.0}],
                               None,
                               'HIV1B-seed')
        expected_resistance_consensus_csv = """\
seed,region,coord_region,version,offset,sequence
HIV1B-seed,RT,RT,1.0,0,A-MS
"""
        write_consensus(writer, amino_list, self.version)

        self.assertEqual(expected_resistance_consensus_csv,
                         resistance_consensus_csv.getvalue())

    def test_offset(self):
        resistance_consensus_csv = StringIO()
        writer = create_consensus_writer(resistance_consensus_csv)
        amino_list = AminoList('RT',
                               [{}, {}, {'M': 1.0}, {'S': 1.0}],
                               None,
                               'HIV1B-seed')
        expected_resistance_consensus_csv = """\
seed,region,coord_region,version,offset,sequence
HIV1B-seed,RT,RT,1.0,2,MS
"""
        write_consensus(writer, amino_list, self.version)

        self.assertEqual(expected_resistance_consensus_csv,
                         resistance_consensus_csv.getvalue())

    def test_tail(self):
        resistance_consensus_csv = StringIO()
        writer = create_consensus_writer(resistance_consensus_csv)
        amino_list = AminoList('RT',
                               [{'A': 1.0}, {'R': 1.0}, {}, {}],
                               None,
                               'HIV1B-seed')
        expected_resistance_consensus_csv = """\
seed,region,coord_region,version,offset,sequence
HIV1B-seed,RT,RT,1.0,0,AR
"""
        write_consensus(writer, amino_list, self.version)

        self.assertEqual(expected_resistance_consensus_csv,
                         resistance_consensus_csv.getvalue())

    def test_mixture(self):
        resistance_consensus_csv = StringIO()
        writer = create_consensus_writer(resistance_consensus_csv)
        amino_list = AminoList('RT',
                               [{'A': 1.0}, {'R': 1.0}, {'T': 0.1, 'M': 0.9}, {'S': 1.0}],
                               None,
                               'HIV1B-seed')
        expected_resistance_consensus_csv = """\
seed,region,coord_region,version,offset,sequence
HIV1B-seed,RT,RT,1.0,0,AR[MT]S
"""
        write_consensus(writer, amino_list, self.version)

        self.assertEqual(expected_resistance_consensus_csv,
                         resistance_consensus_csv.getvalue())

    def test_unmapped_region(self):
        resistance_consensus_csv = StringIO()
        writer = create_consensus_writer(resistance_consensus_csv)
        amino_list = AminoList('RT',
                               [{}, {}, {}, {}],
                               None,
                               'HIV1B-seed')
        expected_resistance_consensus_csv = """\
seed,region,coord_region,version,offset,sequence
"""
        write_consensus(writer, amino_list, self.version)

        self.assertEqual(expected_resistance_consensus_csv,
                         resistance_consensus_csv.getvalue())


# Configuration validation tests for resistance.py
# These catch inconsistencies and ensure resistance processing won't fail

def test_reported_regions_have_wild_types(asi_algorithms):
    """All regions in REPORTED_REGIONS must have wild type sequences.
    
    Prevents: KeyError when trying to get reference sequence for region.
    """
    from micall.resistance.resistance import REPORTED_REGIONS
    from micall.data.landmark_reader import LandmarkReader
    from micall.core.project_config import ProjectConfig
    
    projects = ProjectConfig.loadDefault()
    landmark_reader = LandmarkReader.load()
    
    # Filter to HIV regions
    hiv_regions = {r for r in REPORTED_REGIONS if r not in ('NS3', 'NS5a', 'NS5b')}
    
    # Map region names: IN in REPORTED_REGIONS → INT in references
    region_name_map = {'IN': 'INT'}
    
    missing = []
    for region in hiv_regions:
        ref_name = region_name_map.get(region, region)
        try:
            # Try to get reference
            ref = projects.getReference(ref_name)
            if not ref:
                missing.append(region)
        except KeyError:
            missing.append(region)
    
    assert not missing, (
        f"Regions in REPORTED_REGIONS lack reference sequences: {missing}\n"
        f"Fix: Add wild type sequences to micall/data/landmark_references.yaml"
    )


def test_reported_regions_in_algorithm(asi_algorithms):
    """Regions in REPORTED_REGIONS should have algorithm support.
    
    Prevents: Processing regions that can't generate resistance scores.
    """
    from micall.resistance.resistance import REPORTED_REGIONS
    
    hiv_algorithm = asi_algorithms[None]
    algorithm_regions = set(hiv_algorithm.gene_def.keys())
    
    # Map INT to IN for comparison
    if 'INT' in REPORTED_REGIONS:
        algorithm_regions.add('INT')
    
    hiv_regions = {r for r in REPORTED_REGIONS if r not in ('NS3', 'NS5a', 'NS5b')}
    
    unsupported = hiv_regions - algorithm_regions
    assert not unsupported, (
        f"REPORTED_REGIONS has regions without algorithm support: {unsupported}\n"
        f"These regions will be processed but can't generate resistance scores.\n"
        f"Fix: Either add algorithm support or remove from REPORTED_REGIONS"
    )


def test_algorithm_version_matches_config(asi_algorithms):
    """Algorithm version should match HIVDB_VERSION constant.
    
    Prevents: Version mismatch between loaded algorithm and expected version.
    """
    from micall.resistance.resistance import HIVDB_VERSION
    
    hiv_algorithm = asi_algorithms[None]
    algorithm_version = hiv_algorithm.alg_version
    
    assert algorithm_version == HIVDB_VERSION, (
        f"Algorithm version mismatch: loaded {algorithm_version}, "
        f"expected {HIVDB_VERSION}\n"
        f"Fix: Update HIVDB_VERSION in resistance.py to match algorithm file"
    )


def test_get_reported_region_mapping():
    """Test region name mapping for all known regions.
    
    Prevents: Incorrect region names in output causing downstream errors.
    """
    from micall.resistance.resistance import get_reported_region
    
    # Test known mappings
    assert get_reported_region('INT') == 'IN'
    assert get_reported_region('PR') == 'PR'
    assert get_reported_region('RT') == 'RT'
    assert get_reported_region('CA') == 'CA'
    
    # Test HCV regions
    assert get_reported_region('HCV1A-H77-NS3') == 'NS3'
    assert get_reported_region('HCV1A-H77-NS5a') == 'NS5a'
    assert get_reported_region('HCV1A-H77-NS5b') == 'NS5b'


def test_filter_aminos_adds_all_algorithm_regions(asi_algorithms):
    """filter_aminos should add empty entries for all algorithm regions.
    
    Prevents: Missing regions in resistance output causing incomplete reports.
    """
    from micall.resistance.resistance import filter_aminos, get_algorithm_regions
    
    # Start with single region
    input_aminos = [
        AminoList('PR', [{'K': 1.0}] * 30, None, 'HIV1B-seed', True)
    ]
    
    filtered = filter_aminos(input_aminos, asi_algorithms)
    
    # Should have all algorithm regions
    hiv_algorithm = asi_algorithms[None]
    expected_regions = set(get_algorithm_regions(hiv_algorithm))
    actual_regions = {amino.region for amino in filtered}
    
    missing = expected_regions - actual_regions
    assert not missing, (
        f"filter_aminos didn't add all algorithm regions. Missing: {missing}\n"
        f"Fix: Update filter_aminos logic to handle all algorithm regions"
    )


def test_write_resistance_handles_all_regions(asi_algorithms):
    """write_resistance should handle all regions without errors.
    
    Prevents: Runtime errors when processing specific regions.
    """
    from micall.resistance.resistance import write_resistance, get_algorithm_regions
    
    hiv_algorithm = asi_algorithms[None]
    regions = get_algorithm_regions(hiv_algorithm)
    
    # Create amino lists for all regions
    amino_lists = []
    for region in regions:
        # Map INT→IN for stds lookup (algorithm uses INT, stds uses IN)
        std_name = 'IN' if region == 'INT' else region
        
        # Get the proper length from algorithm's wild type sequence
        wild_type = hiv_algorithm.stds.get(std_name, '')
        if not wild_type:
            pytest.skip(f"No wild type sequence for region {region}")
        
        # Create aminos matching the wild type length (all consensus A)
        aminos = [{'A': 1.0}] * len(wild_type)
        amino_lists.append(AminoList(region, aminos, None, 'HIV1B-seed', True))
    
    resistance_csv = StringIO()
    mutations_csv = StringIO()
    
    # Should not raise any errors
    try:
        write_resistance(amino_lists, resistance_csv, mutations_csv,
                        algorithms=asi_algorithms)
    except Exception as e:
        pytest.fail(f"write_resistance failed for algorithm regions: {e}")
    
    # Verify output was generated
    resistance_csv.seek(0)
    lines = resistance_csv.readlines()
    assert len(lines) > 1, "Should have generated resistance data"


def test_hcv_genotype_extraction():
    """Test HCV genotype extraction from seed names.
    
    Prevents: Incorrect genotype assignment causing wrong algorithm selection.
    """
    from micall.resistance.resistance import get_genotype
    
    # Test various HCV seeds
    assert get_genotype('HCV-1a') == '1A'
    assert get_genotype('HCV-1b') == '1B'
    assert get_genotype('HCV-2') == '2'
    assert get_genotype('HCV-6e') == '6E'
    
    # HIV should return None
    assert get_genotype('HIV1-B-FR-K03455-seed') is None
    assert get_genotype(None) is None


def test_reported_regions_constant_completeness():
    """REPORTED_REGIONS should be complete and not have typos.
    
    Prevents: Typos or missing regions in REPORTED_REGIONS constant.
    """
    from micall.resistance.resistance import REPORTED_REGIONS
    
    # Verify expected regions are present
    expected_hiv = {'PR', 'RT', 'IN', 'CA'}
    expected_hcv = {'NS3', 'NS5a', 'NS5b'}
    
    assert expected_hiv.issubset(REPORTED_REGIONS), (
        f"Missing expected HIV regions in REPORTED_REGIONS. "
        f"Missing: {expected_hiv - REPORTED_REGIONS}"
    )
    
    assert expected_hcv.issubset(REPORTED_REGIONS), (
        f"Missing expected HCV regions in REPORTED_REGIONS. "
        f"Missing: {expected_hcv - REPORTED_REGIONS}"
    )
