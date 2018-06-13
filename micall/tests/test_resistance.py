from csv import DictReader, DictWriter
from io import StringIO
from unittest import TestCase

from micall.resistance.resistance import read_aminos, write_resistance, select_reported_regions, AminoList, \
    filter_aminos, load_asi, get_genotype, combine_aminos, create_fail_writer


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


class CombineAminosTest(TestCase):
    def check_combination(self,
                          amino_csv,
                          midi_amino_csv,
                          expected_csv,
                          expected_fail_csv='seed,region,reason\n'):
        """ Combine two amino count CSV files, and check result.

        :param amino_csv: open file with main amino counts
        :param midi_amino_csv: open file with MIDI amino counts, or the same
            file as amino_csv, in which case it will be ignored.
        :param expected_csv: text expected to find in combined CSV
        :param expected_fail_csv: text of expected failure messages
        """
        self.maxDiff = 2000
        fail_csv = StringIO()
        fail_writer = create_fail_writer(fail_csv)
        combined_rows = list(combine_aminos(amino_csv, midi_amino_csv, fail_writer))

        expected_reader = DictReader(StringIO(expected_csv))
        expected_rows = list(expected_reader)
        self.assertEqual(format_rows(expected_reader.fieldnames, expected_rows),
                         format_rows(expected_reader.fieldnames, combined_rows))
        self.assertEqual(expected_fail_csv, fail_csv.getvalue())

    def test_simple(self):
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

        self.check_combination(amino_csv, amino_csv, expected_csv)

    def test_low_average(self):
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
        expected_fail_csv = """\
seed,region,reason
HCV-1a,HCV1A-H77-NS5b,low average coverage
"""

        self.check_combination(amino_csv, amino_csv, expected_csv, expected_fail_csv)

    def test_ns5b_window(self):
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
        expected_fail_csv = """\
seed,region,reason
HCV-1a,HCV1A-H77-NS5b,low average coverage
"""

        self.check_combination(amino_csv, amino_csv, expected_csv, expected_fail_csv)

    def test_ns3_and_ns5a_windows(self):
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
        expected_fail_csv = """\
seed,region,reason
HCV-1a,HCV1A-H77-NS3,low average coverage
HCV-1a,HCV1A-H77-NS5a,low average coverage
"""

        self.check_combination(amino_csv, amino_csv, expected_csv, expected_fail_csv)

    def test_ns3_and_ns5_missing_tails(self):
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
        expected_fail_csv = """\
seed,region,reason
HCV-1a,HCV1A-H77-NS3,not enough high-coverage amino acids
HCV-1a,HCV1A-H77-NS5a,not enough high-coverage amino acids
"""

        self.check_combination(amino_csv, amino_csv, expected_csv, expected_fail_csv)

    def test_ns3_and_ns5a_missing_heads(self):
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
        expected_fail_csv = """\
seed,region,reason
HCV-1a,HCV1A-H77-NS3,not enough high-coverage amino acids
HCV-1a,HCV1A-H77-NS5a,not enough high-coverage amino acids
"""

        self.check_combination(amino_csv, amino_csv, expected_csv, expected_fail_csv)

    def test_combine_ns5b(self):
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

        self.check_combination(amino_csv, midi_amino_csv, expected_csv)

    def test_combine_ns5b_multiple_seeds(self):
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
HCV-1b,HCV1B-H77-NS5b,15,1,558,0,0,0,0,0,0,0,0,20000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,20000
HCV-1b,HCV1B-H77-NS5b,15,1,559,0,0,0,0,0,0,0,0,20000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,20000
HCV-1b,HCV1B-H77-NS5b,15,4,560,0,0,0,0,20000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,20000
HCV-1b,HCV1B-H77-NS5b,15,7,561,20000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,20000
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

        self.check_combination(amino_csv, midi_amino_csv, expected_csv)

    def test_low_coverage_in_midi(self):
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
        expected_fail_csv = """\
seed,region,reason
HCV-1a,HCV1A-H77-NS5b,MIDI: low average coverage
"""

        self.check_combination(amino_csv, midi_amino_csv, expected_csv, expected_fail_csv)

    def test_midi_other_regions(self):
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

        self.check_combination(amino_csv, midi_amino_csv, expected_csv)

    def test_combine_ignores_early_midi(self):
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

        self.check_combination(amino_csv, midi_amino_csv, expected_csv)

    def test_combine_takes_higher_coverage(self):
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

        self.check_combination(amino_csv, midi_amino_csv, expected_csv)

    def test_combine_ignores_main_tail(self):
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

        self.check_combination(amino_csv, midi_amino_csv, expected_csv)

    def test_combine_ignores_main_tail_with_midi(self):
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

        self.check_combination(amino_csv, midi_amino_csv, expected_csv)

    def test_combine_ignores_main_tail_without_midi(self):
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

        self.check_combination(amino_csv, midi_amino_csv, expected_csv)

    def test_combine_midi_only(self):
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
        expected_fail_csv = """\
seed,region,reason
HCV-1a,HCV1A-H77-NS5b,low average coverage
"""

        self.check_combination(amino_csv, midi_amino_csv, expected_csv, expected_fail_csv)


class ReadAminosTest(TestCase):
    def test_simple(self):
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
                                     None)]

        aminos = list(read_aminos(amino_csv, min_fraction))

        self.assertEqual(expected_aminos, aminos)

    def test_hcv(self):
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
                                     '1A')]

        aminos = list(read_aminos(amino_csv, min_fraction))

        self.assertEqual(expected_aminos, aminos)

    def test_mixtures(self):
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
                                     None)]

        aminos = list(read_aminos(amino_csv, min_fraction))

        self.assertEqual(expected_aminos, aminos)

    def test_no_coverage(self):
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
                                     None)]

        aminos = list(read_aminos(amino_csv, min_fraction))

        self.assertEqual(expected_aminos, aminos)

    def test_insertion_before_coverage(self):
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
                                     None)]

        aminos = list(read_aminos(amino_csv, min_fraction))

        self.assertEqual(expected_aminos, aminos)

    def test_multiple_regions(self):
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
                                     None),
                           AminoList('R2',
                                     [{'K': 1.0}, {'F': 1.0}, {'C': 1.0}],
                                     None)]

        aminos = list(read_aminos(amino_csv, min_fraction))

        self.assertEqual(expected_aminos, aminos)

    def test_reported_regions(self):
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
                                     None),
                           AminoList('R2',
                                     [{'K': 1.0}, {'F': 1.0}, {'C': 1.0}],
                                     None)]

        aminos = list(read_aminos(amino_csv,
                                  min_fraction,
                                  reported_regions=reported_regions))

        self.assertEqual(expected_aminos, aminos)

    def test_low_coverage(self):
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
                                     None),
                           AminoList('R2',
                                     [{'K': 1.0}, {'F': 1.0}, {'C': 1.0}],
                                     None)]

        aminos = list(read_aminos(amino_csv,
                                  min_fraction,
                                  reported_regions=reported_regions,
                                  min_coverage=min_coverage))

        self.assertEqual(expected_aminos, aminos)

    def test_good_average_coverage(self):
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
                                     '1A')]

        aminos = list(read_aminos(amino_csv,
                                  min_fraction,
                                  min_coverage=min_coverage))

        self.assertEqual(expected_aminos, aminos)

    def test_low_average_coverage(self):
        amino_csv = DictReader(StringIO("""\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,g2p_overlap,coverage
HCV-1a,HCV-Foo-NS5a,15,1,100,0,0,0,0,0,0,0,0,505,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,505
HCV-1a,HCV-Foo-NS5a,15,4,101,0,0,0,0,403,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,403
"""))
        min_fraction = 0.2
        min_coverage = 9
        expected_aminos = []

        aminos = list(read_aminos(amino_csv,
                                  min_fraction,
                                  min_coverage=min_coverage))

        self.assertEqual(expected_aminos, aminos)

    def test_coverage_outside_window(self):
        amino_csv = DictReader(StringIO("""\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,g2p_overlap,coverage
HCV-1a,HCV-Foo-NS5a,15,1,101,0,0,0,0,0,0,0,0,500,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,500
HCV-1a,HCV-Foo-NS5a,15,4,102,0,0,0,0,500,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,500
"""))
        min_fraction = 0.2
        min_coverage = 9
        expected_aminos = []

        aminos = list(read_aminos(amino_csv,
                                  min_fraction,
                                  min_coverage=min_coverage))

        self.assertEqual(expected_aminos, aminos)

    def test_coverage_missing_at_end(self):
        amino_csv = DictReader(StringIO("""\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,g2p_overlap,coverage
HCV-1a,HCV-Foo-NS5a,15,1,99,0,0,0,0,0,0,0,0,505,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,505
HCV-1a,HCV-Foo-NS5a,15,4,100,0,0,0,0,404,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,404
"""))
        min_fraction = 0.2
        min_coverage = 9
        expected_aminos = []

        aminos = list(read_aminos(amino_csv,
                                  min_fraction,
                                  min_coverage=min_coverage))

        self.assertEqual(expected_aminos, aminos)

    def test_resistant_even_with_missing_midi(self):
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
                                     '1A')]

        aminos = list(read_aminos(amino_csv,
                                  min_fraction,
                                  min_coverage=min_coverage))

        self.assertEqual(expected_aminos, aminos)

    def test_missing_region(self):
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
                                     None)]

        aminos = list(read_aminos(amino_csv,
                                  min_fraction,
                                  reported_regions=reported_regions))

        self.assertEqual(expected_aminos, aminos)

    def test_deletions(self):
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
                                     None)]

        aminos = list(read_aminos(amino_csv, min_fraction))

        self.assertEqual(expected_aminos, aminos)

    def test_insertions(self):
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
                                     None)]

        aminos = list(read_aminos(amino_csv, min_fraction))

        self.assertEqual(expected_aminos, aminos)

    def test_stop_codons(self):
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
                                     None)]

        aminos = list(read_aminos(amino_csv, min_fraction))

        self.assertEqual(expected_aminos, aminos)

    def test_missing_position(self):
        amino_csv = DictReader(StringIO("""\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,g2p_overlap,coverage
R1-seed,R1,15,4,2,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
R1-seed,R1,15,7,3,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
"""))
        min_fraction = 0.2
        expected_aminos = [AminoList('R1',
                                     [{}, {'F': 1.0}, {'A': 1.0}],
                                     None)]

        aminos = list(read_aminos(amino_csv, min_fraction))

        self.assertEqual(expected_aminos, aminos)


class FilterAminosTest(TestCase):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.algorithms = load_asi()

    def test_all(self):
        all_aminos = [AminoList('PR',
                                [{'K': 1.0}, {'F': 1.0}, {'A': 1.0}],
                                None),
                      AminoList('RT',
                                [{'C': 1.0}, {'A': 1.0}, {'N': 1.0}],
                                None),
                      AminoList('INT',
                                [{'E': 1.0}, {'A': 1.0}, {'T': 1.0}],
                                None)]
        expected_aminos = all_aminos

        aminos = filter_aminos(all_aminos, self.algorithms)

        self.assertEqual(expected_aminos, aminos)

    def test_missing(self):
        all_aminos = [AminoList('RT',
                                [{'K': 1.0}, {'F': 1.0}, {'A': 1.0}],
                                None)]
        expected_aminos = [AminoList('RT',
                                     [{'K': 1.0}, {'F': 1.0}, {'A': 1.0}],
                                     None),
                           AminoList('INT', [{}]*288, None),
                           AminoList('PR', [{}]*99, None)]

        aminos = filter_aminos(all_aminos, self.algorithms)

        self.assertEqual(expected_aminos, aminos)

    def test_missing_with_genotypes(self):
        all_aminos = [AminoList('HCV1A-H77-NS5b',
                                [{'K': 1.0}, {'F': 1.0}, {'A': 1.0}],
                                '1A')]
        expected_aminos = [AminoList('HCV1A-H77-NS5b',
                                     [{'K': 1.0}, {'F': 1.0}, {'A': 1.0}],
                                     '1A'),
                           AminoList('HCV1A-H77-NS3', [{}]*631, '1A'),
                           AminoList('HCV1A-H77-NS5a', [{}]*448, '1A')]

        aminos = filter_aminos(all_aminos, self.algorithms)

        self.assertEqual(expected_aminos, aminos)

    def test_exclude(self):
        all_aminos = [AminoList('PR',
                                [{'K': 1.0}, {'F': 1.0}, {'A': 1.0}],
                                None),
                      AminoList('RT',
                                [{'K': 1.0}, {'F': 1.0}, {'A': 1.0}],
                                None),
                      AminoList('INT',
                                [{'K': 1.0}, {'F': 1.0}, {'A': 1.0}],
                                None),
                      AminoList('NS3',
                                [{}, {}, {}],
                                '1A')]
        expected_aminos = [AminoList('PR',
                                     [{'K': 1.0}, {'F': 1.0}, {'A': 1.0}],
                                     None),
                           AminoList('RT',
                                     [{'K': 1.0}, {'F': 1.0}, {'A': 1.0}],
                                     None),
                           AminoList('INT',
                                     [{'K': 1.0}, {'F': 1.0}, {'A': 1.0}],
                                     None)]

        aminos = filter_aminos(all_aminos, self.algorithms)

        self.assertEqual(expected_aminos, aminos)

    def test_no_good(self):
        all_aminos = [AminoList('PR', [{}, {}, {}], None),
                      AminoList('RT', [{}, {}, {}], None),
                      AminoList('HCV1A-H77-NS3', [{}, {}, {}], '1A')]
        expected_aminos = []

        aminos = filter_aminos(all_aminos, self.algorithms)

        self.assertEqual(expected_aminos, aminos)


class WriteResistanceTest(TestCase):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.algorithms = load_asi()

    def test_simple(self):
        rt_aminos = [{c: 1.0} for c in self.algorithms[None].stds['RT']]
        rt_aminos[40] = {'L': 1.0}  # pos 41
        aminos = [AminoList('RT', rt_aminos, None)]
        resistance_csv = StringIO()
        mutations_csv = StringIO()
        expected_resistance = """\
region,drug_class,drug,drug_name,level,level_name,score,genotype
RT,NRTI,3TC,lamivudine,1,Susceptible,0.0,
RT,NRTI,ABC,abacavir,1,Susceptible,5.0,
RT,NRTI,AZT,zidovudine,3,Low-Level Resistance,15.0,
RT,NRTI,D4T,stavudine,3,Low-Level Resistance,15.0,
RT,NRTI,DDI,didanosine,2,Susceptible,10.0,
RT,NRTI,FTC,emtricitabine,1,Susceptible,0.0,
RT,NRTI,TDF,tenofovir,1,Susceptible,5.0,
RT,NNRTI,EFV,efavirenz,1,Susceptible,0.0,
RT,NNRTI,ETR,etravirine,1,Susceptible,0.0,
RT,NNRTI,NVP,nevirapine,1,Susceptible,0.0,
RT,NNRTI,RPV,rilpivirine,1,Susceptible,0.0,
"""
        expected_mutations = """\
drug_class,mutation,prevalence,genotype
NRTI,M41L,1.0,
"""

        write_resistance(aminos, resistance_csv, mutations_csv, self.algorithms)

        self.assertEqual(expected_resistance, resistance_csv.getvalue())
        self.assertEqual(expected_mutations, mutations_csv.getvalue())

    def test_low_coverage(self):
        rt_aminos = [{c: 1.0} for c in self.algorithms[None].stds['RT']]
        rt_aminos[40] = {'L': 1.0}  # pos 41
        aminos = [AminoList('PR', [{}]*99, None),
                  AminoList('RT', rt_aminos, None)]
        resistance_csv = StringIO()
        mutations_csv = StringIO()
        expected_resistance = """\
region,drug_class,drug,drug_name,level,level_name,score,genotype
PR,PI,ATV/r,atazanavir/r,0,Sequence does not meet quality-control standards,0.0,
PR,PI,DRV/r,darunavir/r,0,Sequence does not meet quality-control standards,0.0,
PR,PI,FPV/r,fosamprenavir/r,0,Sequence does not meet quality-control standards,0.0,
PR,PI,IDV/r,indinavir/r,0,Sequence does not meet quality-control standards,0.0,
PR,PI,LPV/r,lopinavir/r,0,Sequence does not meet quality-control standards,0.0,
PR,PI,NFV,nelfinavir,0,Sequence does not meet quality-control standards,0.0,
PR,PI,SQV/r,saquinavir/r,0,Sequence does not meet quality-control standards,0.0,
PR,PI,TPV/r,tipranavir/r,0,Sequence does not meet quality-control standards,0.0,
RT,NRTI,3TC,lamivudine,1,Susceptible,0.0,
RT,NRTI,ABC,abacavir,1,Susceptible,5.0,
RT,NRTI,AZT,zidovudine,3,Low-Level Resistance,15.0,
RT,NRTI,D4T,stavudine,3,Low-Level Resistance,15.0,
RT,NRTI,DDI,didanosine,2,Susceptible,10.0,
RT,NRTI,FTC,emtricitabine,1,Susceptible,0.0,
RT,NRTI,TDF,tenofovir,1,Susceptible,5.0,
RT,NNRTI,EFV,efavirenz,1,Susceptible,0.0,
RT,NNRTI,ETR,etravirine,1,Susceptible,0.0,
RT,NNRTI,NVP,nevirapine,1,Susceptible,0.0,
RT,NNRTI,RPV,rilpivirine,1,Susceptible,0.0,
"""
        expected_mutations = """\
drug_class,mutation,prevalence,genotype
NRTI,M41L,1.0,
"""

        write_resistance(aminos, resistance_csv, mutations_csv, self.algorithms)

        self.assertEqual(expected_resistance, resistance_csv.getvalue())
        self.assertEqual(expected_mutations, mutations_csv.getvalue())

    def test_mixture(self):
        rt_aminos = [{c: 1.0} for c in self.algorithms[None].stds['RT']]
        rt_aminos[40] = {'L': 0.3, 'd': 0.7}  # pos 41
        aminos = [AminoList('RT', rt_aminos, None)]
        resistance_csv = StringIO()
        mutations_csv = StringIO()
        expected_resistance = """\
region,drug_class,drug,drug_name,level,level_name,score,genotype
RT,NRTI,3TC,lamivudine,1,Susceptible,0.0,
RT,NRTI,ABC,abacavir,1,Susceptible,5.0,
RT,NRTI,AZT,zidovudine,3,Low-Level Resistance,15.0,
RT,NRTI,D4T,stavudine,3,Low-Level Resistance,15.0,
RT,NRTI,DDI,didanosine,2,Susceptible,10.0,
RT,NRTI,FTC,emtricitabine,1,Susceptible,0.0,
RT,NRTI,TDF,tenofovir,1,Susceptible,5.0,
RT,NNRTI,EFV,efavirenz,1,Susceptible,0.0,
RT,NNRTI,ETR,etravirine,1,Susceptible,0.0,
RT,NNRTI,NVP,nevirapine,1,Susceptible,0.0,
RT,NNRTI,RPV,rilpivirine,1,Susceptible,0.0,
"""
        expected_mutations = """\
drug_class,mutation,prevalence,genotype
NRTI,M41L,0.3,
"""

        write_resistance(aminos, resistance_csv, mutations_csv, self.algorithms)

        self.assertEqual(expected_resistance, resistance_csv.getvalue())
        self.assertEqual(expected_mutations, mutations_csv.getvalue())

    def test_hcv(self):
        aminos = [AminoList('HCV6-EUHK2-NS5b',
                            [{'T': 1.0}] * 591,
                            '6')]
        resistance_csv = StringIO()
        mutations_csv = StringIO()
        expected_resistance = """\
region,drug_class,drug,drug_name,level,level_name,score,genotype
NS5b,NS5b,DSV,Dasabuvir,2,Not Indicated,0.0,6
NS5b,NS5b,SOF-EPC,Sofosbuvir in Epclusa,5,Resistance Likely,8.0,6
NS5b,NS5b,SOF-HAR,Sofosbuvir in Harvoni,2,Not Indicated,0.0,6
"""
        expected_mutations = """\
drug_class,mutation,prevalence,genotype
NS5b,L159T,1.0,6
NS5b,S282T,1.0,6
NS5b,C316T,1.0,6
NS5b,L320T,1.0,6
NS5b,V321T,1.0,6
"""

        write_resistance(aminos, resistance_csv, mutations_csv, self.algorithms)

        self.assertEqual(expected_resistance, resistance_csv.getvalue())
        self.assertEqual(expected_mutations, mutations_csv.getvalue())

    def test_hcv_missing_midi(self):
        ns5b_aminos = [{c: 1.0}
                       for c in self.algorithms[None].stds['HCV1B-Con1-NS5b']]
        ns5b_aminos[315] = {'Y': 1.0}  # C316Y mutation => Resistance Likely
        ns5b_aminos = ns5b_aminos[:336]  # Missing MIDI
        aminos = [AminoList('HCV1B-Con1-NS5b', ns5b_aminos, '1B')]
        resistance_csv = StringIO()
        mutations_csv = StringIO()
        expected_resistance = """\
region,drug_class,drug,drug_name,level,level_name,score,genotype
NS5b,NS5b,DSV,Dasabuvir,5,Resistance Likely,8.0,1B
NS5b,NS5b,SOF-EPC,Sofosbuvir in Epclusa,0,Sequence does not meet quality-control standards,0.0,1B
NS5b,NS5b,SOF-HAR,Sofosbuvir in Harvoni,0,Sequence does not meet quality-control standards,0.0,1B
"""
        expected_mutations = """\
drug_class,mutation,prevalence,genotype
NS5b,C316Y,1.0,1B
NS5b,V321V,1.0,1B
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
        aminos = [AminoList('HCV1A-H77-NS3', ns3_aminos, '1A')]
        resistance_csv = StringIO()
        mutations_csv = StringIO()
        expected_resistance = """\
region,drug_class,drug,drug_name,level,level_name,score,genotype
NS3,NS3,ASV,Asunaprevir,0,Sequence does not meet quality-control standards,0.0,1A
NS3,NS3,GLE,Glecaprevir,1,Likely Susceptible,0.0,1A
NS3,NS3,GZR,Grazoprevir,1,Likely Susceptible,0.0,1A
NS3,NS3,PTV,Paritaprevir,1,Likely Susceptible,0.0,1A
NS3,NS3,SPV,Simeprevir,0,Sequence does not meet quality-control standards,0.0,1A
NS3,NS3,VOX,Voxilaprevir,1,Likely Susceptible,0.0,1A
"""
        expected_mutations = """\
drug_class,mutation,prevalence,genotype
"""

        write_resistance(aminos, resistance_csv, mutations_csv, self.algorithms)

        self.assertEqual(expected_mutations, mutations_csv.getvalue())
        self.assertEqual(expected_resistance, resistance_csv.getvalue())

    def test_hcv_low_coverage(self):
        aminos = []
        resistance_csv = StringIO()
        mutations_csv = StringIO()
        expected_resistance = """\
region,drug_class,drug,drug_name,level,level_name,score,genotype
"""
        expected_mutations = """\
drug_class,mutation,prevalence,genotype
"""

        write_resistance(aminos, resistance_csv, mutations_csv, self.algorithms)

        self.assertEqual(expected_resistance, resistance_csv.getvalue())
        self.assertEqual(expected_mutations, mutations_csv.getvalue())

    def test_hcv_mostly_low_coverage(self):
        aminos = [AminoList('HCV1B-Con1-NS5b',
                            [{'T': 1.0}] + [{}] * 319 + [{'I': 1.0}] + [{}] * 200,
                            '1B')]
        resistance_csv = StringIO()
        mutations_csv = StringIO()
        expected_resistance = """\
region,drug_class,drug,drug_name,level,level_name,score,genotype
"""
        expected_mutations = """\
drug_class,mutation,prevalence,genotype
"""

        write_resistance(aminos, resistance_csv, mutations_csv, self.algorithms)

        self.assertEqual(expected_resistance, resistance_csv.getvalue())
        self.assertEqual(expected_mutations, mutations_csv.getvalue())


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
