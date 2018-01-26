from io import StringIO
from unittest import TestCase

from micall.hivdb.hivdb import read_aminos, write_resistance, select_reported_regions, AminoList, \
    filter_aminos, load_asi, get_genotype


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


class ReadAminosTest(TestCase):
    def test_simple(self):
        amino_csv = StringIO("""\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,g2p_overlap,coverage
R1-seed,R1,15,1,1,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
R1-seed,R1,15,4,2,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
R1-seed,R1,15,7,3,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
""")
        min_fraction = 0.2
        expected_aminos = [AminoList('R1',
                                     [{'K': 1.0}, {'F': 1.0}, {'A': 1.0}],
                                     None)]

        aminos = list(read_aminos(amino_csv, min_fraction))

        self.assertEqual(expected_aminos, aminos)

    def test_hcv(self):
        amino_csv = StringIO("""\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,g2p_overlap,coverage
HCV-1a,HCV1A-H77-NS5b,15,1,1,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
HCV-1a,HCV1A-H77-NS5b,15,4,2,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
HCV-1a,HCV1A-H77-NS5b,15,7,3,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
""")
        min_fraction = 0.2
        expected_aminos = [AminoList('HCV1A-H77-NS5b',
                                     [{'K': 1.0}, {'F': 1.0}, {'A': 1.0}],
                                     '1A')]

        aminos = list(read_aminos(amino_csv, min_fraction))

        self.assertEqual(expected_aminos, aminos)

    def test_mixtures(self):
        amino_csv = StringIO("""\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,g2p_overlap,coverage
R1-seed,R1,15,1,1,1,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10
R1-seed,R1,15,4,2,2,0,0,0,8,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10
R1-seed,R1,15,,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
""")
        min_fraction = 0.2
        expected_aminos = [AminoList('R1',
                                     [{'K': 0.9}, {'F': 0.8, 'A': 0.2}, {}],
                                     None)]

        aminos = list(read_aminos(amino_csv, min_fraction))

        self.assertEqual(expected_aminos, aminos)

    def test_no_coverage(self):
        amino_csv = StringIO("""\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,g2p_overlap,coverage
R1-seed,R1,15,1,1,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
R1-seed,R1,15,4,2,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
R1-seed,R1,15,,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
""")
        min_fraction = 0.2
        expected_aminos = [AminoList('R1',
                                     [{'K': 1.0}, {'F': 1.0}, {}],
                                     None)]

        aminos = list(read_aminos(amino_csv, min_fraction))

        self.assertEqual(expected_aminos, aminos)

    def test_insertion_before_coverage(self):
        amino_csv = StringIO("""\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,g2p_overlap,coverage
R1-seed,R1,15,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9,0,0,0
R1-seed,R1,15,4,2,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
R1-seed,R1,15,7,3,0,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
""")
        min_fraction = 0.2
        expected_aminos = [AminoList('R1',
                                     [{}, {'F': 1.0}, {'L': 1.0}],
                                     None)]

        aminos = list(read_aminos(amino_csv, min_fraction))

        self.assertEqual(expected_aminos, aminos)

    def test_multiple_regions(self):
        amino_csv = StringIO("""\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,g2p_overlap,coverage
R1-seed,R1,15,1,1,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
R1-seed,R1,15,4,2,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
R1-seed,R1,15,,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
R2-seed,R2,15,1,1,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
R2-seed,R2,15,4,2,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
R2-seed,R2,15,7,3,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
""")
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
        amino_csv = StringIO("""\
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
""")
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
        amino_csv = StringIO("""\
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
""")
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
        amino_csv = StringIO("""\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,g2p_overlap,coverage
HCV-1a,HCV-Foo-NS5a,15,1,100,0,0,0,0,0,0,0,0,505,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,505
HCV-1a,HCV-Foo-NS5a,15,4,101,0,0,0,0,404,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,404
""")
        min_fraction = 0.2
        min_coverage = 9
        expected_aminos = [AminoList('HCV-Foo-NS5a',
                                     [{'K': 1.0}, {'F': 1.0}],
                                     '1A')]

        aminos = list(read_aminos(amino_csv,
                                  min_fraction,
                                  min_coverage=min_coverage))

        self.assertEqual(expected_aminos, aminos)

    def test_low_average_coverage(self):
        amino_csv = StringIO("""\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,g2p_overlap,coverage
HCV-1a,HCV-Foo-NS5a,15,1,100,0,0,0,0,0,0,0,0,505,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,505
HCV-1a,HCV-Foo-NS5a,15,4,101,0,0,0,0,403,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,403
""")
        min_fraction = 0.2
        min_coverage = 9
        expected_aminos = []

        aminos = list(read_aminos(amino_csv,
                                  min_fraction,
                                  min_coverage=min_coverage))

        self.assertEqual(expected_aminos, aminos)

    def test_coverage_outside_window(self):
        amino_csv = StringIO("""\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,g2p_overlap,coverage
HCV-1a,HCV-Foo-NS5a,15,1,101,0,0,0,0,0,0,0,0,500,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,500
HCV-1a,HCV-Foo-NS5a,15,4,102,0,0,0,0,500,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,500
""")
        min_fraction = 0.2
        min_coverage = 9
        expected_aminos = []

        aminos = list(read_aminos(amino_csv,
                                  min_fraction,
                                  min_coverage=min_coverage))

        self.assertEqual(expected_aminos, aminos)

    def test_coverage_missing_at_end(self):
        amino_csv = StringIO("""\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,g2p_overlap,coverage
HCV-1a,HCV-Foo-NS5a,15,1,99,0,0,0,0,0,0,0,0,505,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,505
HCV-1a,HCV-Foo-NS5a,15,4,100,0,0,0,0,404,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,404
""")
        min_fraction = 0.2
        min_coverage = 9
        expected_aminos = []

        aminos = list(read_aminos(amino_csv,
                                  min_fraction,
                                  min_coverage=min_coverage))

        self.assertEqual(expected_aminos, aminos)

    def test_missing_region(self):
        amino_csv = StringIO("""\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,g2p_overlap,coverage
R2-seed,R2,15,1,1,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
R2-seed,R2,15,4,2,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
R2-seed,R2,15,7,3,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9
""")
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
        amino_csv = StringIO("""\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,g2p_overlap,coverage
R1-seed,R1,15,1,1,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,10
R1-seed,R1,15,4,2,0,0,0,0,8,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,10
R1-seed,R1,15,,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
""")
        min_fraction = 0.2
        expected_aminos = [AminoList('R1',
                                     [{'K': 0.9}, {'F': 0.8, 'd': 0.2}, {}],
                                     None)]

        aminos = list(read_aminos(amino_csv, min_fraction))

        self.assertEqual(expected_aminos, aminos)

    def test_insertions(self):
        amino_csv = StringIO("""\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,g2p_overlap,coverage
R1-seed,R1,15,1,1,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,9
R1-seed,R1,15,4,2,0,0,0,0,8,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,8
R1-seed,R1,15,,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
""")
        min_fraction = 0.2
        expected_aminos = [AminoList('R1',
                                     [{'K': 1.0}, {'F': 1.0, 'i': 0.25}, {}],
                                     None)]

        aminos = list(read_aminos(amino_csv, min_fraction))

        self.assertEqual(expected_aminos, aminos)

    def test_stop_codons(self):
        amino_csv = StringIO("""\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,g2p_overlap,coverage
R1-seed,R1,15,1,1,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,10
R1-seed,R1,15,4,2,0,0,0,0,8,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,10
R1-seed,R1,15,,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
""")
        min_fraction = 0.2
        expected_aminos = [AminoList('R1',
                                     [{'K': 0.9}, {'F': 0.8}, {}],
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
RT,NRTI,DDI,didanosine,2,Potential Low-Level Resistance,10.0,
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
RT,NRTI,DDI,didanosine,2,Potential Low-Level Resistance,10.0,
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
RT,NRTI,DDI,didanosine,2,Potential Low-Level Resistance,10.0,
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
NS5b,NS5b,DSV,DSV,2,Not Indicated,0.0,6
NS5b,NS5b,SOF-EPC,SOF-EPC,5,Resistance Likely,8.0,6
NS5b,NS5b,SOF-HAR,SOF-HAR,2,Not Indicated,0.0,6
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
        aminos = [AminoList('HCV1A-H77-NS5b',
                            [{'T': 1.0}] + [{}] * 590,
                            '1A')]
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
