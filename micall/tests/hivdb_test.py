from io import StringIO
from unittest import TestCase

from micall.hivdb.hivdb import read_aminos, write_resistance, find_good_regions


class FindGoodRegionsTest(TestCase):
    def test_good_coverage(self):
        coverage_scores_csv = StringIO("""\
project,region,on.score
R1,R1,4
R2,R2,4
R3,R3,4
""")
        original_regions = {'R1': 'Region 1',
                            'R2': 'Region 2',
                            'R3': 'Region 3'}
        expected_regions = {'R1': ['Region 1', True],
                            'R2': ['Region 2', True],
                            'R3': ['Region 3', True]}

        good_regions = find_good_regions(original_regions, coverage_scores_csv)

        self.assertEqual(expected_regions, good_regions)

    def test_bad_coverage(self):
        coverage_scores_csv = StringIO("""\
project,region,on.score
R1,R1,4
R2,R2,3
R3,R3,4
""")
        original_regions = {'R1': 'Region 1',
                            'R2': 'Region 2',
                            'R3': 'Region 3'}
        expected_regions = {'R1': ['Region 1', True],
                            'R2': ['Region 2', False],
                            'R3': ['Region 3', True]}

        good_regions = find_good_regions(original_regions, coverage_scores_csv)

        self.assertEqual(expected_regions, good_regions)

    def test_mixed_coverage(self):
        coverage_scores_csv = StringIO("""\
project,region,on.score
R1,R1,4
R2,R2,4
R2,R2,3
R3,R3,4
""")
        original_regions = {'R1': 'Region 1',
                            'R2': 'Region 2',
                            'R3': 'Region 3'}
        expected_regions = {'R1': ['Region 1', True],
                            'R2': ['Region 2', True],
                            'R3': ['Region 3', True]}

        good_regions = find_good_regions(original_regions, coverage_scores_csv)

        self.assertEqual(expected_regions, good_regions)

    def test_other_regions(self):
        coverage_scores_csv = StringIO("""\
project,region,on.score
R1,R1,4
R2,R2,4
R3,R3,4
R4,R4,4
""")
        original_regions = {'R1': 'Region 1',
                            'R2': 'Region 2',
                            'R3': 'Region 3'}
        expected_regions = {'R1': ['Region 1', True],
                            'R2': ['Region 2', True],
                            'R3': ['Region 3', True]}

        good_regions = find_good_regions(original_regions, coverage_scores_csv)

        self.assertEqual(expected_regions, good_regions)

    def test_missing_regions(self):
        coverage_scores_csv = StringIO("""\
project,region,on.score
R2,R2,4
R3,R3,4
""")
        original_regions = {'R1': 'Region 1',
                            'R2': 'Region 2',
                            'R3': 'Region 3'}
        expected_regions = {'R1': ['Region 1', False],
                            'R2': ['Region 2', True],
                            'R3': ['Region 3', True]}

        good_regions = find_good_regions(original_regions, coverage_scores_csv)

        self.assertEqual(expected_regions, good_regions)


class ReadAminosTest(TestCase):
    def test_simple(self):
        amino_csv = StringIO("""\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,g2p_overlap
R1-seed,R1,15,1,1,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
R1-seed,R1,15,4,2,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
R1-seed,R1,15,7,3,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
""")
        min_fraction = 0.2
        expected_aminos = [('R1', [['K'], ['F'], ['A']])]

        aminos = list(read_aminos(amino_csv, min_fraction))

        self.assertEqual(expected_aminos, aminos)

    def test_mixtures(self):
        amino_csv = StringIO("""\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,g2p_overlap
R1-seed,R1,15,1,1,1,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
R1-seed,R1,15,4,2,2,0,0,0,8,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
R1-seed,R1,15,,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
""")
        min_fraction = 0.2
        expected_aminos = [('R1', [['K'], ['A', 'F'], []])]

        aminos = list(read_aminos(amino_csv, min_fraction))

        self.assertEqual(expected_aminos, aminos)

    def test_no_coverage(self):
        amino_csv = StringIO("""\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,g2p_overlap
R1-seed,R1,15,1,1,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
R1-seed,R1,15,4,2,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
R1-seed,R1,15,,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
""")
        min_fraction = 0.2
        expected_aminos = [('R1', [['K'], ['F'], []])]

        aminos = list(read_aminos(amino_csv, min_fraction))

        self.assertEqual(expected_aminos, aminos)

    def test_multiple_regions(self):
        amino_csv = StringIO("""\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,g2p_overlap
R1-seed,R1,15,1,1,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
R1-seed,R1,15,4,2,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
R1-seed,R1,15,,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
R2-seed,R2,15,1,1,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
R2-seed,R2,15,4,2,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
R2-seed,R2,15,7,3,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
""")
        min_fraction = 0.2
        expected_aminos = [('R1', [['K'], ['F'], []]),
                           ('R2', [['K'], ['F'], ['C']])]

        aminos = list(read_aminos(amino_csv, min_fraction))

        self.assertEqual(expected_aminos, aminos)

    def test_reported_regions(self):
        amino_csv = StringIO("""\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,g2p_overlap
R1-seed,R1,15,1,1,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
R1-seed,R1,15,4,2,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
R1-seed,R1,15,,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
R2-seed,R2,15,1,1,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
R2-seed,R2,15,4,2,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
R2-seed,R2,15,7,3,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
R3-seed,R3,15,1,1,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
R3-seed,R3,15,4,2,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
R3-seed,R3,15,7,3,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
""")
        min_fraction = 0.2
        reported_regions = {'R1': ['Region1', False],  # No data
                            'R2': ['Region2', True]}   # Others are skipped
        expected_aminos = [('Region1', None),
                           ('Region2', [['K'], ['F'], ['C']])]

        aminos = list(read_aminos(amino_csv,
                                  min_fraction,
                                  reported_regions=reported_regions))

        self.assertEqual(expected_aminos, aminos)

    def test_missing_region(self):
        amino_csv = StringIO("""\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,g2p_overlap
R2-seed,R2,15,1,1,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
R2-seed,R2,15,4,2,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
R2-seed,R2,15,7,3,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
""")
        min_fraction = 0.2
        reported_regions = {'R1': ['Region1', False],  # No data
                            'R2': ['Region2', True]}   # Others are skipped
        expected_aminos = [('Region2', [['K'], ['F'], ['C']]),
                           ('Region1', None)]

        aminos = list(read_aminos(amino_csv,
                                  min_fraction,
                                  reported_regions=reported_regions))

        self.assertEqual(expected_aminos, aminos)

    def test_deletions(self):
        amino_csv = StringIO("""\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,g2p_overlap
R1-seed,R1,15,1,1,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0
R1-seed,R1,15,4,2,0,0,0,0,8,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0
R1-seed,R1,15,,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
""")
        min_fraction = 0.2
        expected_aminos = [('R1', [['K'], ['F', 'd'], []])]

        aminos = list(read_aminos(amino_csv, min_fraction))

        self.assertEqual(expected_aminos, aminos)

    def test_insertions(self):
        amino_csv = StringIO("""\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,g2p_overlap
R1-seed,R1,15,1,1,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0
R1-seed,R1,15,4,2,0,0,0,0,8,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0
R1-seed,R1,15,,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
""")
        min_fraction = 0.2
        expected_aminos = [('R1', [['K'], ['F', 'i'], []])]

        aminos = list(read_aminos(amino_csv, min_fraction))

        self.assertEqual(expected_aminos, aminos)

    def test_stop_codons(self):
        amino_csv = StringIO("""\
seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,g2p_overlap
R1-seed,R1,15,1,1,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0
R1-seed,R1,15,4,2,0,0,0,0,8,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0
R1-seed,R1,15,,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
""")
        min_fraction = 0.2
        expected_aminos = [('R1', [['K'], ['F'], []])]

        aminos = list(read_aminos(amino_csv, min_fraction))

        self.assertEqual(expected_aminos, aminos)


class WriteResistanceTest(TestCase):
    def test_simple(self):
        self.maxDiff = None
        aminos = [('RT', [['A']] * 40 + [['L']])]
        resistance_csv = StringIO()
        mutations_csv = StringIO()
        expected_resistance = """\
region,drug_class,drug,drug_name,level,level_name,score
RT,NRTI,3TC,lamivudine,1,Susceptible,0.0
RT,NRTI,ABC,abacavir,1,Susceptible,5.0
RT,NRTI,AZT,zidovudine,3,Low-Level Resistance,15.0
RT,NRTI,D4T,stavudine,3,Low-Level Resistance,15.0
RT,NRTI,DDI,didanosine,2,Potential Low-Level Resistance,10.0
RT,NRTI,FTC,emtricitabine,1,Susceptible,0.0
RT,NRTI,TDF,tenofovir,1,Susceptible,5.0
RT,NNRTI,EFV,efavirenz,1,Susceptible,0.0
RT,NNRTI,ETR,etravirine,1,Susceptible,0.0
RT,NNRTI,NVP,nevirapine,1,Susceptible,0.0
RT,NNRTI,RPV,rilpivirine,1,Susceptible,0.0
"""
        expected_mutations = """\
drug_class,mutation
NRTI,M41L
"""

        write_resistance(aminos, resistance_csv, mutations_csv)

        self.assertEqual(expected_resistance, resistance_csv.getvalue())
        self.assertEqual(expected_mutations, mutations_csv.getvalue())

    def test_low_coverage(self):
        aminos = [('PR', None),
                  ('RT', [['A']] * 40 + [['L']])]
        resistance_csv = StringIO()
        mutations_csv = StringIO()
        expected_resistance = """\
region,drug_class,drug,drug_name,level,level_name,score
PR,PI,ATV/r,atazanavir/r,0,Insufficient data available,0.0
PR,PI,DRV/r,darunavir/r,0,Insufficient data available,0.0
PR,PI,FPV/r,fosamprenavir/r,0,Insufficient data available,0.0
PR,PI,IDV/r,indinavir/r,0,Insufficient data available,0.0
PR,PI,LPV/r,lopinavir/r,0,Insufficient data available,0.0
PR,PI,NFV,nelfinavir,0,Insufficient data available,0.0
PR,PI,SQV/r,saquinavir/r,0,Insufficient data available,0.0
PR,PI,TPV/r,tipranavir/r,0,Insufficient data available,0.0
RT,NRTI,3TC,lamivudine,1,Susceptible,0.0
RT,NRTI,ABC,abacavir,1,Susceptible,5.0
RT,NRTI,AZT,zidovudine,3,Low-Level Resistance,15.0
RT,NRTI,D4T,stavudine,3,Low-Level Resistance,15.0
RT,NRTI,DDI,didanosine,2,Potential Low-Level Resistance,10.0
RT,NRTI,FTC,emtricitabine,1,Susceptible,0.0
RT,NRTI,TDF,tenofovir,1,Susceptible,5.0
RT,NNRTI,EFV,efavirenz,1,Susceptible,0.0
RT,NNRTI,ETR,etravirine,1,Susceptible,0.0
RT,NNRTI,NVP,nevirapine,1,Susceptible,0.0
RT,NNRTI,RPV,rilpivirine,1,Susceptible,0.0
"""
        expected_mutations = """\
drug_class,mutation
NRTI,M41L
"""

        write_resistance(aminos, resistance_csv, mutations_csv)

        self.assertEqual(expected_resistance, resistance_csv.getvalue())
        self.assertEqual(expected_mutations, mutations_csv.getvalue())
