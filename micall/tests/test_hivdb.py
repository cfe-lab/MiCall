from io import StringIO
from unittest import TestCase

from micall.hivdb.hivdb import read_aminos, write_resistance, find_good_regions, select_reported_regions, AminoList


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


class FindGoodRegionsTest(TestCase):
    def test_good_coverage(self):
        coverage_scores_csv = StringIO("""\
project,region,on.score
PR,PR,4
RT,RT,4
INT,INT,4
""")
        original_regions = {'PR', 'RT', 'IN'}
        expected_regions = {'PR': ['PR', True],
                            'RT': ['RT', True],
                            'INT': ['IN', True]}

        good_regions = find_good_regions(original_regions, coverage_scores_csv)

        self.assertEqual(expected_regions, good_regions)

    def test_bad_coverage(self):
        coverage_scores_csv = StringIO("""\
project,region,on.score
PR,PR,4
RT,RT,3
INT,INT,4
""")
        original_regions = {'PR', 'RT', 'IN'}
        expected_regions = {'PR': ['PR', True],
                            'RT': ['RT', False],
                            'INT': ['IN', True]}

        good_regions = find_good_regions(original_regions, coverage_scores_csv)

        self.assertEqual(expected_regions, good_regions)

    def test_mixed_coverage(self):
        coverage_scores_csv = StringIO("""\
project,region,on.score
PR,PR,4
RT,RT,4
RT,RT,3
INT,INT,4
""")
        original_regions = {'PR', 'RT', 'IN'}
        expected_regions = {'PR': ['PR', True],
                            'RT': ['RT', True],
                            'INT': ['IN', True]}

        good_regions = find_good_regions(original_regions, coverage_scores_csv)

        self.assertEqual(expected_regions, good_regions)

    def test_other_regions(self):
        coverage_scores_csv = StringIO("""\
project,region,on.score
PR,PR,4
RT,RT,4
INT,INT,4
R4,R4,4
""")
        original_regions = {'PR', 'RT', 'IN'}
        expected_regions = {'PR': ['PR', True],
                            'RT': ['RT', True],
                            'INT': ['IN', True]}

        good_regions = find_good_regions(original_regions, coverage_scores_csv)

        self.assertEqual(expected_regions, good_regions)

    def test_missing_regions(self):
        coverage_scores_csv = StringIO("""\
project,region,on.score
RT,RT,4
INT,INT,4
""")
        original_regions = {'PR', 'RT', 'IN'}
        expected_regions = {'RT': ['RT', True],
                            'INT': ['IN', True]}

        good_regions = find_good_regions(original_regions, coverage_scores_csv)

        self.assertEqual(expected_regions, good_regions)


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
                                     'R1-seed')]

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
                                     'R1-seed')]

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
                                     'R1-seed')]

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
                                     'R1-seed')]

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
                                     'R1-seed'),
                           AminoList('R2',
                                     [{'K': 1.0}, {'F': 1.0}, {'C': 1.0}],
                                     'R2-seed')]

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
        reported_regions = {'R1': ['Region1', False],  # No data
                            'R2': ['Region2', True]}   # Others are skipped
        expected_aminos = [AminoList('R1', None, None),
                           AminoList('R2',
                                     [{'K': 1.0}, {'F': 1.0}, {'C': 1.0}],
                                     'R2-seed')]

        aminos = list(read_aminos(amino_csv,
                                  min_fraction,
                                  reported_regions=reported_regions))

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
        reported_regions = {'R1': ['Region1', False],  # No data
                            'R2': ['Region2', True]}   # Others are skipped
        expected_aminos = [AminoList('R2',
                                     [{'K': 1.0}, {'F': 1.0}, {'C': 1.0}],
                                     'R2-seed'),
                           AminoList('R1', None, None)]

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
                                     'R1-seed')]

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
                                     'R1-seed')]

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
                                     'R1-seed')]

        aminos = list(read_aminos(amino_csv, min_fraction))

        self.assertEqual(expected_aminos, aminos)


class WriteResistanceTest(TestCase):
    def test_simple(self):
        self.maxDiff = None
        aminos = [AminoList('RT',
                            [{'A': 1.0}] * 40 + [{'L': 1.0}] + [{}] * 399,
                            'HIV1-A1-CD-AM000053-seed')]
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

        write_resistance(aminos, resistance_csv, mutations_csv)

        self.assertEqual(expected_resistance, resistance_csv.getvalue())
        self.assertEqual(expected_mutations, mutations_csv.getvalue())

    def test_low_coverage(self):
        aminos = [AminoList('PR', None, None),
                  AminoList('RT',
                            [{'A': 1.0}] * 40 + [{'L': 1.0}] + [{}] * 399,
                            'HIV1-A1-CD-AM000053-seed')]
        resistance_csv = StringIO()
        mutations_csv = StringIO()
        expected_resistance = """\
region,drug_class,drug,drug_name,level,level_name,score,genotype
PR,PI,ATV/r,atazanavir/r,0,Insufficient data available,0.0,
PR,PI,DRV/r,darunavir/r,0,Insufficient data available,0.0,
PR,PI,FPV/r,fosamprenavir/r,0,Insufficient data available,0.0,
PR,PI,IDV/r,indinavir/r,0,Insufficient data available,0.0,
PR,PI,LPV/r,lopinavir/r,0,Insufficient data available,0.0,
PR,PI,NFV,nelfinavir,0,Insufficient data available,0.0,
PR,PI,SQV/r,saquinavir/r,0,Insufficient data available,0.0,
PR,PI,TPV/r,tipranavir/r,0,Insufficient data available,0.0,
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

        write_resistance(aminos, resistance_csv, mutations_csv)

        self.assertEqual(expected_resistance, resistance_csv.getvalue())
        self.assertEqual(expected_mutations, mutations_csv.getvalue())

    def test_mixture(self):
        self.maxDiff = None
        aminos = [AminoList('RT',
                            [{'A': 1.0}] * 40 + [{'L': 0.3}] + [{}] * 399,
                            'HIV1-A1-CD-AM000053-seed')]
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

        write_resistance(aminos, resistance_csv, mutations_csv)

        self.assertEqual(expected_resistance, resistance_csv.getvalue())
        self.assertEqual(expected_mutations, mutations_csv.getvalue())

    def test_hcv(self):
        self.maxDiff = None
        aminos = [AminoList('HCV6-EUHK2-NS5b',
                            [{'T': 1.0}] * 591,
                            'HCV-6f')]
        resistance_csv = StringIO()
        mutations_csv = StringIO()
        expected_resistance = """\
region,drug_class,drug,drug_name,level,level_name,score,genotype
NS5b,NS5b,SOF-EPC,SOF-EPC,3,Resistance Likely,40.0,6
"""
        expected_mutations = """\
drug_class,mutation,prevalence,genotype
NS5b,S282T,1.0,6
"""

        write_resistance(aminos, resistance_csv, mutations_csv)

        self.assertEqual(expected_resistance, resistance_csv.getvalue())
        self.assertEqual(expected_mutations, mutations_csv.getvalue())

    def test_multiple_mutations(self):
        self.maxDiff = None
        aminos = [AminoList('HCV1A-H77-NS5b',
                            [{'R': 1.0}] * 591,
                            'HCV-1a')]
        resistance_csv = StringIO()
        mutations_csv = StringIO()
        expected_resistance = """\
region,drug_class,drug,drug_name,level,level_name,score,genotype
NS5b,NS5b,DSV,DSV,3,Resistance Likely,16.0,1
NS5b,NS5b,SOF-EPC,SOF-EPC,3,Resistance Likely,8.0,1
NS5b,NS5b,SOF-HAR,SOF-HAR,3,Resistance Likely,8.0,1
"""
        expected_mutations = """\
drug_class,mutation,prevalence,genotype
NS5b,G307R,1.0,1
NS5b,S556R,1.0,1
NS5b,G558R,1.0,1
"""

        write_resistance(aminos, resistance_csv, mutations_csv)

        self.assertEqual(expected_resistance, resistance_csv.getvalue())
        self.assertEqual(expected_mutations, mutations_csv.getvalue())
