from io import StringIO
from unittest import TestCase

from micall.hivdb.hivdb import read_aminos, write_resistance


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
""")
        min_fraction = 0.2
        reported_regions = {'R2': 'Region2'}  # Others are skipped
        expected_aminos = [('Region2', [['K'], ['F'], ['C']])]

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
        aminos = [('RT', [['A']] * 40 + [['L']])]
        resistance_csv = StringIO()
        mutations_csv = StringIO()
        expected_resistance = """\
region,drug,drug_name,level,level_name,score
RT,3TC,lamivudine,1,Susceptible,0.0
RT,ABC,abacavir,1,Susceptible,5.0
RT,AZT,zidovudine,3,Low-Level Resistance,15.0
RT,D4T,stavudine,3,Low-Level Resistance,15.0
RT,DDI,didanosine,2,Potential Low-Level Resistance,10.0
RT,FTC,emtricitabine,1,Susceptible,0.0
RT,TDF,tenofovir,1,Susceptible,5.0
RT,EFV,efavirenz,1,Susceptible,0.0
RT,ETR,etravirine,1,Susceptible,0.0
RT,NVP,nevirapine,1,Susceptible,0.0
RT,RPV,rilpivirine,1,Susceptible,0.0
"""
        expected_mutations = """\
region,mutation
RT,M41L
"""

        write_resistance(aminos, resistance_csv, mutations_csv)

        self.assertEqual(expected_resistance, resistance_csv.getvalue())
        self.assertEqual(expected_mutations, mutations_csv.getvalue())
