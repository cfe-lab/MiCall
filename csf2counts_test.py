import StringIO
import unittest

import csf2counts

class Csf2CountsTest(unittest.TestCase):
    def setUp(self):
        self.writer = csf2counts.AminoFrequencyWriter(aafile=StringIO.StringIO(),
                                                      sample_name = 'E1234_S1',
                                                      refseqs = {'R1': 'XXX',
                                                                 'R2': 'YYYY'})
        self.amino_counts = {}
        self.qindex_to_refcoord = {}
        self.inserts = []
        self.region = 'R1'
        self.qcut = 15
        
    def testUnmappedRegion(self):
        expected_text = ''
        
        self.writer.write(region='R1',
                          qcut=15,
                          qindex_to_refcoord={},
                          amino_counts={},
                          inserts=[])
        
        self.assertEqual(self.writer.aafile.getvalue(), expected_text)
        
    def testWriteAminosFullRegion(self):
        expected_text = """\
sample,region,q-cutoff,query.aa.pos,refseq.aa.pos,A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*
E1234_S1,R1,15,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0
E1234_S1,R1,15,1,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0
E1234_S1,R1,15,2,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0
"""
        
        self.writer.write(region='R1',
                          qcut=15,
                          qindex_to_refcoord={0:0, 1:1, 2:2},
                          amino_counts={0: {'Q': 1}, 1: {'R': 1}, 2: {'S': 1}},
                          inserts=[])
        
        self.assertEqual(self.writer.aafile.getvalue(), expected_text)
        
    def testWriteAminosEndOfRegion(self):
        expected_text = """\
sample,region,q-cutoff,query.aa.pos,refseq.aa.pos,A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*
E1234_S1,R1,15,,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
E1234_S1,R1,15,1,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0
E1234_S1,R1,15,2,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0
"""

        self.writer.write(region='R1',
                          qcut=15,
                          qindex_to_refcoord={0:1, 1:2},
                          amino_counts={1: {'R': 1}, 2: {'S': 1}},
                          inserts=[])
        
        self.assertEqual(self.writer.aafile.getvalue(), expected_text)
        
    def testWriteAminosStartOfRegion(self):
        expected_text = """\
sample,region,q-cutoff,query.aa.pos,refseq.aa.pos,A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*
E1234_S1,R1,15,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0
E1234_S1,R1,15,1,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0
E1234_S1,R1,15,,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
"""

        self.writer.write(region='R1',
                          qcut=15,
                          qindex_to_refcoord={0:0, 1:1},
                          amino_counts={0: {'Q': 1}, 1: {'R': 1}},
                          inserts=[])
        
        self.assertEqual(self.writer.aafile.getvalue(), expected_text)
        
    def testWriteAminosWithInsert(self):
        expected_text = """\
sample,region,q-cutoff,query.aa.pos,refseq.aa.pos,A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*
E1234_S1,R1,15,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0
E1234_S1,R1,15,2,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0
E1234_S1,R1,15,3,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0
"""

        self.writer.write(region='R1',
                          qcut=15,
                          qindex_to_refcoord={0:0, 2:1, 3:2},
                          amino_counts={0: {'Q': 1},
                                        1: {'F': 1}, # This is the insert
                                        2: {'R': 1},
                                        3: {'S': 1}},
                          inserts=[1])
        
        self.assertEqual(self.writer.aafile.getvalue(), expected_text)
        
    def testWriteAminosWithDeletion(self):
        expected_text = """\
sample,region,q-cutoff,query.aa.pos,refseq.aa.pos,A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*
E1234_S1,R1,15,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0
E1234_S1,R1,15,,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
E1234_S1,R1,15,1,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0
"""

        self.writer.write(region='R1',
                          qcut=15,
                          qindex_to_refcoord={0:0, 1:2},
                          amino_counts={0: {'Q': 1}, 1: {'S': 1}},
                          inserts=[])
        
        self.assertEqual(self.writer.aafile.getvalue(), expected_text)
        
    def testWriteAminosWithGap(self):
        expected_text = """\
sample,region,q-cutoff,query.aa.pos,refseq.aa.pos,A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*
E1234_S1,R1,15,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0
E1234_S1,R1,15,,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
E1234_S1,R1,15,2,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0
"""
        
        self.writer.write(region='R1',
                          qcut=15,
                          qindex_to_refcoord={0:0, 2:2},
                          amino_counts={0: {'Q': 1}, 2: {'S': 1}},
                          inserts=[])
        
        self.assertEqual(self.writer.aafile.getvalue(), expected_text)

        
    def testWriteAminosTwoRegions(self):
        expected_text = """\
sample,region,q-cutoff,query.aa.pos,refseq.aa.pos,A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*
E1234_S1,R1,15,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0
E1234_S1,R1,15,1,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0
E1234_S1,R1,15,2,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0
E1234_S1,R2,15,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0
E1234_S1,R2,15,1,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0
E1234_S1,R2,15,2,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0
E1234_S1,R2,15,3,4,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0
"""
        
        self.writer.write(region='R1',
                          qcut=15,
                          qindex_to_refcoord={0:0, 1:1, 2:2},
                          amino_counts={0: {'Q': 1}, 1: {'R': 1}, 2: {'S': 1}},
                          inserts=[])
        self.writer.write(region='R2',
                          qcut=15,
                          qindex_to_refcoord={0:0, 1:1, 2:2, 3:3},
                          amino_counts={0: {'T': 1}, 1: {'S': 1}, 2: {'R': 1}, 3: {'Q': 1}},
                          inserts=[])
        
        self.assertEqual(self.writer.aafile.getvalue(), expected_text)
