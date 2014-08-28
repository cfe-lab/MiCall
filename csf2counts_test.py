import StringIO
import unittest

import csf2counts

class AminoFrequencyWriterTest(unittest.TestCase):
    def setUp(self):
        self.writer = csf2counts.AminoFrequencyWriter(aafile=StringIO.StringIO(),
                                                      refseqs = {'R1': 'XXX',
                                                                 'R2': 'YYYY'})
        
    def testUnmappedRegion(self):
        expected_text = """\
sample,region,q-cutoff,query.aa.pos,refseq.aa.pos,A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*
"""
        
        self.writer.write(sample_name = 'E1234_S1',
                          region='R1',
                          qcut=15,
                          qindex_to_refcoord={},
                          amino_counts={},
                          inserts=[])
        
        self.assertMultiLineEqual(expected_text, self.writer.aafile.getvalue())
        
    def testWriteAminosFullRegion(self):
        expected_text = """\
sample,region,q-cutoff,query.aa.pos,refseq.aa.pos,A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*
E1234_S1,R1,15,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0
E1234_S1,R1,15,1,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0
E1234_S1,R1,15,2,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0
"""
        
        self.writer.write(sample_name = 'E1234_S1',
                          region='R1',
                          qcut=15,
                          qindex_to_refcoord={0:0, 1:1, 2:2},
                          amino_counts={0: {'Q': 1}, 1: {'R': 1}, 2: {'S': 1}},
                          inserts=[])
        
        self.assertMultiLineEqual(expected_text, self.writer.aafile.getvalue())
        
    def testWriteAminosEndOfRegion(self):
        expected_text = """\
sample,region,q-cutoff,query.aa.pos,refseq.aa.pos,A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*
E1234_S1,R1,15,,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
E1234_S1,R1,15,1,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0
E1234_S1,R1,15,2,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0
"""

        self.writer.write(sample_name = 'E1234_S1',
                          region='R1',
                          qcut=15,
                          qindex_to_refcoord={0:1, 1:2},
                          amino_counts={1: {'R': 1}, 2: {'S': 1}},
                          inserts=[])
        
        self.assertMultiLineEqual(expected_text, self.writer.aafile.getvalue())
        
    def testWriteAminosStartOfRegion(self):
        expected_text = """\
sample,region,q-cutoff,query.aa.pos,refseq.aa.pos,A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*
E1234_S1,R1,15,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0
E1234_S1,R1,15,1,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0
E1234_S1,R1,15,,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
"""

        self.writer.write(sample_name = 'E1234_S1',
                          region='R1',
                          qcut=15,
                          qindex_to_refcoord={0:0, 1:1},
                          amino_counts={0: {'Q': 1}, 1: {'R': 1}},
                          inserts=[])
        
        self.assertMultiLineEqual(expected_text, self.writer.aafile.getvalue())
        
    def testWriteAminosWithInsert(self):
        expected_text = """\
sample,region,q-cutoff,query.aa.pos,refseq.aa.pos,A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*
E1234_S1,R1,15,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0
E1234_S1,R1,15,2,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0
E1234_S1,R1,15,3,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0
"""

        self.writer.write(sample_name = 'E1234_S1',
                          region='R1',
                          qcut=15,
                          qindex_to_refcoord={0:0, 2:1, 3:2},
                          amino_counts={0: {'Q': 1},
                                        1: {'F': 1}, # This is the insert
                                        2: {'R': 1},
                                        3: {'S': 1}},
                          inserts=[1])
        
        self.assertMultiLineEqual(expected_text, self.writer.aafile.getvalue())
        
    def testWriteAminosWithDeletion(self):
        expected_text = """\
sample,region,q-cutoff,query.aa.pos,refseq.aa.pos,A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*
E1234_S1,R1,15,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0
E1234_S1,R1,15,,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
E1234_S1,R1,15,1,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0
"""

        self.writer.write(sample_name = 'E1234_S1',
                          region='R1',
                          qcut=15,
                          qindex_to_refcoord={0:0, 1:2},
                          amino_counts={0: {'Q': 1}, 1: {'S': 1}},
                          inserts=[])
        
        self.assertMultiLineEqual(expected_text, self.writer.aafile.getvalue())
        
    def testWriteAminosWithGap(self):
        expected_text = """\
sample,region,q-cutoff,query.aa.pos,refseq.aa.pos,A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*
E1234_S1,R1,15,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0
E1234_S1,R1,15,,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
E1234_S1,R1,15,2,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0
"""
        
        self.writer.write(sample_name = 'E1234_S1',
                          region='R1',
                          qcut=15,
                          qindex_to_refcoord={0:0, 2:2},
                          amino_counts={0: {'Q': 1}, 2: {'S': 1}},
                          inserts=[])
        
        self.assertMultiLineEqual(expected_text, self.writer.aafile.getvalue())

        
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
        
        self.writer.write(sample_name = 'E1234_S1',
                          region='R1',
                          qcut=15,
                          qindex_to_refcoord={0:0, 1:1, 2:2},
                          amino_counts={0: {'Q': 1}, 1: {'R': 1}, 2: {'S': 1}},
                          inserts=[])
        self.writer.write(sample_name = 'E1234_S1',
                          region='R2',
                          qcut=15,
                          qindex_to_refcoord={0:0, 1:1, 2:2, 3:3},
                          amino_counts={0: {'T': 1}, 1: {'S': 1}, 2: {'R': 1}, 3: {'Q': 1}},
                          inserts=[])
        
        self.assertMultiLineEqual(expected_text, self.writer.aafile.getvalue())

class NucleotideFrequencyWriterTest(unittest.TestCase):
    def setUp(self):
        self.writer = csf2counts.NucleotideFrequencyWriter(
            nucfile=StringIO.StringIO(),
            amino_ref_seqs = {'R1': 'XXX',
                              'R2': 'YYYY'})
        
    def testNoWrites(self):
        expected_text = """\
sample,region,q-cutoff,query.nuc.pos,refseq.nuc.pos,A,C,G,T
"""
        
        self.assertMultiLineEqual(expected_text, self.writer.nucfile.getvalue())
        
    def testWithoutReference(self):
        # A region with no amino reference, like HLA.
        expected_text = """\
sample,region,q-cutoff,query.nuc.pos,refseq.nuc.pos,A,C,G,T
E1234_S1,R3,15,1,,1,2,0,0
E1234_S1,R3,15,2,,0,0,0,3
"""
        
        self.writer.write(sample_name = 'E1234_S1',
                          region='R3',
                          qcut=15,
                          nuc_counts={0: {'A': 1, 'C': 2},
                                      1: {'T': 3}})
        
        self.assertMultiLineEqual(expected_text, self.writer.nucfile.getvalue())
        
    def testWithReference(self):
        expected_text = """\
sample,region,q-cutoff,query.nuc.pos,refseq.nuc.pos,A,C,G,T
E1234_S1,R1,15,1,1,1,2,0,0
E1234_S1,R1,15,2,2,0,0,0,3
"""
        
        self.writer.write(sample_name = 'E1234_S1',
                          region='R1',
                          qcut=15,
                          nuc_counts={0: {'A': 1, 'C': 2},
                                      1: {'T': 3}},
                          qindex_to_refcoord={0: 0},
                          min_offset=0)
        
        self.assertMultiLineEqual(expected_text, self.writer.nucfile.getvalue())
        
    def testOffsetWithoutReference(self):
        # A region with no amino reference, like HLA.
        expected_text = """\
sample,region,q-cutoff,query.nuc.pos,refseq.nuc.pos,A,C,G,T
E1234_S1,R3,15,6,,1,2,0,0
E1234_S1,R3,15,7,,0,0,0,3
"""
        
        self.writer.write(sample_name = 'E1234_S1',
                          region='R3',
                          qcut=15,
                          nuc_counts={5: {'A': 1, 'C': 2},
                                      6: {'T': 3}})
        
        self.assertMultiLineEqual(expected_text, self.writer.nucfile.getvalue())
        
    def testOffsetWithReference(self):
        # No coverage at the start of the reference, and first codon starts
        # with its second nucleotide. The partial codon is output.
        expected_text = """\
sample,region,q-cutoff,query.nuc.pos,refseq.nuc.pos,A,C,G,T
E1234_S1,R1,15,7,7,0,0,1,0
E1234_S1,R1,15,8,8,0,0,1,0
E1234_S1,R1,15,9,9,1,0,0,0
"""
        
        self.writer.write(sample_name = 'E1234_S1',
                          region='R1',
                          qcut=15,
                          nuc_counts={4: {'T': 1},
                                      5: {'G': 1},
                                      6: {'G': 1},
                                      7: {'G': 1},
                                      8: {'A': 1}},
                          qindex_to_refcoord={1: 2},
                          min_offset=4)
        
        self.assertMultiLineEqual(expected_text, self.writer.nucfile.getvalue())
        
    def testGapWithoutReference(self):
        expected_text = """\
sample,region,q-cutoff,query.nuc.pos,refseq.nuc.pos,A,C,G,T
E1234_S1,R3,15,1,,1,2,0,0
E1234_S1,R3,15,2,,0,0,0,3
E1234_S1,R3,15,3,,0,0,0,0
E1234_S1,R3,15,4,,0,0,3,0
"""
        
        self.writer.write(sample_name = 'E1234_S1',
                          region='R3',
                          qcut=15,
                          nuc_counts={0: {'A': 1, 'C': 2},
                                      1: {'T': 3},
                                      3: {'G': 3}})
        
        self.assertMultiLineEqual(expected_text, self.writer.nucfile.getvalue())
        
    def testGapWithReference(self):
        expected_text = """\
sample,region,q-cutoff,query.nuc.pos,refseq.nuc.pos,A,C,G,T
E1234_S1,R1,15,1,1,1,2,0,0
E1234_S1,R1,15,2,2,0,0,0,3
E1234_S1,R1,15,4,4,0,0,3,0
"""
        
        self.writer.write(sample_name = 'E1234_S1',
                          region='R1',
                          qcut=15,
                          nuc_counts={0: {'A': 1, 'C': 2},
                                      1: {'T': 3},
                                      3: {'G': 3}},
                          qindex_to_refcoord={0: 0, 1: 1},
                          min_offset=0)
        
        self.assertMultiLineEqual(expected_text, self.writer.nucfile.getvalue())
        
    def testInsertionWithReference(self):
        expected_text = """\
sample,region,q-cutoff,query.nuc.pos,refseq.nuc.pos,A,C,G,T
E1234_S1,R1,15,1,1,1,2,0,0
E1234_S1,R1,15,2,2,0,0,0,3
E1234_S1,R1,15,3,3,0,2,0,0
E1234_S1,R1,15,4,7,0,0,3,0
"""
        
        self.writer.write(sample_name = 'E1234_S1',
                          region='R1',
                          qcut=15,
                          nuc_counts={0: {'A': 1, 'C': 2},
                                      1: {'T': 3},
                                      2: {'C': 2},
                                      3: {'G': 3}},
                          qindex_to_refcoord={0: 0, 1: 2},
                          min_offset=0)
        
        self.assertMultiLineEqual(expected_text, self.writer.nucfile.getvalue())

class CoordinateMapTest(unittest.TestCase):
    def testStraightMapping(self):
        query_sequence =     'CTRPNNN'
        reference_sequence = 'CTRPNNN'
        expected_inserts = []
        expected_mapping = {0:0, 1:1, 2:2, 3:3, 4:4, 5:5, 6:6}
        
        qindex_to_refcoord, inserts = csf2counts.coordinate_map(
            query_sequence,
            reference_sequence)
        
        self.assertEqual(qindex_to_refcoord, expected_mapping)
        self.assertEqual(inserts, expected_inserts)
        
    def testInsertion(self):
        query_sequence =     'CTNPRPNNN'
        reference_sequence = 'CT--RPNNN'
        expected_inserts = [2, 3]
        expected_mapping = {0:0, 1:1, 4:2, 5:3, 6:4, 7:5, 8:6}
        
        qindex_to_refcoord, inserts = csf2counts.coordinate_map(
            query_sequence,
            reference_sequence)
        
        self.assertEqual(expected_mapping, qindex_to_refcoord)
        self.assertEqual(expected_inserts, inserts)
        
    def testDeletion(self):
        query_sequence =     'CT--NNN'
        reference_sequence = 'CTRPNNN'
        expected_inserts = []
        expected_mapping = {0:0, 1:1, 2:4, 3:5, 4:6}
        
        qindex_to_refcoord, inserts = csf2counts.coordinate_map(
            query_sequence,
            reference_sequence)
        
        self.assertEqual(expected_mapping, qindex_to_refcoord)
        self.assertEqual(expected_inserts, inserts)
        
    def testDeletionAndInsertion(self):
        query_sequence =     'CT--NNCPN'
        reference_sequence = 'CTRPNN--N'
        expected_inserts = [4, 5] # note that these are the non-blank indexes
        expected_mapping = {0:0, 1:1, 2:4, 3:5, 6:6}
        
        qindex_to_refcoord, inserts = csf2counts.coordinate_map(
            query_sequence,
            reference_sequence)
        
        self.assertEqual(qindex_to_refcoord, expected_mapping)
        self.assertEqual(inserts, expected_inserts)
        
    def testQueryStartsBeforeReference(self):
        query_sequence =     'NPCTRPNNN'
        reference_sequence = '--CTRPNNN'
        expected_inserts = []
        expected_mapping = {2:0, 3:1, 4:2, 5:3, 6:4, 7:5, 8:6}
        
        qindex_to_refcoord, inserts = csf2counts.coordinate_map(
            query_sequence,
            reference_sequence)
        
        self.assertEqual(expected_mapping, qindex_to_refcoord)
        self.assertEqual(expected_inserts, inserts)
        
    def testQueryEndsAfterReference(self):
        query_sequence =     'CTRPNNNNP'
        reference_sequence = 'CTRPNNN--'
        expected_inserts = []
        expected_mapping = {0:0, 1:1, 2:2, 3:3, 4:4, 5:5, 6:6}
        
        qindex_to_refcoord, inserts = csf2counts.coordinate_map(
            query_sequence,
            reference_sequence)
        
        self.assertEqual(expected_mapping, qindex_to_refcoord)
        self.assertEqual(expected_inserts, inserts)
        
    def testAmbiguous(self):
        query_sequence =     'CT?PNNN'
        reference_sequence = 'CTRPNNN'
        expected_inserts = []
        expected_mapping = {0:0, 1:1, 3:3, 4:4, 5:5, 6:6}
        
        qindex_to_refcoord, inserts = csf2counts.coordinate_map(
            query_sequence,
            reference_sequence)
        
        self.assertEqual(qindex_to_refcoord, expected_mapping)
        self.assertEqual(inserts, expected_inserts)

class IndelWriterTest(unittest.TestCase):
    def setUp(self):
        self.writer = csf2counts.IndelWriter(indelfile=StringIO.StringIO())
        self.writer.start_group(sample_name='E1234_S1', region='R1', qcut=15)
        
    def testNoInserts(self):
        expected_text = """\
sample,region,qcut,left,insert,count
"""
        
        self.writer.add_read(offset_sequence='ACDEF', count=1)
        self.writer.write(inserts=[])
        
        self.assertEqual(expected_text, self.writer.indelfile.getvalue())
        
    def testInsert(self):
        expected_text = """\
sample,region,qcut,left,insert,count
E1234_S1,R1,15,2,D,1
"""
        
        self.writer.add_read(offset_sequence='ACDEF', count=1)
        self.writer.write(inserts=[2])
        
        self.assertMultiLineEqual(expected_text, self.writer.indelfile.getvalue())
        
    def testInsertWithOffset(self):
        expected_text = """\
sample,region,qcut,left,insert,count
E1234_S1,R1,15,2,D,1
"""
        
        self.writer.add_read(offset_sequence='-CDEFG', count=1)
        self.writer.write(inserts=[1], min_offset=1)
        
        self.assertMultiLineEqual(expected_text, self.writer.indelfile.getvalue())
        
    def testTwoInsertsWithOffset(self):
        expected_text = """\
sample,region,qcut,left,insert,count
E1234_S1,R1,15,2,D,1
E1234_S1,R1,15,4,F,1
"""
        
        self.writer.add_read(offset_sequence='-CDEFG', count=1)
        self.writer.write(inserts=[1, 3], min_offset=1)
        
        self.assertMultiLineEqual(expected_text, self.writer.indelfile.getvalue())

    def testInsertsWithVariants(self):
        expected_text = """\
sample,region,qcut,left,insert,count
E1234_S1,R1,15,2,D,2
"""
        
        self.writer.add_read(offset_sequence='ACDEF', count=1)
        self.writer.add_read(offset_sequence='AFDEF', count=1)
        self.writer.write(inserts=[2])
        
        self.assertMultiLineEqual(expected_text, self.writer.indelfile.getvalue())

    def testDifferentInserts(self):
        expected_text = """\
sample,region,qcut,left,insert,count
E1234_S1,R1,15,2,D,2
E1234_S1,R1,15,2,F,3
"""
        
        self.writer.add_read(offset_sequence='ACDEF', count=2)
        self.writer.add_read(offset_sequence='ACFEF', count=3)
        self.writer.write(inserts=[2])
        
        self.assertMultiLineEqual(expected_text, self.writer.indelfile.getvalue())

    def testMulticharacterInsert(self):
        expected_text = """\
sample,region,qcut,left,insert,count
E1234_S1,R1,15,2,DE,1
"""
        
        self.writer.add_read(offset_sequence='ACDEF', count=1)
        self.writer.write(inserts=[2,3])
        
        self.assertMultiLineEqual(expected_text, self.writer.indelfile.getvalue())
