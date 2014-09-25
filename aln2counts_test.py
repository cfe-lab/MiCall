import StringIO
import unittest

import aln2counts
from testfixtures.logcapture import LogCapture

class AminoFrequencyWriterTest(unittest.TestCase):
    def setUp(self):
        self.writer = aln2counts.AminoFrequencyWriter(aafile=StringIO.StringIO(),
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
        self.writer = aln2counts.NucleotideFrequencyWriter(
            nucfile=StringIO.StringIO(),
            amino_ref_seqs = {'R1': 'EE', # contents irrelevant, length matters
                              'R2': 'EEE'},
            nuc_ref_seqs = {'R3': 'AAA',
                            'R4': 'AAAAAAAAA'}) # contents irrelevant, length matters
        
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
E1234_S1,R3,15,3,,0,0,0,0
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
E1234_S1,R1,15,,3,0,0,0,0
E1234_S1,R1,15,,4,0,0,0,0
E1234_S1,R1,15,,5,0,0,0,0
E1234_S1,R1,15,,6,0,0,0,0
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
E1234_S1,R4,15,1,,0,0,0,0
E1234_S1,R4,15,2,,0,0,0,0
E1234_S1,R4,15,3,,0,0,0,0
E1234_S1,R4,15,4,,0,0,0,0
E1234_S1,R4,15,5,,0,0,0,0
E1234_S1,R4,15,6,,1,2,0,0
E1234_S1,R4,15,7,,0,0,0,3
E1234_S1,R4,15,8,,0,0,0,0
E1234_S1,R4,15,9,,0,0,0,0
"""
        
        self.writer.write(sample_name = 'E1234_S1',
                          region='R4',
                          qcut=15,
                          nuc_counts={5: {'A': 1, 'C': 2},
                                      6: {'T': 3}})
        
        self.assertMultiLineEqual(expected_text, self.writer.nucfile.getvalue())
        
    def testOffsetWithReference(self):
        # No coverage at the start of the reference, and first codon starts
        # with its second nucleotide. The partial codon is output.
        expected_text = """\
sample,region,q-cutoff,query.nuc.pos,refseq.nuc.pos,A,C,G,T
E1234_S1,R2,15,,1,0,0,0,0
E1234_S1,R2,15,,2,0,0,0,0
E1234_S1,R2,15,,3,0,0,0,0
E1234_S1,R2,15,,4,0,0,0,0
E1234_S1,R2,15,,5,0,0,0,0
E1234_S1,R2,15,,6,0,0,0,0
E1234_S1,R2,15,7,7,0,0,1,0
E1234_S1,R2,15,8,8,0,0,1,0
E1234_S1,R2,15,9,9,1,0,0,0
"""
        
        with LogCapture() as log:
            self.writer.write(sample_name = 'E1234_S1',
                              region='R2',
                              qcut=15,
                              nuc_counts={4: {'T': 1},
                                          5: {'G': 1},
                                          6: {'G': 1},
                                          7: {'G': 1},
                                          8: {'A': 1}},
                              qindex_to_refcoord={1: 2},
                              min_offset=4)
        
        self.assertMultiLineEqual(expected_text, self.writer.nucfile.getvalue())
        log.check(('root',
                   'WARNING',
                   "No coordinate mapping for query nuc 4 (amino 0) in E1234_S1"),
                  ('root',
                   'WARNING',
                   "No coordinate mapping for query nuc 5 (amino 0) in E1234_S1"))
        
    def testGapWithoutReference(self):
        expected_text = """\
sample,region,q-cutoff,query.nuc.pos,refseq.nuc.pos,A,C,G,T
E1234_S1,R4,15,1,,1,2,0,0
E1234_S1,R4,15,2,,0,0,0,3
E1234_S1,R4,15,3,,0,0,0,0
E1234_S1,R4,15,4,,0,0,3,0
E1234_S1,R4,15,5,,0,0,0,0
E1234_S1,R4,15,6,,0,0,0,0
E1234_S1,R4,15,7,,0,0,0,0
E1234_S1,R4,15,8,,0,0,0,0
E1234_S1,R4,15,9,,0,0,0,0
"""
        
        self.writer.write(sample_name = 'E1234_S1',
                          region='R4',
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
E1234_S1,R1,15,,3,0,0,0,0
E1234_S1,R1,15,4,4,0,0,3,0
E1234_S1,R1,15,,5,0,0,0,0
E1234_S1,R1,15,,6,0,0,0,0
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
E1234_S1,R2,15,1,1,1,2,0,0
E1234_S1,R2,15,2,2,0,0,0,3
E1234_S1,R2,15,3,3,0,2,0,0
E1234_S1,R2,15,,4,0,0,0,0
E1234_S1,R2,15,,5,0,0,0,0
E1234_S1,R2,15,,6,0,0,0,0
E1234_S1,R2,15,4,7,0,0,3,0
E1234_S1,R2,15,,8,0,0,0,0
E1234_S1,R2,15,,9,0,0,0,0
"""
        
        self.writer.write(sample_name = 'E1234_S1',
                          region='R2',
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
        
        qindex_to_refcoord, inserts = aln2counts.coordinate_map(
            query_sequence,
            reference_sequence)
        
        self.assertEqual(qindex_to_refcoord, expected_mapping)
        self.assertEqual(inserts, expected_inserts)
        
    def testInsertion(self):
        query_sequence =     'CTNPRPNNN'
        reference_sequence = 'CT--RPNNN'
        expected_inserts = [2, 3]
        expected_mapping = {0:0, 1:1, 4:2, 5:3, 6:4, 7:5, 8:6}
        
        qindex_to_refcoord, inserts = aln2counts.coordinate_map(
            query_sequence,
            reference_sequence)
        
        self.assertEqual(expected_mapping, qindex_to_refcoord)
        self.assertEqual(expected_inserts, inserts)
        
    def testDeletion(self):
        query_sequence =     'CT--NNN'
        reference_sequence = 'CTRPNNN'
        expected_inserts = []
        expected_mapping = {0:0, 1:1, 2:4, 3:5, 4:6}
        
        qindex_to_refcoord, inserts = aln2counts.coordinate_map(
            query_sequence,
            reference_sequence)
        
        self.assertEqual(expected_mapping, qindex_to_refcoord)
        self.assertEqual(expected_inserts, inserts)
        
    def testDeletionAndInsertion(self):
        query_sequence =     'CT--NNCPN'
        reference_sequence = 'CTRPNN--N'
        expected_inserts = [4, 5] # note that these are the non-blank indexes
        expected_mapping = {0:0, 1:1, 2:4, 3:5, 6:6}
        
        qindex_to_refcoord, inserts = aln2counts.coordinate_map(
            query_sequence,
            reference_sequence)
        
        self.assertEqual(qindex_to_refcoord, expected_mapping)
        self.assertEqual(inserts, expected_inserts)
        
    def testQueryStartsBeforeReference(self):
        query_sequence =     'NPCTRPNNN'
        reference_sequence = '--CTRPNNN'
        expected_inserts = []
        expected_mapping = {2:0, 3:1, 4:2, 5:3, 6:4, 7:5, 8:6}
        
        qindex_to_refcoord, inserts = aln2counts.coordinate_map(
            query_sequence,
            reference_sequence)
        
        self.assertEqual(expected_mapping, qindex_to_refcoord)
        self.assertEqual(expected_inserts, inserts)
        
    def testQueryEndsAfterReference(self):
        query_sequence =     'CTRPNNNNP'
        reference_sequence = 'CTRPNNN--'
        expected_inserts = []
        expected_mapping = {0:0, 1:1, 2:2, 3:3, 4:4, 5:5, 6:6}
        
        qindex_to_refcoord, inserts = aln2counts.coordinate_map(
            query_sequence,
            reference_sequence)
        
        self.assertEqual(expected_mapping, qindex_to_refcoord)
        self.assertEqual(expected_inserts, inserts)
        
    def testAmbiguous(self):
        query_sequence =     'CT?PNNN'
        reference_sequence = 'CTRPNNN'
        expected_inserts = []
        expected_mapping = {0:0, 1:1, 3:3, 4:4, 5:5, 6:6}
        
        qindex_to_refcoord, inserts = aln2counts.coordinate_map(
            query_sequence,
            reference_sequence)
        
        self.assertEqual(qindex_to_refcoord, expected_mapping)
        self.assertEqual(inserts, expected_inserts)

class IndelWriterTest(unittest.TestCase):
    def setUp(self):
        self.writer = aln2counts.IndelWriter(indelfile=StringIO.StringIO())
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

class TranslateTest(unittest.TestCase):
    def testSingleCodon(self):
        nucs = 'TTT'
        offset = 0
        expected_aminos = 'F'

        aminos = aln2counts.translate(nucs, offset)
        
        self.assertEqual(expected_aminos, aminos)
        
    def testPartialCodon(self):
        nucs = 'TTTC'
        offset = 0
        expected_aminos = 'F'

        aminos = aln2counts.translate(nucs, offset)
        
        self.assertEqual(expected_aminos, aminos)
        
    def testTwoCodons(self):
        nucs = 'TTTCCT'
        offset = 0
        expected_aminos = 'FP'

        aminos = aln2counts.translate(nucs, offset)
        
        self.assertEqual(expected_aminos, aminos)
        
    def testOffset(self):
        nucs = "TTTCCT"
        offset = 3
        expected_aminos = "-FP"
        
        aminos = aln2counts.translate(nucs, offset)
        
        self.assertEqual(expected_aminos, aminos)
        
    def testSingleDashAmbiguous(self):
        nucs = '-TT'
        offset = 0
        expected_aminos = '?'

        aminos = aln2counts.translate(nucs, offset)
        
        self.assertEqual(expected_aminos, aminos)
        
    def testSingleDashUnambiguous(self):
        nucs = 'CG-' # CGA, CGC, CGG, CGT all map to R
        offset = 0
        expected_aminos = 'R'

        aminos = aln2counts.translate(nucs, offset)
        
        self.assertEqual(expected_aminos, aminos)
        
    def testTwoDashes(self):
        nucs = '--T'
        offset = 0
        expected_aminos = '?'

        aminos = aln2counts.translate(nucs, offset)
        
        self.assertEqual(expected_aminos, aminos)
        
    def testThreeDashes(self):
        nucs = '---'
        offset = 0
        expected_aminos = '-'

        aminos = aln2counts.translate(nucs, offset)
        
        self.assertEqual(expected_aminos, aminos)
        
    def testAmbiguousBasesThatAreSynonyms(self):
        nucs = 'TTY' # TTC or TTT: both map to F
        offset = 0
        expected_aminos = 'F'

        aminos = aln2counts.translate(nucs, offset)
        
        self.assertEqual(expected_aminos, aminos)
        
    def testTwoAmbiguousBasesThatAreSynonyms(self):
        nucs = 'MGR' # CGA, CGG, AGA, or AGG: all map to R
        offset = 0
        expected_aminos = '?'

        aminos = aln2counts.translate(nucs, offset)
        
        self.assertEqual(expected_aminos, aminos)

class MakeConsensusTest(unittest.TestCase):
    def testNoMixes(self):
        nuc_counts = {0: {'A': 1},
                      1: {'C': 1},
                      2: {'T': 1},
                      3: {'G': 1}}
        conseq_mixture_cutoffs = [0.1]
        expected_conseqs = {'MAX': 'ACTG',
                            '0.100': 'ACTG'}
        
        conseqs = aln2counts.make_consensus(nuc_counts, conseq_mixture_cutoffs)
        
        self.assertDictEqual(expected_conseqs, conseqs)

    def testMixed(self):
        nuc_counts = {0: {'A': 3},
                      1: {'C': 2, 'T': 1},
                      2: {'T': 3},
                      3: {'G': 3}}
        conseq_mixture_cutoffs = [0.1]
        expected_conseqs = {'MAX': 'ACTG',
                            '0.100': 'AYTG'} # Y is a mix of C and T
        
        conseqs = aln2counts.make_consensus(nuc_counts, conseq_mixture_cutoffs)
        
        self.assertDictEqual(expected_conseqs, conseqs)

    def testMixedThree(self):
        nuc_counts = {0: {'A': 4},
                      1: {'C': 2, 'T': 1, 'G': 1},
                      2: {'T': 4},
                      3: {'G': 4}}
        conseq_mixture_cutoffs = [0.1]
        expected_conseqs = {'MAX': 'ACTG',
                            '0.100': 'ABTG'} # B is a mix of T, G, and C
        
        conseqs = aln2counts.make_consensus(nuc_counts, conseq_mixture_cutoffs)
        
        self.assertDictEqual(expected_conseqs, conseqs)

    def testMixedAll(self):
        nuc_counts = {0: {'A': 5},
                      1: {'C': 2, 'T': 1, 'G': 1, 'A': 1},
                      2: {'T': 5},
                      3: {'G': 5}}
        conseq_mixture_cutoffs = [0.1]
        expected_conseqs = {'MAX': 'ACTG',
                            '0.100': 'ANTG'} # All four are reported as N
        
        conseqs = aln2counts.make_consensus(nuc_counts, conseq_mixture_cutoffs)
        
        self.assertDictEqual(expected_conseqs, conseqs)

    def testMixedMax(self):
        nuc_counts = {0: {'A': 5},
                      1: {'C': 2, 'T': 2, 'G': 1},
                      2: {'T': 5},
                      3: {'G': 5}}
        conseq_mixture_cutoffs = [0.1]
        expected_conseqs = {'MAX': 'AYTG', # C and T tie for max, mix is Y
                            '0.100': 'ABTG'} # C, T, and G mix is B
        
        conseqs = aln2counts.make_consensus(nuc_counts, conseq_mixture_cutoffs)
        
        self.assertDictEqual(expected_conseqs, conseqs)

    def testCutoff(self):
        nuc_counts = {0: {'A': 3},
                      1: {'C': 2, 'T': 1},
                      2: {'T': 3},
                      3: {'G': 3}}
        conseq_mixture_cutoffs = [0.1, 0.5]
        expected_conseqs = {'MAX': 'ACTG',
                            '0.100': 'AYTG', # Y is a mix of C and T
                            '0.500': 'ACTG'} # T was below the cutoff
        
        conseqs = aln2counts.make_consensus(nuc_counts, conseq_mixture_cutoffs)
        
        self.assertDictEqual(expected_conseqs, conseqs)

    def testCutoffBoundary(self):
        # T is at cutoff in position 1, then below cutoff in position 2
        nuc_counts = {0: {'A': 9950},
                      1: {'C': 9000, 'T': 1000},
                      2: {'C': 9001, 'T':  999},
                      3: {'G': 9799}}
        conseq_mixture_cutoffs = [0.1]
        expected_conseqs = {'MAX': 'ACCG',
                            '0.100': 'AYCG'}
        
        conseqs = aln2counts.make_consensus(nuc_counts, conseq_mixture_cutoffs)
        
        self.assertDictEqual(expected_conseqs, conseqs)

    def testMixedWithPoorQuality(self):
        nuc_counts = {0: {'A': 3},
                      1: {'N': 2, 'T': 1},
                      2: {'T': 3},
                      3: {'G': 3}}
        conseq_mixture_cutoffs = [0.1]
        expected_conseqs = {'MAX': 'ATTG',
                            '0.100': 'ATTG'} # N always overruled
        
        conseqs = aln2counts.make_consensus(nuc_counts, conseq_mixture_cutoffs)
        
        self.assertDictEqual(expected_conseqs, conseqs)

    def testMixedWithGap(self):
        nuc_counts = {0: {'A': 3},
                      1: {'-': 2, 'T': 1},
                      2: {'T': 3},
                      3: {'G': 3}}
        conseq_mixture_cutoffs = [0.1]
        expected_conseqs = {'MAX': 'ATTG',
                            '0.100': 'ATTG'} # dash always overruled
        
        conseqs = aln2counts.make_consensus(nuc_counts, conseq_mixture_cutoffs)
        
        self.assertDictEqual(expected_conseqs, conseqs)

    def testMixedWithGapAndPoorQuality(self):
        nuc_counts = {0: {'A': 6},
                      1: {'N': 3, '-': 2, 'T': 1},
                      2: {'T': 6},
                      3: {'G': 6}}
        conseq_mixture_cutoffs = [0.1]
        expected_conseqs = {'MAX': 'ATTG',
                            '0.100': 'ATTG'}
        
        conseqs = aln2counts.make_consensus(nuc_counts, conseq_mixture_cutoffs)
        
        self.assertDictEqual(expected_conseqs, conseqs)

    def testPoorQualityOnly(self):
        nuc_counts = {0: {'A': 1},
                      1: {'N': 1},
                      2: {'T': 1},
                      3: {'G': 1}}
        conseq_mixture_cutoffs = [0.1]
        expected_conseqs = {'MAX': 'ANTG',
                            '0.100': 'ANTG'}
        
        conseqs = aln2counts.make_consensus(nuc_counts, conseq_mixture_cutoffs)
        
        self.assertDictEqual(expected_conseqs, conseqs)

    def testMixedGapAndPoorQualityOnly(self):
        nuc_counts = {0: {'A': 5},
                      1: {'N': 3, '-': 2},
                      2: {'T': 5},
                      3: {'G': 5}}
        conseq_mixture_cutoffs = [0.1]
        expected_conseqs = {'MAX': 'ANTG',
                            '0.100': 'ANTG'}
        
        conseqs = aln2counts.make_consensus(nuc_counts, conseq_mixture_cutoffs)
        
        self.assertDictEqual(expected_conseqs, conseqs)

    def testAllBelowCutoff(self):
        nuc_counts = {0: {'A': 300},
                      1: {'C': 101, 'T': 100, 'G': 99},
                      2: {'T': 300},
                      3: {'G': 300}}
        conseq_mixture_cutoffs = [0.5]
        expected_conseqs = {'MAX': 'ACTG',
                            '0.500': 'ANTG'}
        
        conseqs = aln2counts.make_consensus(nuc_counts, conseq_mixture_cutoffs)
        
        self.assertDictEqual(expected_conseqs, conseqs)

    def testMissingPositions(self):
        """ Missing positions are ignored in the consensus sequence. """
        nuc_counts = {0: {'A': 1},
                      1: {'C': 1},
                      2: {'T': 1},
                      4: {'G': 1}}
        conseq_mixture_cutoffs = [0.1]
        expected_conseqs = {'MAX': 'ACTG',
                            '0.100': 'ACTG'}
        
        conseqs = aln2counts.make_consensus(nuc_counts, conseq_mixture_cutoffs)
        
        self.assertDictEqual(expected_conseqs, conseqs)

    def testOrdering(self):
        """ Position keys are sorted. """
        nuc_counts = {2: {'A': 1},
                      3: {'C': 1},
                      4: {'T': 1},
                      100: {'G': 1}}
        conseq_mixture_cutoffs = [0.1]
        expected_conseqs = {'MAX': 'ACTG',
                            '0.100': 'ACTG'}
        
        conseqs = aln2counts.make_consensus(nuc_counts, conseq_mixture_cutoffs)
        
        self.assertDictEqual(expected_conseqs, conseqs)
