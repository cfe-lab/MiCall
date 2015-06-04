import StringIO
import unittest

from micall.core import remap
from micall.core.remap import is_first_read, is_primer

class RemapTest(unittest.TestCase):
    def testSampleName(self):
        filename = '/some/path/2020A-V3LOOP_S4_L001_R1_001.fastq'
        expected_sample_name = '2020A-V3LOOP_S4'
        
        sample_name = remap.calculate_sample_name(filename)
        
        self.assertEqual(sample_name, expected_sample_name)
        
    def assertCigarIsPrimer(self, cigar, is_primer_expected):
        row = {'cigar': cigar}
        max_primer_length = 29
        self.assertEqual(is_primer_expected, is_primer(row, max_primer_length))
        
    def testIsPrimerForLongRead(self):
        self.assertCigarIsPrimer('45M', False)
        
    def testIsPrimerForShortRead(self):
        self.assertCigarIsPrimer('10M', True)
        
    def testIsPrimerForShortReadWithClipping(self):
        self.assertCigarIsPrimer('45S10M', True)
        
    def testIsPrimerForReadWithMultipleMatches(self):
        self.assertCigarIsPrimer('10M3D45M', False)
        

class IsFirstReadTest(unittest.TestCase):
    def testFirstRead(self):
        flag = '99'
        isFirstExpected = True
        
        isFirst = is_first_read(flag)
        
        self.assertEqual(isFirstExpected, isFirst)

    def testSecondRead(self):
        flag = '147'
        isFirstExpected = False
        
        isFirst = is_first_read(flag)
        
        self.assertEqual(isFirstExpected, isFirst)        
        
    def testSmallFlag(self):
        flag = '3'
        isFirstExpected = False
        
        isFirst = is_first_read(flag)
        
        self.assertEqual(isFirstExpected, isFirst)

class PileupToConseqTest(unittest.TestCase):
    def testInsertion(self):
        pileupIO = StringIO.StringIO(
            "GP41-seed\t1\tN\t2\t^MG^Mg\tAA\n" +
            "GP41-seed\t2\tN\t2\tCc\tAA\n" +
            "GP41-seed\t3\tN\t2\tC+3ATAc+3ata\tAA\n" +
            "GP41-seed\t4\tN\t2\tGg\tAA\n" +
            "GP41-seed\t2\tN\t2\tCc\tAA\n" +
            "GP41-seed\t2\tN\t2\tCc\tAA\n")
        qCutoff = 20
        expected_conseq = "GCCATAGCC"
        
        conseq = remap.pileup_to_conseq(pileupIO, qCutoff)
        
        self.assertSequenceEqual(expected_conseq, conseq)
        
    def testInsertionShiftsFrame(self):
        pileupIO = StringIO.StringIO(
            "GP41-seed\t1\tN\t2\t^MG^Mg\tAA\n" +
            "GP41-seed\t2\tN\t2\tCc\tAA\n" +
            "GP41-seed\t3\tN\t2\tC+2ATc+2at\tAA\n" +
            "GP41-seed\t4\tN\t2\tGg\tAA\n" +
            "GP41-seed\t5\tN\t2\tCc\tAA\n" +
            "GP41-seed\t6\tN\t2\tCc\tAA\n")
        qCutoff = 20
        expected_conseq = "GCCGCC"
        
        conseq = remap.pileup_to_conseq(pileupIO, qCutoff)
        
        self.assertSequenceEqual(expected_conseq, conseq)
        
    def testLongInsertion(self):
        pileupIO = StringIO.StringIO(
            "GP41-seed\t1\tN\t2\t^MG^Mg\tAA\n" +
            "GP41-seed\t2\tN\t2\tCc\tAA\n" +
            "GP41-seed\t3\tN\t2\tC+12AAATTTAAATTTc+12aaatttaaattt\tAA\n" +
            "GP41-seed\t4\tN\t2\tGg\tAA\n" +
            "GP41-seed\t5\tN\t2\tCc\tAA\n" +
            "GP41-seed\t6\tN\t2\tCc\tAA\n")
        qCutoff = 20
        expected_conseq = "GCCAAATTTAAATTTGCC"
        
        conseq = remap.pileup_to_conseq(pileupIO, qCutoff)
        
        self.assertSequenceEqual(expected_conseq, conseq)
        
    def testDeletion(self):
        pileupIO = StringIO.StringIO(
            "GP41-seed\t1\tN\t2\t^MG^Mg\tAA\n" +
            "GP41-seed\t2\tN\t2\tCc\tAA\n" +
            "GP41-seed\t3\tN\t2\tC-3NNNc-3nnn\tAA\n" +
            "GP41-seed\t4\tN\t2\t**\tAA\n" +
            "GP41-seed\t5\tN\t2\t**\tAA\n" +
            "GP41-seed\t6\tN\t2\t**\tAA\n" +
            "GP41-seed\t7\tN\t2\tGg\tAA\n" +
            "GP41-seed\t8\tN\t2\tCc\tAA\n" +
            "GP41-seed\t9\tN\t2\tCc\tAA\n")
        qCutoff = 20
        expected_conseq = "GCCGCC"
        
        conseq = remap.pileup_to_conseq(pileupIO, qCutoff)
        
        self.assertSequenceEqual(expected_conseq, conseq)
        
    def testOffset(self):
        pileupIO = StringIO.StringIO(
            "GP41-seed\t2\tN\t2\t^MC^Mc\tAA\n" +
            "GP41-seed\t3\tN\t2\tCc\tAA\n" +
            "GP41-seed\t4\tN\t2\tGg\tAA\n" +
            "GP41-seed\t5\tN\t2\tCc\tAA\n" +
            "GP41-seed\t6\tN\t2\tCc\tAA\n")
        qCutoff = 20
        expected_conseq = "NCCGCC"
        
        conseq = remap.pileup_to_conseq(pileupIO, qCutoff)
        
        self.assertSequenceEqual(expected_conseq, conseq)

class SamToPileupTest(unittest.TestCase):
    def testOffset(self):
        samIO = StringIO.StringIO(
            "qname,flag,rname,pos,mapq,cigar,rnext,pnext,tlen,seq,qual\n"
            "test1,99,test,3,44,12M,=,3,12,ACAAGACCCAAC,JJJJJJJJJJJJ\n"
            "test1,147,test,3,44,12M,=,3,-12,ACAAGACCCAAC,JJJJJJJJJJJJ\n"
        )
        expected_pileup = {'test': {3: {'s': '^MA^Ma', 'q': 'JJ'},
                                    4: {'s': 'Cc', 'q': 'JJ'},
                                    5: {'s': 'Aa', 'q': 'JJ'},
                                    6: {'s': 'Aa', 'q': 'JJ'},
                                    7: {'s': 'Gg', 'q': 'JJ'},
                                    8: {'s': 'Aa', 'q': 'JJ'},
                                    9: {'s': 'Cc', 'q': 'JJ'},
                                    10: {'s': 'Cc', 'q': 'JJ'},
                                    11: {'s': 'Cc', 'q': 'JJ'},
                                    12: {'s': 'Aa', 'q': 'JJ'},
                                    13: {'s': 'Aa', 'q': 'JJ'},
                                    14: {'s': 'C$c$', 'q': 'JJ'}
                                    }}
        expected_counts = {'test': 2}
        pileup, counts = remap.csv_to_pileup(samIO)
        self.maxDiff = None
        self.assertEqual(pileup, expected_pileup)
        self.assertEqual(counts, expected_counts)

    def testSimpleInsertion(self):
        samIO = StringIO.StringIO(
            "qname,flag,rname,pos,mapq,cigar,rnext,pnext,tlen,seq,qual\n"
            "test1,99,test,3,44,3M3I9M,=,3,12,ACAGGGAGACCCAAC,JJJJJJJJJJJJJJJ\n"
            "test1,147,test,3,44,3M3I9M,=,3,-12,ACAGGGAGACCCAAC,JJJJJJJJJJJJJJJ\n"
        )
        expected_pileup = {'test': {3: {'s': '^MA^Ma', 'q': 'JJ'},
                                    4: {'s': 'Cc', 'q': 'JJ'},
                                    5: {'s': 'A+3GGGa+3ggg', 'q': 'JJ'},
                                    6: {'s': 'Aa', 'q': 'JJ'},
                                    7: {'s': 'Gg', 'q': 'JJ'},
                                    8: {'s': 'Aa', 'q': 'JJ'},
                                    9: {'s': 'Cc', 'q': 'JJ'},
                                    10: {'s': 'Cc', 'q': 'JJ'},
                                    11: {'s': 'Cc', 'q': 'JJ'},
                                    12: {'s': 'Aa', 'q': 'JJ'},
                                    13: {'s': 'Aa', 'q': 'JJ'},
                                    14: {'s': 'C$c$', 'q': 'JJ'}
                                    }}
        expected_counts = {'test': 2}
        pileup, counts = remap.csv_to_pileup(samIO)
        self.maxDiff = None
        self.assertEqual(pileup, expected_pileup)
        self.assertEqual(counts, expected_counts)

    def testComplexInsertion(self):
        samIO = StringIO.StringIO(
            "qname,flag,rname,pos,mapq,cigar,rnext,pnext,tlen,seq,qual\n"
            "test1,99,test,3,44,3M1I3M2I6M,=,3,12,ACAGAGAGGCCCAAC,JJJJJJJJJJJJJJJ\n"
            "test1,147,test,3,44,3M1I3M2I6M,=,3,-12,ACAGAGAGGCCCAAC,JJJJJJJJJJJJJJJ\n"
        )
        expected_pileup = {'test': {3: {'s': '^MA^Ma', 'q': 'JJ'},
                                    4: {'s': 'Cc', 'q': 'JJ'},
                                    5: {'s': 'A+1Ga+1g', 'q': 'JJ'},
                                    6: {'s': 'Aa', 'q': 'JJ'},
                                    7: {'s': 'Gg', 'q': 'JJ'},
                                    8: {'s': 'A+2GGa+2gg', 'q': 'JJ'},
                                    9: {'s': 'Cc', 'q': 'JJ'},
                                    10: {'s': 'Cc', 'q': 'JJ'},
                                    11: {'s': 'Cc', 'q': 'JJ'},
                                    12: {'s': 'Aa', 'q': 'JJ'},
                                    13: {'s': 'Aa', 'q': 'JJ'},
                                    14: {'s': 'C$c$', 'q': 'JJ'}
                                    }}
        expected_counts = {'test': 2}
        pileup, counts = remap.csv_to_pileup(samIO)
        self.maxDiff = None
        self.assertEqual(pileup, expected_pileup)
        self.assertEqual(counts, expected_counts)

    def testStaggeredPair(self):
        samIO = StringIO.StringIO(
            "qname,flag,rname,pos,mapq,cigar,rnext,pnext,tlen,seq,qual\n"
            "test1,99,test,3,44,12M,=,3,12,ACAAGACCCAAC,JJJJJJJJJJJJ\n"
            "test1,147,test,9,44,12M,=,9,-12,CCCAACAACAAT,JJJJJJJJJJJJ\n"
        )
        expected_pileup = {'test': {3: {'s': '^MA', 'q': 'J'},
                                    4: {'s': 'C', 'q': 'J'},
                                    5: {'s': 'A', 'q': 'J'},
                                    6: {'s': 'A', 'q': 'J'},
                                    7: {'s': 'G', 'q': 'J'},
                                    8: {'s': 'A', 'q': 'J'},
                                    9: {'s': 'C^Mc', 'q': 'JJ'},
                                    10: {'s': 'Cc', 'q': 'JJ'},
                                    11: {'s': 'Cc', 'q': 'JJ'},
                                    12: {'s': 'Aa', 'q': 'JJ'},
                                    13: {'s': 'Aa', 'q': 'JJ'},
                                    14: {'s': 'C$c', 'q': 'JJ'},
                                    15: {'s': 'a', 'q': 'J'},
                                    16: {'s': 'a', 'q': 'J'},
                                    17: {'s': 'c', 'q': 'J'},
                                    18: {'s': 'a', 'q': 'J'},
                                    19: {'s': 'a', 'q': 'J'},
                                    20: {'s': 't$', 'q': 'J'}
                                    }}
        expected_counts = {'test': 2}
        pileup, counts = remap.csv_to_pileup(samIO)
        self.maxDiff = None
        self.assertEqual(pileup, expected_pileup)
        self.assertEqual(counts, expected_counts)
