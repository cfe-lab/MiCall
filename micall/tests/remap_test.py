import StringIO
import unittest

from micall.core import remap
from micall.core.remap import is_first_read, is_short_read

class RemapTest(unittest.TestCase):
    def assertCigarIsPrimer(self, cigar, is_primer_expected):
        row = {'cigar': cigar}
        max_primer_length = 29
        self.assertEqual(is_primer_expected, is_short_read(row, max_primer_length))
        
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
            "GP41-seed\t5\tN\t2\tCc\tAA\n" +
            "GP41-seed\t6\tN\t2\tCc\tAA\n")
        qCutoff = 20
        expected_conseqs = {'GP41-seed': "GCCATAGCC"}
        
        conseqs = remap.pileup_to_conseq(pileupIO, qCutoff)
        
        self.assertSequenceEqual(expected_conseqs, conseqs)
        
    def testInsertionShiftsFrame(self):
        pileupIO = StringIO.StringIO(
            "GP41-seed\t1\tN\t2\t^MG^Mg\tAA\n" +
            "GP41-seed\t2\tN\t2\tCc\tAA\n" +
            "GP41-seed\t3\tN\t2\tC+2ATc+2at\tAA\n" +
            "GP41-seed\t4\tN\t2\tGg\tAA\n" +
            "GP41-seed\t5\tN\t2\tCc\tAA\n" +
            "GP41-seed\t6\tN\t2\tCc\tAA\n")
        qCutoff = 20
        expected_conseqs = {'GP41-seed': "GCCGCC"}
        
        conseqs = remap.pileup_to_conseq(pileupIO, qCutoff)

        self.assertSequenceEqual(expected_conseqs, conseqs)
        
    def testLongInsertion(self):
        pileupIO = StringIO.StringIO(
            "GP41-seed\t1\tN\t2\t^MG^Mg\tAA\n" +
            "GP41-seed\t2\tN\t2\tCc\tAA\n" +
            "GP41-seed\t3\tN\t2\tC+12AAATTTAAATTTc+12aaatttaaattt\tAA\n" +
            "GP41-seed\t4\tN\t2\tGg\tAA\n" +
            "GP41-seed\t5\tN\t2\tCc\tAA\n" +
            "GP41-seed\t6\tN\t2\tCc\tAA\n")
        qCutoff = 20
        expected_conseqs = {'GP41-seed': "GCCAAATTTAAATTTGCC"}
        
        conseqs = remap.pileup_to_conseq(pileupIO, qCutoff)
        
        self.assertSequenceEqual(expected_conseqs, conseqs)
        
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
        expected_conseqs = {'GP41-seed': "GCCGCC"}
        
        conseqs= remap.pileup_to_conseq(pileupIO, qCutoff)
        
        self.assertSequenceEqual(expected_conseqs, conseqs)
        
    def testOffset(self):
        pileupIO = StringIO.StringIO(
            "GP41-seed\t2\tN\t2\t^MC^Mc\tAA\n" +
            "GP41-seed\t3\tN\t2\tCc\tAA\n" +
            "GP41-seed\t4\tN\t2\tGg\tAA\n" +
            "GP41-seed\t5\tN\t2\tCc\tAA\n" +
            "GP41-seed\t6\tN\t2\tCc\tAA\n")
        qCutoff = 20
        expected_conseqs= {'GP41-seed': "NCCGCC"}
        
        conseqs = remap.pileup_to_conseq(pileupIO, qCutoff)
        
        self.assertSequenceEqual(expected_conseqs, conseqs)


class MakeConsensusTest(unittest.TestCase):
    def testInsertion(self):
        pileup = {'GP41-seed': {
            1: {'s': '^MG^Mg', 'q': 'AA'},
            2: {'s': 'Cc', 'q': 'AA'},
            3: {'s': 'C+3ATAc+3ata', 'q': 'AA'},
            4: {'s': 'Gg', 'q': 'AA'},
            5: {'s': 'Cc', 'q': 'AA'},
            6: {'s': 'Cc', 'q': 'AA'}
        }}
        last_conseqs = {'GP41-seed': 'GCCGCC'}
        qCutoff = 20
        expected_conseqs = {'GP41-seed': 'GCCATAGCC'}

        conseqs = remap.make_consensus(pileup=pileup, last_conseqs=last_conseqs, qCutoff=qCutoff)

        self.assertEqual(conseqs, expected_conseqs)

    def testInsertionShiftsFrame(self):
        pileup = {'GP41-seed': {
            1: {'s': '^MG^Mg', 'q': 'AA'},
            2: {'s': 'Cc', 'q': 'AA'},
            3: {'s': 'C+2ATc+2at', 'q': 'AA'},
            4: {'s': 'Gg', 'q': 'AA'},
            5: {'s': 'Cc', 'q': 'AA'},
            6: {'s': 'Cc', 'q': 'AA'}
        }}
        last_conseqs = {'GP41-seed': 'GCCGCC'}
        qCutoff = 20
        expected_conseqs = {'GP41-seed': 'GCCGCC'}

        conseqs = remap.make_consensus(pileup=pileup, last_conseqs=last_conseqs, qCutoff=qCutoff)

        self.assertEqual(conseqs, expected_conseqs)

    def testLongInsertion(self):
        pileup = {'GP41-seed': {
            1: {'s': '^MG^Mg', 'q': 'AA'},
            2: {'s': 'Cc', 'q': 'AA'},
            3: {'s': 'C+12AAATTTAAATTTc+12aaatttaaattt', 'q': 'AA'},
            4: {'s': 'Gg', 'q': 'AA'},
            5: {'s': 'Cc', 'q': 'AA'},
            6: {'s': 'Cc', 'q': 'AA'}
        }}
        last_conseqs = {'GP41-seed': 'GCCGCC'}
        qCutoff = 20
        expected_conseqs = {'GP41-seed': 'GCCAAATTTAAATTTGCC'}

        conseqs = remap.make_consensus(pileup=pileup, last_conseqs=last_conseqs, qCutoff=qCutoff)

        self.assertEqual(conseqs, expected_conseqs)

    def testDeletion(self):
        pileup = {'GP41-seed': {
            1: {'s': '^MG^Mg', 'q': 'AA'},
            2: {'s': 'Cc', 'q': 'AA'},
            3: {'s': 'C-3NNNc-3nnn', 'q': 'AA'},
            4: {'s': '**', 'q': 'AA'},
            5: {'s': '**', 'q': 'AA'},
            6: {'s': '**', 'q': 'AA'},
            7: {'s': 'Gg', 'q': 'AA'},
            8: {'s': 'Cc', 'q': 'AA'},
            9: {'s': 'Cc', 'q': 'AA'}
        }}
        last_conseqs = {'GP41-seed': 'GCCAAAGCC'}
        qCutoff = 20
        expected_conseqs = {'GP41-seed': 'GCCGCC'}

        conseqs = remap.make_consensus(pileup=pileup, last_conseqs=last_conseqs, qCutoff=qCutoff)

        self.assertEqual(conseqs, expected_conseqs)

    def testOffset(self):
        pileup = {'GP41-seed': {
            2: {'s': '^MC^Mc', 'q': 'AA'},
            3: {'s': 'Cc', 'q': 'AA'},
            4: {'s': 'Gg', 'q': 'AA'},
            5: {'s': 'Cc', 'q': 'AA'},
            6: {'s': 'Cc', 'q': 'AA'}
        }}
        last_conseqs = {'GP41-seed': 'GCCGCC'}
        qCutoff = 20
        expected_conseqs = {'GP41-seed': 'NCCGCC'}

        conseqs = remap.make_consensus(pileup=pileup, last_conseqs=last_conseqs, qCutoff=qCutoff)

        self.assertEqual(conseqs, expected_conseqs)


class SamToPileupTest(unittest.TestCase):
    def testOffset(self):
        samIO = StringIO.StringIO(
            "test1\t99\ttest\t3\t44\t12M\t=\t3\t12\tACAAGACCCAAC\tJJJJJJJJJJJJ\n"
            "test1\t147\ttest\t3\t44\t12M\t=\t3\t-12\tACAAGACCCAAC\tJJJJJJJJJJJJ\n"
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
        pileup, _counts = remap.sam_to_pileup(samIO, max_primer_length=0)
        self.maxDiff = None
        self.assertEqual(pileup, expected_pileup)

    def testSimpleInsertion(self):
        samIO = StringIO.StringIO(
            "test1\t99\ttest\t3\t44\t3M3I9M\t=\t3\t12\tACAGGGAGACCCAAC\tJJJJJJJJJJJJJJJ\n"
            "test1\t147\ttest\t3\t44\t3M3I9M\t=\t3\t-12\tACAGGGAGACCCAAC\tJJJJJJJJJJJJJJJ\n"
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
        pileup, _counts = remap.sam_to_pileup(samIO, max_primer_length=0)
        self.maxDiff = None
        self.assertEqual(pileup, expected_pileup)

    def testComplexInsertion(self):
        samIO = StringIO.StringIO(
            "test1\t99\ttest\t3\t44\t3M1I3M2I6M\t=\t3\t12\tACAGAGAGGCCCAAC\tJJJJJJJJJJJJJJJ\n"
            "test1\t147\ttest\t3\t44\t3M1I3M2I6M\t=\t3\t-12\tACAGAGAGGCCCAAC\tJJJJJJJJJJJJJJJ\n"
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
        pileup, _counts = remap.sam_to_pileup(samIO, max_primer_length=0)
        self.maxDiff = None
        self.assertEqual(pileup, expected_pileup)

    def testStaggeredPair(self):
        samIO = StringIO.StringIO(
            "test1\t99\ttest\t3\t44\t12M\t=\t3\t12\tACAAGACCCAAC\tJJJJJJJJJJJJ\n"
            "test1\t147\ttest\t9\t44\t12M\t=\t9\t-12\tCCCAACAACAAT\tJJJJJJJJJJJJ\n"
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
        pileup, _counts = remap.sam_to_pileup(samIO, max_primer_length=0)
        self.maxDiff = None
        self.assertEqual(pileup, expected_pileup)