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
        
        self.assertDictEqual(expected_conseqs, conseqs)
        
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

        self.assertDictEqual(expected_conseqs, conseqs)
        
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
        
        self.assertDictEqual(expected_conseqs, conseqs)
        
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
        
        self.assertDictEqual(expected_conseqs, conseqs)
        
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
        
        self.assertDictEqual(expected_conseqs, conseqs)
        
    def testOffsetInSecondRegion(self):
        pileupIO = StringIO.StringIO(
            "V3LOOP-seed\t1\tN\t1\t^MC^Mc\tAA\n" +
            "GP41-seed\t2\tN\t2\t^MC^Mc\tAA\n" +
            "GP41-seed\t3\tN\t2\tCc\tAA\n" +
            "GP41-seed\t4\tN\t2\tGg\tAA\n" +
            "GP41-seed\t5\tN\t2\tCc\tAA\n" +
            "GP41-seed\t6\tN\t2\tCc\tAA\n")
        qCutoff = 20
        expected_conseqs= {'V3LOOP-seed': 'C', 'GP41-seed': "NCCGCC"}
        
        conseqs = remap.pileup_to_conseq(pileupIO, qCutoff)
        
        self.assertDictEqual(expected_conseqs, conseqs)
        
    def testGapBetweenReads(self):
        pileupIO = StringIO.StringIO(
            "GP41-seed\t1\tN\t1\t^MC^MC\tAA\n" +
            "GP41-seed\t2\tN\t2\tG$G$\tAA\n" +
            "GP41-seed\t5\tN\t2\t^MT^MT\tB3\n" +
            "GP41-seed\t6\tN\t2\tA$A$\tAA\n")
        qCutoff = 20
        expected_conseqs= {'GP41-seed': "CGNNTA"}
        
        conseqs = remap.pileup_to_conseq(pileupIO, qCutoff)
        
        self.assertDictEqual(expected_conseqs, conseqs)
        
    def testMajorityLowQuality(self):
        pileupIO = StringIO.StringIO(
            "GP41-seed\t1\tN\t1\t^MC^MC\tAA\n" +
            "GP41-seed\t2\tN\t2\tTTT\tA33\n" +
            "GP41-seed\t3\tN\t2\tG$G$\tAA\n")
        qCutoff = 20
        expected_conseqs= {'GP41-seed': "CTG"}
        
        conseqs = remap.pileup_to_conseq(pileupIO, qCutoff)
        
        self.assertDictEqual(expected_conseqs, conseqs)


class MakeConsensusTest(unittest.TestCase):
    def testSimple(self):
        pileup = {'GP41-seed': {
            1: {'s': '^MC^Mc', 'q': 'AA'},
            2: {'s': 'Cc', 'q': 'AA'},
            3: {'s': 'Gg', 'q': 'AA'},
            4: {'s': 'Cc', 'q': 'AA'},
            5: {'s': 'Cc', 'q': 'AA'}
        }}
        qCutoff = 20
        expected_conseqs = {'GP41-seed': 'CCGCC'}

        conseqs = remap.make_consensus(pileup=pileup, qCutoff=qCutoff)

        self.assertEqual(conseqs, expected_conseqs)

    def testInsertion(self):
        pileup = {'GP41-seed': {
            1: {'s': '^MG^Mg', 'q': 'AA'},
            2: {'s': 'Cc', 'q': 'AA'},
            3: {'s': 'C+3ATAc+3ata', 'q': 'AA'},
            4: {'s': 'Gg', 'q': 'AA'},
            5: {'s': 'Cc', 'q': 'AA'},
            6: {'s': 'Cc', 'q': 'AA'}
        }}
        qCutoff = 20
        expected_conseqs = {'GP41-seed': 'GCCATAGCC'}

        conseqs = remap.make_consensus(pileup=pileup, qCutoff=qCutoff)

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
        qCutoff = 20
        expected_conseqs = {'GP41-seed': 'GCCGCC'}

        conseqs = remap.make_consensus(pileup=pileup, qCutoff=qCutoff)

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
        qCutoff = 20
        expected_conseqs = {'GP41-seed': 'GCCAAATTTAAATTTGCC'}

        conseqs = remap.make_consensus(pileup=pileup, qCutoff=qCutoff)

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
        qCutoff = 20
        expected_conseqs = {'GP41-seed': 'GCCGCC'}

        conseqs = remap.make_consensus(pileup=pileup, qCutoff=qCutoff)

        self.assertEqual(conseqs, expected_conseqs)

    def testOffset(self):
        pileup = {'GP41-seed': {
            2: {'s': '^MC^Mc', 'q': 'AA'},
            3: {'s': 'Cc', 'q': 'AA'},
            4: {'s': 'Gg', 'q': 'AA'},
            5: {'s': 'Cc', 'q': 'AA'},
            6: {'s': 'Cc', 'q': 'AA'}
        }}
        qCutoff = 20
        expected_conseqs = {'GP41-seed': 'NCCGCC'}

        conseqs = remap.make_consensus(pileup=pileup, qCutoff=qCutoff)

        self.assertEqual(conseqs, expected_conseqs)

    def testLowQuality(self):
        pileup = {'GP41-seed': {
            1: {'s': 'CCC', 'q': 'A//'}
        }}
        qCutoff = 20
        expected_conseqs = {'GP41-seed': 'C'}

        conseqs = remap.make_consensus(pileup=pileup, qCutoff=qCutoff)

        self.assertEqual(conseqs, expected_conseqs)


class SamToPileupTest(unittest.TestCase):
    def testOffset(self):
        #SAM:qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual
        samIO = StringIO.StringIO(
            "test1\t99\ttest\t3\t44\t12M\t=\t3\t12\tACAAGACCCAAC\tJJJJJJJJJJJJ\n"
            "test1\t147\ttest\t3\t44\t12M\t=\t3\t-12\tACAAGACCCAAC\tJJJJJJJJJJJJ\n"
        )
        expected_pileup = {'test': {3: {'s': '^MA', 'q': 'J'},
                                    4: {'s': 'C', 'q': 'J'},
                                    5: {'s': 'A', 'q': 'J'},
                                    6: {'s': 'A', 'q': 'J'},
                                    7: {'s': 'G', 'q': 'J'},
                                    8: {'s': 'A', 'q': 'J'},
                                    9: {'s': 'C', 'q': 'J'},
                                    10: {'s': 'C', 'q': 'J'},
                                    11: {'s': 'C', 'q': 'J'},
                                    12: {'s': 'A', 'q': 'J'},
                                    13: {'s': 'A', 'q': 'J'},
                                    14: {'s': 'C$', 'q': 'J'}
                                    }}
        pileup, _counts = remap.sam_to_pileup(samIO, max_primer_length=0)
        self.maxDiff = None
        self.assertEqual(expected_pileup, pileup)

    def testSimpleInsertion(self):
        samIO = StringIO.StringIO(
            "test1\t99\ttest\t3\t44\t3M3I9M\t=\t3\t12\tACAGGGAGACCCAAC\tJJJJJJJJJJJJJJJ\n"
            "test1\t147\ttest\t3\t44\t3M3I9M\t=\t3\t-12\tACAGGGAGACCCAAC\tJJJJJJJJJJJJJJJ\n"
        )
        expected_pileup = {'test': {3: {'s': '^MA', 'q': 'J'},
                                    4: {'s': 'C', 'q': 'J'},
                                    5: {'s': 'A+3GGG', 'q': 'J'},
                                    6: {'s': 'A', 'q': 'J'},
                                    7: {'s': 'G', 'q': 'J'},
                                    8: {'s': 'A', 'q': 'J'},
                                    9: {'s': 'C', 'q': 'J'},
                                    10: {'s': 'C', 'q': 'J'},
                                    11: {'s': 'C', 'q': 'J'},
                                    12: {'s': 'A', 'q': 'J'},
                                    13: {'s': 'A', 'q': 'J'},
                                    14: {'s': 'C$', 'q': 'J'}
                                    }}
        pileup, _counts = remap.sam_to_pileup(samIO, max_primer_length=0)
        self.maxDiff = None
        self.assertEqual(expected_pileup, pileup)

    def testComplexInsertion(self):
        samIO = StringIO.StringIO(
            "test1\t99\ttest\t3\t44\t3M1I3M2I6M\t=\t3\t12\tACAGAGAGGCCCAAC\tJJJJJJJJJJJJJJJ\n"
            "test1\t147\ttest\t3\t44\t3M1I3M2I6M\t=\t3\t-12\tACAGAGAGGCCCAAC\tJJJJJJJJJJJJJJJ\n"
        )
        expected_pileup = {'test': {3: {'s': '^MA', 'q': 'J'},
                                    4: {'s': 'C', 'q': 'J'},
                                    5: {'s': 'A+1G', 'q': 'J'},
                                    6: {'s': 'A', 'q': 'J'},
                                    7: {'s': 'G', 'q': 'J'},
                                    8: {'s': 'A+2GG', 'q': 'J'},
                                    9: {'s': 'C', 'q': 'J'},
                                    10: {'s': 'C', 'q': 'J'},
                                    11: {'s': 'C', 'q': 'J'},
                                    12: {'s': 'A', 'q': 'J'},
                                    13: {'s': 'A', 'q': 'J'},
                                    14: {'s': 'C$', 'q': 'J'}
                                    }}
        pileup, _counts = remap.sam_to_pileup(samIO, max_primer_length=0)
        self.maxDiff = None
        self.assertEqual(expected_pileup, pileup)

    def testDeletion(self):
        samIO = StringIO.StringIO(
            "test1\t99\ttest\t3\t44\t3M3D3M\t=\t3\t6\tACAGGG\tJJJJJJ\n"
            "test1\t147\ttest\t3\t44\t3M3D3M\t=\t3\t-6\tACAGGG\tJJJJJJ\n"
        )
        expected_pileup = {'test': {3: {'s': '^MA', 'q': 'J'},
                                    4: {'s': 'C', 'q': 'J'},
                                    5: {'s': 'A-3NNN', 'q': 'J'},
                                    6: {'s': '*', 'q': ''},
                                    7: {'s': '*', 'q': ''},
                                    8: {'s': '*', 'q': ''},
                                    9: {'s': 'G', 'q': 'J'},
                                    10: {'s': 'G', 'q': 'J'},
                                    11: {'s': 'G$', 'q': 'J'}
                                    }}
        pileup, _counts = remap.sam_to_pileup(samIO, max_primer_length=0)
        self.maxDiff = None
        self.assertEqual(expected_pileup, pileup)

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
                                    9: {'s': 'C', 'q': 'J'},
                                    10: {'s': 'C', 'q': 'J'},
                                    11: {'s': 'C', 'q': 'J'},
                                    12: {'s': 'A', 'q': 'J'},
                                    13: {'s': 'A', 'q': 'J'},
                                    14: {'s': 'C$', 'q': 'J'},
                                    15: {'s': 'a', 'q': 'J'},
                                    16: {'s': 'a', 'q': 'J'},
                                    17: {'s': 'c', 'q': 'J'},
                                    18: {'s': 'a', 'q': 'J'},
                                    19: {'s': 'a', 'q': 'J'},
                                    20: {'s': 't$', 'q': 'J'}
                                    }}
        pileup, _counts = remap.sam_to_pileup(samIO, max_primer_length=0)
        self.maxDiff = None
        self.assertEqual(expected_pileup, pileup)

    def testDeletionNearOverlap(self):
        samIO = StringIO.StringIO(
            "test1\t99\ttest\t3\t44\t6M\t=\t3\t6\tACAGGG\tJJJJJJ\n"
            "test1\t147\ttest\t6\t44\t3M3D3M\t=\t3\t-6\tGGGCAT\tJJJJJJ\n"
        )
        expected_pileup = {'test': {3: {'s': '^MA', 'q': 'J'},
                                    4: {'s': 'C', 'q': 'J'},
                                    5: {'s': 'A', 'q': 'J'},
                                    6: {'s': 'G', 'q': 'J'},
                                    7: {'s': 'G', 'q': 'J'},
                                    8: {'s': 'G-3nnn$', 'q': 'J'},
                                    9: {'s': '*', 'q': ''},
                                    10: {'s': '*', 'q': ''},
                                    11: {'s': '*', 'q': ''},
                                    12: {'s': 'c', 'q': 'J'},
                                    13: {'s': 'a', 'q': 'J'},
                                    14: {'s': 't$', 'q': 'J'}
                                    }}
        pileup, _counts = remap.sam_to_pileup(samIO, max_primer_length=0)
        self.maxDiff = None
        self.assertEqual(expected_pileup, pileup)

    def testPairMapsToTwoReferences(self):
        samIO = StringIO.StringIO(
            "test1\t99\ttestX\t3\t44\t12M\t=\t3\t12\tACAAGACCCAAC\tJJJJJJJJJJJJ\n"
            "test1\t147\ttestY\t9\t44\t12M\t=\t9\t-12\tCCCAACAACAAT\tJJJJJJJJJJJJ\n"
        )
        expected_pileup = {}
        pileup, _counts = remap.sam_to_pileup(samIO, max_primer_length=0)
        self.maxDiff = None
        self.assertEqual(pileup, expected_pileup)

class SamToConseqsTest(unittest.TestCase):
    def testSimple(self):
        #SAM:qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual
        samIO = StringIO.StringIO(
            "test1\t99\ttest\t1\t44\t12M\t=\t1\t12\tACAAGACCCAAC\tJJJJJJJJJJJJ\n"
            "test1\t147\ttest\t1\t44\t12M\t=\t1\t-12\tACAAGACCCAAC\tJJJJJJJJJJJJ\n"
        )
        expected_conseqs = {'test': 'ACAAGACCCAAC'}
        conseqs = remap.sam_to_conseqs(samIO)
        self.assertDictEqual(expected_conseqs, conseqs)
        
    def testOffset(self):
        samIO = StringIO.StringIO(
            "test1\t99\ttest\t4\t44\t12M\t=\t3\t12\tACAAGACCCAAC\tJJJJJJJJJJJJ\n"
            "test1\t147\ttest\t4\t44\t12M\t=\t3\t-12\tACAAGACCCAAC\tJJJJJJJJJJJJ\n"
        )
        expected_conseqs = {'test': 'NNNACAAGACCCAAC'}
        conseqs = remap.sam_to_conseqs(samIO)
        self.assertDictEqual(expected_conseqs, conseqs)
        
    def testHeaders(self):
        samIO = StringIO.StringIO(
            "@SH some header\n"
            "@AHI all headers are ignored\n"
            "test1\t99\ttest\t1\t44\t12M\t=\t1\t3\tACA\tJJJ\n"
            "test1\t147\ttest\t1\t44\t12M\t=\t1\t-3\tACA\tJJJ\n"
        )
        expected_conseqs = {'test': 'ACA'}
        conseqs = remap.sam_to_conseqs(samIO)
        self.assertDictEqual(expected_conseqs, conseqs)
        
    def testExtraFields(self):
        samIO = StringIO.StringIO(
            "test1\t99\ttest\t1\t44\t12M\t=\t1\t3\tACA\tJJJ\tAS:i:236\tNM:i:12\n"
            "test1\t147\ttest\t1\t44\t12M\t=\t1\t-3\tACA\tJJJ\tAS:i:236\tNM:i:12\n"
        )
        expected_conseqs = {'test': 'ACA'}
        conseqs = remap.sam_to_conseqs(samIO)
        self.assertDictEqual(expected_conseqs, conseqs)
        
    def testMaxConsensus(self):
        samIO = StringIO.StringIO(
            "test1\t99\ttest\t1\t44\t12M\t=\t1\t3\tACA\tJJJ\n"
            "test1\t147\ttest\t1\t44\t12M\t=\t1\t-3\tACA\tJJJ\n"
            "test2\t99\ttest\t1\t44\t12M\t=\t1\t3\tACA\tJJJ\n"
            "test2\t147\ttest\t1\t44\t12M\t=\t1\t-3\tACA\tJJJ\n"
            "test3\t99\ttest\t1\t44\t12M\t=\t1\t3\tTCA\tJJJ\n"
            "test3\t147\ttest\t1\t44\t12M\t=\t1\t-3\tTCA\tJJJ\n"
        )
        expected_conseqs = {'test': 'ACA'}
        conseqs = remap.sam_to_conseqs(samIO)
        self.assertDictEqual(expected_conseqs, conseqs)
 
    def testSoftClip(self):
        samIO = StringIO.StringIO(
            "test1\t99\ttest\t1\t44\t3S5M1S\t=\t1\t9\tACAGGGAGA\tJJJJJJJJJJJJ\n"
            "test1\t147\ttest\t1\t44\t3S5M1S\t=\t1\t-9\tACAGGGAGA\tJJJJJJJJJJJJ\n"
        )
        expected_conseqs = {'test': 'GGGAG'}
        conseqs = remap.sam_to_conseqs(samIO)
        self.assertDictEqual(expected_conseqs, conseqs)
 
    def testSimpleInsertion(self):
        samIO = StringIO.StringIO(
            "test1\t99\ttest\t1\t44\t3M3I3M\t=\t1\t9\tACAGGGAGA\tJJJJJJJJJJJJ\n"
            "test1\t147\ttest\t1\t44\t3M3I3M\t=\t1\t-9\tACAGGGAGA\tJJJJJJJJJJJJ\n"
        )
        expected_conseqs = {'test': 'ACAGGGAGA'}
        conseqs = remap.sam_to_conseqs(samIO)
        self.assertDictEqual(expected_conseqs, conseqs)
  
    def testInsertionAndOffset(self):
        samIO = StringIO.StringIO(
            "test1\t99\ttest\t1\t44\t3M3I3M\t=\t1\t9\tACAGGGAGA\tJJJJJJJJJJJJ\n"
            "test1\t147\ttest\t1\t44\t3M3I3M\t=\t1\t-9\tACAGGGAGA\tJJJJJJJJJJJJ\n"
            "test2\t99\ttest\t5\t44\t5M\t=\t1\t5\tGACCC\tJJJJJ\n"
            "test2\t147\ttest\t5\t44\t5M\t=\t1\t-5\tGACCC\tJJJJJ\n"
        )
        expected_conseqs = {'test': 'ACAGGGAGACCC'}
        conseqs = remap.sam_to_conseqs(samIO)
        self.assertDictEqual(expected_conseqs, conseqs)
 
    def testComplexInsertion(self):
        # Insertions are ignored if not a multiple of three
        samIO = StringIO.StringIO(
            "test1\t99\ttest\t1\t44\t3M1I3M2I6M\t=\t1\t12\tACAGAGAGGCCCAAC\tJJJJJJJJJJJJJJJ\n"
            "test1\t147\ttest\t1\t44\t3M1I3M2I6M\t=\t1\t-12\tACAGAGAGGCCCAAC\tJJJJJJJJJJJJJJJ\n"
        )
        expected_conseqs = {'test': 'ACAAGACCCAAC'}
        conseqs = remap.sam_to_conseqs(samIO)
        self.assertDictEqual(expected_conseqs, conseqs)
 
    def testDeletion(self):
        samIO = StringIO.StringIO(
            "test1\t99\ttest\t1\t44\t3M3D3M\t=\t3\t6\tACAGGG\tJJJJJJ\n"
            "test1\t147\ttest\t1\t44\t3M3D3M\t=\t3\t-6\tACAGGG\tJJJJJJ\n"
        )
        expected_conseqs = {'test': 'ACAGGG'}
        conseqs = remap.sam_to_conseqs(samIO)
        self.assertDictEqual(expected_conseqs, conseqs)
  
    def testOverlapsCountOnce(self):
        samIO = StringIO.StringIO(
            "test1\t99\ttest\t1\t44\t3M\t=\t1\t3\tACG\tJJJ\n"
            "test1\t147\ttest\t1\t44\t3M\t=\t1\t-3\tACG\tJJJ\n"
            "test2\t99\ttest\t1\t44\t3M\t=\t1\t3\tACG\tJJJ\n"
            "test2\t147\ttest\t1\t44\t3M\t=\t1\t-3\tACG\tJJJ\n"
            "test3\t99\ttest\t1\t44\t3M\t=\t3\t3\tATG\tJJJ\n"
            "test3\t147\ttest\t3\t44\t3M\t=\t1\t-3\tGCC\tJJJ\n"
            "test4\t99\ttest\t1\t44\t3M\t=\t3\t3\tATG\tJJJ\n"
            "test4\t147\ttest\t3\t44\t3M\t=\t1\t-3\tGCC\tJJJ\n"
            "test5\t99\ttest\t1\t44\t3M\t=\t3\t3\tATG\tJJJ\n"
            "test5\t147\ttest\t3\t44\t3M\t=\t1\t-3\tGCC\tJJJ\n"
        )
        expected_conseqs = {'test': 'ATGCC'}
        conseqs = remap.sam_to_conseqs(samIO)
        self.assertDictEqual(expected_conseqs, conseqs)
 
    def testPairMapsToTwoReferences(self):
        samIO = StringIO.StringIO(
            "test1\t99\ttestX\t1\t44\t3M\t=\t1\t3\tACG\tJJJ\n"
            "test1\t147\ttestY\t1\t44\t3M\t=\t1\t-3\tACG\tJJJ\n"
        )
        expected_conseqs = {}
        conseqs = remap.sam_to_conseqs(samIO)
        self.assertDictEqual(expected_conseqs, conseqs)
 
    def testLowQuality(self):
        # Note that we ignore the overlapped portion of the reverse read,
        # even if it has higher quality.
        samIO = StringIO.StringIO(
            "test1\t99\ttest\t1\t44\t3M\t=\t1\t3\tACG\tJ/J\n"
            "test1\t147\ttest\t1\t44\t3M\t=\t1\t-3\tACG\tJJJ\n"
        )
        expected_conseqs = {'test': 'ANG'}
        conseqs = remap.sam_to_conseqs(samIO, quality_cutoff=32)
        self.assertDictEqual(expected_conseqs, conseqs)
 
    def testDebugReports(self):
        samIO = StringIO.StringIO(
            "test1\t99\ttest\t1\t44\t3M3I9M\t=\t1\t12\tACTGGGAGACCCAAC\tJIJJJJJJJJJJJJJ\n"
            "test1\t147\ttest\t1\t44\t3M3I9M\t=\t1\t-12\tACTGGGAGACCCAAC\tJKJJJJJJJJJJJJJ\n"
            "test1\t99\ttest\t1\t44\t3M3I9M\t=\t1\t12\tATTGGGAGACCCAAC\tJHJJJJJJJJJJJJJ\n"
            "test1\t147\ttest\t1\t44\t3M3I9M\t=\t1\t-12\tATTGGGAGACCCAAC\tJFJJJJJJJJJJJJJ\n"
        )
        reports = {('test', 2): None}
        expected_reports = {('test', 2): 'H{C: 1, T: 1}, I{C: 1}'}
        
        remap.sam_to_conseqs(samIO, debug_reports=reports)
        
        self.assertDictEqual(expected_reports, reports)
