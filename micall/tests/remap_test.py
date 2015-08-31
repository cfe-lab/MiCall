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
        
    def testDeletionWithoutStars(self):
        pileupIO = StringIO.StringIO(
            "GP41-seed\t1\tN\t1\t^Kg-1n\tC\n" +
            "GP41-seed\t2\tN\t0\n" +
            "GP41-seed\t3\tN\t0\n")
        qCutoff = 20
        expected_conseqs = {'GP41-seed': "G-N"}
        
        conseqs= remap.pileup_to_conseq(pileupIO, qCutoff)
        
        self.assertDictEqual(expected_conseqs, conseqs)
        
    def testDeletionAfterLowQuality(self):
        pileupIO = StringIO.StringIO(
            "GP41-seed\t1\tN\t1\t^Kg-3nnn\t#\n" +
            "GP41-seed\t2\tN\t1\t*\tC\n" +
            "GP41-seed\t3\tN\t1\t*\tC\n" +
            "GP41-seed\t4\tN\t1\t*\tC\n" +
            "GP41-seed\t5\tN\t1\ta$\tC\n")
        qCutoff = 20
        expected_conseqs = {'GP41-seed': "NA"}
        
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
            "test1\t99\ttest\t1\t44\t3S5M1S\t=\t1\t9\tACAGGGAGA\tJJJJJJJJJ\n"
        )
        expected_conseqs = {'test': 'GGGAG'}
        conseqs = remap.sam_to_conseqs(samIO)
        self.assertDictEqual(expected_conseqs, conseqs)
 
    def testSimpleInsertion(self):
        samIO = StringIO.StringIO(
            "test1\t99\ttest\t1\t44\t3M3I3M\t=\t1\t9\tACAGGGAGA\tJJJJJJJJJ\n"
        )
        expected_conseqs = {'test': 'ACAGGGAGA'}
        conseqs = remap.sam_to_conseqs(samIO)
        self.assertDictEqual(expected_conseqs, conseqs)
 
    def testLowQualityInsertion(self):
        samIO = StringIO.StringIO(
            "test1\t99\ttest\t1\t44\t3M3I3M\t=\t1\t9\tACAGGGAGA\tJJJJ/JJJJ\n"
        )
        expected_conseqs = {'test': 'ACAAGA'}
        conseqs = remap.sam_to_conseqs(samIO, quality_cutoff=32)
        self.assertDictEqual(expected_conseqs, conseqs)
 
    def testInsertionAfterLowQuality(self):
        samIO = StringIO.StringIO(
            "test1\t99\ttest\t1\t44\t3M3I3M\t=\t1\t9\tACAGGGAGA\tJJ/JJJJJJ\n"
        )
        expected_conseqs = {'test': 'ACNAGA'}
        conseqs = remap.sam_to_conseqs(samIO, quality_cutoff=32)
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
 
    def testDeletionInSomeReads(self):
        samIO = StringIO.StringIO(
            "test1\t99\ttest\t1\t44\t3M3D3M\t=\t3\t6\tACAGGG\tJJJJJJ\n"
            "test1\t147\ttest\t1\t44\t3M3D3M\t=\t3\t-6\tACAGGG\tJJJJJJ\n"
            "test2\t99\ttest\t1\t44\t3M3D3M\t=\t3\t6\tACAGGG\tJJJJJJ\n"
            "test2\t147\ttest\t1\t44\t3M3D3M\t=\t3\t-6\tACAGGG\tJJJJJJ\n"
            "test3\t99\ttest\t1\t44\t9M\t=\t3\t9\tACATTTGGG\tJJJJJJJJJ\n"
            "test3\t147\ttest\t1\t44\t9M\t=\t3\t-9\tACATTTGGG\tJJJJJJJJJ\n"
        )
        expected_conseqs = {'test': 'ACATTTGGG'}
        conseqs = remap.sam_to_conseqs(samIO)
        self.assertDictEqual(expected_conseqs, conseqs)
  
    def testDeletionWithFrameShift(self):
        samIO = StringIO.StringIO(
            "test1\t99\ttest\t1\t44\t3M1D3M\t=\t3\t6\tACAGGG\tJJJJJJ\n"
            "test1\t147\ttest\t1\t44\t3M1D3M\t=\t3\t-6\tACAGGG\tJJJJJJ\n"
        )
        expected_conseqs = {'test': 'ACA-GGG'}
        conseqs = remap.sam_to_conseqs(samIO)
        self.assertDictEqual(expected_conseqs, conseqs)
  
#     def testOverlapsCountOnce(self):
#         samIO = StringIO.StringIO(
#             "test1\t99\ttest\t1\t44\t3M\t=\t1\t3\tACG\tJJJ\n"
#             "test1\t147\ttest\t1\t44\t3M\t=\t1\t-3\tACG\tJJJ\n"
#             "test2\t99\ttest\t1\t44\t3M\t=\t1\t3\tACG\tJJJ\n"
#             "test2\t147\ttest\t1\t44\t3M\t=\t1\t-3\tACG\tJJJ\n"
#             "test3\t99\ttest\t1\t44\t3M\t=\t3\t3\tATG\tJJJ\n"
#             "test3\t147\ttest\t3\t44\t3M\t=\t1\t-3\tGCC\tJJJ\n"
#             "test4\t99\ttest\t1\t44\t3M\t=\t3\t3\tATG\tJJJ\n"
#             "test4\t147\ttest\t3\t44\t3M\t=\t1\t-3\tGCC\tJJJ\n"
#             "test5\t99\ttest\t1\t44\t3M\t=\t3\t3\tATG\tJJJ\n"
#             "test5\t147\ttest\t3\t44\t3M\t=\t1\t-3\tGCC\tJJJ\n"
#         )
#         expected_conseqs = {'test': 'ATGCC'}
#         conseqs = remap.sam_to_conseqs(samIO)
#         self.assertDictEqual(expected_conseqs, conseqs)
#   
    def testOverlapsCountTwice(self):
        #TODO: Switch to counting once after we remove samtools
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
        expected_conseqs = {'test': 'ACGCC'}
        conseqs = remap.sam_to_conseqs(samIO)
        self.assertDictEqual(expected_conseqs, conseqs)
  
    def testReverseLeftOfForward(self):
        samIO = StringIO.StringIO(
            "test1\t99\ttest\t2\t44\t1M\t=\t1\t1\tC\tJ\n"
            "test1\t147\ttest\t1\t44\t1M\t=\t2\t-1\tA\tJ\n"
        )
        expected_conseqs = {'test': 'AC'}
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
            "test1\t147\ttest\t1\t44\t3M\t=\t1\t-3\tACG\tJ/J\n"
        )
        expected_conseqs = {'test': 'ANG'}
        conseqs = remap.sam_to_conseqs(samIO, quality_cutoff=32)
        self.assertDictEqual(expected_conseqs, conseqs)
 
    def testLowQualityAtEnd(self):
        samIO = StringIO.StringIO(
            "test1\t99\ttest\t1\t44\t3M\t=\t1\t3\tACG\tJJ/\n"
            "test1\t147\ttest\t1\t44\t3M\t=\t1\t-3\tACG\tJJ/\n"
        )
        expected_conseqs = {'test': 'ACN'}
        conseqs = remap.sam_to_conseqs(samIO, quality_cutoff=32)
        self.assertDictEqual(expected_conseqs, conseqs)
 
    def testAllLowQuality(self):
        #SAM:qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual
        samIO = StringIO.StringIO(
            "test1\t147\tRT-seed\t1\t24\t1M\t=\t1\t-1\tT\t#\n"
        )
        expected_conseqs = {}
        
        conseqs = remap.sam_to_conseqs(samIO, quality_cutoff=32)
        
        self.assertDictEqual(expected_conseqs, conseqs)
 
    def testBadAlignmentFlag(self):
        #SAM:qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual
        samIO = StringIO.StringIO(
            "test1\t145\tRT-seed\t1\t24\t1M\t=\t1\t-1\tT\tF\n"
        )
        expected_conseqs = {}
        
        conseqs = remap.sam_to_conseqs(samIO, quality_cutoff=32)
        
        self.assertDictEqual(expected_conseqs, conseqs)
 
    def testLowQualityAndDeletion(self):
        samIO = StringIO.StringIO(
            "test1\t99\ttest\t1\t44\t3M3D3M\t=\t3\t6\tACAGGG\tJJJJJJ\n"
            "test1\t147\ttest\t1\t44\t3M3D3M\t=\t3\t-6\tACAGGG\tJJJJJJ\n"
            "test2\t99\ttest\t1\t44\t9M\t=\t3\t9\tACATTTGGG\tJJJ///JJJ\n"
            "test2\t147\ttest\t1\t44\t9M\t=\t3\t-9\tACATTTGGG\tJJJ///JJJ\n"
        )
        expected_conseqs = {'test': 'ACANNNGGG'}

        conseqs = remap.sam_to_conseqs(samIO, quality_cutoff=32)

        self.assertDictEqual(expected_conseqs, conseqs)
 
    def testDebugReports(self):
        samIO = StringIO.StringIO(
            "test1\t99\ttest\t1\t44\t3M3I9M\t=\t1\t12\tACTGGGAGACCCAAC\tJIJJJJJJJJJJJJJ\n"
            "test1\t147\ttest\t1\t44\t3M3I9M\t=\t1\t-12\tACTGGGAGACCCAAC\tJIJJJJJJJJJJJJJ\n"
            "test1\t99\ttest\t1\t44\t3M3I9M\t=\t1\t12\tATTGGGAGACCCAAC\tJHJJJJJJJJJJJJJ\n"
            "test1\t147\ttest\t1\t44\t3M3I9M\t=\t1\t-12\tATTGGGAGACCCAAC\tJHJJJJJJJJJJJJJ\n"
        )
        reports = {('test', 2): None}
        expected_reports = {('test', 2): 'H{C: 2, T: 2}, I{C: 2}'}
        
        remap.sam_to_conseqs(samIO, debug_reports=reports)
        
        self.assertDictEqual(expected_reports, reports)
