from io import StringIO
import unittest

from micall.core import remap
from micall.core.remap import is_first_read, is_short_read, \
    MixedReferenceSplitter


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


class SamToConseqsTest(unittest.TestCase):
    def testSimple(self):
        # SAM:qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual
        samIO = StringIO(
            "@SQ\tSN:test\n"
            "test1\t99\ttest\t1\t44\t12M\t=\t1\t12\tACAAGACCCAAC\tJJJJJJJJJJJJ\n"
        )
        expected_conseqs = {'test': 'ACAAGACCCAAC'}
        conseqs = remap.sam_to_conseqs(samIO)
        self.assertDictEqual(expected_conseqs, conseqs)

    def testOffset(self):
        samIO = StringIO(
            "@SQ\tSN:test\n"
            "test1\t147\ttest\t4\t44\t12M\t=\t3\t-12\tACAAGACCCAAC\tJJJJJJJJJJJJ\n"
        )
        expected_conseqs = {'test': 'NNNACAAGACCCAAC'}
        conseqs = remap.sam_to_conseqs(samIO)
        self.assertDictEqual(expected_conseqs, conseqs)

    def testHeaders(self):
        samIO = StringIO(
            "@SH\tsome header\n"
            "@MHI\tmost headers are ignored, except SQ for sequence reference\n"
            "@SQ\tSN:test\n"
            "test1\t99\ttest\t1\t44\t3M\t=\t1\t3\tACA\tJJJ\n"
        )
        expected_conseqs = {'test': 'ACA'}
        conseqs = remap.sam_to_conseqs(samIO)
        self.assertDictEqual(expected_conseqs, conseqs)

    def testUnknownReferenceName(self):
        samIO = StringIO(
            "@SQ\tSN:testX\n"
            "test1\t99\ttestY\t1\t44\t12M\t=\t1\t3\tACA\tJJJ\n"
        )
        expected_conseqs = {}
        conseqs = remap.sam_to_conseqs(samIO)
        self.assertDictEqual(expected_conseqs, conseqs)

    def testHeaderFields(self):
        samIO = StringIO(
            "@SQ\tOF:other field: ignored\tSN:test\n"
            "test1\t99\ttest\t1\t44\t3M\t=\t1\t3\tACA\tJJJ\n"
        )
        expected_conseqs = {'test': 'ACA'}
        conseqs = remap.sam_to_conseqs(samIO)
        self.assertDictEqual(expected_conseqs, conseqs)

    def testExtraFields(self):
        samIO = StringIO(
            "@SQ\tSN:test\n"
            "test1\t99\ttest\t1\t44\t3M\t=\t1\t3\tACA\tJJJ\tAS:i:236\tNM:i:12\n"
        )
        expected_conseqs = {'test': 'ACA'}
        conseqs = remap.sam_to_conseqs(samIO)
        self.assertDictEqual(expected_conseqs, conseqs)

    def testMaxConsensus(self):
        samIO = StringIO(
            "@SQ\tSN:test\n"
            "test1\t99\ttest\t1\t44\t3M\t=\t1\t3\tACA\tJJJ\n"
            "test2\t147\ttest\t1\t44\t3M\t=\t1\t-3\tACA\tJJJ\n"
            "test3\t99\ttest\t1\t44\t3M\t=\t1\t3\tTCA\tJJJ\n"
        )
        expected_conseqs = {'test': 'ACA'}
        conseqs = remap.sam_to_conseqs(samIO)
        self.assertDictEqual(expected_conseqs, conseqs)

    def testTie(self):
        samIO = StringIO(
            "@SQ\tSN:test\n"
            "test1\t99\ttest\t1\t44\t3M\t=\t1\t3\tGCA\tJJJ\n"
            "test2\t147\ttest\t1\t44\t3M\t=\t1\t-3\tTCA\tJJJ\n"
        )
        expected_conseqs = {'test': 'GCA'}
        conseqs = remap.sam_to_conseqs(samIO)
        self.assertDictEqual(expected_conseqs, conseqs)

    def testSoftClip(self):
        samIO = StringIO(
            "@SQ\tSN:test\n"
            "test1\t99\ttest\t1\t44\t3S5M1S\t=\t1\t9\tACAGGGAGA\tJJJJJJJJJ\n"
        )
        expected_conseqs = {'test': 'GGGAG'}
        conseqs = remap.sam_to_conseqs(samIO)
        self.assertDictEqual(expected_conseqs, conseqs)

    def testSimpleInsertion(self):
        samIO = StringIO(
            "@SQ\tSN:test\n"
            "test1\t99\ttest\t1\t44\t3M3I3M\t=\t1\t9\tACAGGGAGA\tJJJJJJJJJ\n"
        )
        expected_conseqs = {'test': 'ACAGGGAGA'}
        conseqs = remap.sam_to_conseqs(samIO)
        self.assertDictEqual(expected_conseqs, conseqs)

    def testLowQualityInsertion(self):
        samIO = StringIO(
            "@SQ\tSN:test\n"
            "test1\t99\ttest\t1\t44\t3M3I3M\t=\t1\t9\tACAGGGAGA\tJJJJ/JJJJ\n"
        )
        expected_conseqs = {'test': 'ACAAGA'}
        conseqs = remap.sam_to_conseqs(samIO, quality_cutoff=32)
        self.assertDictEqual(expected_conseqs, conseqs)

    def testInsertionAfterLowQuality(self):
        samIO = StringIO(
            "@SQ\tSN:test\n"
            "test1\t99\ttest\t1\t44\t3M3I3M\t=\t1\t9\tACAGGGAGA\tJJ/JJJJJJ\n"
        )
        expected_conseqs = {'test': 'ACNAGA'}
        conseqs = remap.sam_to_conseqs(samIO, quality_cutoff=32)
        self.assertDictEqual(expected_conseqs, conseqs)

    def testInsertionAndOffset(self):
        samIO = StringIO(
            "@SQ\tSN:test\n"
            "test1\t99\ttest\t1\t44\t3M3I3M\t=\t1\t9\tACAGGGAGA\tJJJJJJJJJJJJ\n"
            "test2\t99\ttest\t5\t44\t5M\t=\t1\t5\tGACCC\tJJJJJ\n"
        )
        expected_conseqs = {'test': 'ACAGGGAGACCC'}
        conseqs = remap.sam_to_conseqs(samIO)
        self.assertDictEqual(expected_conseqs, conseqs)

    def testComplexInsertion(self):
        # Insertions are ignored if not a multiple of three
        samIO = StringIO(
            "@SQ\tSN:test\n"
            "test1\t99\ttest\t1\t44\t3M1I3M2I6M\t=\t1\t12\tACAGAGAGGCCCAAC\tJJJJJJJJJJJJJJJ\n"
        )
        expected_conseqs = {'test': 'ACAAGACCCAAC'}
        conseqs = remap.sam_to_conseqs(samIO)
        self.assertDictEqual(expected_conseqs, conseqs)

    def testDeletion(self):
        samIO = StringIO(
            "@SQ\tSN:test\n"
            "test1\t99\ttest\t1\t44\t3M3D3M\t=\t3\t6\tACAGGG\tJJJJJJ\n"
        )
        expected_conseqs = {'test': 'ACAGGG'}
        conseqs = remap.sam_to_conseqs(samIO)
        self.assertDictEqual(expected_conseqs, conseqs)

    def testDeletionInSomeReads(self):
        samIO = StringIO(
            "@SQ\tSN:test\n"
            "test1\t99\ttest\t1\t44\t3M3D3M\t=\t3\t6\tACAGGG\tJJJJJJ\n"
            "test2\t99\ttest\t1\t44\t3M3D3M\t=\t3\t6\tACAGGG\tJJJJJJ\n"
            "test3\t99\ttest\t1\t44\t9M\t=\t3\t9\tACATTTGGG\tJJJJJJJJJ\n"
        )
        expected_conseqs = {'test': 'ACATTTGGG'}
        conseqs = remap.sam_to_conseqs(samIO)
        self.assertDictEqual(expected_conseqs, conseqs)

    def testDeletionWithFrameShift(self):
        samIO = StringIO(
            "@SQ\tSN:test\n"
            "test1\t99\ttest\t1\t44\t3M1D3M\t=\t3\t6\tACAGGG\tJJJJJJ\n"
        )
        expected_conseqs = {'test': 'ACA-GGG'}
        conseqs = remap.sam_to_conseqs(samIO)
        self.assertDictEqual(expected_conseqs, conseqs)

    def testBigDeletionWithFrameShift(self):
        samIO = StringIO(
            "@SQ\tSN:test\n"
            "test1\t99\ttest\t1\t44\t3M4D3M\t=\t3\t6\tACAGGG\tJJJJJJ\n"
        )
        expected_conseqs = {'test': 'ACA----GGG'}
        conseqs = remap.sam_to_conseqs(samIO)
        self.assertDictEqual(expected_conseqs, conseqs)

    def testOverlapsCountOnce(self):
        samIO = StringIO(
            "@SQ\tSN:test\n"
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

    def testReverseLeftOfForward(self):
        samIO = StringIO(
            "@SQ\tSN:test\n"
            "test1\t99\ttest\t2\t44\t1M\t=\t1\t1\tC\tJ\n"
            "test1\t147\ttest\t1\t44\t1M\t=\t2\t-1\tA\tJ\n"
        )
        expected_conseqs = {'test': 'AC'}
        conseqs = remap.sam_to_conseqs(samIO)
        self.assertDictEqual(expected_conseqs, conseqs)

    def testPairMapsToTwoReferences(self):
        samIO = StringIO(
            "@SQ\tSN:testX\n"
            "@SQ\tSN:testY\n"
            "test1\t99\ttestX\t1\t44\t3M\t=\t1\t3\tACG\tJJJ\n"
            "test1\t147\ttestY\t1\t44\t3M\t=\t1\t-3\tACG\tJJJ\n"
        )
        expected_conseqs = {}
        conseqs = remap.sam_to_conseqs(samIO)
        self.assertDictEqual(expected_conseqs, conseqs)

    def testLowQuality(self):
        # Note that we ignore the overlapped portion of the reverse read,
        # even if it has higher quality.
        samIO = StringIO(
            "@SQ\tSN:test\n"
            "test1\t99\ttest\t1\t44\t3M\t=\t1\t3\tACG\tJ/J\n"
        )
        expected_conseqs = {'test': 'ANG'}
        conseqs = remap.sam_to_conseqs(samIO, quality_cutoff=32)
        self.assertDictEqual(expected_conseqs, conseqs)

    def testLowQualityAtEnd(self):
        samIO = StringIO(
            "@SQ\tSN:test\n"
            "test1\t99\ttest\t1\t44\t3M\t=\t1\t3\tACG\tJJ/\n"
        )
        expected_conseqs = {'test': 'ACN'}
        conseqs = remap.sam_to_conseqs(samIO, quality_cutoff=32)
        self.assertDictEqual(expected_conseqs, conseqs)

    def testLowQualityForward(self):
        samIO = StringIO(
            "@SQ\tSN:test\n"
            "test1\t99\ttest\t1\t44\t3M\t=\t3\t3\tATA\tJJA\n"
            "test1\t147\ttest\t3\t44\t3M\t=\t1\t-3\tGCC\tJJJ\n"
        )
        expected_conseqs = {'test': 'ATGCC'}
        conseqs = remap.sam_to_conseqs(samIO)
        self.assertDictEqual(expected_conseqs, conseqs)

    def testAllLowQuality(self):
        # SAM:qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual
        samIO = StringIO(
            "@SQ\tSN:test\n"
            "test1\t147\ttest\t1\t24\t1M\t=\t1\t-1\tT\t#\n"
        )
        expected_conseqs = {}

        conseqs = remap.sam_to_conseqs(samIO, quality_cutoff=32)

        self.assertDictEqual(expected_conseqs, conseqs)

    def testBadPairFlag(self):
        """ Even if the pair isn't concordant, still include in consensus.

        SAM flag 145 does not have bit 2 for properly aligned.
        """
        # SAM:qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual
        samIO = StringIO(
            "@SQ\tSN:test\n"
            "test1\t145\ttest\t1\t24\t1M\t=\t1\t-1\tT\tF\n"
        )
        expected_conseqs = {'test': 'T'}

        conseqs = remap.sam_to_conseqs(samIO, quality_cutoff=32)

        self.assertEqual(expected_conseqs, conseqs)

    def testUnmappedFlag(self):
        """ If the read is unmapped, don't include in consensus.

        SAM flag 149 has bit 4 for unmapped. Region is irrelevant.
        """
        # SAM:qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual
        samIO = StringIO(
            "@SQ\tSN:test\n"
            "test1\t149\ttest\t1\t24\t1M\t=\t1\t-1\tT\tF\n"
        )
        expected_conseqs = {}

        conseqs = remap.sam_to_conseqs(samIO, quality_cutoff=32)

        self.assertEqual(expected_conseqs, conseqs)

    def testLowQualityAndDeletion(self):
        samIO = StringIO(
            "@SQ\tSN:test\n"
            "test1\t99\ttest\t1\t44\t3M3D3M\t=\t3\t6\tACAGGG\tJJJJJJ\n"
            "test2\t99\ttest\t1\t44\t9M\t=\t3\t9\tACATTTGGG\tJJJ///JJJ\n"
        )
        expected_conseqs = {'test': 'ACANNNGGG'}

        conseqs = remap.sam_to_conseqs(samIO, quality_cutoff=32)

        self.assertDictEqual(expected_conseqs, conseqs)

    def testSeeds(self):
        samIO = StringIO(
            "@SQ\tSN:test\n"
            "test1\t99\ttest\t4\t44\t3M\t=\t10\t3\tTAT\tJJJ\n"
            "test2\t99\ttest\t10\t44\t3M\t=\t4\t-3\tCAC\tJJJ\n"
        )
        seeds = {'test': 'ACATTTGGGCAC'}
        expected_conseqs = {'test': 'ACATATGGGCAC'}

        conseqs = remap.sam_to_conseqs(samIO, seeds=seeds)

        self.assertDictEqual(expected_conseqs, conseqs)

    def testSeedsNeedSomeReads(self):
        samIO = StringIO(
            "@SQ\tSN:test\n"
            "test1\t99\ttest\t4\t44\t3M\t=\t10\t3\tTAT\tJJJ\n"
        )
        seeds = {'test': 'ACATTTGGGCAC',
                 'other': 'TATGCACCC'}
        expected_conseqs = {'test': 'ACATATGGGCAC'}

        conseqs = remap.sam_to_conseqs(samIO, seeds=seeds)

        self.assertDictEqual(expected_conseqs, conseqs)

    def testSeedsWithLowQuality(self):
        samIO = StringIO(
            "@SQ\tSN:test\n"
            "test1\t99\ttest\t4\t44\t3M\t=\t10\t3\tTAT\tJJ/\n"
        )
        seeds = {'test': 'ACATTTGGGCAC'}
        expected_conseqs = {'test': 'ACATATGGGCAC'}

        conseqs = remap.sam_to_conseqs(samIO,
                                       seeds=seeds,
                                       quality_cutoff=32)

        self.assertDictEqual(expected_conseqs, conseqs)

    def testDebugReports(self):
        samIO = StringIO(
            "@SQ\tSN:test\n"
            "test1\t99\ttest\t1\t44\t3M3I9M\t=\t1\t12\tACTGGGAGACCCAAC\tJIJJJJJJJJJJJJJ\n"
            "test1\t147\ttest\t1\t44\t3M3I9M\t=\t1\t-12\tACTGGGAGACCCAAC\tJIJJJJJJJJJJJJJ\n"
            "test1\t99\ttest\t1\t44\t3M3I9M\t=\t1\t12\tATTGGGAGACCCAAC\tJHJJJJJJJJJJJJJ\n"
            "test1\t147\ttest\t1\t44\t3M3I9M\t=\t1\t-12\tATTGGGAGACCCAAC\tJHJJJJJJJJJJJJJ\n"
        )
        reports = {('test', 2): None}
        expected_reports = {('test', 2): 'H{C: 1, T: 1}, I{C: 1}'}

        remap.sam_to_conseqs(samIO, debug_reports=reports)

        self.assertDictEqual(expected_reports, reports)

    def testDebugReportsOnReverseRead(self):
        samIO = StringIO(
            "@SQ\tSN:test\n"
            "test1\t99\ttest\t1\t44\t3M3I2M\t=\t1\t8\tACTGGGAG\tJJJJJJJJ\n"
            "test1\t147\ttest\t5\t44\t8M\t=\t1\t-8\tGACCCAAC\tJJJJJIJJ\n"
            "test1\t99\ttest\t1\t44\t3M3I2M\t=\t1\t12\tATTGGGAG\tJJJJJJJJ\n"
            "test1\t147\ttest\t5\t44\t8M\t=\t1\t-12\tGACCCAAC\tJJJJJHJJ\n"
        )
        reports = {('test', 10): None}
        expected_reports = {('test', 10): 'H{A: 2}, I{A: 1}'}

        remap.sam_to_conseqs(samIO, debug_reports=reports)

        self.assertDictEqual(expected_reports, reports)

    def testSeedsConverged(self):
        # SAM:qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual
        samIO = StringIO(
            "@SQ\tSN:test\tSN:other\tSN:wayoff\n"
            "test1\t99\ttest\t1\t44\t10M\t=\t1\t10\tATGAGGAGTA\tJJJJJJJJJJJJ\n"
            "other1\t99\tother\t1\t44\t10M\t=\t1\t10\tATGACCAGTA\tJJJJJJJJJJJJ\n"
            "wayoff1\t99\twayoff\t1\t44\t10M\t=\t1\t10\tATGAGGGTAC\tJJJJJJJJJJJJ\n"
        )
        seeds = {'test': 'ATGAAGTA',
                 'other': 'AAGCCGAA',
                 'wayoff': 'TCATGTAC'}
        expected_conseqs = {'test': 'ATGAGGAGTA'}
        expected_distances = {'test': dict(seed_dist=2,
                                           other_dist=5,
                                           other_seed='other'),
                              'other': dict(seed_dist=4,
                                            other_dist=2,
                                            other_seed='test'),
                              'wayoff': dict(seed_dist=4,
                                             other_dist=3,
                                             other_seed='test')}
        distances = {}

        conseqs = remap.sam_to_conseqs(samIO,
                                       seeds=seeds,
                                       is_filtered=True,
                                       distance_report=distances)

        self.maxDiff = 1000
        self.assertEqual(expected_conseqs, conseqs)
        self.assertEqual(expected_distances, distances)

    def testSeedsConvergedWithDifferentAlignment(self):
        """ Seeds have similar regions, but at different positions.
        """
        samIO = StringIO(
            "@SQ\tSN:test\tSN:other\n"
            "test1\t99\ttest\t1\t44\t10M\t=\t1\t10\tATGAGGAGTA\tJJJJJJJJJJJJ\n"
            "other1\t99\tother\t11\t44\t10M\t=\t1\t10\tATGACCAGTA\tJJJJJJJJJJJJ\n"
        )
        seeds = {'test': 'ATGAAGTA',
                 'other': 'TCTCTCTCTCAAGCCGAA'}
        expected_conseqs = {'test': 'ATGAGGAGTA'}

        conseqs = remap.sam_to_conseqs(samIO,
                                       seeds=seeds,
                                       is_filtered=True)

        self.assertDictEqual(expected_conseqs, conseqs)

    def testSeedsConvergedWithDifferentAlignmentAndGap(self):
        """ Gaps between areas with coverage.
        """
        samIO = StringIO(
            "@SQ\tSN:test\tSN:other\n"
            "test1\t99\ttest\t1\t44\t10M\t=\t1\t10\tATGAGGAGTA\tJJJJJJJJJJJJ\n"
            "other1\t99\tother\t11\t44\t5M\t=\t1\t5\tATGAC\tJJJJJJJ\n"
            "other2\t99\tother\t26\t44\t5M\t=\t1\t5\tCAGTA\tJJJJJJJ\n"
        )
        seeds = {'test': 'ATGAAGTA',
                 'other': 'TCTCTCTCTCAAGCTATATATATACGAA'}
        expected_conseqs = {'test': 'ATGAGGAGTA'}

        conseqs = remap.sam_to_conseqs(samIO,
                                       seeds=seeds,
                                       is_filtered=True)

        self.assertDictEqual(expected_conseqs, conseqs)

    def testSeedsConvergedWithConfusingGap(self):
        """ Reads match other seed well, but with a big gap.
        """
        samIO = StringIO(
            "@SQ\tSN:test\tSN:other\n"
            "test1\t99\ttest\t1\t44\t8M\t=\t1\t8\tATGTCGTA\tJJJJJJJJ\n"
            "other1\t99\tother\t14\t44\t9M\t=\t1\t9\tAAGCTATAT\tJJJJJJJJJ\n"
        )
        seeds = {'test': 'ATGAAGTA',
                 'other': 'ATGTCTCTCTCTCAAGCTATATATATACGAAGTA'}
        expected_conseqs = {'test': 'ATGTCGTA',
                            'other': 'ATGTCTCTCTCTCAAGCTATATATATACGAAGTA'}

        conseqs = remap.sam_to_conseqs(samIO,
                                       seeds=seeds,
                                       is_filtered=True)

        self.assertDictEqual(expected_conseqs, conseqs)

    def testSeedsConvergedPlusOtherLowCoverage(self):
        """ Portion with decent coverage has converged, other hasn't.
        """
        samIO = StringIO(
            "@SQ\tSN:test\tSN:other\n"
            "test1\t99\ttest\t1\t44\t10M\t=\t1\t10\tATGAGGAGTA\tJJJJJJJJJJJJ\n"
            "test2\t99\ttest\t1\t44\t10M\t=\t1\t10\tATGAGGAGTA\tJJJJJJJJJJJJ\n"
            "other1\t99\tother\t1\t44\t10M\t=\t1\t10\tATGACCAGTA\tJJJJJJJJJJJJ\n"
            "other2\t99\tother\t1\t44\t10M\t=\t1\t10\tATGACCAGTA\tJJJJJJJJJJJJ\n"
            "other3\t99\tother\t11\t44\t6M\t=\t1\t16\tGTGTGT\tJJJJJJ\n"
        )
        seeds = {'test': 'ATGAAGTACTCTCT',
                 'other': 'AAGCCGAAGTGTGT'}
        expected_conseqs = {'test': 'ATGAGGAGTACTCT'}

        conseqs = remap.sam_to_conseqs(samIO,
                                       seeds=seeds,
                                       is_filtered=True,
                                       filter_coverage=2)

        self.assertDictEqual(expected_conseqs, conseqs)

    def testSeedsBothConverged(self):
        """ Both references are now closer to the other seed than the start.

        Don't drop both. Keep test because it has more reads.
        """
        # SAM:qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual
        samIO = StringIO(
            "@SQ\tSN:test\tSN:other\tSN:unrelated\n"
            "test1\t99\ttest\t1\t44\t8M\t=\t1\t10\tAAGCCGTA\tJJJJJJJJJJ\n"
            "test2\t99\ttest\t1\t44\t8M\t=\t1\t10\tAAGCCGTA\tJJJJJJJJJJ\n"
            "other1\t99\tother\t1\t44\t8M\t=\t1\t10\tATGAAGTA\tJJJJJJJJJJ\n"
            "unrelated1\t99\tunrelated\t1\t44\t9M\t=\t1\t10\tGGGTTTGGG\tJJJJJJJJJ\n"
        )
        seeds = {'test': 'ATGAAGTA',
                 'other': 'AAGCCGAA',
                 'unrelated': 'GGGTTTGGG'}
        expected_conseqs = {'test': 'AAGCCGTA',
                            'unrelated': 'GGGTTTGGG'}
        expected_distances = {'test': dict(seed_dist=3,
                                           other_dist=1,
                                           other_seed='other'),
                              'other': dict(seed_dist=5,
                                            other_dist=0,
                                            other_seed='test'),
                              'unrelated': dict(seed_dist=0,
                                                other_dist=7,
                                                other_seed='test')}
        distances = {}

        conseqs = remap.sam_to_conseqs(samIO,
                                       seeds=seeds,
                                       is_filtered=True,
                                       distance_report=distances)

        self.maxDiff = 1000
        self.assertEqual(expected_distances, distances)
        self.assertEqual(expected_conseqs, conseqs)

    def testAllSeedsLowCoverage(self):
        "Multiple seeds mapped, but none have good coverage. Choose most reads."

        samIO = StringIO(
            "@SQ\tSN:test\tSN:other\n"
            "test1\t99\ttest\t1\t44\t10M\t=\t1\t10\tATGAGGAGTA\tJJJJJJJJJJJJ\n"
            "test2\t99\ttest\t11\t44\t6M\t=\t1\t10\tCTCTCT\tJJJJJJ\n"
            "other1\t99\tother\t1\t44\t10M\t=\t1\t10\tATGACCAGTA\tJJJJJJJJJJJJ\n"
        )
        seeds = {'test': 'ATGAAGTACTCTCT',
                 'other': 'AAGCCGAAGTGTGT'}
        expected_conseqs = {'test': 'ATGAGGAGTACTCTCT'}

        conseqs = remap.sam_to_conseqs(samIO,
                                       seeds=seeds,
                                       is_filtered=True,
                                       filter_coverage=2)

        self.assertDictEqual(expected_conseqs, conseqs)

    def testExtractRelevantSeeds(self):
        expectations = [  # (aligned_conseq, aligned_seed, expected_seed)
            ('ACTG',
             'ATTG',
             'ATTG'),
            ('-ACTG-',
             'CATTGT',
             'ATTG'),
            ('-AC-TG--',
             'CATATGT',
             'ATATG'),
            ('-AC-TG-AT-',
             'CATATGTATC',
             'ATATGTAT'),
            ('--T--',
             'CATAT',
             'T'),
            ('TACG----',
             '----GGCC',
             '')]
        for aligned_conseq, aligned_seed, expected_seed in expectations:
            relevant = remap.extract_relevant_seed(aligned_conseq, aligned_seed)
            self.assertEqual((aligned_conseq, aligned_seed, expected_seed),
                             (aligned_conseq, aligned_seed, relevant))

    def testNothingMapped(self):
        samIO = StringIO(
            "@SQ\tSN:test\tSN:other\n"
        )
        seeds = {'test': 'ATGAAGTACTCTCT',
                 'other': 'AAGCCGAAGTGTGT'}
        expected_conseqs = {}

        conseqs = remap.sam_to_conseqs(samIO,
                                       seeds=seeds,
                                       is_filtered=True,
                                       filter_coverage=2)

        self.assertDictEqual(expected_conseqs, conseqs)


class MixedReferenceMemorySplitter(MixedReferenceSplitter):
    """ Dummy class to hold split reads in memory. Useful for testing. """
    def create_split_file(self, refname, direction):
        self.is_closed = False
        return StringIO()

    def close_split_file(self, split_file):
        self.is_closed = True


class MixedReferenceSplitterTest(unittest.TestCase):
    def setUp(self):
        super(MixedReferenceSplitterTest, self).setUp()
        self.addTypeEqualityFunc(str, self.assertMultiLineEqual)

    def testSimple(self):
        samIO = StringIO(
            "@SQ\tSN:r\n"
            "r1\t99\tr\t1\t44\t3M\t=\t1\t3\tACG\tJJJ\n"
            "r1\t147\tr\t1\t44\t3M\t=\t1\t-3\tACG\tJJJ\n")
        expected_rows = [
            ["r1", "99", "r", "1", "44", "3M", "=", "1", "3", "ACG", "JJJ"],
            ["r1", "147", "r", "1", "44", "3M", "=", "1", "-3", "ACG", "JJJ"]]

        splitter = MixedReferenceSplitter()
        rows = list(splitter.split(samIO))

        self.assertEqual(expected_rows, rows)

    def testTrimOptionalFields(self):
        samIO = StringIO(
            "@SQ\tSN:r\n"
            "r1\t99\tr\t1\t44\t3M\t=\t1\t3\tACG\tJJJ\tAS:i:100\n"
            "r1\t147\tr\t1\t44\t3M\t=\t1\t-3\tACG\tJJJ\tYS:Z:UP\n")
        expected_rows = [
            ["r1", "99", "r", "1", "44", "3M", "=", "1", "3", "ACG", "JJJ"],
            ["r1", "147", "r", "1", "44", "3M", "=", "1", "-3", "ACG", "JJJ"]]

        splitter = MixedReferenceSplitter()
        rows = list(splitter.split(samIO))

        self.assertEqual(expected_rows, rows)

    def testUnmapped(self):
        samIO = StringIO(
            "@SQ\tSN:r\n"
            "r1\t107\tr\t1\t44\t3M\t*\t1\t3\tACG\tJJJ\n"
            "r1\t149\t*\t*\t*\t*\tr\t*\t*\tACG\tJJJ\n")
        expected_rows = [
            ["r1", "107", "r", "1", "44", "3M", "*", "1", "3", "ACG", "JJJ"],
            ["r1", "149", "*", "*", "*", "*", "r", "*", "*", "ACG", "JJJ"]]

        splitter = MixedReferenceSplitter()
        rows = list(splitter.split(samIO))

        self.assertEqual(expected_rows, rows)

    def testSplit(self):
        """ If a pair is split over two references, choose one.

        Use the reference that gave the higher mapq score.
        """
        samIO = StringIO(
            "@SQ\tSN:r\n"
            "r1\t99\tRX\t1\t44\t3M\t=\t1\t3\tACG\tJJJ\tAS:i:100\n"
            "r1\t147\tRX\t1\t44\t3M\t=\t1\t-3\tACG\tJJJ\tYS:Z:UP\n"
            "r2\t99\tRX\t1\t44\t3M\tRY\t1\t3\tACG\tJJJ\tAS:i:100\n"
            "r2\t147\tRY\t1\t11\t3M\tRX\t1\t-3\tACC\tJJK\tAS:i:200\n")
        expected_rows = [
            ["r1", "99", "RX", "1", "44", "3M", "=", "1", "3", "ACG", "JJJ"],
            ["r1", "147", "RX", "1", "44", "3M", "=", "1", "-3", "ACG", "JJJ"]]
        expected_fastq1 = """\
@r2
ACG
+
JJJ
"""
        expected_fastq2 = """\
@r2
GGT
+
KJJ
"""

        splitter = MixedReferenceMemorySplitter()
        rows = list(splitter.split(samIO))
        is_closed = splitter.is_closed
        splits = splitter.splits

        self.assertEqual(expected_rows, rows)
        self.assertEqual(['RX'], list(splits.keys()))
        fastq1, fastq2 = splits['RX']
        self.assertEqual(expected_fastq1, fastq1.getvalue())
        self.assertEqual(expected_fastq2, fastq2.getvalue())
        self.assertTrue(is_closed)

    def testTiedMapQ(self):
        """ If both mates have the same mapq, choose higher alignment score.

        Use the reference that gave the higher mapq score.
        """
        samIO = StringIO(
            "@SQ\tSN:r\n"
            "r1\t99\tRX\t1\t44\t3M\t=\t1\t3\tACG\tJJJ\tAS:i:100\n"
            "r1\t147\tRX\t1\t44\t3M\t=\t1\t-3\tACG\tJJJ\tYS:Z:UP\n"
            "r2\t99\tRX\t1\t11\t3M\tRY\t1\t3\tACG\tJJJ\tAS:i:100\n"
            "r2\t147\tRY\t1\t11\t3M\tRX\t1\t-3\tACC\tJJK\tAS:i:200\n")

        splitter = MixedReferenceMemorySplitter()
        list(splitter.split(samIO))
        splits = splitter.splits

        self.assertEqual(['RY'], list(splits.keys()))

    def testWalk(self):
        samIO = StringIO(
            "@SQ\tSN:r\n"
            "r1\t99\tRX\t1\t44\t3M\t=\t1\t3\tACG\tJJJ\tAS:i:100\n"
            "r1\t147\tRX\t1\t44\t3M\t=\t1\t-3\tACG\tJJJ\tYS:Z:UP\n"
            "r2\t99\tRX\t1\t44\t3M\tRY\t1\t3\tACG\tJJJ\tAS:i:100\n"
            "r2\t147\tRY\t1\t44\t3M\tRX\t1\t-3\tACT\tKKK\tAS:i:200\n")
        expected_rows = [
            ["r1", "99", "RX", "1", "44", "3M", "=", "1", "3", "ACG", "JJJ", "AS:i:100"],
            ["r1", "147", "RX", "1", "44", "3M", "=", "1", "-3", "ACG", "JJJ", "YS:Z:UP"],
            ["r2", "99", "RX", "1", "44", "3M", "RY", "1", "3", "ACG", "JJJ", "AS:i:100"],
            ["r2", "147", "RY", "1", "44", "3M", "RX", "1", "-3", "ACT", "KKK", "AS:i:200"]]

        splitter = MixedReferenceMemorySplitter()
        rows = list(splitter.walk(samIO))
        splits = splitter.splits

        self.assertEqual(expected_rows, rows)
        self.assertEqual({}, splits)
