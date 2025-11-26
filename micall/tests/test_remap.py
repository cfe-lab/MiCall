from csv import DictWriter
from io import StringIO
import os
import unittest
from pathlib import Path
from unittest.mock import patch, Mock, DEFAULT

from pytest import fixture

from micall.core import remap
from micall.core.project_config import ProjectConfig
from micall.core.remap import is_first_read, is_short_read, \
    MixedReferenceSplitter, write_remap_counts, convert_prelim, read_contigs
from micall.utils.externals import Bowtie2, Bowtie2Build

HXB2_NAME = "HIV1-B-FR-K03455-seed"


@fixture(name='projects', scope="session")
def load_projects():
    yield ProjectConfig.loadDefault()


class IsShortReadTest(unittest.TestCase):
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
        is_first_expected = True

        is_first = is_first_read(flag)

        self.assertEqual(is_first_expected, is_first)

    def testSecondRead(self):
        flag = '147'
        is_first_expected = False

        is_first = is_first_read(flag)

        self.assertEqual(is_first_expected, is_first)

    def testSmallFlag(self):
        flag = '3'
        is_first_expected = False

        is_first = is_first_read(flag)

        self.assertEqual(is_first_expected, is_first)


class SamToConseqsTest(unittest.TestCase):
    def testSimple(self):
        # SAM:qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual
        sam_file = StringIO(
            "@SQ\tSN:test\n"
            "test1\t99\ttest\t1\t44\t12M\t=\t1\t12\tACAAGACCCAAC\tJJJJJJJJJJJJ\n"
        )
        expected_conseqs = {'test': 'ACAAGACCCAAC'}
        conseqs = remap.sam_to_conseqs(sam_file)
        self.assertDictEqual(expected_conseqs, conseqs)

    def testOffset(self):
        sam_file = StringIO(
            "@SQ\tSN:test\n"
            "test1\t147\ttest\t4\t44\t12M\t=\t3\t-12\tACAAGACCCAAC\tJJJJJJJJJJJJ\n"
        )
        expected_conseqs = {'test': 'NNNACAAGACCCAAC'}
        conseqs = remap.sam_to_conseqs(sam_file)
        self.assertDictEqual(expected_conseqs, conseqs)

    def testHeaders(self):
        sam_file = StringIO(
            "@SH\tsome header\n"
            "@MHI\tmost headers are ignored, except SQ for sequence reference\n"
            "@SQ\tSN:test\n"
            "test1\t99\ttest\t1\t44\t3M\t=\t1\t3\tACA\tJJJ\n"
        )
        expected_conseqs = {'test': 'ACA'}
        conseqs = remap.sam_to_conseqs(sam_file)
        self.assertDictEqual(expected_conseqs, conseqs)

    def testUnknownReferenceName(self):
        sam_file = StringIO(
            "@SQ\tSN:testX\n"
            "test1\t99\ttestY\t1\t44\t12M\t=\t1\t3\tACA\tJJJ\n"
        )
        expected_conseqs = {}
        conseqs = remap.sam_to_conseqs(sam_file)
        self.assertDictEqual(expected_conseqs, conseqs)

    def testHeaderFields(self):
        sam_file = StringIO(
            "@SQ\tOF:other field: ignored\tSN:test\n"
            "test1\t99\ttest\t1\t44\t3M\t=\t1\t3\tACA\tJJJ\n"
        )
        expected_conseqs = {'test': 'ACA'}
        conseqs = remap.sam_to_conseqs(sam_file)
        self.assertDictEqual(expected_conseqs, conseqs)

    def testExtraFields(self):
        sam_file = StringIO(
            "@SQ\tSN:test\n"
            "test1\t99\ttest\t1\t44\t3M\t=\t1\t3\tACA\tJJJ\tAS:i:236\tNM:i:12\n"
        )
        expected_conseqs = {'test': 'ACA'}
        conseqs = remap.sam_to_conseqs(sam_file)
        self.assertDictEqual(expected_conseqs, conseqs)

    def testMaxConsensus(self):
        sam_file = StringIO(
            "@SQ\tSN:test\n"
            "test1\t99\ttest\t1\t44\t3M\t=\t1\t3\tACA\tJJJ\n"
            "test2\t147\ttest\t1\t44\t3M\t=\t1\t-3\tACA\tJJJ\n"
            "test3\t99\ttest\t1\t44\t3M\t=\t1\t3\tTCA\tJJJ\n"
        )
        expected_conseqs = {'test': 'ACA'}
        conseqs = remap.sam_to_conseqs(sam_file)
        self.assertDictEqual(expected_conseqs, conseqs)

    def testTie(self):
        sam_file = StringIO(
            "@SQ\tSN:test\n"
            "test1\t99\ttest\t1\t44\t3M\t=\t1\t3\tGCA\tJJJ\n"
            "test2\t147\ttest\t1\t44\t3M\t=\t1\t-3\tTCA\tJJJ\n"
        )
        expected_conseqs = {'test': 'GCA'}
        conseqs = remap.sam_to_conseqs(sam_file)
        self.assertDictEqual(expected_conseqs, conseqs)

    def testSoftClip(self):
        sam_file = StringIO(
            "@SQ\tSN:test\n"
            "test1\t99\ttest\t1\t44\t3S5M1S\t=\t1\t9\tACAGGGAGA\tJJJJJJJJJ\n"
        )
        expected_conseqs = {'test': 'GGGAG'}
        conseqs = remap.sam_to_conseqs(sam_file)
        self.assertDictEqual(expected_conseqs, conseqs)

    def testSimpleInsertion(self):
        sam_file = StringIO(
            "@SQ\tSN:test\n"
            "test1\t99\ttest\t1\t44\t3M3I3M\t=\t1\t9\tACAGGGAGA\tJJJJJJJJJ\n"
        )
        expected_conseqs = {'test': 'ACAGGGAGA'}
        conseqs = remap.sam_to_conseqs(sam_file)
        self.assertDictEqual(expected_conseqs, conseqs)

    def testLowQualityInsertion(self):
        sam_file = StringIO(
            "@SQ\tSN:test\n"
            "test1\t99\ttest\t1\t44\t3M3I3M\t=\t1\t9\tACAGGGAGA\tJJJJ/JJJJ\n"
        )
        expected_conseqs = {'test': 'ACAAGA'}
        conseqs = remap.sam_to_conseqs(sam_file, quality_cutoff=32)
        self.assertDictEqual(expected_conseqs, conseqs)

    def testInsertionAfterLowQuality(self):
        sam_file = StringIO(
            "@SQ\tSN:test\n"
            "test1\t99\ttest\t1\t44\t3M3I3M\t=\t1\t9\tACAGGGAGA\tJJ/JJJJJJ\n"
        )
        expected_conseqs = {'test': 'ACNAGA'}
        conseqs = remap.sam_to_conseqs(sam_file, quality_cutoff=32)
        self.assertDictEqual(expected_conseqs, conseqs)

    def testInsertionAndOffset(self):
        sam_file = StringIO(
            "@SQ\tSN:test\n"
            "test1\t99\ttest\t1\t44\t3M3I3M\t=\t1\t9\tACAGGGAGA\tJJJJJJJJJJJJ\n"
            "test2\t99\ttest\t5\t44\t5M\t=\t1\t5\tGACCC\tJJJJJ\n"
        )
        expected_conseqs = {'test': 'ACAGGGAGACCC'}
        conseqs = remap.sam_to_conseqs(sam_file)
        self.assertDictEqual(expected_conseqs, conseqs)

    def testComplexInsertion(self):
        # Insertions are ignored if not a multiple of three
        sam_file = StringIO(
            "@SQ\tSN:test\n"
            "test1\t99\ttest\t1\t44\t3M1I3M2I6M\t=\t1\t12\tACAGAGAGGCCCAAC\tJJJJJJJJJJJJJJJ\n"
        )
        expected_conseqs = {'test': 'ACAAGACCCAAC'}
        conseqs = remap.sam_to_conseqs(sam_file)
        self.assertDictEqual(expected_conseqs, conseqs)

    def testDeletion(self):
        sam_file = StringIO(
            "@SQ\tSN:test\n"
            "test1\t99\ttest\t1\t44\t3M3D3M\t=\t3\t6\tACAGGG\tJJJJJJ\n"
        )
        expected_conseqs = {'test': 'ACAGGG'}
        conseqs = remap.sam_to_conseqs(sam_file)
        self.assertDictEqual(expected_conseqs, conseqs)

    def testDeletionInSomeReads(self):
        sam_file = StringIO(
            "@SQ\tSN:test\n"
            "test1\t99\ttest\t1\t44\t3M3D3M\t=\t3\t6\tACAGGG\tJJJJJJ\n"
            "test2\t99\ttest\t1\t44\t3M3D3M\t=\t3\t6\tACAGGG\tJJJJJJ\n"
            "test3\t99\ttest\t1\t44\t9M\t=\t3\t9\tACATTTGGG\tJJJJJJJJJ\n"
        )
        expected_conseqs = {'test': 'ACATTTGGG'}
        conseqs = remap.sam_to_conseqs(sam_file)
        self.assertDictEqual(expected_conseqs, conseqs)

    def testDeletionWithFrameShift(self):
        sam_file = StringIO(
            "@SQ\tSN:test\n"
            "test1\t99\ttest\t1\t44\t3M1D3M\t=\t3\t6\tACAGGG\tJJJJJJ\n"
        )
        expected_conseqs = {'test': 'ACA-GGG'}
        conseqs = remap.sam_to_conseqs(sam_file)
        self.assertDictEqual(expected_conseqs, conseqs)

    def testBigDeletionWithFrameShift(self):
        sam_file = StringIO(
            "@SQ\tSN:test\n"
            "test1\t99\ttest\t1\t44\t3M4D3M\t=\t3\t6\tACAGGG\tJJJJJJ\n"
        )
        expected_conseqs = {'test': 'ACA----GGG'}
        conseqs = remap.sam_to_conseqs(sam_file)
        self.assertDictEqual(expected_conseqs, conseqs)

    def testOverlapsCountOnce(self):
        sam_file = StringIO(
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
        conseqs = remap.sam_to_conseqs(sam_file)
        self.assertDictEqual(expected_conseqs, conseqs)

    def testReverseLeftOfForward(self):
        sam_file = StringIO(
            "@SQ\tSN:test\n"
            "test1\t99\ttest\t2\t44\t1M\t=\t1\t1\tC\tJ\n"
            "test1\t147\ttest\t1\t44\t1M\t=\t2\t-1\tA\tJ\n"
        )
        expected_conseqs = {'test': 'AC'}
        conseqs = remap.sam_to_conseqs(sam_file)
        self.assertDictEqual(expected_conseqs, conseqs)

    def testPairMapsToTwoReferences(self):
        sam_file = StringIO(
            "@SQ\tSN:testX\n"
            "@SQ\tSN:testY\n"
            "test1\t99\ttestX\t1\t44\t3M\t=\t1\t3\tACG\tJJJ\n"
            "test1\t147\ttestY\t1\t44\t3M\t=\t1\t-3\tACG\tJJJ\n"
        )
        expected_conseqs = {}
        conseqs = remap.sam_to_conseqs(sam_file)
        self.assertDictEqual(expected_conseqs, conseqs)

    def testLowQuality(self):
        # Note that we ignore the overlapped portion of the reverse read,
        # even if it has higher quality.
        sam_file = StringIO(
            "@SQ\tSN:test\n"
            "test1\t99\ttest\t1\t44\t3M\t=\t1\t3\tACG\tJ/J\n"
        )
        expected_conseqs = {'test': 'ANG'}
        conseqs = remap.sam_to_conseqs(sam_file, quality_cutoff=32)
        self.assertDictEqual(expected_conseqs, conseqs)

    def testLowQualityAtEnd(self):
        sam_file = StringIO(
            "@SQ\tSN:test\n"
            "test1\t99\ttest\t1\t44\t3M\t=\t1\t3\tACG\tJJ/\n"
        )
        expected_conseqs = {'test': 'ACN'}
        conseqs = remap.sam_to_conseqs(sam_file, quality_cutoff=32)
        self.assertDictEqual(expected_conseqs, conseqs)

    def testLowQualityForward(self):
        sam_file = StringIO(
            "@SQ\tSN:test\n"
            "test1\t99\ttest\t1\t44\t3M\t=\t3\t3\tATA\tJJA\n"
            "test1\t147\ttest\t3\t44\t3M\t=\t1\t-3\tGCC\tJJJ\n"
        )
        expected_conseqs = {'test': 'ATGCC'}
        conseqs = remap.sam_to_conseqs(sam_file)
        self.assertDictEqual(expected_conseqs, conseqs)

    def testAllLowQuality(self):
        # SAM:qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual
        sam_file = StringIO(
            "@SQ\tSN:test\n"
            "test1\t147\ttest\t1\t24\t1M\t=\t1\t-1\tT\t#\n"
        )
        expected_conseqs = {}

        conseqs = remap.sam_to_conseqs(sam_file, quality_cutoff=32)

        self.assertDictEqual(expected_conseqs, conseqs)

    def testBadPairFlag(self):
        """ Even if the pair isn't concordant, still include in consensus.

        SAM flag 145 does not have bit 2 for properly aligned.
        """
        # SAM:qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual
        sam_file = StringIO(
            "@SQ\tSN:test\n"
            "test1\t145\ttest\t1\t24\t1M\t=\t1\t-1\tT\tF\n"
        )
        expected_conseqs = {'test': 'T'}

        conseqs = remap.sam_to_conseqs(sam_file, quality_cutoff=32)

        self.assertEqual(expected_conseqs, conseqs)

    def testUnmappedFlag(self):
        """ If the read is unmapped, don't include in consensus.

        SAM flag 149 has bit 4 for unmapped. Region is irrelevant.
        """
        # SAM:qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual
        sam_file = StringIO(
            "@SQ\tSN:test\n"
            "test1\t149\ttest\t1\t24\t1M\t=\t1\t-1\tT\tF\n"
        )
        expected_conseqs = {}

        conseqs = remap.sam_to_conseqs(sam_file, quality_cutoff=32)

        self.assertEqual(expected_conseqs, conseqs)

    def testLowQualityAndDeletion(self):
        sam_file = StringIO(
            "@SQ\tSN:test\n"
            "test1\t99\ttest\t1\t44\t3M3D3M\t=\t3\t6\tACAGGG\tJJJJJJ\n"
            "test2\t99\ttest\t1\t44\t9M\t=\t3\t9\tACATTTGGG\tJJJ///JJJ\n"
        )
        expected_conseqs = {'test': 'ACANNNGGG'}

        conseqs = remap.sam_to_conseqs(sam_file, quality_cutoff=32)

        self.assertDictEqual(expected_conseqs, conseqs)

    def testSeeds(self):
        sam_file = StringIO(
            "@SQ\tSN:test\n"
            "test1\t99\ttest\t4\t44\t3M\t=\t10\t3\tTAT\tJJJ\n"
            "test2\t99\ttest\t10\t44\t3M\t=\t4\t-3\tCAC\tJJJ\n"
        )
        seeds = {'test': 'ACATTTGGGCAC'}
        expected_conseqs = {'test': 'ACATATGGGCAC'}

        conseqs = remap.sam_to_conseqs(sam_file, seeds=seeds)

        self.assertDictEqual(expected_conseqs, conseqs)

    def testSeedsNeedSomeReads(self):
        sam_file = StringIO(
            "@SQ\tSN:test\n"
            "test1\t99\ttest\t4\t44\t3M\t=\t10\t3\tTAT\tJJJ\n"
        )
        seeds = {'test': 'ACATTTGGGCAC',
                 'other': 'TATGCACCC'}
        expected_conseqs = {'test': 'ACATATGGGCAC'}

        conseqs = remap.sam_to_conseqs(sam_file, seeds=seeds)

        self.assertDictEqual(expected_conseqs, conseqs)

    def testSeedsWithLowQuality(self):
        sam_file = StringIO(
            "@SQ\tSN:test\n"
            "test1\t99\ttest\t4\t44\t3M\t=\t10\t3\tTAT\tJJ/\n"
        )
        seeds = {'test': 'ACATTTGGGCAC'}
        expected_conseqs = {'test': 'ACATATGGGCAC'}

        conseqs = remap.sam_to_conseqs(sam_file,
                                       seeds=seeds,
                                       quality_cutoff=32)

        self.assertDictEqual(expected_conseqs, conseqs)

    def testSeedsWithPartialDeletion(self):
        sam_file = StringIO(
            "@SQ\tSN:test\n"
            "test1\t99\ttest\t4\t44\t1M1D1M\t=\t10\t3\tTT\tJJ\n"
        )
        seeds = {'test': 'ACATTTGGGCAC'}
        expected_conseqs = {'test': 'ACATTTGGGCAC'}

        conseqs = remap.sam_to_conseqs(sam_file,
                                       seeds=seeds,
                                       quality_cutoff=32)

        self.assertDictEqual(expected_conseqs, conseqs)

    def testSeedsWithCodonDeletion(self):
        sam_file = StringIO(
            "@SQ\tSN:test\n"
            "test1\t99\ttest\t1\t44\t3M3D3M\t=\t10\t6\tACAGGG\tJJJJJJ\n"
        )
        seeds = {'test': 'ACATTTGGGCAC'}
        expected_conseqs = {'test': 'ACATTTGGGCAC'}

        conseqs = remap.sam_to_conseqs(sam_file,
                                       seeds=seeds,
                                       quality_cutoff=32)

        self.assertDictEqual(expected_conseqs, conseqs)

    def testDebugReports(self):
        sam_file = StringIO(
            "@SQ\tSN:test\n"
            "test1\t99\ttest\t1\t44\t3M3I9M\t=\t1\t12\tACTGGGAGACCCAAC\tJIJJJJJJJJJJJJJ\n"
            "test1\t147\ttest\t1\t44\t3M3I9M\t=\t1\t-12\tACTGGGAGACCCAAC\tJIJJJJJJJJJJJJJ\n"
            "test1\t99\ttest\t1\t44\t3M3I9M\t=\t1\t12\tATTGGGAGACCCAAC\tJHJJJJJJJJJJJJJ\n"
            "test1\t147\ttest\t1\t44\t3M3I9M\t=\t1\t-12\tATTGGGAGACCCAAC\tJHJJJJJJJJJJJJJ\n"
        )
        reports = {('test', 2): None}
        expected_reports = {('test', 2): 'H{C: 1, T: 1}, I{C: 1}'}

        remap.sam_to_conseqs(sam_file, debug_reports=reports)

        self.assertDictEqual(expected_reports, reports)

    def testDebugReportsOnReverseRead(self):
        sam_file = StringIO(
            "@SQ\tSN:test\n"
            "test1\t99\ttest\t1\t44\t3M3I2M\t=\t1\t8\tACTGGGAG\tJJJJJJJJ\n"
            "test1\t147\ttest\t5\t44\t8M\t=\t1\t-8\tGACCCAAC\tJJJJJIJJ\n"
            "test1\t99\ttest\t1\t44\t3M3I2M\t=\t1\t12\tATTGGGAG\tJJJJJJJJ\n"
            "test1\t147\ttest\t5\t44\t8M\t=\t1\t-12\tGACCCAAC\tJJJJJHJJ\n"
        )
        reports = {('test', 10): None}
        expected_reports = {('test', 10): 'H{A: 2}, I{A: 1}'}

        remap.sam_to_conseqs(sam_file, debug_reports=reports)

        self.assertDictEqual(expected_reports, reports)


def test_drop_drifters_seeds_converged():
    relevant_conseqs = dict(test='ATGAGGAGTA',
                            other='ATGACCAGTA',
                            wayoff='ATGAGGGTAC')
    original_seeds = dict(test='ATGAAGTA',
                          other='AAGCCGAA',
                          wayoff='TCATGTAC')
    read_counts = dict(test=1, other=1, wayoff=1)
    distance_report = {}

    expected_seed_names = {'test'}
    expected_distances = dict(test=dict(seed_dist=2,
                                        other_dist=5,
                                        other_seed='other'),
                              other=dict(seed_dist=4,
                                         other_dist=2,
                                         other_seed='test'),
                              wayoff=dict(seed_dist=4,
                                          other_dist=3,
                                          other_seed='test'))

    remap.drop_drifters(relevant_conseqs,
                        original_seeds,
                        distance_report,
                        read_counts)

    assert distance_report == expected_distances
    assert set(relevant_conseqs) == expected_seed_names


def test_drop_drifters_seeds_converged_with_different_alignment():
    """ Seeds have similar regions, but at different positions.
    """
    relevant_conseqs = dict(test='ATGAGGAGTA',
                            other='ATGACCAGTA')
    original_seeds = dict(test='ATGAAGTA',
                          other='TCTCTCTCTCAAGCCGAA')
    read_counts = dict(test=1, other=1)
    distance_report = {}
    expected_seed_names = {'test'}

    remap.drop_drifters(relevant_conseqs,
                        original_seeds,
                        distance_report,
                        read_counts)

    assert set(relevant_conseqs) == expected_seed_names


def test_sam_to_conseqs_seeds_converged_with_different_alignment_and_gap():
    """ Gaps between areas with coverage.
    """
    sam_file = StringIO(
        f"@SQ\tSN:test\tSN:other\n"
        f"test1\t99\ttest\t1\t44\t10M\t=\t1\t10\tATGAGGAGTA\tJJJJJJJJJJJJ\n"
        f"other1\t99\tother\t11\t44\t5M\t=\t1\t5\tATGAC\tJJJJJJJ\n"
        f"other2\t99\tother\t26\t44\t5M\t=\t1\t5\tCAGTA\tJJJJJJJ\n"
    )
    seeds = {'test': 'ATGAAGTA',
             'other': 'TCTCTCTCTCAAGCTATATATATACGAA'}
    expected_conseqs = {'test': 'ATGAGGAGTA'}

    conseqs = remap.sam_to_conseqs(sam_file,
                                   seeds=seeds,
                                   original_seeds=seeds,
                                   is_filtered=True)

    assert conseqs == expected_conseqs


def test_drop_drifters_seeds_converged_with_confusing_gap():
    """ Reads match other seed well, but with a big gap.
    """
    #                             vvvvvvvv
    relevant_conseqs = dict(test='ATGTCGTA',
                            #      |||||||||
                            other='AAGCTATAT')
    original_seeds = dict(
        #     vvv??vvv
        test='ATGAAGTA',
        #      vvvvv        |||||||||         vvv
        other='ATGTCTCTCTCTCAAGCTATATATATACGAAGTA')
    read_counts = dict(test=1, other=1)
    distance_report = {}
    expected_seed_names = {'test', 'other'}

    remap.drop_drifters(relevant_conseqs,
                        original_seeds,
                        distance_report,
                        read_counts)

    assert set(relevant_conseqs) == expected_seed_names


def test_drop_drifters_seeds_converged_plus_other_low_coverage():
    """ Portion with decent coverage has converged, other hasn't.
    """
    relevant_conseqs = dict(test='ATGAGGAGTA', other='ATGACCAGTA')
    original_seeds = dict(test='ATGAAGTACTCTCT', other='AAGCCGAAGTGTGT')
    read_counts = dict(test=2, other=3)
    distance_report = {}

    expected_seed_names = {'test'}

    remap.drop_drifters(relevant_conseqs,
                        original_seeds,
                        distance_report,
                        read_counts)

    assert set(relevant_conseqs) == expected_seed_names


def test_drop_drifters_seeds_both_converged(projects):
    """ Both references are now closer to the other seed than the start.

    Don't drop both. Keep test because it has more reads.
    """
    hxb2_end = projects.getReference(HXB2_NAME)[-200:]
    relevant_conseqs = dict(test='AAGCCGTA' + hxb2_end,
                            #      ^ ^^
                            other='ATGAAGTA' + hxb2_end,
                            #       ^ ^^ ^
                            unrelated='GGGTTTGGG' + hxb2_end)
    original_seeds = dict(test='ATGAAGTA' + hxb2_end,
                          other='AAGCCGAA' + hxb2_end,
                          unrelated='GGGTTTGGG' + hxb2_end)
    read_counts = dict(test=2, other=1, unrelated=1)
    distance_report = {}

    expected_seed_names = {'test', 'unrelated'}
    expected_distances = dict(test=dict(seed_dist=3,
                                        other_dist=1,
                                        other_seed='other'),
                              other=dict(seed_dist=4,
                                         other_dist=0,
                                         other_seed='test'),
                              unrelated=dict(seed_dist=0,
                                             other_dist=7,
                                             other_seed='other'))

    remap.drop_drifters(relevant_conseqs,
                        original_seeds,
                        distance_report,
                        read_counts)

    assert distance_report == expected_distances
    assert set(relevant_conseqs) == expected_seed_names


def test_sam_to_conseqs_all_seeds_low_coverage():
    """ Multiple seeds mapped, but none have good coverage.

    Choose most reads.
    """

    sam_file = StringIO(
        "@SQ\tSN:test\tSN:other\n"
        "test1\t99\ttest\t1\t44\t10M\t=\t1\t10\tATGAGGAGTA\tJJJJJJJJJJJJ\n"
        "test2\t99\ttest\t11\t44\t6M\t=\t1\t10\tCTCTCT\tJJJJJJ\n"
        "other1\t99\tother\t1\t44\t10M\t=\t1\t10\tATGACCAGTA\tJJJJJJJJJJJJ\n"
    )
    seeds = {'test': 'ATGAAGTACTCTCT',
             'other': 'AAGCCGAAGTGTGT'}
    expected_conseqs = {'test': 'ATGAGGAGTACTCTCT'}

    conseqs = remap.sam_to_conseqs(sam_file,
                                   seeds=seeds,
                                   original_seeds=seeds,
                                   is_filtered=True,
                                   filter_coverage=2)

    assert conseqs == expected_conseqs


def test_extract_relevant_seeds():
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
        assert (aligned_conseq,
                aligned_seed,
                relevant) == (aligned_conseq, aligned_seed, expected_seed)


def test_sam_to_conseqs_nothing_mapped():
    sam_file = StringIO(
        "@SQ\tSN:test\tSN:other\n"
    )
    seeds = {'test': 'ATGAAGTACTCTCT',
             'other': 'AAGCCGAAGTGTGT'}
    expected_conseqs = {}

    conseqs = remap.sam_to_conseqs(sam_file,
                                   seeds=seeds,
                                   original_seeds=seeds,
                                   is_filtered=True,
                                   filter_coverage=2)

    assert conseqs == expected_conseqs


class MixedReferenceMemorySplitter(MixedReferenceSplitter):
    """ Dummy class to hold split reads in memory. Useful for testing. """
    def __init__(self, work_path):
        super().__init__(work_path)
        self.is_closed = True

    def create_split_file(self, refname, direction):
        self.is_closed = False
        return StringIO()

    def close_split_file(self, split_file):
        self.is_closed = True


# noinspection DuplicatedCode
class MixedReferenceSplitterTest(unittest.TestCase):
    def setUp(self):
        super(MixedReferenceSplitterTest, self).setUp()
        self.addTypeEqualityFunc(str, self.assertMultiLineEqual)
        self.work_path = os.path.dirname(__file__)

    def testSimple(self):
        sam_file = StringIO(
            "@SQ\tSN:r\n"
            "r1\t99\tr\t1\t44\t3M\t=\t1\t3\tACG\tJJJ\n"
            "r1\t147\tr\t1\t44\t3M\t=\t1\t-3\tACG\tJJJ\n")
        expected_rows = [
            ["r1", "99", "r", "1", "44", "3M", "=", "1", "3", "ACG", "JJJ"],
            ["r1", "147", "r", "1", "44", "3M", "=", "1", "-3", "ACG", "JJJ"]]

        splitter = MixedReferenceSplitter(self.work_path)
        rows = list(splitter.split(sam_file))

        self.assertEqual(expected_rows, rows)

    def testTrimOptionalFields(self):
        sam_file = StringIO(
            "@SQ\tSN:r\n"
            "r1\t99\tr\t1\t44\t3M\t=\t1\t3\tACG\tJJJ\tAS:i:100\n"
            "r1\t147\tr\t1\t44\t3M\t=\t1\t-3\tACG\tJJJ\tYS:Z:UP\n")
        expected_rows = [
            ["r1", "99", "r", "1", "44", "3M", "=", "1", "3", "ACG", "JJJ"],
            ["r1", "147", "r", "1", "44", "3M", "=", "1", "-3", "ACG", "JJJ"]]

        splitter = MixedReferenceSplitter(self.work_path)
        rows = list(splitter.split(sam_file))

        self.assertEqual(expected_rows, rows)

    def testUnmapped(self):
        sam_file = StringIO(
            "@SQ\tSN:r\n"
            "r1\t107\tr\t1\t44\t3M\t*\t1\t3\tACG\tJJJ\n"
            "r1\t149\t*\t*\t*\t*\tr\t*\t*\tACG\tJJJ\n")
        expected_rows = []

        splitter = MixedReferenceSplitter(self.work_path)
        rows = list(splitter.split(sam_file))

        self.assertEqual(expected_rows, rows)

    def testSplit(self):
        """ If a pair is split over two references, choose one.

        Use the reference that gave the higher mapq score.
        """
        sam_file = StringIO(
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

        splitter = MixedReferenceMemorySplitter(self.work_path)
        rows = list(splitter.split(sam_file))
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
        sam_file = StringIO(
            "@SQ\tSN:r\n"
            "r1\t99\tRX\t1\t44\t3M\t=\t1\t3\tACG\tJJJ\tAS:i:100\n"
            "r1\t147\tRX\t1\t44\t3M\t=\t1\t-3\tACG\tJJJ\tYS:Z:UP\n"
            "r2\t99\tRX\t1\t11\t3M\tRY\t1\t3\tACG\tJJJ\tAS:i:100\n"
            "r2\t147\tRY\t1\t11\t3M\tRX\t1\t-3\tACC\tJJK\tAS:i:200\n")

        splitter = MixedReferenceMemorySplitter(self.work_path)
        list(splitter.split(sam_file))
        splits = splitter.splits

        self.assertEqual(['RY'], list(splits.keys()))

    def testWalk(self):
        sam_file = StringIO(
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

        splitter = MixedReferenceMemorySplitter(self.work_path)
        rows = list(splitter.walk(sam_file))
        splits = splitter.splits

        self.assertEqual(expected_rows, rows)
        self.assertEqual({}, splits)


class RemapCountsTest(unittest.TestCase):
    def test_simple(self):
        report = StringIO()
        writer = DictWriter(
            report,
            ['type', 'count', 'seed_dist', 'other_dist', 'other_seed'],
            lineterminator=os.linesep)
        counts = {'r1': 100, 'r2': 200}
        expected_report = """\
prelim r1,100,,,
prelim r2,200,,,
"""

        write_remap_counts(writer, counts, 'prelim')

        self.assertEqual(expected_report, report.getvalue())

    def test_distance(self):
        report = StringIO()
        writer = DictWriter(
            report,
            ['type', 'count', 'seed_dist', 'other_dist', 'other_seed'],
            lineterminator=os.linesep)
        counts = {'r1': 100, 'r2': 200}
        distance_report = {'r1': {'seed_dist': 1,
                                  'other_dist': 10,
                                  'other_seed': 'r2'},
                           'r2': {'seed_dist': 2,
                                  'other_dist': 20,
                                  'other_seed': 'r1'}}
        expected_report = """\
remap r1,100,1,10,r2
remap r2,200,2,20,r1
"""

        write_remap_counts(writer, counts, 'remap', distance_report)

        self.assertEqual(expected_report, report.getvalue())


# noinspection DuplicatedCode
class ConvertPrelimTest(unittest.TestCase):
    def setUp(self):
        self.projects = ProjectConfig()
        self.projects.load(StringIO("""\
            {
              "regions": {
                "R1-seed": {
                  "seed_group": "main",
                  "reference": ["ACTAAAGGG"]
                },
                "R2-seed": {
                  "seed_group": "main",
                  "reference": ["ACTAAAGGGAAA"]
                }
              }
            }
            """))
        self.sam_file = StringIO()
        self.remap_counts = StringIO()
        self.remap_counts_writer = DictWriter(
            self.remap_counts,
            ['type', 'filtered_count', 'count'],
            lineterminator=os.linesep)
        self.remap_counts_writer.writeheader()

    def test_simple(self):
        prelim_csv = StringIO("""\
qname,flag,rname,pos,mapq,cigar,rnext,pnext,tlen,seq,qual
example1,89,R1-seed,1,0,9M,=,1,0,AAACCCTTT,BBBBBBBBB
""")
        count_threshold = 2
        expected_sam_file = """\
@HD	VN:1.0	SO:unsorted
@SQ	SN:R1-seed	LN:9
@SQ	SN:R2-seed	LN:12
@PG	ID:bowtie2	PN:bowtie2	VN:2.2.3	CL:""
example1\t89\tR1-seed\t1\t0\t9M\t=\t1\t0\tAAACCCTTT\tBBBBBBBBB
"""
        expected_remap_counts = """\
type,filtered_count,count
prelim R1-seed,0,1
"""
        expected_seed_counts = {}

        seed_counts = convert_prelim(prelim_csv,
                                     self.sam_file,
                                     self.remap_counts_writer,
                                     count_threshold,
                                     self.projects)

        self.assertEqual(expected_sam_file, self.sam_file.getvalue())
        self.assertEqual(expected_remap_counts, self.remap_counts.getvalue())
        self.assertEqual(expected_seed_counts, seed_counts)

    def test_two_regions(self):
        prelim_csv = StringIO("""\
qname,flag,rname,pos,mapq,cigar,rnext,pnext,tlen,seq,qual
example1,89,R1-seed,1,0,9M,=,1,0,AAACCCTTT,BBBBBBBBB
example2,89,R2-seed,1,0,9M,=,1,0,AAAACCTTT,BBBBBBBBB
example3,89,R2-seed,1,0,9M,=,1,0,AAAAACTTT,BBBBBBBBB
""")
        count_threshold = 2
        expected_sam_file = """\
@HD	VN:1.0	SO:unsorted
@SQ	SN:R1-seed	LN:9
@SQ	SN:R2-seed	LN:12
@PG	ID:bowtie2	PN:bowtie2	VN:2.2.3	CL:""
example1\t89\tR1-seed\t1\t0\t9M\t=\t1\t0\tAAACCCTTT\tBBBBBBBBB
example2\t89\tR2-seed\t1\t0\t9M\t=\t1\t0\tAAAACCTTT\tBBBBBBBBB
example3\t89\tR2-seed\t1\t0\t9M\t=\t1\t0\tAAAAACTTT\tBBBBBBBBB
"""
        expected_remap_counts = """\
type,filtered_count,count
prelim R1-seed,0,1
prelim R2-seed,0,2
"""
        expected_seed_counts = {}

        seed_counts = convert_prelim(prelim_csv,
                                     self.sam_file,
                                     self.remap_counts_writer,
                                     count_threshold,
                                     self.projects)

        self.assertEqual(expected_sam_file, self.sam_file.getvalue())
        self.assertEqual(expected_remap_counts, self.remap_counts.getvalue())
        self.assertEqual(expected_seed_counts, seed_counts)

    def test_long_reads(self):
        self.maxDiff = None
        prelim_csv = StringIO("""\
qname,flag,rname,pos,mapq,cigar,rnext,pnext,tlen,seq,qual
example1,89,R1-seed,1,0,54M,=,1,0,\
AAACCCTTTAAACCCTTTAAACCCTTTAAACCCTTTAAACCCTTTAAACCCTTT,\
BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
example2,89,R1-seed,1,0,54M,=,1,0,\
AAAACCTTTAAAACCTTTAAAACCTTTAAAACCTTTAAAACCTTTAAAACCTTT,\
BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
""")
        count_threshold = 2
        expected_sam_file = """\
@HD	VN:1.0	SO:unsorted
@SQ	SN:R1-seed	LN:9
@SQ	SN:R2-seed	LN:12
@PG	ID:bowtie2	PN:bowtie2	VN:2.2.3	CL:""
example1\t89\tR1-seed\t1\t0\t54M\t=\t1\t0\t\
AAACCCTTTAAACCCTTTAAACCCTTTAAACCCTTTAAACCCTTTAAACCCTTT\t\
BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
example2\t89\tR1-seed\t1\t0\t54M\t=\t1\t0\t\
AAAACCTTTAAAACCTTTAAAACCTTTAAAACCTTTAAAACCTTTAAAACCTTT\t\
BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
"""
        expected_remap_counts = """\
type,filtered_count,count
prelim R1-seed,2,2
"""
        expected_seed_counts = {'R1-seed': 2}

        seed_counts = convert_prelim(prelim_csv,
                                     self.sam_file,
                                     self.remap_counts_writer,
                                     count_threshold,
                                     self.projects)

        self.assertEqual(expected_sam_file, self.sam_file.getvalue())
        self.assertEqual(expected_remap_counts, self.remap_counts.getvalue())
        self.assertEqual(expected_seed_counts, seed_counts)

    def test_star_region(self):
        self.maxDiff = None
        prelim_csv = StringIO("""\
qname,flag,rname,pos,mapq,cigar,rnext,pnext,tlen,seq,qual
example1,89,R1-seed,1,0,54M,=,1,0,\
AAACCCTTTAAACCCTTTAAACCCTTTAAACCCTTTAAACCCTTTAAACCCTTT,\
BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
example2,89,R1-seed,1,0,54M,=,1,0,\
AAAACCTTTAAAACCTTTAAAACCTTTAAAACCTTTAAAACCTTTAAAACCTTT,\
BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
example3,93,*,*,*,*,*,*,*,\
AAAACCTTTAAAACCTTTAAAACCTTTAAAACCTTTAAAACCTTTAAAACCTTT,\
BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
""")
        count_threshold = 2
        expected_sam_file = """\
@HD	VN:1.0	SO:unsorted
@SQ	SN:R1-seed	LN:9
@SQ	SN:R2-seed	LN:12
@PG	ID:bowtie2	PN:bowtie2	VN:2.2.3	CL:""
example1\t89\tR1-seed\t1\t0\t54M\t=\t1\t0\t\
AAACCCTTTAAACCCTTTAAACCCTTTAAACCCTTTAAACCCTTTAAACCCTTT\t\
BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
example2\t89\tR1-seed\t1\t0\t54M\t=\t1\t0\t\
AAAACCTTTAAAACCTTTAAAACCTTTAAAACCTTTAAAACCTTTAAAACCTTT\t\
BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
example3\t93\t*\t*\t*\t*\t*\t*\t*\t\
AAAACCTTTAAAACCTTTAAAACCTTTAAAACCTTTAAAACCTTTAAAACCTTT\t\
BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
"""
        expected_remap_counts = """\
type,filtered_count,count
prelim *,0,1
prelim R1-seed,2,2
"""
        expected_seed_counts = {'R1-seed': 2}

        seed_counts = convert_prelim(prelim_csv,
                                     self.sam_file,
                                     self.remap_counts_writer,
                                     count_threshold,
                                     self.projects)

        self.assertEqual(expected_sam_file, self.sam_file.getvalue())
        self.assertEqual(expected_remap_counts, self.remap_counts.getvalue())
        self.assertEqual(expected_seed_counts, seed_counts)

    def test_best_in_group(self):
        self.maxDiff = None
        prelim_csv = StringIO("""\
qname,flag,rname,pos,mapq,cigar,rnext,pnext,tlen,seq,qual
example1,89,R1-seed,1,0,54M,=,1,0,\
AAACCCTTTAAACCCTTTAAACCCTTTAAACCCTTTAAACCCTTTAAACCCTTT,\
BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
example2,89,R2-seed,1,0,54M,=,1,0,\
AAAACCTTTAAAACCTTTAAAACCTTTAAAACCTTTAAAACCTTTAAAACCTTT,\
BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
example3,89,R1-seed,1,0,54M,=,1,0,\
AAAAAATTTAAAACCTTTAAAACCTTTAAAACCTTTAAAACCTTTAAAACCTTT,\
BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
example4,89,R2-seed,1,0,54M,=,1,0,\
AAAAAAAATAAAACCTTTAAAACCTTTAAAACCTTTAAAACCTTTAAAACCTTT,\
BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
example5,89,R2-seed,1,0,54M,=,1,0,\
AAAAAAAAAAAAACCTTTAAAACCTTTAAAACCTTTAAAACCTTTAAAACCTTT,\
BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
""")
        count_threshold = 2
        expected_sam_file = """\
@HD	VN:1.0	SO:unsorted
@SQ	SN:R1-seed	LN:9
@SQ	SN:R2-seed	LN:12
@PG	ID:bowtie2	PN:bowtie2	VN:2.2.3	CL:""
example1\t89\tR1-seed\t1\t0\t54M\t=\t1\t0\t\
AAACCCTTTAAACCCTTTAAACCCTTTAAACCCTTTAAACCCTTTAAACCCTTT\t\
BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
example2\t89\tR2-seed\t1\t0\t54M\t=\t1\t0\t\
AAAACCTTTAAAACCTTTAAAACCTTTAAAACCTTTAAAACCTTTAAAACCTTT\t\
BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
example3\t89\tR1-seed\t1\t0\t54M\t=\t1\t0\t\
AAAAAATTTAAAACCTTTAAAACCTTTAAAACCTTTAAAACCTTTAAAACCTTT\t\
BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
example4\t89\tR2-seed\t1\t0\t54M\t=\t1\t0\t\
AAAAAAAATAAAACCTTTAAAACCTTTAAAACCTTTAAAACCTTTAAAACCTTT\t\
BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
example5\t89\tR2-seed\t1\t0\t54M\t=\t1\t0\t\
AAAAAAAAAAAAACCTTTAAAACCTTTAAAACCTTTAAAACCTTTAAAACCTTT\t\
BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
"""
        expected_remap_counts = """\
type,filtered_count,count
prelim R1-seed,2,2
prelim R2-seed,3,3
"""
        expected_seed_counts = {'R2-seed': 3}

        seed_counts = convert_prelim(prelim_csv,
                                     self.sam_file,
                                     self.remap_counts_writer,
                                     count_threshold,
                                     self.projects)

        self.assertEqual(expected_sam_file, self.sam_file.getvalue())
        self.assertEqual(expected_remap_counts, self.remap_counts.getvalue())
        self.assertEqual(expected_seed_counts, seed_counts)

    def test_unmapped_read(self):
        self.maxDiff = None
        prelim_csv = StringIO("""\
qname,flag,rname,pos,mapq,cigar,rnext,pnext,tlen,seq,qual
example1,89,R1-seed,1,0,54M,=,1,0,\
AAACCCTTTAAACCCTTTAAACCCTTTAAACCCTTTAAACCCTTTAAACCCTTT,\
BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
example2,93,R1-seed,1,0,54M,=,1,0,\
AAAACCTTTAAAACCTTTAAAACCTTTAAAACCTTTAAAACCTTTAAAACCTTT,\
BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
""")
        count_threshold = 2
        expected_sam_file = """\
@HD	VN:1.0	SO:unsorted
@SQ	SN:R1-seed	LN:9
@SQ	SN:R2-seed	LN:12
@PG	ID:bowtie2	PN:bowtie2	VN:2.2.3	CL:""
example1\t89\tR1-seed\t1\t0\t54M\t=\t1\t0\t\
AAACCCTTTAAACCCTTTAAACCCTTTAAACCCTTTAAACCCTTTAAACCCTTT\t\
BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
example2\t93\tR1-seed\t1\t0\t54M\t=\t1\t0\t\
AAAACCTTTAAAACCTTTAAAACCTTTAAAACCTTTAAAACCTTTAAAACCTTT\t\
BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
"""
        expected_remap_counts = """\
type,filtered_count,count
prelim *,0,1
prelim R1-seed,1,1
"""
        expected_seed_counts = {}

        seed_counts = convert_prelim(prelim_csv,
                                     self.sam_file,
                                     self.remap_counts_writer,
                                     count_threshold,
                                     self.projects)

        self.assertEqual(expected_sam_file, self.sam_file.getvalue())
        self.assertEqual(expected_remap_counts, self.remap_counts.getvalue())
        self.assertEqual(expected_seed_counts, seed_counts)


class RemapTest(unittest.TestCase):
    def setUp(self):
        patcher = patch.multiple(Bowtie2, __init__=Mock(return_value=None), yield_output=DEFAULT)
        self.bowtie2_output = []
        mocks = patcher.start()
        mocks['yield_output'].return_value = self.bowtie2_output

        self.addCleanup(patcher.stop)
        patcher = patch.multiple(Bowtie2Build,
                                 __init__=Mock(return_value=None),
                                 build=DEFAULT)
        patcher.start()
        self.addCleanup(patcher.stop)

        mock_refs = {'R1': 'GTGGG',
                     'R2': 'ACAAA'}
        patcher = patch.object(ProjectConfig, 'loadDefault')
        mock_projects = patcher.start()
        self.addCleanup(patcher.stop)
        mock_projects.return_value.getAllReferences.return_value = mock_refs
        mock_projects.return_value.getReference.side_effect = mock_refs.__getitem__
        patcher = patch('micall.core.remap.is_short_read', Mock(return_value=False))
        patcher.start()
        self.addCleanup(patcher.stop)

    def test_good_contig(self):
        contigs_csv = StringIO("""\
ref,match,group_ref,contig
R1,1.0,R1,TCACCAGGACAGCGGGTTGAATTCCTCGTGCAAGCGTGGAA
""")
        self.bowtie2_output.extend([
            "read1\t99\t1-R1\t1\t44\t5M\t=\t1\t-81\tGTGGG\tAAAAA\n",
            "read1\t147\t1-R1\t1\t44\t5M\t=\t1\t-81\tGTGGG\tAAAAA\n",
            "read2\t99\t1-R1\t1\t44\t5M\t=\t1\t-81\tGTGGG\tAAAAA\n",
            "read2\t147\t1-R1\t1\t44\t5M\t=\t1\t-81\tGTGGG\tAAAAA\n",
            "read3\t99\t1-R1\t1\t44\t5M\t=\t1\t-81\tGTGGG\tAAAAA\n",
            "read3\t147\t1-R1\t1\t44\t5M\t=\t1\t-81\tGTGGG\tAAAAA\n",
            "read4\t99\t1-R1\t1\t44\t5M\t=\t1\t-81\tGTGGG\tAAAAA\n",
            "read4\t147\t1-R1\t1\t44\t5M\t=\t1\t-81\tGTGGG\tAAAAA\n",
            "read5\t77\t*\t0\t0\t*\t*\t0\t0\tGTAAA\tAAAAA\n",
            "read5\t141\t*\t0\t0\t*\t*\t0\t0\tGTAAA\tAAAAA\n"])
        expected_remap_counts_csv = """\
type,count,filtered_count,seed_dist,other_dist,other_seed
raw,20,,,,
remap 1-R1,8,,,,
remap-final 1-R1,8,,,,
unmapped,2,,,,
"""
        expected_remap_csv = """\
qname,flag,rname,pos,mapq,cigar,rnext,pnext,tlen,seq,qual
read1,99,1-R1,1,44,5M,=,1,-81,GTGGG,AAAAA
read1,147,1-R1,1,44,5M,=,1,-81,GTGGG,AAAAA
read2,99,1-R1,1,44,5M,=,1,-81,GTGGG,AAAAA
read2,147,1-R1,1,44,5M,=,1,-81,GTGGG,AAAAA
read3,99,1-R1,1,44,5M,=,1,-81,GTGGG,AAAAA
read3,147,1-R1,1,44,5M,=,1,-81,GTGGG,AAAAA
read4,99,1-R1,1,44,5M,=,1,-81,GTGGG,AAAAA
read4,147,1-R1,1,44,5M,=,1,-81,GTGGG,AAAAA
"""
        self.assertMapsToContigs(contigs_csv, expected_remap_csv, expected_remap_counts_csv)

    def assertMapsToContigs(self, contigs_csv, expected_remap_csv, expected_remap_counts_csv):
        test_path = os.path.dirname(__file__)
        remap_counts_csv = StringIO()
        remap_csv = StringIO()
        remap.map_to_contigs(
            os.path.join(test_path,
                         'microtest',
                         '1234A-V3LOOP_S1_L001_R1_001.fastq'),
            os.path.join(test_path,
                         'microtest',
                         '1234A-V3LOOP_S1_L001_R2_001.fastq'),
            contigs_csv,
            remap_csv,
            remap_counts_csv,
            StringIO(),
            StringIO(),
            StringIO(),
            work_path=os.path.join(test_path, 'working'))

        self.assertEqual(expected_remap_counts_csv, remap_counts_csv.getvalue())
        self.assertEqual(expected_remap_csv, remap_csv.getvalue())

    def test_bad_contig(self):
        contigs_csv = StringIO("""\
ref,match,group_ref,contig
R1,1.0,R1,TCACCAGGACAGCGGGTTGAATTCCTCGTGCAAGCGTGGAA
R2,0.2,R2,TCACCAGGACAGCGGGTTGAATTCCTCGTGCAAGCGTGGAA
""")
        self.bowtie2_output.extend([
            "read1\t99\t1-R1\t1\t44\t5M\t=\t1\t-81\tGTGGG\tAAAAA\n",
            "read1\t147\t1-R1\t1\t44\t5M\t=\t1\t-81\tGTGGG\tAAAAA\n",
            "read2\t99\t1-R1\t1\t44\t5M\t=\t1\t-81\tGTGGG\tAAAAA\n",
            "read2\t147\t1-R1\t1\t44\t5M\t=\t1\t-81\tGTGGG\tAAAAA\n",
            "read3\t99\t1-R1\t1\t44\t5M\t=\t1\t-81\tGTGGG\tAAAAA\n",
            "read3\t147\t1-R1\t1\t44\t5M\t=\t1\t-81\tGTGGG\tAAAAA\n",
            "read4\t99\t2-R2-partial\t1\t44\t5M\t=\t1\t-81\tGTGGG\tAAAAA\n",
            "read4\t147\t2-R2-partial\t1\t44\t5M\t=\t1\t-81\tGTGGG\tAAAAA\n",
            "read5\t77\t*\t0\t0\t*\t*\t0\t0\tGTAAA\tAAAAA\n",
            "read5\t141\t*\t0\t0\t*\t*\t0\t0\tGTAAA\tAAAAA\n"])
        expected_remap_counts_csv = """\
type,count,filtered_count,seed_dist,other_dist,other_seed
raw,20,,,,
remap 1-R1,6,,,,
remap 2-R2-partial,2,,,,
remap-final 1-R1,6,,,,
remap-final 2-R2-partial,2,,,,
unmapped,2,,,,
"""
        expected_remap_csv = """\
qname,flag,rname,pos,mapq,cigar,rnext,pnext,tlen,seq,qual
read1,99,1-R1,1,44,5M,=,1,-81,GTGGG,AAAAA
read1,147,1-R1,1,44,5M,=,1,-81,GTGGG,AAAAA
read2,99,1-R1,1,44,5M,=,1,-81,GTGGG,AAAAA
read2,147,1-R1,1,44,5M,=,1,-81,GTGGG,AAAAA
read3,99,1-R1,1,44,5M,=,1,-81,GTGGG,AAAAA
read3,147,1-R1,1,44,5M,=,1,-81,GTGGG,AAAAA
read4,99,2-R2-partial,1,44,5M,=,1,-81,GTGGG,AAAAA
read4,147,2-R2-partial,1,44,5M,=,1,-81,GTGGG,AAAAA
"""
        self.assertMapsToContigs(contigs_csv, expected_remap_csv, expected_remap_counts_csv)

    def test_excluded_contig(self):
        test_path = os.path.dirname(__file__)
        contigs_csv = StringIO("""\
ref,match,group_ref,contig
R1,1.0,R1,GTGGG
R2,1.0,R2,ACAAA
""")
        self.bowtie2_output.extend([
            "read1\t99\t1-R1\t1\t44\t5M\t=\t1\t-81\tGTGGG\tAAAAA\n",
            "read1\t147\t1-R1\t1\t44\t5M\t=\t1\t-81\tGTGGG\tAAAAA\n",
            "read2\t99\t1-R1\t1\t44\t5M\t=\t1\t-81\tGTGGG\tAAAAA\n",
            "read2\t147\t1-R1\t1\t44\t5M\t=\t1\t-81\tGTGGG\tAAAAA\n",
            "read3\t99\t1-R1\t1\t44\t5M\t=\t1\t-81\tGTGGG\tAAAAA\n",
            "read3\t147\t1-R1\t1\t44\t5M\t=\t1\t-81\tGTGGG\tAAAAA\n",
            "read4\t99\t2-R2-excluded\t1\t44\t5M\t=\t1\t-81\tGTGGG\tAAAAA\n",
            "read4\t147\t2-R2-excluded\t1\t44\t5M\t=\t1\t-81\tGTGGG\tAAAAA\n",
            "read5\t77\t*\t0\t0\t*\t*\t0\t0\tGTAAA\tAAAAA\n",
            "read5\t141\t*\t0\t0\t*\t*\t0\t0\tGTAAA\tAAAAA\n"])
        excluded_seeds = ['R2']
        expected_remap_counts_csv = """\
type,count,filtered_count,seed_dist,other_dist,other_seed
raw,20,,,,
remap 1-R1,6,,,,
remap 2-R2-excluded,2,,,,
remap-final 1-R1,6,,,,
remap-final 2-R2-excluded,2,,,,
unmapped,2,,,,
"""
        expected_remap_csv = """\
qname,flag,rname,pos,mapq,cigar,rnext,pnext,tlen,seq,qual
read1,99,1-R1,1,44,5M,=,1,-81,GTGGG,AAAAA
read1,147,1-R1,1,44,5M,=,1,-81,GTGGG,AAAAA
read2,99,1-R1,1,44,5M,=,1,-81,GTGGG,AAAAA
read2,147,1-R1,1,44,5M,=,1,-81,GTGGG,AAAAA
read3,99,1-R1,1,44,5M,=,1,-81,GTGGG,AAAAA
read3,147,1-R1,1,44,5M,=,1,-81,GTGGG,AAAAA
"""
        expected_remap_conseq_csv = """\
region,sequence
1-R1,GTGGG
"""
        remap_counts_csv = StringIO()
        remap_csv = StringIO()
        remap_conseq_csv = StringIO()

        remap.map_to_contigs(
            os.path.join(test_path,
                         'microtest',
                         '1234A-V3LOOP_S1_L001_R1_001.fastq'),
            os.path.join(test_path,
                         'microtest',
                         '1234A-V3LOOP_S1_L001_R2_001.fastq'),
            contigs_csv,
            remap_csv,
            remap_counts_csv,
            remap_conseq_csv,
            StringIO(),
            StringIO(),
            work_path=os.path.join(test_path, 'working'),
            excluded_seeds=excluded_seeds)

        self.assertEqual(expected_remap_counts_csv, remap_counts_csv.getvalue())
        self.assertEqual(expected_remap_csv, remap_csv.getvalue())
        self.assertEqual(expected_remap_conseq_csv, remap_conseq_csv.getvalue())


def test_read_contigs(projects):
    contigs_csv = StringIO("""\
ref,match,group_ref,contig
HCV-1a,1.0,HCV-1a,TCACCAGGACAGCGGGTTGAATTCCTCGTGCAAGCGTGGAA
HCV-2b,1.0,HCV-2b,GCCCGCCCCCTGATGGGGGCGACACTCCGCCA
""")
    expected_conseqs = {
        '1-HCV-1a': 'TCACCAGGACAGCGGGTTGAATTCCTCGTGCAAGCGTGGAA',
        '2-HCV-2b': 'GCCCGCCCCCTGATGGGGGCGACACTCCGCCA'}

    conseqs = read_contigs(contigs_csv)

    assert expected_conseqs == conseqs


def test_read_contigs_filter():
    contigs_csv = StringIO("""\
ref,match,group_ref,contig
HCV-1a,0.24,HCV-1a,TCACCAGGACAGCGGGTTGAATTCCTCGTGCAAGCGTGGAA
HCV-2b,0.25,HCV-2b,GCCCGCCCCCTGATGGGGGCGACACTCCGCCA
""")
    expected_conseqs = {
        '1-HCV-1a-partial': 'TCACCAGGACAGCGGGTTGAATTCCTCGTGCAAGCGTGGAA',
        '2-HCV-2b': 'GCCCGCCCCCTGATGGGGGCGACACTCCGCCA'}

    conseqs = read_contigs(contigs_csv)

    assert expected_conseqs == conseqs


def test_read_contigs_reversed(projects):
    contigs_csv = StringIO("""\
ref,match,group_ref,contig
HCV-1a,-1.0,HCV-1a,TCACCAGGACAGCGGGTTGAATTCCTCGTGCAAGCGTGGAA
HCV-2b,-0.1,HCV-2b,GCCCGCCCCCTGATGGGGGCGACACTCCGCCA
""")
    expected_conseqs = {
        '1-HCV-1a-reversed': 'TCACCAGGACAGCGGGTTGAATTCCTCGTGCAAGCGTGGAA',
        '2-HCV-2b-reversed': 'GCCCGCCCCCTGATGGGGGCGACACTCCGCCA'}

    conseqs = read_contigs(contigs_csv)

    assert expected_conseqs == conseqs


def test_read_contigs_excluded():
    contigs_csv = StringIO("""\
ref,match,group_ref,contig
HCV-1a,1.0,HCV-1a,TCACCAGGACAGCGGGTTGAATTCCTCGTGCAAGCGTGGAA
HLA-B-seed,0.02,HLA-B-seed,ATGCGGGTCACGGCACCCCGAACCGT
""")
    excluded_seeds = ['HLA-B-seed']
    expected_conseqs = {
        '1-HCV-1a': 'TCACCAGGACAGCGGGTTGAATTCCTCGTGCAAGCGTGGAA',
        '2-HLA-B-seed-excluded': 'ATGCGGGTCACGGCACCCCGAACCGT'}

    conseqs = read_contigs(contigs_csv, excluded_seeds)

    assert expected_conseqs == conseqs


def test_read_contigs_mutations():
    contigs_csv = StringIO("""\
ref,match,group_ref,contig
HCV-1a,1.0,HCV-1a,TCACCAGGACAGCGGGTTGAATTCCTCGTGCAAGCGTGGAA
HCV-2b,1.0,HCV-2b,GCCCGACCCATGATGGGGGCGACACTCCGCCA
""")
    expected_conseqs = {
        '1-HCV-1a': 'TCACCAGGACAGCGGGTTGAATTCCTCGTGCAAGCGTGGAA',
        '2-HCV-2b': 'GCCCGACCCATGATGGGGGCGACACTCCGCCA'}

    conseqs = read_contigs(contigs_csv)

    assert expected_conseqs == conseqs


def test_read_contigs_merged():
    contigs_csv = StringIO("""\
ref,match,group_ref,contig
HCV-1a,1.0,HCV-1a,TCACCAGGACAGCGGGTTGAATTCCTCGTGCAAGCGTGGAA
HCV-2b,1.0,HCV-2b,GCCCGCCCCCTGATGGGGGCGACACTCCTCCA
HCV-2a,1.0,HCV-2b,CTCCACCATGAATCACTCCCCTG
""")
    # Changes:        ^G->A                   ^G->T
    # TODO: Merge contigs 2 and 3, because they overlap.
    #   'GCCCGCCCCCTGATGGGGGCGACACTCCTCCATGAATCACTCCCCTG'
    # Change:                        ^
    expected_conseqs = {
        '1-HCV-1a': 'TCACCAGGACAGCGGGTTGAATTCCTCGTGCAAGCGTGGAA',
        '2-HCV-2b': 'GCCCGCCCCCTGATGGGGGCGACACTCCTCCA',
        '3-HCV-2a': 'CTCCACCATGAATCACTCCCCTG'}

    conseqs = read_contigs(contigs_csv)

    assert expected_conseqs == conseqs


def test_read_contigs_untrimmed_left():
    """ Don't trim contigs that extend past reference start. """
    contigs_csv = StringIO("""\
ref,match,group_ref,contig
HCV-1a,1.0,HCV-1a,TCACCAGGACAGCGGGTTGAATTCCTCGTGCAAGCGTGGAA
HCV-2b,1.0,HCV-2b,TAGACATATTACCGCCCGCCCCCTGATGGGGGCGACACTCCGCCATGAATCACTCCCCTGT
""")
    expected_conseqs = {
        '1-HCV-1a': 'TCACCAGGACAGCGGGTTGAATTCCTCGTGCAAGCGTGGAA',
        '2-HCV-2b': 'TAGACATATTACCGCCCGCCCCCTGATGGGGGCGACACTCCGCCATGAATCACTCCCCTGT'}

    conseqs = read_contigs(contigs_csv)

    assert expected_conseqs == conseqs


def test_read_contigs_untrimmed_right():
    """ Don't trim contigs that extend past reference end. """
    contigs_csv = StringIO("""\
ref,match,group_ref,contig
HCV-1a,1.0,HCV-1a,TCACCAGGACAGCGGGTTGAATTCCTCGTGCAAGCGTGGAA
HCV-2b,1.0,HCV-2b,TAGTTTCCGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTAGACATATTACC
""")
    expected_conseqs = {
        '1-HCV-1a': 'TCACCAGGACAGCGGGTTGAATTCCTCGTGCAAGCGTGGAA',
        '2-HCV-2b': 'TAGTTTCCGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTAGACATATTACC'}

    conseqs = read_contigs(contigs_csv)

    assert expected_conseqs == conseqs


def test_full_remap(tmp_path):
    """ Test the full process of the remapping step. """
    microtest_path = Path(__file__).parent / 'microtest'
    fastq1 = microtest_path / '1234A-V3LOOP_S1_L001_R1_001.fastq'
    fastq2 = microtest_path / '1234A-V3LOOP_S1_L001_R2_001.fastq'
    prelim_csv = tmp_path / 'prelim.csv'
    prelim_csv.write_text("""\
qname,flag,rname,pos,mapq,cigar,rnext,pnext,tlen,seq,qual
M01234:01:000000000-AAAAA:1:1101:01234:0001,99,HIV1-C-BR-JX140663-seed,6535,36,51M,=,6535,-51,\
TGCACAAGACCCAACAACAATACAAGAAAAAGTATAAGGATAGGACCAGGA,AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
M01234:01:000000000-AAAAA:1:1101:01234:0001,147,HIV1-C-BR-JX140663-seed,6535,36,51M,=,6535,-51,\
TGCACAAGACCCAACAACAATACAAGAAAAAGTATAAGGATAGGACCAGGA,AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
""")
    remap_csv = tmp_path / 'remap.csv'
    remap_counts_csv = tmp_path / 'remap_counts.csv'
    remap_conseq_csv = tmp_path / 'remap_conseq.csv'
    unmapped1_fastq = tmp_path / 'unmapped1.fastq'
    unmapped2_fastq = tmp_path / 'unmapped2.fastq'
    expected_remap_counts = """\
type,count,filtered_count,seed_dist,other_dist,other_seed
raw,20,,,,
prelim HIV1-C-BR-JX140663-seed,2,2,,,
remap-1 HIV1-C-BR-JX140663-seed,20,,,,
remap-final HIV1-C-BR-JX140663-seed,20,,,,
unmapped,0,,,,
"""

    from micall.utils.work_dir import WorkDir
    with WorkDir.using(tmp_path):
        remap.remap(fastq1,
                    fastq2,
                    prelim_csv,
                    remap_csv,
                    remap_counts_csv,
                    remap_conseq_csv,
                    unmapped1_fastq,
                    unmapped2_fastq,
                    count_threshold=1)

    assert remap_counts_csv.read_text() == expected_remap_counts
