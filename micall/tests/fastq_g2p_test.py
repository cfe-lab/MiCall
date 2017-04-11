import os
from io import StringIO
import unittest

from micall.g2p.pssm_lib import Pssm
from micall.g2p.fastq_g2p import fastq_g2p, FastqReader, FastqError, merge_reads, \
    trim_reads, count_reads, write_rows


class DummyFile(StringIO):
    def __repr__(self):
        s = self.getvalue()
        if len(s) > 50:
            s = '...' + s[-47:]
        return 'DummyFile({!r})'.format(s)


def prepare_g2p(test_case):
    if os.path.exists('../g2p/g2p_fpr.txt'):
        test_case.pssm = Pssm(path_to_lookup='../g2p/g2p_fpr.txt',
                              path_to_matrix='../g2p/g2p.matrix')
    else:
        test_case.pssm = Pssm(path_to_lookup='micall/g2p/g2p_fpr.txt',
                              path_to_matrix='micall/g2p/g2p.matrix')

    test_case.g2p_csv = DummyFile()
    test_case.g2p_summary_csv = DummyFile()


class WriteRowsTest(unittest.TestCase):
    def setUp(self):
        super().setUp()
        prepare_g2p(self)

    def testSimple(self):
        counts = [("TGTACAAGA", 1)]
        expected_g2p_csv = """\
rank,count,g2p,fpr,call,seq,aligned,error,comment
1,1,,,,CTR,,cysteines,
"""

        write_rows(self.pssm, counts, self.g2p_csv)

        self.assertEqual(expected_g2p_csv, self.g2p_csv.getvalue())

    def testSummarySuccess(self):
        counts = [("TGTACAAGACCCAACAACAATACAAGAAAAAGTATACATATAGGACCAGGGAGAGCAT"
                   "TTTATGCAACAGGAGAAATAATAGGAGATATAAGACAAGCACATTGT",
                   1)]
        expected_g2p_csv = """\
rank,count,g2p,fpr,call,seq,aligned,error,comment
1,1,0.067754,42.3,R5,CTRPNNNTRKSIHIGPGRAFYATGEIIGDIRQAHC,CTRPN-NNT--RKSIHI---GPGR---AFYAT----GEIIGDI--RQAHC,,
"""
        expected_summary_csv = """\
mapped,valid,X4calls,X4pct,final
1,1,0,0.00,R5
"""

        write_rows(self.pssm, counts, self.g2p_csv, self.g2p_summary_csv)

        self.assertEqual(expected_g2p_csv, self.g2p_csv.getvalue())
        self.assertEqual(expected_summary_csv, self.g2p_summary_csv.getvalue())

    def testSummaryFailed(self):
        counts = [("TGTACAAGA", 1)]
        expected_g2p_csv = """\
rank,count,g2p,fpr,call,seq,aligned,error,comment
1,1,,,,CTR,,cysteines,
"""
        expected_summary_csv = """\
mapped,valid,X4calls,X4pct,final
1,0,0,,
"""

        write_rows(self.pssm, counts, self.g2p_csv, self.g2p_summary_csv)

        self.assertEqual(expected_g2p_csv, self.g2p_csv.getvalue())
        self.assertEqual(expected_summary_csv, self.g2p_summary_csv.getvalue())

    def testSummaryX4(self):
        counts = [("TGTATGAGACCCAACAACAATACAAGAAAAAGTATACATATAGGACCAGGGAGAGC"
                   "ATTTTATGCAACAGGAGAAATAATAGGAGATATAAGACGAGCACATTGT",
                   2),
                  ("TGTACAAGACCCAACAACAATACAAGAAAAAGTATACATATAGGACCAGGGAGAGC"
                   "ATTTTATGCAACAGGAGAAATAATAGGAGATATAAGACAAGCACATTGT",
                   1)]
        expected_g2p_csv = """\
rank,count,g2p,fpr,call,seq,aligned,error,comment
1,2,0.454349,2.6,X4,CMRPNNNTRKSIHIGPGRAFYATGEIIGDIRRAHC,CMRPN-NNT--RKSIHI---GPGR---AFYAT----GEIIGDI--RRAHC,,
2,1,0.067754,42.3,R5,CTRPNNNTRKSIHIGPGRAFYATGEIIGDIRQAHC,CTRPN-NNT--RKSIHI---GPGR---AFYAT----GEIIGDI--RQAHC,,
"""
        expected_summary_csv = """\
mapped,valid,X4calls,X4pct,final
3,3,2,66.67,X4
"""

        write_rows(self.pssm, counts, self.g2p_csv, self.g2p_summary_csv)

        self.assertEqual(expected_g2p_csv, self.g2p_csv.getvalue())
        self.assertEqual(expected_summary_csv, self.g2p_summary_csv.getvalue())

    def testMinCount(self):
        counts = [("TGTACAAGA", 3),
                  ("TGTACAGGG", 2),
                  ("TGTACAGAA", 2)]
        expected_g2p_csv = """\
rank,count,g2p,fpr,call,seq,aligned,error,comment
1,3,,,,CTR,,cysteines,
2,4,,,,,,count < 3,
"""
        expected_summary_csv = """\
mapped,valid,X4calls,X4pct,final
7,0,0,,
"""

        write_rows(self.pssm, counts, self.g2p_csv, self.g2p_summary_csv, min_count=3)

        self.assertEqual(expected_g2p_csv, self.g2p_csv.getvalue())
        self.assertEqual(expected_summary_csv, self.g2p_summary_csv.getvalue())

    def testSynonymMixture(self):
        """ Marking position 12 as low quality means codon 4 has to be P.
        """
        counts = [("TGTACAAGACCNAACAACAATACAAGAAAAAGTATACATATAGGACCAGGGAGAGC"
                   "ATTTTATGCAACAGGAGAAATAATAGGAGATATAAGACAAGCACATTGT",
                   1)]
        expected_g2p_csv = """\
rank,count,g2p,fpr,call,seq,aligned,error,comment
1,1,0.067754,42.3,R5,CTRPNNNTRKSIHIGPGRAFYATGEIIGDIRQAHC,CTRPN-NNT--RKSIHI---GPGR---AFYAT----GEIIGDI--RQAHC,,
"""

        write_rows(self.pssm, counts, self.g2p_csv, self.g2p_summary_csv)

        self.assertEqual(expected_g2p_csv, self.g2p_csv.getvalue())

    def testAmbiguousMixture(self):
        """ Marking position 9 as low quality means codon 3 could be S or R.
        """
        counts = [("TGTACAAGNCCCAACAACAATACAAGAAAAAGTATACATATAGGACCAGGGAGAGC"
                   "ATTTTATGCAACAGGAGAAATAATAGGAGATATAAGACAAGCACATTGT",
                   1)]
        expected_g2p_csv = """\
rank,count,g2p,fpr,call,seq,aligned,error,comment
1,1,0.066305,43.0,R5,CT[RS]PNNNTRKSIHIGPGRAFYATGEIIGDIRQAHC,CT[RS]PN-NNT--RKSIHI---GPGR---AFYAT----GEIIGDI--RQAHC,,ambiguous
"""

        write_rows(self.pssm, counts, self.g2p_csv, self.g2p_summary_csv)

        self.assertEqual(expected_g2p_csv, self.g2p_csv.getvalue())

    def testAmbiguousAtTwoPositions(self):
        """ Same thing with codons 9 and 18 - rejected. """
        counts = [("TGTACAAGACCCAACAACAATACAAGNAAAAGTATACATATAGGACCAGGGAGNGC"
                   "ATTTTATGCAACAGGAGAAATAATAGGAGATATAAGACAAGCACATTGT",
                   1)]
        expected_g2p_csv = """\
rank,count,g2p,fpr,call,seq,aligned,error,comment
1,1,,,,CTRPNNNTXKSIHIGPGXAFYATGEIIGDIRQAHC,,> 2 ambiguous,
"""

        write_rows(self.pssm, counts, self.g2p_csv, self.g2p_summary_csv)

        self.assertEqual(expected_g2p_csv, self.g2p_csv.getvalue())

    def testAmbiguousMixtureThreeChoices(self):
        """ Marking position 14 as low quality means codon 5 could be L, S, or *.
        """
        counts = [("TGTACAAGACCCTNAAACTGT", 1)]
        expected_g2p_csv = """\
rank,count,g2p,fpr,call,seq,aligned,error,comment
1,1,,,,CTRPXNC,,> 2 ambiguous,
"""

        write_rows(self.pssm, counts, self.g2p_csv, self.g2p_summary_csv)

        self.assertEqual(expected_g2p_csv, self.g2p_csv.getvalue())

    def testLowQuality(self):
        counts = [("TNTNNNGGN", 1)]
        expected_g2p_csv = """\
rank,count,g2p,fpr,call,seq,aligned,error,comment
1,1,,,,,,low quality,
"""

        write_rows(self.pssm, counts, self.g2p_csv, self.g2p_summary_csv)

        self.assertEqual(expected_g2p_csv, self.g2p_csv.getvalue())

    def testPartialCodon(self):
        counts = [("TGTACAGG", 1)]
        expected_g2p_csv = """\
rank,count,g2p,fpr,call,seq,aligned,error,comment
1,1,,,,CT,,notdiv3,
"""

        write_rows(self.pssm, counts, self.g2p_csv, self.g2p_summary_csv)

        self.assertEqual(expected_g2p_csv, self.g2p_csv.getvalue())

    def testStopCodon(self):
        counts = [("TGTTAGTGT", 1)]
        expected_g2p_csv = """\
rank,count,g2p,fpr,call,seq,aligned,error,comment
1,1,,,,C*C,,stop codons,
"""

        write_rows(self.pssm, counts, self.g2p_csv, self.g2p_summary_csv)

        self.assertEqual(expected_g2p_csv, self.g2p_csv.getvalue())

    def testLengthMinimum(self):
        counts = [("TGTGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGAAA"
                   "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGTGT",
                   1)]
        expected_g2p_csv = """\
rank,count,g2p,fpr,call,seq,aligned,error,comment
1,1,0.806327,1.5,X4,CGGGGGGGGGGGGGGGKGGGGGGGGGGGGGGC,---CG-GGG--GGGGGG---GGGG---GKGGG----GGGGGGG--GGGGC,,
"""

        write_rows(self.pssm, counts, self.g2p_csv, self.g2p_summary_csv)

        self.assertEqual(expected_g2p_csv, self.g2p_csv.getvalue())

    def testLengthTooShort(self):
        counts = [("TGTGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGAAA"
                   "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGTGT",
                   1)]
        expected_g2p_csv = """\
rank,count,g2p,fpr,call,seq,aligned,error,comment
1,1,,,,CGGGGGGGGGGGGGGGKGGGGGGGGGGGGGC,,length,
"""

        write_rows(self.pssm, counts, self.g2p_csv, self.g2p_summary_csv)

        self.assertEqual(expected_g2p_csv, self.g2p_csv.getvalue())


class FastqG2PTest(unittest.TestCase):
    def setUp(self):
        super().setUp()
        prepare_g2p(self)

    def testSummarySuccess(self):
        fastq1 = StringIO("""\
@A:B:C X:Y
TGTACAAGACCCAACAACAATACAAGAAAAAGTATACATATAGGACCAGGGAGAGC
+
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
""")
        fastq2 = StringIO("""\
@A:B:C Q:R
ACAATGTGCTTGTCTTATATCTCCTATTATTTCTCCTGTTGCATAAAATGCTCTCC
+
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
""")
        expected_g2p_csv = """\
rank,count,g2p,fpr,call,seq,aligned,error,comment
1,1,0.067754,42.3,R5,CTRPNNNTRKSIHIGPGRAFYATGEIIGDIRQAHC,CTRPN-NNT--RKSIHI---GPGR---AFYAT----GEIIGDI--RQAHC,,
"""
        expected_summary_csv = """\
mapped,valid,X4calls,X4pct,final
1,1,0,0.00,R5
"""

        fastq_g2p(self.pssm, fastq1, fastq2, self.g2p_csv, self.g2p_summary_csv)

        self.assertEqual(expected_g2p_csv, self.g2p_csv.getvalue())
        self.assertEqual(expected_summary_csv, self.g2p_summary_csv.getvalue())


class FastqReaderTest(unittest.TestCase):
    def test_one_pair(self):
        self.fastq1 = StringIO("""\
@A:B:C X:Y
ACGT
+
QUAL
""")
        self.fastq2 = StringIO("""\
@A:B:C Q:R
TTGG
+
LAUQ
""")
        expected_reads = [("A:B:C",
                           ("X:Y", "ACGT", "QUAL"),
                           ("Q:R", "TTGG", "LAUQ"))]
        reader = FastqReader(self.fastq1, self.fastq2)

        reads = list(reader)

        self.assertEqual(expected_reads, reads)

    def test_two_pairs(self):
        self.fastq1 = StringIO("""\
@A:B:C X:Y
ACGT
+
QUAL
@A:B:E X:Y
AAGT
+
QUAL
""")
        self.fastq2 = StringIO("""\
@A:B:C Q:R
TTGG
+
LAUQ
@A:B:E Q:R
TAGA
+
LAUQ
""")
        expected_reads = [("A:B:C",
                           ("X:Y", "ACGT", "QUAL"),
                           ("Q:R", "TTGG", "LAUQ")),
                          ("A:B:E",
                           ("X:Y", "AAGT", "QUAL"),
                           ("Q:R", "TAGA", "LAUQ"))]
        reader = FastqReader(self.fastq1, self.fastq2)

        reads = list(reader)

        self.assertEqual(expected_reads, reads)

    def test_two_pairs_unordered(self):
        self.fastq1 = StringIO("""\
@A:B:C X:Y
ACGT
+
QUAL
@A:B:E X:Y
AAGT
+
QUAL
""")
        self.fastq2 = StringIO("""\
@A:B:E Q:R
TAGA
+
LAUQ
@A:B:C Q:R
TTGG
+
LAUQ
""")
        expected_reads = [("A:B:C",
                           ("X:Y", "ACGT", "QUAL"),
                           ("Q:R", "TTGG", "LAUQ")),
                          ("A:B:E",
                           ("X:Y", "AAGT", "QUAL"),
                           ("Q:R", "TAGA", "LAUQ"))]
        reader = FastqReader(self.fastq1, self.fastq2)

        reads = list(reader)

        self.assertEqual(expected_reads, reads)

    def test_unmatched_pair(self):
        self.fastq1 = StringIO("""\
@A:B:C X:Y
ACGT
+
QUAL
""")
        self.fastq2 = StringIO("""\
@A:B:X Q:R
TTGG
+
LAUQ
""")
        expected_reads = [("A:B:C",
                           ("X:Y", "ACGT", "QUAL"),
                           ("Q:R", "TTGG", "LAUQ"))]
        reader = FastqReader(self.fastq1, self.fastq2)

        with self.assertRaisesRegex(FastqError, 'No match for read A:B:C.'):
            list(reader)


class MergeReadsTest(unittest.TestCase):
    def test_overlap(self):
        reads = [("A:B:C",
                  ("X:Y", "AAACCCTTTGGGAAA", "BBBBBBBBBBBBBBB"),
                  ("Q:R", "GGGTTTCCCAAA", "BBBBBBBBBBBB"))]
        expected_merged_reads = [("A:B:C", "AAACCCTTTGGGAAACCC")]

        merged_reads = list(merge_reads(reads))

        self.assertEqual(expected_merged_reads, merged_reads)

    def test_multiple_reads(self):
        reads = [("A:B:C",
                  ("X:Y", "AAACCCTTTGGGAAA", "BBBBBBBBBBBBBBB"),
                  ("Q:R", "GGGTTTCCCAAA", "BBBBBBBBBBBB")),
                 ("A:B:E",
                  ("X:Y", "TTTCCCTTTGGGAAA", "BBBBBBBBBBBBBBB"),
                  ("Q:R", "GGGTTTCCCAAA", "BBBBBBBBBBBB"))]
        expected_merged_reads = [("A:B:C", "AAACCCTTTGGGAAACCC"),
                                 ("A:B:E", "TTTCCCTTTGGGAAACCC")]

        merged_reads = list(merge_reads(reads))

        self.assertEqual(expected_merged_reads, merged_reads)

    def test_disagreement(self):
        reads = [("A:B:C",
                  ("X:Y", "AAACCCTTTGGGAAA", "BBBBBBBBBBBBBBB"),
                  ("Q:R", "GGGTTTCACAAA", "@@@@@@@Y@@@@"))]
        expected_merged_reads = [("A:B:C", "AAACCCTTTGTGAAACCC")]

        merged_reads = list(merge_reads(reads))

        self.assertEqual(expected_merged_reads, merged_reads)

    def test_low_quality(self):
        reads = [("A:B:C",
                  ("X:Y", "AAACCCTTTGGGAAA", "B!BBBBBBBBBBBBB"),
                  ("Q:R", "GGGTTTCCCAAA", "BBBBBBBBBBBB"))]
        expected_merged_reads = [("A:B:C", "ANACCCTTTGGGAAACCC")]

        merged_reads = list(merge_reads(reads))

        self.assertEqual(expected_merged_reads, merged_reads)

    def test_no_overlap(self):
        reads = [("A:B:C",
                  ("X:Y", "AAACCCTTTGGGAAA", "BBBBBBBBBBBBBBB"),
                  ("Q:R", "ACACACACACAC", "BBBBBBBBBBBB"))]
        expected_merged_reads = []

        merged_reads = list(merge_reads(reads))

        self.assertEqual(expected_merged_reads, merged_reads)

    def test_reverse_overlap(self):
        """ The start of read 1 aligns with the end of read 2. Reject! """
        reads = [("A:B:C",
                  ("Q:R", "TTTGGGAAACCC", "BBBBBBBBBBBB"),
                  ("X:Y", "TTTCCCAAAGGGTTT", "BBBBBBBBBBBBBBB"))]
        expected_merged_reads = []

        merged_reads = list(merge_reads(reads))

        self.assertEqual(expected_merged_reads, merged_reads)


class TrimReadsTest(unittest.TestCase):
    def test_untrimmed(self):
        reads = [("A:B:C", "TGTACAAGACC")]
        expected_reads = [("A:B:C", "TGTACAAGACC")]

        trimmed_reads = list(trim_reads(reads))

        self.assertEqual(expected_reads, trimmed_reads)

    def test_trimmed(self):
        reads = [("A:B:C", "AAGTGTACAAGACC")]
        expected_reads = [("A:B:C", "TGTACAAGACC")]

        trimmed_reads = list(trim_reads(reads))

        self.assertEqual(expected_reads, trimmed_reads)

    def test_multiple(self):
        reads = [("A:B:C", "TGTACAAGACC"),
                 ("A:B:E", "AAGTGTACAAGACC")]
        expected_reads = [("A:B:C", "TGTACAAGACC"),
                          ("A:B:E", "TGTACAAGACC")]

        trimmed_reads = list(trim_reads(reads))

        self.assertEqual(expected_reads, trimmed_reads)

    def test_not_v3loop(self):
        reads = [("A:B:C", "ATATATATATAT")]
        expected_reads = []

        trimmed_reads = list(trim_reads(reads))

        self.assertEqual(expected_reads, trimmed_reads)


class CountReadsTest(unittest.TestCase):
    def test_counts(self):
        reads = [("A:B:C", "TGTACAAGA"),
                 ("A:B:D", "AGAACAAGA"),
                 ("A:B:E", "TGTACAAGA")]
        expected_counts = [("TGTACAAGA", 2),
                           ("AGAACAAGA", 1)]

        counts = count_reads(reads)

        self.assertEqual(expected_counts, counts)
