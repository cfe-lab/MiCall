import os
from io import StringIO
import unittest

from micall.g2p.pssm_lib import Pssm
from micall.g2p.fastq_g2p import fastq_g2p, FastqReader, FastqError, merge_reads, \
    trim_reads, count_reads, write_rows, write_unmapped_reads, write_aligned_reads, \
    extract_target, get_top_counts

TEMP_PREFIX = os.path.join(os.path.dirname(__file__), 'g2p_temp')


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
        self.pssm = self.g2p_csv = self.g2p_summary_csv = None
        prepare_g2p(self)

    def testSimple(self):
        counts = [((("TGTACAAGACCCAACAACAATACAAGAAAAAGAATCCGTATCCAGAGAGGACCAGGGA"
                    "GAGCATTTGTTACAATAGGAAAAATAGGAAATATGAGACAAGCACATTGT"), ("TGTACAAGA-------------------------------------------------"
                    "--------------------------------------------------"), ''), 1)]
        expected_g2p_csv = """\
rank,count,g2p,fpr,call,seq,aligned,error,comment
1,1,,,,CTR,,cysteines,
"""

        write_rows(self.pssm, counts, self.g2p_csv)

        self.assertEqual(expected_g2p_csv, self.g2p_csv.getvalue())

    def testSummarySuccess(self):
        counts = [((("TGTACAAGACCCAACAACAATACAAGAAAAAGAATCCGTATCCAGAGAGGACCAGGGA"
                    "GAGCATTT---GTTACAATAGGAAAAATAGGAAATATGAGACAAGCACATTGT"), ("TGTACAAGACCCAACAACAATACAAGAAAAA------GTATACATATAGGACCAGGGA"
                    "GAGCATTTTATGCAACAGGAGAAATAATAGGAGATATAAGACAAGCACATTGT"), ''), 1)]
        expected_g2p_csv = """\
rank,count,g2p,fpr,call,seq,aligned,error,comment
1,1,0.067754,42.3,R5,CTRPNNNTRKSIHIGPGRAFYATGEIIGDIRQAHC,CTRPN-NNT--RKSIHI---GPGR---AFYAT----GEIIGDI--RQAHC,,
"""
        expected_summary_csv = """\
mapped,valid,X4calls,X4pct,final,validpct
1,1,0,0.00,R5,100.00
"""

        write_rows(self.pssm, counts, self.g2p_csv, self.g2p_summary_csv)

        self.assertEqual(expected_g2p_csv, self.g2p_csv.getvalue())
        self.assertEqual(expected_summary_csv, self.g2p_summary_csv.getvalue())

    def testSummaryFailed(self):
        counts = [((("TGTACAAGACCCAACAACAATACAAGAAAAAGAATCCGTATCCAGAGAGGACCAGGGA"
                    "GAGCATTTGTTACAATAGGAAAAATAGGAAATATGAGACAAGCACATTGT"), ("TGTACAAGA-------------------------------------------------"
                    "--------------------------------------------------"), ''), 1)]
        expected_g2p_csv = """\
rank,count,g2p,fpr,call,seq,aligned,error,comment
1,1,,,,CTR,,cysteines,
"""
        expected_summary_csv = """\
mapped,valid,X4calls,X4pct,final,validpct
1,0,0,,,0.00
"""

        write_rows(self.pssm, counts, self.g2p_csv, self.g2p_summary_csv)

        self.assertEqual(expected_g2p_csv, self.g2p_csv.getvalue())
        self.assertEqual(expected_summary_csv, self.g2p_summary_csv.getvalue())

    def testSummaryX4(self):
        counts = [((("TGTACAAGACCCAACAACAATACAAGAAAAAGAATCCGTATCCAGAGAGGACCAGGGA"
                    "GAGCATTT---GTTACAATAGGAAAAATAGGAAATATGAGACAAGCACATTGT"), ("TGTATGAGACCCAACAACAATACAAGAAAAAGTATACATAT------AGGACCAGGGA"
                    "GAGCATTTTATGCAACAGGAGAAATAATAGGAGATATAAGACGAGCACATTGT"), ''), 2),
                  ((("TGTACAAGACCCAACAACAATACAAGAAAAAGAATCCGTATCCAGAGAGGACCAGGGA"
                    "GAGCATTT---GTTACAATAGGAAAAATAGGAAATATGAGACAAGCACATTGT"), ("TGTACAAGACCCAACAACAATACAAGAAAAAGTATACATAT------AGGACCAGGGA"
                    "GAGCATTTTATGCAACAGGAGAAATAATAGGAGATATAAGACAAGCACATTGT"), ''), 1)]
        expected_g2p_csv = """\
rank,count,g2p,fpr,call,seq,aligned,error,comment
1,2,0.454349,2.6,X4,CMRPNNNTRKSIHIGPGRAFYATGEIIGDIRRAHC,CMRPN-NNT--RKSIHI---GPGR---AFYAT----GEIIGDI--RRAHC,,
2,1,0.067754,42.3,R5,CTRPNNNTRKSIHIGPGRAFYATGEIIGDIRQAHC,CTRPN-NNT--RKSIHI---GPGR---AFYAT----GEIIGDI--RQAHC,,
"""
        expected_summary_csv = """\
mapped,valid,X4calls,X4pct,final,validpct
3,3,2,66.67,X4,100.00
"""

        write_rows(self.pssm, counts, self.g2p_csv, self.g2p_summary_csv)

        self.assertEqual(expected_g2p_csv, self.g2p_csv.getvalue())
        self.assertEqual(expected_summary_csv, self.g2p_summary_csv.getvalue())

    def testSummaryThresholdsPassed(self):
        counts = [((("TGTACAAGACCCAACAACAATACAAGAAAAAGAATCCGTATCCAGAGAGGACCAGGGA"
                    "GAGCATTT---GTTACAATAGGAAAAATAGGAAATATGAGACAAGCACATTGT"), ("TGTACAAGACCCAACAACAATACAAGAAAAA------GTATACATATAGGACCAGGGA"
                    "GAGCATTTTATGCAACAGGAGAAATAATAGGAGATATAAGACAAGCACATTGT"), ''), 300),
                  ((("TGTACAAGACCCAACAACAATACAAGAAAAAGAATCCGTATCCAGAGAGGACCAGGGA"
                    "GAGCATTTGTTACAATAGGAAAAATAGGAAATATGAGACAAGCACATTGT"), ("TGTACAAGA-------------------------------------------------"
                    "--------------------------------------------------"), ''), 100)]
        expected_g2p_csv = """\
rank,count,g2p,fpr,call,seq,aligned,error,comment
1,300,0.067754,42.3,R5,CTRPNNNTRKSIHIGPGRAFYATGEIIGDIRQAHC,CTRPN-NNT--RKSIHI---GPGR---AFYAT----GEIIGDI--RQAHC,,
2,100,,,,CTR,,cysteines,
"""
        expected_summary_csv = """\
mapped,valid,X4calls,X4pct,final,validpct
400,300,0,0.00,R5,75.00
"""

        write_rows(self.pssm,
                   counts,
                   self.g2p_csv,
                   self.g2p_summary_csv,
                   min_valid=300,
                   min_valid_percent=75.0)

        self.assertEqual(expected_g2p_csv, self.g2p_csv.getvalue())
        self.assertEqual(expected_summary_csv, self.g2p_summary_csv.getvalue())

    def testSummaryValidCountThresholdFailed(self):
        counts = [((("TGTACAAGACCCAACAACAATACAAGAAAAAGAATCCGTATCCAGAGAGGACCAGGGA"
                    "GAGCATTT---GTTACAATAGGAAAAATAGGAAATATGAGACAAGCACATTGT"), ("TGTACAAGACCCAACAACAATACAAGAAAAA------GTATACATATAGGACCAGGGA"
                    "GAGCATTTTATGCAACAGGAGAAATAATAGGAGATATAAGACAAGCACATTGT"), ''), 300),
                  ((("TGTACAAGACCCAACAACAATACAAGAAAAAGAATCCGTATCCAGAGAGGACCAGGGA"
                    "GAGCATTTGTTACAATAGGAAAAATAGGAAATATGAGACAAGCACATTGT"), ("TGTACAAGA-------------------------------------------------"
                    "--------------------------------------------------"), ''), 100)]
        expected_g2p_csv = """\
rank,count,g2p,fpr,call,seq,aligned,error,comment
1,300,0.067754,42.3,R5,CTRPNNNTRKSIHIGPGRAFYATGEIIGDIRQAHC,CTRPN-NNT--RKSIHI---GPGR---AFYAT----GEIIGDI--RQAHC,,
2,100,,,,CTR,,cysteines,
"""
        expected_summary_csv = """\
mapped,valid,X4calls,X4pct,final,validpct
400,300,0,0.00,,75.00
"""

        write_rows(self.pssm,
                   counts,
                   self.g2p_csv,
                   self.g2p_summary_csv,
                   min_valid=301,
                   min_valid_percent=75.0)

        self.assertEqual(expected_g2p_csv, self.g2p_csv.getvalue())
        self.assertEqual(expected_summary_csv, self.g2p_summary_csv.getvalue())

    def testSummaryValidPercentageThresholdFailed(self):
        counts = [((("TGTACAAGACCCAACAACAATACAAGAAAAAGAATCCGTATCCAGAGAGGACCAGGGA"
                    "GAGCATTT---GTTACAATAGGAAAAATAGGAAATATGAGACAAGCACATTGT"), ("TGTACAAGACCCAACAACAATACAAGAAAAA------GTATACATATAGGACCAGGGA"
                    "GAGCATTTTATGCAACAGGAGAAATAATAGGAGATATAAGACAAGCACATTGT"), ''), 300),
                  ((("TGTACAAGACCCAACAACAATACAAGAAAAAGAATCCGTATCCAGAGAGGACCAGGGA"
                    "GAGCATTTGTTACAATAGGAAAAATAGGAAATATGAGACAAGCACATTGT"), ("TGTACAAGA-------------------------------------------------"
                    "--------------------------------------------------"), ''), 100)]
        expected_g2p_csv = """\
rank,count,g2p,fpr,call,seq,aligned,error,comment
1,300,0.067754,42.3,R5,CTRPNNNTRKSIHIGPGRAFYATGEIIGDIRQAHC,CTRPN-NNT--RKSIHI---GPGR---AFYAT----GEIIGDI--RQAHC,,
2,100,,,,CTR,,cysteines,
"""
        expected_summary_csv = """\
mapped,valid,X4calls,X4pct,final,validpct
400,300,0,0.00,,75.00
"""

        write_rows(self.pssm,
                   counts,
                   self.g2p_csv,
                   self.g2p_summary_csv,
                   min_valid=300,
                   min_valid_percent=75.1)

        self.assertEqual(expected_g2p_csv, self.g2p_csv.getvalue())
        self.assertEqual(expected_summary_csv, self.g2p_summary_csv.getvalue())

    def testMinCount(self):
        counts = [((("TGTACAAGACCCAACAACAATACAAGAAAAAGAATCCGTATCCAGAGAGGACCAGGGA"
                    "GAGCATTTGTTACAATAGGAAAAATAGGAAATATGAGACAAGCACATTGT"), ("TGTACAAGA-------------------------------------------------"
                    "--------------------------------------------------"), ''), 3),
                  ((("TGTACAAGACCCAACAACAATACAAGAAAAAGAATCCGTATCCAGAGAGGACCAGGGA"
                    "GAGCATTTGTTACAATAGGAAAAATAGGAAATATGAGACAAGCACATTGT"), ("TGTACAGGG-------------------------------------------------"
                    "--------------------------------------------------"), ''), 2),
                  ((("TGTACAAGACCCAACAACAATACAAGAAAAAGAATCCGTATCCAGAGAGGACCAGGGA"
                    "GAGCATTTGTTACAATAGGAAAAATAGGAAATATGAGACAAGCACATTGT"), ("TGTACAGAA-------------------------------------------------"
                    "--------------------------------------------------"), ''), 2)]
        expected_g2p_csv = """\
rank,count,g2p,fpr,call,seq,aligned,error,comment
1,3,,,,CTR,,cysteines,
2,4,,,,,,count < 3,
"""
        expected_summary_csv = """\
mapped,valid,X4calls,X4pct,final,validpct
7,0,0,,,0.00
"""

        write_rows(self.pssm, counts, self.g2p_csv, self.g2p_summary_csv, min_count=3)

        self.assertEqual(expected_g2p_csv, self.g2p_csv.getvalue())
        self.assertEqual(expected_summary_csv, self.g2p_summary_csv.getvalue())

    def testSynonymMixture(self):
        """ Marking position 12 as low quality means codon 4 has to be P.
        """
        counts = [((("TGTACAAGACCCAACAACAATACAAGAAAAAGAATCCGTATCCAGAGAGGACCAGGGA"
                    "GAGCATTT---GTTACAATAGGAAAAATAGGAAATATGAGACAAGCACATTGT"), ("TGTACAAGACCNAACAACAATACAAGAAAAAG------TATACATATAGGACCAGGGA"
                    "GAGCATTTTATGCAACAGGAGAAATAATAGGAGATATAAGACAAGCACATTGT"), ''), 1)]
        expected_g2p_csv = """\
rank,count,g2p,fpr,call,seq,aligned,error,comment
1,1,0.067754,42.3,R5,CTRPNNNTRKSIHIGPGRAFYATGEIIGDIRQAHC,CTRPN-NNT--RKSIHI---GPGR---AFYAT----GEIIGDI--RQAHC,,
"""

        write_rows(self.pssm, counts, self.g2p_csv, self.g2p_summary_csv)

        self.assertEqual(expected_g2p_csv, self.g2p_csv.getvalue())

    def testAmbiguousMixture(self):
        """ Marking position 9 as low quality means codon 3 could be S or R.
        """
        counts = [((("TGTACAAGACCCAACAACAATACAAGAAAAAGAATCCGTATCCAGAGAGGACCAGGGA"
                    "GAGCATTT---GTTACAATAGGAAAAATAGGAAATATGAGACAAGCACATTGT"), ("TGTACAAGNCCCAACAACAATACAAGAAAAAG------TATACATATAGGACCAGGGA"
                    "GAGCATTTTATGCAACAGGAGAAATAATAGGAGATATAAGACAAGCACATTGT"), ''), 1)]
        expected_g2p_csv = """\
rank,count,g2p,fpr,call,seq,aligned,error,comment
1,1,0.066305,43.0,R5,CT[RS]PNNNTRKSIHIGPGRAFYATGEIIGDIRQAHC,CT[RS]PN-NNT--RKSIHI---GPGR---AFYAT----GEIIGDI--RQAHC,,ambiguous
"""

        write_rows(self.pssm, counts, self.g2p_csv, self.g2p_summary_csv)

        self.assertEqual(expected_g2p_csv, self.g2p_csv.getvalue())

    def testAmbiguousAtTwoPositions(self):
        """ Same thing with codons 9 and 18 - rejected. """
        counts = [((("TGTACAAGACCCAACAACAATACAAGAAAAAGAATCCGTATCCAGAGAGGACCAGGGA"
                    "GAGCATTT---GTTACAATAGGAAAAATAGGAAATATGAGACAAGCACATTGT"), ("TGTACAAGACCCAACAACAATACAAGNAAAAG------TATACATATAGGACCAGGGA"
                    "GNGCATTTTATGCAACAGGAGAAATAATAGGAGATATAAGACAAGCACATTGT"), ''), 1)]
        expected_g2p_csv = """\
rank,count,g2p,fpr,call,seq,aligned,error,comment
1,1,,,,CTRPNNNTXKSIHIGPGXAFYATGEIIGDIRQAHC,,> 2 ambiguous,
"""

        write_rows(self.pssm, counts, self.g2p_csv, self.g2p_summary_csv)

        self.assertEqual(expected_g2p_csv, self.g2p_csv.getvalue())

    def testAmbiguousMixtureThreeChoices(self):
        """ Marking position 14 as low quality means codon 5 could be L, S, or *.
        """
        counts = [((("TGTACAAGACCCAACAACAATACAAGAAAAAGAATCCGTATCCAGAGAGGACCAGGGA"
                    "GAGCATTTGTTACAATAGGAAAAATAGGAAATATGAGACAAGCACATTGT"), ("TGTACAAGACCCTNAAACTGT-------------------------------------"
                    "--------------------------------------------------"), ''), 1)]
        expected_g2p_csv = """\
rank,count,g2p,fpr,call,seq,aligned,error,comment
1,1,,,,CTRPXNC,,> 2 ambiguous,
"""

        write_rows(self.pssm, counts, self.g2p_csv, self.g2p_summary_csv)

        self.assertEqual(expected_g2p_csv, self.g2p_csv.getvalue())

    def testLowQuality(self):
        counts = [((("TGTACAAGACCCAACAACAATACAAGAAAAAGAATCCGTATCCAGAGAGGACCAGGGA"
                    "GAGCATTTGTTACAATAGGAAAAATAGGAAATATGAGACAAGCACATTGT"), ("TNTNNNGGN-------------------------------------------------"
                    "--------------------------------------------------"), ''), 1)]
        expected_g2p_csv = """\
rank,count,g2p,fpr,call,seq,aligned,error,comment
1,1,,,,,,low quality,
"""

        write_rows(self.pssm, counts, self.g2p_csv, self.g2p_summary_csv)

        self.assertEqual(expected_g2p_csv, self.g2p_csv.getvalue())

    def testPartialCodon(self):
        counts = [((("TGTACAAGACCCAACAACAATACAAGAAAAAGAATCCGTATCCAGAGAGGACCAGGGA"
                    "GAGCATTTGTTACAATAGGAAAAATAGGAAATATGAGACAAGCACATTGT"), ("TGTACAGG--------------------------------------------------"
                    "--------------------------------------------------"), ''), 1)]
        expected_g2p_csv = """\
rank,count,g2p,fpr,call,seq,aligned,error,comment
1,1,,,,CT,,notdiv3,
"""

        write_rows(self.pssm, counts, self.g2p_csv, self.g2p_summary_csv)

        self.assertEqual(expected_g2p_csv, self.g2p_csv.getvalue())

    def testStopCodon(self):
        counts = [((("TGTACAAGACCCAACAACAATACAAGAAAAAGAATCCGTATCCAGAGAGGACCAGGGA"
                    "GAGCATTTGTTACAATAGGAAAAATAGGAAATATGAGACAAGCACATTGT"), ("TGTTAGTGT-------------------------------------------------"
                    "--------------------------------------------------"), ''), 1)]
        expected_g2p_csv = """\
rank,count,g2p,fpr,call,seq,aligned,error,comment
1,1,,,,C*C,,stop codons,
"""

        write_rows(self.pssm, counts, self.g2p_csv, self.g2p_summary_csv)

        self.assertEqual(expected_g2p_csv, self.g2p_csv.getvalue())

    def testLengthMinimum(self):
        counts = [((("TGTACAAGACCCAACAACAATACAAGAAAAAGAATCCGTATCCAGAGAGGACCAGGGA"
                    "GAGCATTTGTTACAATAGGAAAAATAGGAAATATGAGACAAGCACATTGT"), ("TGTGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGAAA-------"
                    "-----GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGTGT"), ''), 1)]
        expected_g2p_csv = """\
rank,count,g2p,fpr,call,seq,aligned,error,comment
1,1,0.806327,1.5,X4,CGGGGGGGGGGGGGGGKGGGGGGGGGGGGGGC,---CG-GGG--GGGGGG---GGGG---GKGGG----GGGGGGG--GGGGC,,
"""

        write_rows(self.pssm, counts, self.g2p_csv, self.g2p_summary_csv)

        self.assertEqual(expected_g2p_csv, self.g2p_csv.getvalue())

    def testLengthTooShort(self):
        counts = [((("TGTACAAGACCCAACAACAATACAAGAAAAAGAATCCGTATCCAGAGAGGACCAGGGA"
                    "GAGCATTTGTTACAATAGGAAAAATAGGAAATATGAGACAAGCACATTGT"), ("TGTGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGAAA-------"
                    "--------GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGTGT"), ''), 1)]
        expected_g2p_csv = """\
rank,count,g2p,fpr,call,seq,aligned,error,comment
1,1,,,,CGGGGGGGGGGGGGGGKGGGGGGGGGGGGGC,,length,
"""

        write_rows(self.pssm, counts, self.g2p_csv, self.g2p_summary_csv)

        self.assertEqual(expected_g2p_csv, self.g2p_csv.getvalue())


class FastqG2PTest(unittest.TestCase):
    def setUp(self):
        super().setUp()
        self.pssm = self.g2p_csv = self.g2p_summary_csv = None
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
mapped,valid,X4calls,X4pct,final,validpct
1,1,0,0.00,R5,100.00
"""

        fastq_g2p(self.pssm,
                  fastq1,
                  fastq2,
                  self.g2p_csv,
                  self.g2p_summary_csv)

        self.assertEqual(expected_g2p_csv, self.g2p_csv.getvalue())
        self.assertEqual(expected_summary_csv, self.g2p_summary_csv.getvalue())


class ExtractTargetTest(unittest.TestCase):
    def test_extract(self):
        seed = 'GCGGGGCATGCGAGTACA'  # Amino seq: AGHAST
        coord = 'HAS'
        expected_extract = 'CATGCGAGT'  # Amino seq: HAS

        extract = extract_target(seed, coord)

        self.assertEqual(expected_extract, extract)
        
    def test_reading_frame2(self):
        seed = 'AGCGGGGCATGCGAGTACA'  # Amino seq: ?AGHAST (shifted two bases)
        coord = 'HAS'
        expected_extract = 'CATGCGAGT'  # Amino seq: HAS

        extract = extract_target(seed, coord)

        self.assertEqual(expected_extract, extract)
        
    def test_reading_frame1(self):
        seed = 'AAGCGGGGCATGCGAGTACA'  # Amino seq: ?AGHAST (shifted one base)
        coord = 'HAS'
        expected_extract = 'CATGCGAGT'  # Amino seq: HAS

        extract = extract_target(seed, coord)

        self.assertEqual(expected_extract, extract)


# noinspection DuplicatedCode
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
        reader = FastqReader(self.fastq1, self.fastq2)

        with self.assertRaisesRegex(FastqError, 'No match for read A:B:C.'):
            list(reader)

    def test_extra_description(self):
        self.fastq1 = StringIO("""\
@A:B:C X:Y Z
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
                           ("X:Y Z", "ACGT", "QUAL"),
                           ("Q:R", "TTGG", "LAUQ"))]
        reader = FastqReader(self.fastq1, self.fastq2)

        reads = list(reader)

        self.assertEqual(expected_reads, reads)


class MergeReadsTest(unittest.TestCase):
    def test_overlap(self):
        reads = [("A:B:C",
                  ("X:Y", "AAACCCTTTGGGAAA", "BBBBBBBBBBBBBBB"),
                  ("Q:R", "GGGTTTCCCAAA", "BBBBBBBBBBBB"))]
        expected_merged_reads = [("A:B:C",
                                  ("X:Y", "AAACCCTTTGGGAAA", "BBBBBBBBBBBBBBB"),
                                  ("Q:R", "GGGTTTCCCAAA", "BBBBBBBBBBBB"),
                                  "AAACCCTTTGGGAAACCC",
                                  "BBBBBBBBBBBBBBBBBB")]

        merged_reads = list(merge_reads(reads))

        self.assertEqual(expected_merged_reads, merged_reads)

    def test_multiple_reads(self):
        reads = [("A:B:C",
                  ("X:Y", "AAACCCTTTGGGAAA", "BBBBBBBBBBBBBBB"),
                  ("Q:R", "GGGTTTCCCAAA", "BBBBBBBBBBBB")),
                 ("A:B:E",
                  ("X:Y", "TTTCCCTTTGGGAAA", "BBBBBBBBBBBBBBB"),
                  ("Q:R", "GGGTTTCCCAAA", "BBBBBBBBBBBB"))]
        expected_merged_reads = [("A:B:C",
                                  ("X:Y", "AAACCCTTTGGGAAA", "BBBBBBBBBBBBBBB"),
                                  ("Q:R", "GGGTTTCCCAAA", "BBBBBBBBBBBB"),
                                  "AAACCCTTTGGGAAACCC",
                                  "BBBBBBBBBBBBBBBBBB"),
                                 ("A:B:E",
                                  ("X:Y", "TTTCCCTTTGGGAAA", "BBBBBBBBBBBBBBB"),
                                  ("Q:R", "GGGTTTCCCAAA", "BBBBBBBBBBBB"),
                                  "TTTCCCTTTGGGAAACCC",
                                  "BBBBBBBBBBBBBBBBBB")]

        merged_reads = list(merge_reads(reads))

        self.assertEqual(expected_merged_reads, merged_reads)

    def test_disagreement(self):
        reads = [("A:B:C",
                  ("X:Y", "AAACCCTTTGGGAAA", "BBBBBBBBBBBBBBB"),
                  ("Q:R", "GGGTTTCACAAA", "@@@@@@@Y@@@@"))]
        expected_merged_reads = [("A:B:C",
                                  ("X:Y", "AAACCCTTTGGGAAA", "BBBBBBBBBBBBBBB"),
                                  ("Q:R", "GGGTTTCACAAA", "@@@@@@@Y@@@@"),
                                  "AAACCCTTTGTGAAACCC",
                                  "BBBBBBBBBBYBBBB@@@")]

        merged_reads = list(merge_reads(reads))

        self.assertEqual(expected_merged_reads, merged_reads)

    def test_disagreement_ambiguous_quality(self):
        """ Mates disagree with high quality: merged N carries no quality.

        Mate 1 calls A and mate 2 calls G at the same position, both Q40,
        so merge_pairs resolves the base to N. The merged quality must not
        report Q40, otherwise the ambiguous base could later count as
        Q30-qualified insertion support.
        """
        reads = [("A:B:C",
                  ("X:Y", "AAACCCTTTGGGAAA", "IIIIIIIIIIIIIII"),
                  ("Q:R", "GGGTTTCTCAAA", "IIIIIIIIIIII"))]
        expected_merged_reads = [("A:B:C",
                                  ("X:Y", "AAACCCTTTGGGAAA", "IIIIIIIIIIIIIII"),
                                  ("Q:R", "GGGTTTCTCAAA", "IIIIIIIIIIII"),
                                  "AAACCCTTTGNGAAACCC",
                                  "IIIIIIIIII!IIIIIII")]

        merged_reads = list(merge_reads(reads))

        self.assertEqual(expected_merged_reads, merged_reads)

    def test_low_quality(self):
        reads = [("A:B:C",
                  ("X:Y", "AAACCCTTTGGGAAA", "B!BBBBBBBBBBBBB"),
                  ("Q:R", "GGGTTTCCCAAA", "BBBBBBBBBBBB"))]
        expected_merged_reads = [("A:B:C",
                                  ("X:Y", "AAACCCTTTGGGAAA", "B!BBBBBBBBBBBBB"),
                                  ("Q:R", "GGGTTTCCCAAA", "BBBBBBBBBBBB"),
                                  "ANACCCTTTGGGAAACCC",
                                  "B!BBBBBBBBBBBBBBBB")]

        merged_reads = list(merge_reads(reads))

        self.assertEqual(expected_merged_reads, merged_reads)

    def test_no_overlap(self):
        reads = [("A:B:C",
                  ("X:Y", "AAACCCTTTGGGAAA", "BBBBBBBBBBBBBBB"),
                  ("Q:R", "ACACACACACAC", "BBBBBBBBBBBB"))]
        expected_merged_reads = [("A:B:C",
                                  ("X:Y", "AAACCCTTTGGGAAA", "BBBBBBBBBBBBBBB"),
                                  ("Q:R", "ACACACACACAC", "BBBBBBBBBBBB"),
                                  None,
                                  None)]

        merged_reads = list(merge_reads(reads))

        self.assertEqual(expected_merged_reads, merged_reads)

    def test_reverse_overlap(self):
        """ The start of read 1 aligns with the end of read 2. Reject! """
        reads = [("A:B:C",
                  ("Q:R", "TTTGGGAAACCC", "BBBBBBBBBBBB"),
                  ("X:Y", "TTTCCCAAAGGGTTT", "BBBBBBBBBBBBBBB"))]
        expected_merged_reads = [("A:B:C",
                                  ("Q:R", "TTTGGGAAACCC", "BBBBBBBBBBBB"),
                                  ("X:Y", "TTTCCCAAAGGGTTT", "BBBBBBBBBBBBBBB"),
                                  None,
                                  None)]

        merged_reads = list(merge_reads(reads))

        self.assertEqual(expected_merged_reads, merged_reads)


class TrimReadsTest(unittest.TestCase):
    def test_untrimmed(self):
        v3loop_ref = 'TGTACAAGACCCAACAAC'
        reads = [("A:B:C",
                  ("name1", "bases1", "qual1"),
                  ("name2", "bases2", "qual2"),
                  "TGTACAAGACC",
                  "BBBBBBBBBBB")]
        expected_reads = [
            ("A:B:C",
             ("name1", "bases1", "qual1"),
             ("name2", "bases2", "qual2"),
             ('TGTACAAGACCCAACAAC',
              'TGTACAAGACC-------',
              'BBBBBBBBBBB       '))]

        trimmed_reads = list(trim_reads(reads, v3loop_ref))

        self.assertEqual(expected_reads, trimmed_reads)

    def test_trimmed(self):
        v3loop_ref = 'TGTACAAGACCCAACAAC'
        reads = [("A:B:C",
                  ("name1", "bases1", "qual1"),
                  ("name2", "bases2", "qual2"),
                  "AAGTGTACAAGACC",
                  "BBBBBBBBBBBBBB")]
        expected_reads = [
            ("A:B:C",
             ("name1", "bases1", "qual1"),
             ("name2", "bases2", "qual2"),
             ('TGTACAAGACCCAACAAC',
              'TGTACAAGACC-------',
              'BBBBBBBBBBB       '))]

        trimmed_reads = list(trim_reads(reads, v3loop_ref))

        self.assertEqual(expected_reads, trimmed_reads)

    def test_multiple(self):
        v3loop_ref = 'TGTACAAGACCCAACAAC'
        reads = [("A:B:C",
                  ("name1", "bases1", "qual1"),
                  ("name2", "bases2", "qual2"),
                  "TGTACAAGACC",
                  "BBBBBBBBBBB"),
                 ("A:B:E",
                  ("name1", "bases1", "qual1"),
                  ("name2", "bases2", "qual2"),
                  "AAGTGTACAAGACC",
                  "BBBBBBBBBBBBBB")]
        expected_reads = [
            ("A:B:C",
             ("name1", "bases1", "qual1"),
             ("name2", "bases2", "qual2"),
             ('TGTACAAGACCCAACAAC',
              'TGTACAAGACC-------',
              'BBBBBBBBBBB       ')),
            ("A:B:E",
             ("name1", "bases1", "qual1"),
             ("name2", "bases2", "qual2"),
             ('TGTACAAGACCCAACAAC',
              'TGTACAAGACC-------',
              'BBBBBBBBBBB       '))]

        trimmed_reads = list(trim_reads(reads, v3loop_ref))

        self.assertEqual(expected_reads, trimmed_reads)

    def test_not_v3loop(self):
        v3loop_ref = 'TGTACAAGACCCAACAACAATACAAGA'
        reads = [("A:B:C",
                  ("name1", "bases1", "qual1"),
                  ("name2", "bases2", "qual2"),
                  "ATATATATATAT",
                  "BBBBBBBBBBBB")]
        expected_reads = [("A:B:C",
                           ("name1", "bases1", "qual1"),
                           ("name2", "bases2", "qual2"),
                           (None, None, None))]

        trimmed_reads = list(trim_reads(reads, v3loop_ref))

        self.assertEqual(expected_reads, trimmed_reads)

    def test_not_merged(self):
        v3loop_ref = 'TGTACAAGACCCAACAAC'
        reads = [("A:B:C",
                  ("name1", "bases1", "qual1"),
                  ("name2", "bases2", "qual2"),
                  None,
                  None)]
        expected_reads = [("A:B:C",
                           ("name1", "bases1", "qual1"),
                           ("name2", "bases2", "qual2"),
                           (None, None, None))]

        trimmed_reads = list(trim_reads(reads, v3loop_ref))

        self.assertEqual(expected_reads, trimmed_reads)


# noinspection DuplicatedCode
class WriteUnmappedTest(unittest.TestCase):
    def test_unmapped(self):
        reads = [("A:B:C",
                  ("name1", "bases1", "qual1"),
                  ("name2", "bases2", "qual2"),
                  (None, None, None))]
        expected_reads = []
        expected_unmapped1 = """\
@A:B:C name1
bases1
+
qual1
"""
        expected_unmapped2 = """\
@A:B:C name2
bases2
+
qual2
"""
        unmapped1 = StringIO()
        unmapped2 = StringIO()

        filtered_reads = list(write_unmapped_reads(reads, unmapped1, unmapped2))

        self.assertEqual(expected_reads, filtered_reads)
        self.assertEqual(expected_unmapped1, unmapped1.getvalue())
        self.assertEqual(expected_unmapped2, unmapped2.getvalue())

    def test_mapped(self):
        reads = [("A:B:C",
                  ("name1", "bases1", "qual1"),
                  ("name2", "bases2", "qual2"),
                  ("ref", "seq", "qual"))]
        expected_reads = [("ref", "seq", "qual")]
        expected_unmapped1 = ""
        expected_unmapped2 = ""
        unmapped1 = StringIO()
        unmapped2 = StringIO()

        filtered_reads = list(write_unmapped_reads(reads, unmapped1, unmapped2))

        self.assertEqual(expected_reads, filtered_reads)
        self.assertEqual(expected_unmapped1, unmapped1.getvalue())
        self.assertEqual(expected_unmapped2, unmapped2.getvalue())

    def test_no_files(self):
        """ Don't write the files, but still do the filtering. """
        reads = [("A:B:C",
                  ("name1", "bases1", "qual1"),
                  ("name2", "bases2", "qual2"),
                  ("ref", "seq", "qual")),
                 ("A:B:D",
                  ("name1", "bases1", "qual1"),
                  ("name2", "bases2", "qual2"),
                  (None, None, None))]
        expected_reads = [("ref", "seq", "qual")]
        unmapped1 = None
        unmapped2 = None

        filtered_reads = list(write_unmapped_reads(reads, unmapped1, unmapped2))

        self.assertEqual(expected_reads, filtered_reads)


class CountReadsTest(unittest.TestCase):
    def test_counts(self):
        reads = [("TGTACAAGACACACA", "TGTACAAGA------", "BBBBBBBBBBBBBBB"),
                 ("TGTACAAGACACACA", "AGAACAAGA------", "BBBBBBBBBBBBBBB"),
                 ("TGTACAAGACACACA", "TGTACAAGA------", "BBBBBBBBBBBBBBB")]
        expected_counts = [(("TGTACAAGACACACA",
                             "AGAACAAGA------",
                             ()), 1),
                           (("TGTACAAGACACACA",
                             "TGTACAAGA------",
                             ()), 2)]

        counts = list(count_reads(reads, file_prefix=TEMP_PREFIX))

        self.assertEqual(expected_counts, counts)

    def test_without_temp_files(self):
        reads = [("TGTACAAGACACACA", "TGTACAAGA------", "BBBBBBBBBBBBBBB"),
                 ("TGTACAAGACACACA", "AGAACAAGA------", "BBBBBBBBBBBBBBB"),
                 ("TGTACAAGACACACA", "TGTACAAGA------", "BBBBBBBBBBBBBBB")]
        expected_counts = [(("TGTACAAGACACACA",
                             "AGAACAAGA------",
                             ()), 1),
                           (("TGTACAAGACACACA",
                             "TGTACAAGA------",
                             ()), 2)]

        counts = sorted(count_reads(reads, file_prefix=None))

        self.assertEqual(expected_counts, counts)


    def test_insertion_support_counts(self):
        """ One low-quality copy does not erase high-quality copies.

        Two identical alignments carry the same GGG insertion, but only
        one copy reaches Q30 on every inserted base. Total count stays 2
        while qualified insertion support is 1.
        """
        reads = [("TGTACA---AGACCCAAC",
                  "TGTACAGGGAGACCCAAC",
                  "BBBBBBBBBBBBBBBBBB"),
                 ("TGTACA---AGACCCAAC",
                  "TGTACAGGGAGACCCAAC",
                  "BBBBBB>>>BBBBBBBBB")]
        expected_counts = [(("TGTACA---AGACCCAAC",
                             "TGTACAGGGAGACCCAAC",
                             ((6, "GGG", 1),)), 2)]

        counts = list(count_reads(reads, file_prefix=None))

        self.assertEqual(expected_counts, counts)

    def test_ambiguous_insertion_no_support(self):
        """ An ambiguous merged base never counts as insertion support.

        The mates disagreed, so the inserted base became N with '!'
        quality. The alignment still counts, but the insertion gets no
        qualified support.
        """
        reads = [("TGTACA---AGACCCAAC",
                  "TGTACAGGNAGACCCAAC",
                  "BBBBBBBB!BBBBBBBBB")]
        expected_counts = [(("TGTACA---AGACCCAAC",
                             "TGTACAGGNAGACCCAAC",
                             ()), 1)]

        counts = list(count_reads(reads, file_prefix=None))

        self.assertEqual(expected_counts, counts)

class TopReadsTest(unittest.TestCase):
    def test_top(self):
        counts = iter([(("TGTACAAGACACACA", "AGAACAAGA------"), 1),
                       (("TGTACAAGACACACA", "TGTACAAGA------"), 2)])
        expected_top = [(("TGTACAAGACACACA", "TGTACAAGA------"), 2),
                        (("TGTACAAGACACACA", "AGAACAAGA------"), 1)]
        expected_ignored = 0

        top_counts, ignored_count = get_top_counts(counts)

        self.assertEqual(expected_top, top_counts)
        self.assertEqual(expected_ignored, ignored_count)

    def test_min_count(self):
        counts = iter([(("TGTACAAGACACACA", "AGAACAAGA------"), 1),
                       (("TGTACAAGACACACA", "CGAACAAGA------"), 4),
                       (("TGTACAAGACACACA", "GGAACAAGA------"), 2),
                       (("TGTACAAGACACACA", "TGAACAAGA------"), 3)])
        expected_top = [(("TGTACAAGACACACA", "CGAACAAGA------"), 4),
                        (("TGTACAAGACACACA", "TGAACAAGA------"), 3)]
        expected_ignored = 3

        top_counts, ignored_count = get_top_counts(counts, min_count=3)

        self.assertEqual(expected_top, top_counts)
        self.assertEqual(expected_ignored, ignored_count)


class WriteAlignedTest(unittest.TestCase):
    def test_counts(self):
        v3loop_ref = 'TGTACAAGACCCAACAAC'
        counts = [(("TGTACAAGACCCAAC", "TGTACAAGACCCAAC", ()), 2),
                  (("TGTACAAGACCCAAC", "AGAACAAGACCCAAC", ()), 1)]
        seed = "AAAAATGTACAAGACACAACAAC"
        aligned_csv = DummyFile()
        expected_aligned_csv = """\
refname,qcut,rank,count,offset,seq,inserts
HIV1-CON-XX-Consensus-seed,15,0,2,5,TGTACAAGACCCAAC,
HIV1-CON-XX-Consensus-seed,15,1,1,5,AGAACAAGACCCAAC,
"""

        counts2 = list(write_aligned_reads(counts, aligned_csv, seed, v3loop_ref))

        self.assertEqual(expected_aligned_csv, aligned_csv.getvalue())
        self.assertEqual(counts, counts2)

    def test_seed_offset(self):
        v3loop_ref = 'TGTACAAGACCCAACAAC'
        counts = [(("TGTACAAGACCCAAC", "TGTACAAGACCCAAC", ()), 2)]
        hiv_seed = "ATGTACAAGACACAACAAC"
        aligned_csv = DummyFile()
        expected_aligned_csv = """\
refname,qcut,rank,count,offset,seq,inserts
HIV1-CON-XX-Consensus-seed,15,0,2,1,TGTACAAGACCCAAC,
"""

        list(write_aligned_reads(counts, aligned_csv, hiv_seed, v3loop_ref))

        self.assertEqual(expected_aligned_csv, aligned_csv.getvalue())

    def test_seq_offset(self):
        v3loop_ref = 'TGTACAAGACCCAACAAC'
        counts = [(("TGTACAAGACCCAAC", "---ACAAGACCCAAC", ()), 2)]
        hiv_seed = "ATGTACAAGACCCAACAAC"
        aligned_csv = DummyFile()
        expected_aligned_csv = """\
refname,qcut,rank,count,offset,seq,inserts
HIV1-CON-XX-Consensus-seed,15,0,2,4,ACAAGACCCAAC,
"""

        list(write_aligned_reads(counts, aligned_csv, hiv_seed, v3loop_ref))

        self.assertEqual(expected_aligned_csv, aligned_csv.getvalue())

    def test_short_seq(self):
        v3loop_ref = 'TGTACAAGACCCAACAAC'
        counts = [(("TGTACAAGACCCAAC", "TGTACAAGACCC---", ()), 2)]
        hiv_seed = "ATGTACAAGACACAACAAC"
        aligned_csv = DummyFile()
        expected_aligned_csv = """\
refname,qcut,rank,count,offset,seq,inserts
HIV1-CON-XX-Consensus-seed,15,0,2,1,TGTACAAGACCC,
"""

        list(write_aligned_reads(counts, aligned_csv, hiv_seed, v3loop_ref))

        self.assertEqual(expected_aligned_csv, aligned_csv.getvalue())

    def test_seq_deletion(self):
        v3loop_ref = 'TGTACAAGACCCAACAAC'
        counts = [(("TGTACAAGACCCAAC", "TGT---AGACCCAAC", ()), 2)]
        hiv_seed = "ATGTACAAGACCCAACAAC"
        aligned_csv = DummyFile()
        expected_aligned_csv = """\
refname,qcut,rank,count,offset,seq,inserts
HIV1-CON-XX-Consensus-seed,15,0,2,1,TGT---AGACCCAAC,
"""

        list(write_aligned_reads(counts, aligned_csv, hiv_seed, v3loop_ref))

        self.assertEqual(expected_aligned_csv, aligned_csv.getvalue())

    def test_ref_deletion(self):
        v3loop_ref = 'TGTACAAGACCCAACAAC'
        counts = [(("TGTACA---AGACCCAAC",
                     "TGTACAGGGAGACCCAAC",
                     ((6, "GGG", 2),)), 2)]
        hiv_seed = "ATGTACAGGGAGACCCAACAAC"
        aligned_csv = DummyFile()
        expected_aligned_csv = """\
refname,qcut,rank,count,offset,seq,inserts
HIV1-CON-XX-Consensus-seed,15,0,2,1,TGTACA---AGACCCAAC,10:GGG:2
"""

        list(write_aligned_reads(counts, aligned_csv, hiv_seed, v3loop_ref))

        self.assertEqual(expected_aligned_csv, aligned_csv.getvalue())

    def test_ref_and_read_deletion(self):
        v3loop_ref = 'TGTACAAGACCCAACAAC'
        counts = [(("TGTACAAGACCCAAC", "TGTACAAGACCCAAC", ()), 2)]
        # deleted codon    vvv should be reported as dashes
        hiv_seed = "ATGTACAGGGAGACCCAACAACAATAC"
        aligned_csv = DummyFile()
        expected_aligned_csv = """\
refname,qcut,rank,count,offset,seq,inserts
HIV1-CON-XX-Consensus-seed,15,0,2,1,TGTACA---AGACCCAAC,
"""

        list(write_aligned_reads(counts, aligned_csv, hiv_seed, v3loop_ref))

        self.assertEqual(expected_aligned_csv, aligned_csv.getvalue())

    def test_ref_insertion(self):
        v3loop_ref = 'TGTACAAGACCCAACAAC'
        counts = [(("TGTACAAGACCCAAC", "TGTACAAGACCCAAC", ()), 2)]
        # inserted codon   ^^^ shouldn't be included in aligned seq.
        hiv_seed = "ATGTACACCCAACAAC"
        aligned_csv = DummyFile()
        expected_aligned_csv = """\
refname,qcut,rank,count,offset,seq,inserts
HIV1-CON-XX-Consensus-seed,15,0,2,1,TGTACACCCAAC,
"""

        list(write_aligned_reads(counts, aligned_csv, hiv_seed, v3loop_ref))

        self.assertEqual(expected_aligned_csv, aligned_csv.getvalue())

    def test_mixed_quality_insertion_support(self):
        """ Identical reads of mixed quality keep full count, split support.

        Two pairs merge to the same V3LOOP alignment, but only the first
        copy reaches Q30 on every inserted base. The group count stays 2
        while qualified insertion support is 1.
        """
        fastq1 = StringIO("""\
@p1 pair1
AAACCCTTTGGGAAA
+
BBBBBBBBBBBBBBB
@p2 pair1
AAACCCTTTGGGAAA
+
BBBBBB>>>BBBBBB
""")
        fastq2 = StringIO("""\
@p1 pair2
GGGTTTCCCAAA
+
BBBBBBBBBBBB
@p2 pair2
GGGTTTCCCAAA
+
BBBBBBBBB>>>
""")
        v3loop_ref = 'AAACCCTGGGAAACCC'
        hiv_seed = 'AAAA' + v3loop_ref + 'AAAA'
        aligned_csv = DummyFile()
        expected_aligned_csv = """\
refname,qcut,rank,count,offset,seq,inserts
HIV1-CON-XX-Consensus-seed,15,0,2,4,AAACCCTGGGAAACCC,11:TT:1
"""

        reader = FastqReader(fastq1, fastq2)
        merged_reads = merge_reads(reader)
        trimmed_reads = trim_reads(merged_reads, v3loop_ref)
        mapped_reads = write_unmapped_reads(trimmed_reads, None, None)
        read_counts = count_reads(mapped_reads, None)
        list(write_aligned_reads(read_counts, aligned_csv, hiv_seed, v3loop_ref))

        self.assertEqual(expected_aligned_csv, aligned_csv.getvalue())

    def test_read_insertion(self):
        v3loop_ref = 'TGTACAAGACCCAACAAC'
        # gaps in the V3 reference: GGG inserted in the read relative to V3LOOP.
        counts = [(("TGTACA---AGACCCAAC",
                     "TGTACAGGGAGACCCAAC",
                     ((6, "GGG", 2),)), 2)]
        hiv_seed = "ATGTACAAGACCCAACAAC"
        aligned_csv = DummyFile()
        expected_aligned_csv = """\
refname,qcut,rank,count,offset,seq,inserts
HIV1-CON-XX-Consensus-seed,15,0,2,1,TGTACAAGACCCAAC,7:GGG:2
"""

        list(write_aligned_reads(counts, aligned_csv, hiv_seed, v3loop_ref))

        self.assertEqual(expected_aligned_csv, aligned_csv.getvalue())
