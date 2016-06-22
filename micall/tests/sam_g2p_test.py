import os
from StringIO import StringIO
import unittest

from micall.g2p.pssm_lib import Pssm
from micall.g2p.sam_g2p import sam_g2p


class DummyFile(StringIO):
    def __repr__(self):
        s = self.getvalue()
        if len(s) > 50:
            s = '...' + s[-47:]
        return 'DummyFile({!r})'.format(s)


class SamG2PTest(unittest.TestCase):
    def setUp(self):
        super(SamG2PTest, self).setUp()
        if os.path.exists('../g2p/g2p_fpr.txt'):
            self.pssm = Pssm(path_to_lookup='../g2p/g2p_fpr.txt',
                             path_to_matrix='../g2p/g2p.matrix')
        else:
            self.pssm = Pssm(path_to_lookup='micall/g2p/g2p_fpr.txt',
                             path_to_matrix='micall/g2p/g2p.matrix')

        self.nuc_csv = StringIO("""\
seed,region,q-cutoff,query.nuc.pos,refseq.nuc.pos,A,C,G,T
HIV1B-env-seed,V3LOOP,15,877,1,0,0,0,100
HIV1B-env-seed,V3LOOP,15,981,105,0,0,0,100
""")
        self.g2p_csv = DummyFile()
        self.g2p_summary_csv = DummyFile()
        self.addTypeEqualityFunc(str, self.assertMultiLineEqual)

    def testSimple(self):
        remap_csv = StringIO("""\
qname,flag,rname,pos,mapq,cigar,rnext,pnext,tlen,seq,qual
Example_read_1,99,HIV1B-env-seed,877,44,9M,=,877,9,TGTACAAGA,AAAAAAAAA
Example_read_1,147,HIV1B-env-seed,877,44,9M,=,877,-9,TGTACAAGA,AAAAAAAAA
""")
        expected_g2p_csv = """\
rank,count,g2p,fpr,call,seq,aligned,error
1,1,,,,CTR,,cysteines
"""

        sam_g2p(self.pssm, remap_csv, self.nuc_csv, self.g2p_csv)

        self.assertEqual(expected_g2p_csv, self.g2p_csv.getvalue())

    def testSummarySuccess(self):
        remap_csv = StringIO("""\
qname,flag,rname,pos,mapq,cigar,rnext,pnext,tlen,seq,qual
Example_read_1,99,HIV1B-env-seed,877,44,56M,=,926,56,TGTACAAGACCCAACAACAATACAAGAAAAAGTATACATATAGGACCAGGGAGAGC,AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
Example_read_1,147,HIV1B-env-seed,926,44,56M,=,877,-56,GGAGAGCATTTTATGCAACAGGAGAAATAATAGGAGATATAAGACAAGCACATTGT,AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
""")
        expected_g2p_csv = """\
rank,count,g2p,fpr,call,seq,aligned,error
1,1,0.0677537070158,42.3,R5,CTRPNNNTRKSIHIGPGRAFYATGEIIGDIRQAHC,CTRPN-NNT--RKSIHI---GPGR---AFYAT----GEIIGDI--RQAHC,
"""
        expected_summary_csv = """\
mapped,valid,X4calls,X4pct,final
1,1,0,0.00,R5
"""

        sam_g2p(self.pssm, remap_csv, self.nuc_csv, self.g2p_csv, self.g2p_summary_csv)

        self.assertEqual(expected_g2p_csv, self.g2p_csv.getvalue())
        self.assertEqual(expected_summary_csv, self.g2p_summary_csv.getvalue())

    def testSummaryFailed(self):
        remap_csv = StringIO("""\
qname,flag,rname,pos,mapq,cigar,rnext,pnext,tlen,seq,qual
Example_read_1,99,HIV1B-env-seed,877,44,9M,=,877,9,TGTACAAGA,AAAAAAAAA
Example_read_1,147,HIV1B-env-seed,877,44,9M,=,877,-9,TGTACAAGA,AAAAAAAAA
""")
        expected_g2p_csv = """\
rank,count,g2p,fpr,call,seq,aligned,error
1,1,,,,CTR,,cysteines
"""
        expected_summary_csv = """\
mapped,valid,X4calls,X4pct,final
1,0,0,,
"""

        sam_g2p(self.pssm, remap_csv, self.nuc_csv, self.g2p_csv, self.g2p_summary_csv)

        self.assertEqual(expected_g2p_csv, self.g2p_csv.getvalue())
        self.assertEqual(expected_summary_csv, self.g2p_summary_csv.getvalue())

    def testSummaryX4(self):
        remap_csv = StringIO("""\
qname,flag,rname,pos,mapq,cigar,rnext,pnext,tlen,seq,qual
Example_read_1,99,HIV1B-env-seed,877,44,56M,=,926,56,TGTACAAGACCCAACAACAATACAAGAAAAAGTATACATATAGGACCAGGGAGAGC,AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
Example_read_1,147,HIV1B-env-seed,926,44,56M,=,877,-56,GGAGAGCATTTTATGCAACAGGAGAAATAATAGGAGATATAAGACAAGCACATTGT,AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
Example_read_2,99,HIV1B-env-seed,877,44,56M,=,926,56,TGTATGAGACCCAACAACAATACAAGAAAAAGTATACATATAGGACCAGGGAGAGC,AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
Example_read_2,147,HIV1B-env-seed,926,44,56M,=,877,-56,GGAGAGCATTTTATGCAACAGGAGAAATAATAGGAGATATAAGACGAGCACATTGT,AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
Example_read_3,99,HIV1B-env-seed,877,44,56M,=,926,56,TGTATGAGACCCAACAACAATACAAGAAAAAGTATACATATAGGACCAGGGAGAGC,AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
Example_read_3,147,HIV1B-env-seed,926,44,56M,=,877,-56,GGAGAGCATTTTATGCAACAGGAGAAATAATAGGAGATATAAGACGAGCACATTGT,AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
""")
        expected_g2p_csv = """\
rank,count,g2p,fpr,call,seq,aligned,error
1,2,0.454349263704,2.6,X4,CMRPNNNTRKSIHIGPGRAFYATGEIIGDIRRAHC,CMRPN-NNT--RKSIHI---GPGR---AFYAT----GEIIGDI--RRAHC,
2,1,0.0677537070158,42.3,R5,CTRPNNNTRKSIHIGPGRAFYATGEIIGDIRQAHC,CTRPN-NNT--RKSIHI---GPGR---AFYAT----GEIIGDI--RQAHC,
"""
        expected_summary_csv = """\
mapped,valid,X4calls,X4pct,final
3,3,2,66.67,X4
"""

        sam_g2p(self.pssm, remap_csv, self.nuc_csv, self.g2p_csv, self.g2p_summary_csv)

        self.assertEqual(expected_g2p_csv, self.g2p_csv.getvalue())
        self.assertEqual(expected_summary_csv, self.g2p_summary_csv.getvalue())

    def testVariants(self):
        remap_csv = StringIO("""\
qname,flag,rname,pos,mapq,cigar,rnext,pnext,tlen,seq,qual
Example_read_1,99,HIV1B-env-seed,877,44,9M,=,877,9,TGTACAAGA,AAAAAAAAA
Example_read_1,147,HIV1B-env-seed,877,44,9M,=,877,-9,TGTACAAGA,AAAAAAAAA
Example_read_2,99,HIV1B-env-seed,877,44,9M,=,877,9,TGTACAAGA,AAAAAAAAA
Example_read_2,147,HIV1B-env-seed,877,44,9M,=,877,-9,TGTACAAGA,AAAAAAAAA
Example_read_3,99,HIV1B-env-seed,877,44,9M,=,877,9,TGTACAGGG,AAAAAAAAA
Example_read_3,147,HIV1B-env-seed,877,44,9M,=,877,-9,TGTACAGGG,AAAAAAAAA
""")
        expected_g2p_csv = """\
rank,count,g2p,fpr,call,seq,aligned,error
1,2,,,,CTR,,cysteines
2,1,,,,CTG,,cysteines
"""

        sam_g2p(self.pssm, remap_csv, self.nuc_csv, self.g2p_csv)

        self.assertEqual(expected_g2p_csv, self.g2p_csv.getvalue())

    def testMinCount(self):
        remap_csv = StringIO("""\
qname,flag,rname,pos,mapq,cigar,rnext,pnext,tlen,seq,qual
variant1_read1,99,HIV1B-env-seed,877,44,9M,=,877,9,TGTACAAGA,AAAAAAAAA
variant1_read1,147,HIV1B-env-seed,877,44,9M,=,877,-9,TGTACAAGA,AAAAAAAAA
variant1_read2,99,HIV1B-env-seed,877,44,9M,=,877,9,TGTACAAGA,AAAAAAAAA
variant1_read2,147,HIV1B-env-seed,877,44,9M,=,877,-9,TGTACAAGA,AAAAAAAAA
variant1_read3,99,HIV1B-env-seed,877,44,9M,=,877,9,TGTACAAGA,AAAAAAAAA
variant1_read3,147,HIV1B-env-seed,877,44,9M,=,877,-9,TGTACAAGA,AAAAAAAAA
variant2_read1,99,HIV1B-env-seed,877,44,9M,=,877,9,TGTACAGGG,AAAAAAAAA
variant2_read1,147,HIV1B-env-seed,877,44,9M,=,877,-9,TGTACAGGG,AAAAAAAAA
variant2_read2,99,HIV1B-env-seed,877,44,9M,=,877,9,TGTACAGGG,AAAAAAAAA
variant2_read2,147,HIV1B-env-seed,877,44,9M,=,877,-9,TGTACAGGG,AAAAAAAAA
variant3_read1,99,HIV1B-env-seed,877,44,9M,=,877,9,TGTACAGAA,AAAAAAAAA
variant3_read1,147,HIV1B-env-seed,877,44,9M,=,877,-9,TGTACAGAA,AAAAAAAAA
variant3_read2,99,HIV1B-env-seed,877,44,9M,=,877,9,TGTACAGAA,AAAAAAAAA
variant3_read2,147,HIV1B-env-seed,877,44,9M,=,877,-9,TGTACAGAA,AAAAAAAAA
""")
        expected_g2p_csv = """\
rank,count,g2p,fpr,call,seq,aligned,error
1,3,,,,CTR,,cysteines
2,4,,,,,,count < 3
"""
        expected_summary_csv = """\
mapped,valid,X4calls,X4pct,final
7,0,0,,
"""

        sam_g2p(self.pssm,
                remap_csv,
                self.nuc_csv,
                self.g2p_csv,
                self.g2p_summary_csv,
                min_count=3)

        self.assertEqual(expected_g2p_csv, self.g2p_csv.getvalue())
        self.assertEqual(expected_summary_csv, self.g2p_summary_csv.getvalue())

    def testOverlap(self):
        remap_csv = StringIO("""\
qname,flag,rname,pos,mapq,cigar,rnext,pnext,tlen,seq,qual
Example_read_1,99,HIV1B-env-seed,877,44,12M,=,886,12,TGTACAAGACCC,AAAAAAAAAAAA
Example_read_1,147,HIV1B-env-seed,886,44,9M,=,877,-9,CCCAACAAC,AAAAAAAAA
""")
        expected_g2p_csv = """\
rank,count,g2p,fpr,call,seq,aligned,error
1,1,,,,CTRPNN,,cysteines
"""

        sam_g2p(self.pssm, remap_csv, self.nuc_csv, self.g2p_csv)

        self.assertEqual(expected_g2p_csv, self.g2p_csv.getvalue())

    def testSynonymMixture(self):
        """ Marking position 12 as low quality means codon 4 has to be P.
        """
        remap_csv = StringIO("""\
qname,flag,rname,pos,mapq,cigar,rnext,pnext,tlen,seq,qual
Example_read_1,99,HIV1B-env-seed,877,44,56M,=,926,56,TGTACAAGACCCAACAACAATACAAGAAAAAGTATACATATAGGACCAGGGAGAGC,AAAAAAAAAAA#AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
Example_read_1,147,HIV1B-env-seed,926,44,56M,=,877,-56,GGAGAGCATTTTATGCAACAGGAGAAATAATAGGAGATATAAGACAAGCACATTGT,AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
""")
        expected_g2p_csv = """\
rank,count,g2p,fpr,call,seq,aligned,error
1,1,0.0677537070158,42.3,R5,CTRPNNNTRKSIHIGPGRAFYATGEIIGDIRQAHC,CTRPN-NNT--RKSIHI---GPGR---AFYAT----GEIIGDI--RQAHC,
"""

        sam_g2p(self.pssm, remap_csv, self.nuc_csv, self.g2p_csv)

        self.assertEqual(expected_g2p_csv, self.g2p_csv.getvalue())

    def testAmbiguousMixture(self):
        """ Marking position 9 as low quality means codon 3 could be S or R.
        """
        remap_csv = StringIO("""\
qname,flag,rname,pos,mapq,cigar,rnext,pnext,tlen,seq,qual
Example_read_1,99,HIV1B-env-seed,877,44,56M,=,926,56,TGTACAAGACCCAACAACAATACAAGAAAAAGTATACATATAGGACCAGGGAGAGC,AAAAAAAA#AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
Example_read_1,147,HIV1B-env-seed,926,44,56M,=,877,-56,GGAGAGCATTTTATGCAACAGGAGAAATAATAGGAGATATAAGACAAGCACATTGT,AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
""")
        expected_g2p_csv = """\
rank,count,g2p,fpr,call,seq,aligned,error
1,1,0.0663051848427,43.0,R5,CT[RS]PNNNTRKSIHIGPGRAFYATGEIIGDIRQAHC,CT[RS]PN-NNT--RKSIHI---GPGR---AFYAT----GEIIGDI--RQAHC,ambiguous
"""

        sam_g2p(self.pssm, remap_csv, self.nuc_csv, self.g2p_csv)

        self.assertEqual(expected_g2p_csv, self.g2p_csv.getvalue())

    def testAmbiguousAtTwoPositions(self):
        """ Same thing with codons 9 and 18 - rejected. """
        remap_csv = StringIO("""\
qname,flag,rname,pos,mapq,cigar,rnext,pnext,tlen,seq,qual
Example_read_1,99,HIV1B-env-seed,877,44,56M,=,926,56,TGTACAAGACCCAACAACAATACAAGAAAAAGTATACATATAGGACCAGGGAGAGC,AAAAAAAAAAAAAAAAAAAAAAAAAA#AAAAAAAAAAAAAAAAAAAAAAAAAA#AA
Example_read_1,147,HIV1B-env-seed,926,44,56M,=,877,-56,GGAGAGCATTTTATGCAACAGGAGAAATAATAGGAGATATAAGACAAGCACATTGT,AAAA#AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
""")
        expected_g2p_csv = """\
rank,count,g2p,fpr,call,seq,aligned,error
1,1,,,,CTRPNNNTXKSIHIGPGXAFYATGEIIGDIRQAHC,,> 2 ambiguous
"""

        sam_g2p(self.pssm, remap_csv, self.nuc_csv, self.g2p_csv)

        self.assertEqual(expected_g2p_csv, self.g2p_csv.getvalue())

    def testAmbiguousMixtureThreeChoices(self):
        """ Marking position 14 as low quality means codon 5 could be L, S, or *.
        """
        remap_csv = StringIO("""\
qname,flag,rname,pos,mapq,cigar,rnext,pnext,tlen,seq,qual
Example_read_1,99,HIV1B-env-seed,877,44,21M,=,877,56,TGTACAAGACCCTTAAACTGT,AAAAAAAAAAAAA#AAAAAAA
Example_read_1,147,HIV1B-env-seed,877,44,21M,=,877,56,TGTACAAGACCCTTAAACTGT,AAAAAAAAAAAAA#AAAAAAA
""")
        expected_g2p_csv = """\
rank,count,g2p,fpr,call,seq,aligned,error
1,1,,,,CTRPXNC,,> 2 ambiguous
"""

        sam_g2p(self.pssm, remap_csv, self.nuc_csv, self.g2p_csv)

        self.assertEqual(expected_g2p_csv, self.g2p_csv.getvalue())

    def testLowQuality(self):
        remap_csv = StringIO("""\
qname,flag,rname,pos,mapq,cigar,rnext,pnext,tlen,seq,qual
Example_read_1,99,HIV1B-env-seed,877,44,9M,=,877,9,TNTNNNGGN,A#A###AA#
Example_read_1,147,HIV1B-env-seed,877,44,9M,=,877,-9,TNTNNNGGN,A#A###AA#
""")
        expected_g2p_csv = """\
rank,count,g2p,fpr,call,seq,aligned,error
1,1,,,,,,low quality
"""

        sam_g2p(self.pssm, remap_csv, self.nuc_csv, self.g2p_csv)

        self.assertEqual(expected_g2p_csv, self.g2p_csv.getvalue())

    def testPartialCodon(self):
        remap_csv = StringIO("""\
qname,flag,rname,pos,mapq,cigar,rnext,pnext,tlen,seq,qual
Example_read_1,99,HIV1B-env-seed,877,44,8M,=,877,8,TGTACAGG,AAAAAAAA
Example_read_1,147,HIV1B-env-seed,877,44,8M,=,877,-8,TGTACAGG,AAAAAAAA
""")
        expected_g2p_csv = """\
rank,count,g2p,fpr,call,seq,aligned,error
1,1,,,,,,notdiv3
"""

        sam_g2p(self.pssm, remap_csv, self.nuc_csv, self.g2p_csv)

        self.assertEqual(expected_g2p_csv, self.g2p_csv.getvalue())

    def testStopCodon(self):
        remap_csv = StringIO("""\
qname,flag,rname,pos,mapq,cigar,rnext,pnext,tlen,seq,qual
Example_read_1,99,HIV1B-env-seed,877,44,9M,=,877,9,TGTTAGTGT,AAAAAAAAA
Example_read_1,147,HIV1B-env-seed,877,44,9M,=,877,-9,TGTTAGTGT,AAAAAAAAA
""")
        expected_g2p_csv = """\
rank,count,g2p,fpr,call,seq,aligned,error
1,1,,,,C*C,,stop codons
"""

        sam_g2p(self.pssm, remap_csv, self.nuc_csv, self.g2p_csv)

        self.assertEqual(expected_g2p_csv, self.g2p_csv.getvalue())

    def testAllClipped(self):
        """ In this scenario, the reads map outside the clipping region. """
        remap_csv = StringIO("""\
qname,flag,rname,pos,mapq,cigar,rnext,pnext,tlen,seq,qual
Example_read_1,99,HIV1B-env-seed,868,44,9M,=,877,9,TGTACAGGG,AAAAAAAAA
Example_read_1,147,HIV1B-env-seed,868,44,9M,=,877,-9,TGTACAGGG,AAAAAAAAA
""")
        expected_g2p_csv = """\
rank,count,g2p,fpr,call,seq,aligned,error
1,1,,,,,,zerolength
"""

        sam_g2p(self.pssm, remap_csv, self.nuc_csv, self.g2p_csv)

        self.assertEqual(expected_g2p_csv, self.g2p_csv.getvalue())

    def testLengthMinimum(self):
        remap_csv = StringIO("""\
qname,flag,rname,pos,mapq,cigar,rnext,pnext,tlen,seq,qual
Example_read_1,99,HIV1B-env-seed,877,44,51M,=,925,51,TGTGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGAAA,AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
Example_read_1,147,HIV1B-env-seed,925,44,48M,=,877,-48,AAAGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGTGT,AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
""")
        expected_g2p_csv = """\
rank,count,g2p,fpr,call,seq,aligned,error
1,1,0.806326707173,1.5,X4,CGGGGGGGGGGGGGGGKGGGGGGGGGGGGGGC,---CG-GGG--GGGGGG---GGGG---GKGGG----GGGGGGG--GGGGC,
"""

        sam_g2p(self.pssm, remap_csv, self.nuc_csv, self.g2p_csv)

        self.assertEqual(expected_g2p_csv, self.g2p_csv.getvalue())

    def testLengthTooShort(self):
        remap_csv = StringIO("""\
qname,flag,rname,pos,mapq,cigar,rnext,pnext,tlen,seq,qual
Example_read_1,99,HIV1B-env-seed,877,44,51M,=,925,51,TGTGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGAAA,AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
Example_read_1,147,HIV1B-env-seed,925,44,45M,=,877,-45,AAAGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGTGT,AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
""")
        expected_g2p_csv = """\
rank,count,g2p,fpr,call,seq,aligned,error
1,1,,,,CGGGGGGGGGGGGGGGKGGGGGGGGGGGGGC,,length
"""

        sam_g2p(self.pssm, remap_csv, self.nuc_csv, self.g2p_csv)

        self.assertEqual(expected_g2p_csv, self.g2p_csv.getvalue())

    def testDeletion(self):
        remap_csv = StringIO("""\
qname,flag,rname,pos,mapq,cigar,rnext,pnext,tlen,seq,qual
Example_read_1,99,HIV1B-env-seed,877,44,3M3D6M,=,877,9,TGTGGGTGT,AAAAAAAAA
Example_read_1,147,HIV1B-env-seed,877,44,3M3D6M,=,877,-9,TGTGGGTGT,AAAAAAAAA
""")
        expected_g2p_csv = """\
rank,count,g2p,fpr,call,seq,aligned,error
1,1,,,,CGC,,length
"""

        sam_g2p(self.pssm, remap_csv, self.nuc_csv, self.g2p_csv)

        self.assertEqual(expected_g2p_csv, self.g2p_csv.getvalue())

    def testDeletionAtStart(self):
        remap_csv = StringIO("""\
qname,flag,rname,pos,mapq,cigar,rnext,pnext,tlen,seq,qual
Example_read_1,99,HIV1B-env-seed,874,44,3M3D6M,=,874,9,TGTGGGTGT,AAAAAAAAA
Example_read_1,147,HIV1B-env-seed,874,44,3M3D6M,=,874,-9,TGTGGGTGT,AAAAAAAAA
""")
        expected_g2p_csv = """\
rank,count,g2p,fpr,call,seq,aligned,error
1,1,,,,-GC,,cysteines
"""

        sam_g2p(self.pssm, remap_csv, self.nuc_csv, self.g2p_csv)

        self.assertEqual(expected_g2p_csv, self.g2p_csv.getvalue())
