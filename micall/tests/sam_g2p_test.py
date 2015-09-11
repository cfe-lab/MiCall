import os
from StringIO import StringIO
import unittest

from micall.g2p.pssm_lib import Pssm
from micall.g2p.sam_g2p import sam_g2p

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
        self.g2p_csv = StringIO()
        self.addTypeEqualityFunc(str, self.assertMultiLineEqual)
        
    def testSimple(self):
        remap_csv = StringIO("""\
qname,flag,rname,pos,mapq,cigar,rnext,pnext,tlen,seq,qual
Example_read_1,99,HIV1B-env-seed,877,44,9M,=,877,9,TGTACAAGA,AAAAAAAAA
Example_read_1,147,HIV1B-env-seed,877,44,9M,=,877,-9,TGTACAAGA,AAAAAAAAA
""")
        expected_g2p_csv = """\
rank,count,g2p,fpr,aligned,error
1,1,,,CTR,cysteines
"""
        
        sam_g2p(self.pssm, remap_csv, self.nuc_csv, self.g2p_csv)

        self.assertEqual(expected_g2p_csv, self.g2p_csv.getvalue())
        
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
rank,count,g2p,fpr,aligned,error
1,2,,,CTR,cysteines
2,1,,,CTG,cysteines
"""
        
        sam_g2p(self.pssm, remap_csv, self.nuc_csv, self.g2p_csv)

        self.assertEqual(expected_g2p_csv, self.g2p_csv.getvalue())
        
    def testOverlap(self):
        remap_csv = StringIO("""\
qname,flag,rname,pos,mapq,cigar,rnext,pnext,tlen,seq,qual
Example_read_1,99,HIV1B-env-seed,877,44,12M,=,886,12,TGTACAAGACCC,AAAAAAAAAAAA
Example_read_1,147,HIV1B-env-seed,886,44,9M,=,877,-9,CCCAACAAC,AAAAAAAAA
""")
        expected_g2p_csv = """\
rank,count,g2p,fpr,aligned,error
1,1,,,CTRPNN,cysteines
"""
        
        sam_g2p(self.pssm, remap_csv, self.nuc_csv, self.g2p_csv)

        self.assertEqual(expected_g2p_csv, self.g2p_csv.getvalue())
        
    def testSynonymMixture(self):
        remap_csv = StringIO("""\
qname,flag,rname,pos,mapq,cigar,rnext,pnext,tlen,seq,qual
Example_read_1,99,HIV1B-env-seed,877,44,12M,=,877,12,TGTACAGGNTGT,AAAAAAAA#AAA
Example_read_1,147,HIV1B-env-seed,877,44,12M,=,877,-12,TGTACAGGNTGT,AAAAAAAA#AAA
""")
        expected_g2p_csv = """\
rank,count,g2p,fpr,aligned,error
1,1,,,CTXC,length
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
rank,count,g2p,fpr,aligned,error
1,1,,,,low quality
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
rank,count,g2p,fpr,aligned,error
1,1,,,,notdiv3
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
rank,count,g2p,fpr,aligned,error
1,1,,,C*C,stop codons
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
rank,count,g2p,fpr,aligned,error
1,1,,,,zerolength
"""

        sam_g2p(self.pssm, remap_csv, self.nuc_csv, self.g2p_csv)

        self.assertEqual(expected_g2p_csv, self.g2p_csv.getvalue())
        
    def testLength(self):
        remap_csv = StringIO("""\
qname,flag,rname,pos,mapq,cigar,rnext,pnext,tlen,seq,qual
Example_read_1,99,HIV1B-env-seed,877,44,9M,=,877,9,TGTGGGTGT,AAAAAAAAA
Example_read_1,147,HIV1B-env-seed,877,44,9M,=,877,-9,TGTGGGTGT,AAAAAAAAA
""")
        expected_g2p_csv = """\
rank,count,g2p,fpr,aligned,error
1,1,,,CGC,length
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
rank,count,g2p,fpr,aligned,error
1,1,,,CGC,length
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
rank,count,g2p,fpr,aligned,error
1,1,,,-GC,cysteines
"""

        sam_g2p(self.pssm, remap_csv, self.nuc_csv, self.g2p_csv)

        self.assertEqual(expected_g2p_csv, self.g2p_csv.getvalue())
