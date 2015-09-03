import unittest
from StringIO import StringIO
from micall.g2p.sam_g2p import sam_g2p, apply_cigar_and_clip
from micall.g2p.pssm_lib import Pssm

class SamG2PTest(unittest.TestCase):
    def setUp(self):
        super(SamG2PTest, self).setUp()
        self.pssm = Pssm(path_to_lookup='../g2p/g2p_fpr.txt',
                         path_to_matrix='../g2p/g2p.matrix')
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
        

class CigarTest(unittest.TestCase):
    def setUp(self):
        super(CigarTest, self).setUp()
        self.addTypeEqualityFunc(str, self.assertMultiLineEqual)
        
    def testTrivial(self):
        cigar = '9M'
        seq     = 'AAACAACCA'
        quality = 'BBBBBBBBB'
        expected_seq = seq
        expected_quality = quality
        
        clipped_seq, clipped_quality = apply_cigar_and_clip(cigar, seq, quality)
        
        self.assertEqual(expected_seq, clipped_seq)
        self.assertEqual(expected_quality, clipped_quality)
     
    def testDeletion(self):
        cigar = '6M3D3M'
        seq              = 'AAACAACCA'
        quality          = 'BBBDDDEEE'
        expected_seq     = 'AAACAA---CCA'
        expected_quality = 'BBBDDD   EEE'
           
        clipped_seq, clipped_quality = apply_cigar_and_clip(cigar, seq, quality)
           
        self.assertEqual(expected_seq, clipped_seq)
        self.assertEqual(expected_quality, clipped_quality)
 
    def testSoftClip(self):
        cigar = '3S6M'
        seq              = 'AAACAACCA'
        quality          = 'BBBDDDEEE'
        expected_seq     =    'CAACCA'
        expected_quality =    'DDDEEE'
           
        clipped_seq, clipped_quality = apply_cigar_and_clip(cigar, seq, quality)
           
        self.assertEqual(expected_seq, clipped_seq)
        self.assertEqual(expected_quality, clipped_quality)
     
    def testInsertion(self):
        cigar = '3M3I6M'
        seq              = 'AAACAACCACCC'
        quality          = 'BBBDDDEEEFFF'
        expected_seq     = 'AAACAACCACCC'
        expected_quality = 'BBBDDDEEEFFF'
           
        clipped_seq, clipped_quality = apply_cigar_and_clip(cigar, seq, quality)
           
        self.assertEqual(expected_seq, clipped_seq)
        self.assertEqual(expected_quality, clipped_quality)
 
    def testInsertionLowQuality(self):
        cigar = '3M3I6M'
        seq              = 'AAACAACCACCC'
        quality          = 'BBBD*DEEEFFF'
        expected_seq     = 'AAACCACCC'
        expected_quality = 'BBBEEEFFF'
           
        clipped_seq, clipped_quality = apply_cigar_and_clip(cigar, seq, quality)
           
        self.assertEqual(expected_seq, clipped_seq)
        self.assertEqual(expected_quality, clipped_quality)
       
    def testLargeToken(self):
        cigar = '12M'
        seq     = 'AAACAACCACCC'
        quality = 'BBBBBBBBBBBB'
        expected_seq = seq
        expected_quality = quality
           
        clipped_seq, clipped_quality = apply_cigar_and_clip(cigar, seq, quality)
           
        self.assertEqual(expected_seq, clipped_seq)
        self.assertEqual(expected_quality, clipped_quality)
       
    def testPadding(self):
        cigar = '12M'
        seq              = 'AAACAACCACCC'
        quality          = 'BBBDDDEEEFFF'
        pos = 3
        expected_seq     = '---AAACAACCACCC'
        expected_quality = '!!!BBBDDDEEEFFF'
           
        clipped_seq, clipped_quality = apply_cigar_and_clip(cigar, seq, quality, pos)
           
        self.assertEqual(expected_seq, clipped_seq)
        self.assertEqual(expected_quality, clipped_quality)
       
    def testClipping(self):
        cigar = '12M'
        seq              = 'AAACAACCACCC'
        quality          = 'BBBDDDEEEFFF'
        pos = 0
        clip_from = 3
        clip_to = 8
        expected_seq     = 'CAACCA'
        expected_quality = 'DDDEEE'
           
        clipped_seq, clipped_quality = apply_cigar_and_clip(
          cigar,
          seq,
          quality,
          pos,
          clip_from,
          clip_to)
           
        self.assertEqual(expected_seq, clipped_seq)
        self.assertEqual(expected_quality, clipped_quality)
   
    def testClipInsertion(self):
        cigar = '6M3I6M'
        seq              = 'AAACAAGGGCCACCC'
        quality          = 'BBBDDDHHHEEEFFF'
        pos = 0
        clip_from = 3
        clip_to = 8
        expected_seq     = 'CAAGGGCCA'
        expected_quality = 'DDDHHHEEE'
           
        clipped_seq, clipped_quality = apply_cigar_and_clip(
          cigar,
          seq,
          quality,
          pos,
          clip_from,
          clip_to)
           
        self.assertEqual(expected_seq, clipped_seq)
        self.assertEqual(expected_quality, clipped_quality)
       
    def testClipInsertionLowQuality(self):
        cigar = '6M3I6M'
        seq              = 'AAACAAGGGCCACCC'
        quality          = 'BBBDDDHH*EEEFFF'
        pos = 0
        clip_from = 3
        clip_to = 8
        expected_seq     = 'CAACCA'
        expected_quality = 'DDDEEE'
           
        clipped_seq, clipped_quality = apply_cigar_and_clip(
          cigar,
          seq,
          quality,
          pos,
          clip_from,
          clip_to)
           
        self.assertEqual(expected_seq, clipped_seq)
        self.assertEqual(expected_quality, clipped_quality)
       
    def testInsertionBeforeClip(self):
        cigar = '3M3I9M'
        seq              = 'AAAGGGCAACCACCC'
        quality          = 'BBBHHHDDDEEEFFF'
        pos = 0
        clip_from = 3
        clip_to = 8
        expected_seq     = 'CAACCA'
        expected_quality = 'DDDEEE'
           
        clipped_seq, clipped_quality = apply_cigar_and_clip(
          cigar,
          seq,
          quality,
          pos,
          clip_from,
          clip_to)
           
        self.assertEqual(expected_seq, clipped_seq)
        self.assertEqual(expected_quality, clipped_quality)
    
    def testInsertionAfterClipWithOffset(self):
        cigar = '2M1I2M'
        seq     = 'TAGCT'
        quality = 'AABCC'
        pos = 3
        clip_from = 4
        clip_to = 20
        expected_seq     = 'AGCT'
        expected_quality = 'ABCC'
           
        clipped_seq, clipped_quality = apply_cigar_and_clip(
          cigar,
          seq,
          quality,
          pos,
          clip_from,
          clip_to)
           
        self.assertEqual(expected_seq, clipped_seq)
        self.assertEqual(expected_quality, clipped_quality)
       
    def testClippingEverything(self):
        cigar = '12M'
        seq              = 'AAACAACCACCC'
        quality          = 'BBBDDDEEEFFF'
        pos = 0
        clip_from = 100
        clip_to = 108
        expected_seq     = ''
        expected_quality = ''
           
        clipped_seq, clipped_quality = apply_cigar_and_clip(
          cigar,
          seq,
          quality,
          pos,
          clip_from,
          clip_to)
           
        self.assertEqual(expected_seq, clipped_seq)
        self.assertEqual(expected_quality, clipped_quality)
   
    def testInvalidCigar(self):
        cigar = '3M...6M'
        seq     = 'AAACAACCACCC'
        quality = 'BBBDDDEEEFFF'
         
        with self.assertRaises(RuntimeError) as result:
            apply_cigar_and_clip(cigar, seq, quality)
        
        self.assertEqual(
            "Invalid CIGAR string: '3M...6M'.",
            result.exception.message)
    
    def testUnsupportedCigarToken(self):
        cigar = '3M3X6M'
        seq     = 'AAACAACCACCC'
        quality = 'BBBDDDEEEFFF'
         
        with self.assertRaises(RuntimeError) as result:
            apply_cigar_and_clip(cigar, seq, quality)
        
        self.assertEqual(
            "Unsupported CIGAR token: '3X'.",
            result.exception.message)
    
    def testShortCigar(self):
        cigar = '8M'
        seq     = 'AAACAACCA'
        quality = 'BBBDDDEEE'
         
        with self.assertRaises(RuntimeError) as result:
            apply_cigar_and_clip(cigar, seq, quality)
        
        self.assertEqual(
            "CIGAR string '8M' is too short for sequence 'AAACAACCA'.",
            result.exception.message)
    
    def testLongCigar(self):
        cigar = '10M'
        seq     = 'AAACAACCA'
        quality = 'BBBDDDEEE'
         
        with self.assertRaises(RuntimeError) as result:
            apply_cigar_and_clip(cigar, seq, quality)
        
        self.assertEqual(
            "CIGAR string '10M' is too long for sequence 'AAACAACCA'.",
            result.exception.message)
        
    def testInsertionAfterClipping(self):
        cigar = '3M3I3M'
        seq     = "ACTTAGAAA"
        quality = 'AAABBBDDD'
        pos = 0
        clip_from = 0
        clip_to = 2
        expected_seq     = 'ACT'
        expected_quality = 'AAA'
           
        clipped_seq, clipped_quality = apply_cigar_and_clip(
          cigar,
          seq,
          quality,
          pos,
          clip_from,
          clip_to)
           
        self.assertEqual(expected_seq, clipped_seq)
        self.assertEqual(expected_quality, clipped_quality)
        
    def testInsertionAtEndOfClipping(self):
        cigar = '3M3I3M'
        seq     = "ACTTAGAAA"
        quality = 'AAABBBDDD'
        pos = 0
        clip_from = 0
        clip_to = 3
        expected_seq     = 'ACTTAGA'
        expected_quality = 'AAABBBD'
           
        clipped_seq, clipped_quality = apply_cigar_and_clip(
          cigar,
          seq,
          quality,
          pos,
          clip_from,
          clip_to)
           
        self.assertEqual(expected_seq, clipped_seq)
        self.assertEqual(expected_quality, clipped_quality)
