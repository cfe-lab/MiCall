import unittest
from StringIO import StringIO

from micall.core import sam2aln
from micall.core.sam2aln import apply_cigar, merge_pairs


class RemapReaderTest(unittest.TestCase):
    def test_basic(self):
        remap_file = StringIO("""\
qname,flag,rname,pos,mapq,cigar,rnext,pnext,tlen,seq,qual
Example_read_1,99,V3LOOP,1,44,32M,=,1,-32,TGTACAAGACCCAACAACAATACAAGAAAAAG,AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
Example_read_1,147,V3LOOP,1,44,32M,=,1,-32,TGTACAAGACCCAACAACAATACAAGAAAAAG,AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
Example_read_2,99,INT,1,44,32M,=,1,-32,TGTACAAGACCCAACAACAATACAAGAAAAAG,AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
Example_read_2,147,INT,1,44,32M,=,1,-32,TGTACAAGACCCAACAACAATACAAGAAAAAG,AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
""")
        expected_aligned_csv = """\
refname,qcut,rank,count,offset,seq
INT,15,0,1,0,TGTACAAGACCCAACAACAATACAAGAAAAAG
V3LOOP,15,0,1,0,TGTACAAGACCCAACAACAATACAAGAAAAAG
"""
        actual_aligned_csv = StringIO()
        sam2aln.sam2aln(remap_file, actual_aligned_csv, StringIO(), StringIO())

        self.assertMultiLineEqual(expected_aligned_csv,
                                  actual_aligned_csv.getvalue())

    def test_escaping(self):
        remap_file = StringIO("""\
qname,flag,rname,pos,mapq,cigar,rnext,pnext,tlen,seq,qual
Example_read_1,99,V3LOOP,1,44,3M,=,1,3,TGT,"A,A"
Example_read_1,147,V3LOOP,1,44,3M,=,1,-3,TGT,"A""A"
""")
        expected_aligned_csv = """\
refname,qcut,rank,count,offset,seq
V3LOOP,15,0,1,0,TNT
"""
        actual_aligned_csv = StringIO()
        sam2aln.sam2aln(remap_file, actual_aligned_csv, StringIO(), StringIO())

        self.assertMultiLineEqual(expected_aligned_csv,
                                  actual_aligned_csv.getvalue())

    def test_low_mapq(self):
        remap_file = StringIO("""\
qname,flag,rname,pos,mapq,cigar,rnext,pnext,tlen,seq,qual
Example_read_1,99,V3LOOP,1,44,32M,=,1,-32,TGTACAAGACCCAACAACAATACAAGAAAAAG,AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
Example_read_1,147,V3LOOP,1,44,32M,=,1,-32,TGTACAAGACCCAACAACAATACAAGAAAAAG,AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
Example_read_2,99,INT,1,8,32M,=,1,-32,TGTACAAGACCCAACAACAATACAAGAAAAAG,AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
Example_read_2,147,INT,1,44,32M,=,1,-32,TGTACAAGACCCAACAACAATACAAGAAAAAG,AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
""")
        expected_aligned_csv = """\
refname,qcut,rank,count,offset,seq
V3LOOP,15,0,1,0,TGTACAAGACCCAACAACAATACAAGAAAAAG
"""
        expected_failed_csv = """\
qname,cause
Example_read_2,mapq
"""
        actual_aligned_csv = StringIO()
        actual_failed_csv = StringIO()
        sam2aln.sam2aln(remap_file,
                        actual_aligned_csv,
                        StringIO(),
                        actual_failed_csv)

        self.assertMultiLineEqual(expected_aligned_csv,
                                  actual_aligned_csv.getvalue())
        self.assertMultiLineEqual(expected_failed_csv,
                                  actual_failed_csv.getvalue())

    def test_low_read_quality(self):
        remap_file = StringIO("""\
qname,flag,rname,pos,mapq,cigar,rnext,pnext,tlen,seq,qual
Example_read_1,99,V3LOOP,1,44,32M,=,1,-32,TGTACAAGACCCAACAACAATACAAGAAAAAG,000000000000000000AAAAAAAAAAAAAA
Example_read_1,147,V3LOOP,1,44,32M,=,1,-32,TGTACAAGACCCAACAACAATACAAGAAAAAG,000000000000000000AAAAAAAAAAAAAA
Example_read_2,99,INT,1,44,32M,=,1,-32,TGTACAAGACCCAACAACAATACAAGAAAAAG,AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
Example_read_2,147,INT,1,44,32M,=,1,-32,TGTACAAGACCCAACAACAATACAAGAAAAAG,AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
""")
        expected_aligned_csv = """\
refname,qcut,rank,count,offset,seq
INT,15,0,1,0,TGTACAAGACCCAACAACAATACAAGAAAAAG
"""
        expected_failed_csv = """\
qname,cause
Example_read_1,manyNs
"""
        actual_aligned_csv = StringIO()
        actual_failed_csv = StringIO()
        sam2aln.sam2aln(remap_file,
                        actual_aligned_csv,
                        StringIO(),
                        actual_failed_csv)

        self.assertMultiLineEqual(expected_aligned_csv,
                                  actual_aligned_csv.getvalue())
        self.assertMultiLineEqual(expected_failed_csv,
                                  actual_failed_csv.getvalue())

    def test_unmatched_read(self):
        remap_file = StringIO("""\
qname,flag,rname,pos,mapq,cigar,rnext,pnext,tlen,seq,qual
Example_read_1,99,V3LOOP,1,44,32M,=,1,-32,TGTACAAGACCCAACAACAATACAAGAAAAAG,AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
""")
        expected_aligned_csv = """\
refname,qcut,rank,count,offset,seq
"""
        expected_failed_csv = """\
qname,cause
Example_read_1,unmatched
"""
        actual_aligned_csv = StringIO()
        actual_failed_csv = StringIO()
        sam2aln.sam2aln(remap_file,
                        actual_aligned_csv,
                        StringIO(),
                        actual_failed_csv)

        self.assertMultiLineEqual(expected_aligned_csv,
                                  actual_aligned_csv.getvalue())
        self.assertMultiLineEqual(expected_failed_csv,
                                  actual_failed_csv.getvalue())

    def test_bad_cigar(self):
        remap_file = StringIO("""\
qname,flag,rname,pos,mapq,cigar,rnext,pnext,tlen,seq,qual
Example_read_1,99,V3LOOP,1,44,32M,=,1,-32,TGTACAAGACCCAACAACAATACAAGAAAAAG,AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
Example_read_1,147,V3LOOP,1,44,*,=,1,-32,TGTACAAGACCCAACAACAATACAAGAAAAAG,AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
""")
        expected_aligned_csv = """\
refname,qcut,rank,count,offset,seq
"""
        expected_failed_csv = """\
qname,cause
Example_read_1,badCigar
"""
        actual_aligned_csv = StringIO()
        actual_failed_csv = StringIO()
        sam2aln.sam2aln(remap_file,
                        actual_aligned_csv,
                        StringIO(),
                        actual_failed_csv)

        self.assertMultiLineEqual(expected_aligned_csv,
                                  actual_aligned_csv.getvalue())
        self.assertMultiLineEqual(expected_failed_csv,
                                  actual_failed_csv.getvalue())

    def test_different_references(self):
        remap_file = StringIO("""\
qname,flag,rname,pos,mapq,cigar,rnext,pnext,tlen,seq,qual
Example_read_1,99,V3LOOP,1,44,32M,=,1,-32,TGTACAAGACCCAACAACAATACAAGAAAAAG,AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
Example_read_1,147,GP41,1,44,32M,=,1,-32,TGTACAAGACCCAACAACAATACAAGAAAAAG,AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
""")
        expected_aligned_csv = """\
refname,qcut,rank,count,offset,seq
"""
        expected_failed_csv = """\
qname,cause
Example_read_1,2refs
"""
        actual_aligned_csv = StringIO()
        actual_failed_csv = StringIO()
        sam2aln.sam2aln(remap_file,
                        actual_aligned_csv,
                        StringIO(),
                        actual_failed_csv)

        self.assertMultiLineEqual(expected_aligned_csv,
                                  actual_aligned_csv.getvalue())
        self.assertMultiLineEqual(expected_failed_csv,
                                  actual_failed_csv.getvalue())

    def test_insertion(self):
        remap_file = StringIO("""\
qname,flag,rname,pos,mapq,cigar,rnext,pnext,tlen,seq,qual
Example_read_1,99,V3LOOP,1,44,12M6I14M,=,1,-32,TGTACAAGACCCAACAACAATACAAGAAAAAG,AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
Example_read_1,147,V3LOOP,1,44,12M6I14M,=,1,-32,TGTACAAGACCCAACAACAATACAAGAAAAAG,AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
""")
        expected_aligned_csv = """\
refname,qcut,rank,count,offset,seq
V3LOOP,15,0,1,0,TGTACAAGACCCAATACAAGAAAAAG
"""
        expected_insert_csv = """\
qname,fwd_rev,refname,pos,insert,qual
Example_read_1,F,V3LOOP,12,AACAAC,AAAAAA
Example_read_1,R,V3LOOP,12,AACAAC,AAAAAA
"""
        actual_aligned_csv = StringIO()
        actual_insert_csv = StringIO()
        sam2aln.sam2aln(remap_file,
                        actual_aligned_csv,
                        actual_insert_csv,
                        StringIO())

        self.assertMultiLineEqual(expected_aligned_csv,
                                  actual_aligned_csv.getvalue())
        self.assertMultiLineEqual(expected_insert_csv,
                                  actual_insert_csv.getvalue())


class CigarTest(unittest.TestCase):
    def setUp(self):
        super(CigarTest, self).setUp()
        self.addTypeEqualityFunc(str, self.assertMultiLineEqual)

    def testTrivial(self):
        cigar = '9M'
        seq     = 'AAACAACCA'  # @IgnorePep8
        quality = 'BBBBBBBBB'
        expected_seq = seq
        expected_quality = quality

        clipped_seq, clipped_quality, inserts = apply_cigar(cigar, seq, quality)

        self.assertEqual(expected_seq, clipped_seq)
        self.assertEqual(expected_quality, clipped_quality)
        self.assertEqual({}, inserts)

    def testDeletion(self):
        cigar = '6M3D3M'
        seq              = 'AAACAACCA'  # @IgnorePep8
        quality          = 'BBBDDDEEE'  # @IgnorePep8
        expected_seq     = 'AAACAA---CCA'  # @IgnorePep8
        expected_quality = 'BBBDDD   EEE'

        clipped_seq, clipped_quality, inserts = apply_cigar(cigar, seq, quality)

        self.assertEqual(expected_seq, clipped_seq)
        self.assertEqual(expected_quality, clipped_quality)
        self.assertEqual({}, inserts)

    def testSoftClip(self):
        cigar = '3S6M'
        seq              = 'AAACAACCA'  # @IgnorePep8
        quality          = 'BBBDDDEEE'  # @IgnorePep8
        expected_seq     =    'CAACCA'  # @IgnorePep8
        expected_quality =    'DDDEEE'  # @IgnorePep8

        clipped_seq, clipped_quality, inserts = apply_cigar(cigar, seq, quality)

        self.assertEqual(expected_seq, clipped_seq)
        self.assertEqual(expected_quality, clipped_quality)
        self.assertEqual({}, inserts)

    def testInsertion(self):
        cigar = '3M3I6M'
        seq              = 'AAACAACCACCC'  # @IgnorePep8
        quality          = 'BBBDDDEEEFFF'  # @IgnorePep8
        expected_seq     = 'AAACCACCC'  # @IgnorePep8
        expected_quality = 'BBBEEEFFF'

        clipped_seq, clipped_quality, inserts = apply_cigar(cigar, seq, quality)

        self.assertEqual(expected_seq, clipped_seq)
        self.assertEqual(expected_quality, clipped_quality)
        self.assertEqual({3: ('CAA', 'DDD')}, inserts)

    def testInsertionLowQuality(self):
        cigar = '3M3I6M'
        seq              = 'AAACAACCACCC'  # @IgnorePep8
        quality          = 'BBBD*DEEEFFF'  # @IgnorePep8
        expected_seq     = 'AAACCACCC'  # @IgnorePep8
        expected_quality = 'BBBEEEFFF'

        clipped_seq, clipped_quality, inserts = apply_cigar(cigar, seq, quality)

        self.assertEqual(expected_seq, clipped_seq)
        self.assertEqual(expected_quality, clipped_quality)
        self.assertEqual({3: ('CAA', 'D*D')}, inserts)

    def testLargeToken(self):
        cigar = '12M'
        seq     = 'AAACAACCACCC'  # @IgnorePep8
        quality = 'BBBBBBBBBBBB'
        expected_seq = seq
        expected_quality = quality

        clipped_seq, clipped_quality, inserts = apply_cigar(cigar, seq, quality)

        self.assertEqual(expected_seq, clipped_seq)
        self.assertEqual(expected_quality, clipped_quality)
        self.assertEqual({}, inserts)

    def testPadding(self):
        cigar = '12M'
        seq              = 'AAACAACCACCC'  # @IgnorePep8
        quality          = 'BBBDDDEEEFFF'  # @IgnorePep8
        pos = 3
        expected_seq     = '---AAACAACCACCC'  # @IgnorePep8
        expected_quality = '!!!BBBDDDEEEFFF'

        clipped_seq, clipped_quality, inserts = apply_cigar(cigar, seq, quality, pos)

        self.assertEqual(expected_seq, clipped_seq)
        self.assertEqual(expected_quality, clipped_quality)
        self.assertEqual({}, inserts)

    def testClipping(self):
        cigar = '12M'
        seq              = 'AAACAACCACCC'  # @IgnorePep8
        quality          = 'BBBDDDEEEFFF'  # @IgnorePep8
        pos = 0
        clip_from = 3
        clip_to = 8
        expected_seq     = 'CAACCA'  # @IgnorePep8
        expected_quality = 'DDDEEE'

        clipped_seq, clipped_quality, inserts = apply_cigar(
          cigar,
          seq,
          quality,
          pos,
          clip_from,
          clip_to)

        self.assertEqual(expected_seq, clipped_seq)
        self.assertEqual(expected_quality, clipped_quality)
        self.assertEqual({}, inserts)

    def testClipInsertion(self):
        cigar = '6M3I6M'
        seq              = 'AAACAAGGGCCACCC'  # @IgnorePep8
        quality          = 'BBBDDDHHHEEEFFF'  # @IgnorePep8
        pos = 0
        clip_from = 3
        clip_to = 8
        expected_seq     = 'CAACCA'  # @IgnorePep8
        expected_quality = 'DDDEEE'  # @IgnorePep8

        clipped_seq, clipped_quality, inserts = apply_cigar(
          cigar,
          seq,
          quality,
          pos,
          clip_from,
          clip_to)

        self.assertEqual(expected_seq, clipped_seq)
        self.assertEqual(expected_quality, clipped_quality)
        self.assertEqual({6: ('GGG', 'HHH')}, inserts)

    def testClipInsertionLowQuality(self):
        cigar = '6M3I6M'
        seq              = 'AAACAAGGGCCACCC'  # @IgnorePep8
        quality          = 'BBBDDDHH*EEEFFF'  # @IgnorePep8
        pos = 0
        clip_from = 3
        clip_to = 8
        expected_seq     = 'CAACCA'  # @IgnorePep8
        expected_quality = 'DDDEEE'

        clipped_seq, clipped_quality, inserts = apply_cigar(
          cigar,
          seq,
          quality,
          pos,
          clip_from,
          clip_to)

        self.assertEqual(expected_seq, clipped_seq)
        self.assertEqual(expected_quality, clipped_quality)
        self.assertEqual({6: ('GGG', 'HH*')}, inserts)

    def testInsertionBeforeClip(self):
        cigar = '3M3I9M'
        seq              = 'AAAGGGCAACCACCC'  # @IgnorePep8
        quality          = 'BBBHHHDDDEEEFFF'  # @IgnorePep8
        pos = 0
        clip_from = 3
        clip_to = 8
        expected_seq     = 'CAACCA'  # @IgnorePep8
        expected_quality = 'DDDEEE'

        clipped_seq, clipped_quality, inserts = apply_cigar(
          cigar,
          seq,
          quality,
          pos,
          clip_from,
          clip_to)

        self.assertEqual(expected_seq, clipped_seq)
        self.assertEqual(expected_quality, clipped_quality)
        self.assertEqual({3: ('GGG', 'HHH')}, inserts)

    def testInsertionAfterClipWithOffset(self):
        cigar = '2M1I2M'
        seq     = 'TAGCT'  # @IgnorePep8
        quality = 'AABCC'
        pos = 3
        clip_from = 4
        clip_to = 20
        expected_seq     = 'ACT'  # @IgnorePep8
        expected_quality = 'ACC'

        clipped_seq, clipped_quality, inserts = apply_cigar(
          cigar,
          seq,
          quality,
          pos,
          clip_from,
          clip_to)

        self.assertEqual(expected_seq, clipped_seq)
        self.assertEqual(expected_quality, clipped_quality)
        self.assertEqual({5: ('G', 'B')}, inserts)

    def testClippingEverything(self):
        cigar = '12M'
        seq              = 'AAACAACCACCC'  # @IgnorePep8
        quality          = 'BBBDDDEEEFFF'  # @IgnorePep8
        pos = 0
        clip_from = 100
        clip_to = 108
        expected_seq     = ''  # @IgnorePep8
        expected_quality = ''

        clipped_seq, clipped_quality, inserts = apply_cigar(
          cigar,
          seq,
          quality,
          pos,
          clip_from,
          clip_to)

        self.assertEqual(expected_seq, clipped_seq)
        self.assertEqual(expected_quality, clipped_quality)
        self.assertEqual({}, inserts)

    def testInvalidCigar(self):
        cigar = '3M...6M'
        seq     = 'AAACAACCACCC'  # @IgnorePep8
        quality = 'BBBDDDEEEFFF'

        with self.assertRaises(RuntimeError) as result:
            apply_cigar(cigar, seq, quality)

        self.assertEqual(
            "Invalid CIGAR string: '3M...6M'.",
            result.exception.message)

    def testUnsupportedCigarToken(self):
        cigar = '3M3X6M'
        seq     = 'AAACAACCACCC'  # @IgnorePep8
        quality = 'BBBDDDEEEFFF'

        with self.assertRaises(RuntimeError) as result:
            apply_cigar(cigar, seq, quality)

        self.assertEqual(
            "Unsupported CIGAR token: '3X'.",
            result.exception.message)

    def testShortCigar(self):
        cigar = '8M'
        seq     = 'AAACAACCA'  # @IgnorePep8
        quality = 'BBBDDDEEE'

        with self.assertRaises(RuntimeError) as result:
            apply_cigar(cigar, seq, quality)

        self.assertEqual(
            "CIGAR string '8M' is too short for sequence 'AAACAACCA'.",
            result.exception.message)

    def testLongCigar(self):
        cigar = '10M'
        seq     = 'AAACAACCA'  # @IgnorePep8
        quality = 'BBBDDDEEE'

        with self.assertRaises(RuntimeError) as result:
            apply_cigar(cigar, seq, quality)

        self.assertEqual(
            "CIGAR string '10M' is too long for sequence 'AAACAACCA'.",
            result.exception.message)

    def testInsertionAfterClipping(self):
        cigar = '3M3I3M'
        seq     = "ACTTAGAAA"  # @IgnorePep8
        quality = 'AAABBBDDD'
        pos = 0
        clip_from = 0
        clip_to = 2
        expected_seq     = 'ACT'  # @IgnorePep8
        expected_quality = 'AAA'

        clipped_seq, clipped_quality, inserts = apply_cigar(
          cigar,
          seq,
          quality,
          pos,
          clip_from,
          clip_to)

        self.assertEqual(expected_seq, clipped_seq)
        self.assertEqual(expected_quality, clipped_quality)
        self.assertEqual({}, inserts)

    def testInsertionAtEndOfClipping(self):
        cigar = '3M3I3M'
        seq     = "ACTTAGAAA"  # @IgnorePep8
        quality = 'AAABBBDDD'
        pos = 0
        clip_from = 0
        clip_to = 3
        expected_seq     = 'ACTA'  # @IgnorePep8
        expected_quality = 'AAAD'

        clipped_seq, clipped_quality, inserts = apply_cigar(
          cigar,
          seq,
          quality,
          pos,
          clip_from,
          clip_to)

        self.assertEqual(expected_seq, clipped_seq)
        self.assertEqual(expected_quality, clipped_quality)
        self.assertEqual({3: ('TAG', 'BBB')}, inserts)


class MergePairsTest(unittest.TestCase):
    def setUp(self):
        self.addTypeEqualityFunc(str, self.assertMultiLineEqual)

    def testSimple(self):
        seq1          = 'ACTGCA'  # @IgnorePep8
        seq2          = 'ACTGCA'  # @IgnorePep8
        qual1         = 'JJJJJJ'  # @IgnorePep8
        qual2         = 'JJJJJJ'  # @IgnorePep8
        expected_mseq = 'ACTGCA'

        mseq = merge_pairs(seq1, seq2, qual1, qual2)

        self.assertEqual(expected_mseq, mseq)

    def testDifferentLength(self):
        seq1          = 'ACTGCATCT'  # @IgnorePep8
        seq2          = 'ACTGCA'  # @IgnorePep8
        qual1         = 'JJJJJJJJJ'  # @IgnorePep8
        qual2         = 'JJJJJJ'  # @IgnorePep8
        expected_mseq = 'ACTGCATCT'

        mseq = merge_pairs(seq1, seq2, qual1, qual2)

        self.assertEqual(expected_mseq, mseq)

    def testOffset(self):
        seq1          = '-CTGCA'  # @IgnorePep8
        seq2          = '---GCATCT'  # @IgnorePep8
        qual1         = '!JJJJJ'  # @IgnorePep8
        qual2         = '!!!JJJJJJ'  # @IgnorePep8
        expected_mseq = '-CTGCATCT'

        mseq = merge_pairs(seq1, seq2, qual1, qual2)

        self.assertEqual(expected_mseq, mseq)

    def testDisagreementWithDifferentQuality(self):
        seq1          = 'AGTGCA'  # @IgnorePep8
        seq2          = 'ACTGCA'  # @IgnorePep8
        qual1         = 'JAJJJJ'  # @IgnorePep8
        qual2         = 'JJJJJJ'  # @IgnorePep8
        expected_mseq = 'ACTGCA'

        mseq = merge_pairs(seq1, seq2, qual1, qual2)

        self.assertEqual(expected_mseq, mseq)

    def testDisagreementWithCloseQuality(self):
        seq1          = 'AGTGCA'  # @IgnorePep8
        seq2          = 'ACTGCA'  # @IgnorePep8
        qual1         = 'JHJJJJ'  # @IgnorePep8
        qual2         = 'JJJJJJ'  # @IgnorePep8
        expected_mseq = 'ANTGCA'

        mseq = merge_pairs(seq1, seq2, qual1, qual2)

        self.assertEqual(expected_mseq, mseq)

    def testDisagreementWithLowQuality(self):
        seq1          = 'AGTGCA'  # @IgnorePep8
        seq2          = 'ACTGCA'  # @IgnorePep8
        qual1         = 'J!JJJJ'  # @IgnorePep8
        qual2         = 'J*JJJJ'  # @IgnorePep8
        expected_mseq = 'ANTGCA'

        mseq = merge_pairs(seq1, seq2, qual1, qual2)

        self.assertEqual(expected_mseq, mseq)

    def testGap(self):
        seq1          = 'AGT'  # @IgnorePep8
        seq2          = '------GCA'  # @IgnorePep8
        qual1         = 'JJJ'  # @IgnorePep8
        qual2         = '!!!!!!JJJ'  # @IgnorePep8
        expected_mseq = 'AGTnnnGCA'

        mseq = merge_pairs(seq1, seq2, qual1, qual2)

        self.assertEqual(expected_mseq, mseq)

    def testLowQualityInSecondRead(self):
        seq1          = 'AGT'  # @IgnorePep8
        seq2          = '---GCA'  # @IgnorePep8
        qual1         = 'JJJ'  # @IgnorePep8
        qual2         = '!!!J*J'  # @IgnorePep8
        expected_mseq = 'AGTGNA'

        mseq = merge_pairs(seq1, seq2, qual1, qual2)

        self.assertEqual(expected_mseq, mseq)

    def testOneInsertion(self):
        seq1          = 'AGT'  # @IgnorePep8
        seq2          = '---GCA'  # @IgnorePep8
        qual1         = 'JJJ'  # @IgnorePep8
        qual2         = '!!!JJJ'  # @IgnorePep8
        ins1 = {2: ('CCC', 'JJJ')}
        ins2 = {}
        expected_mseq = 'AGCCCTGCA'

        mseq = merge_pairs(seq1, seq2, qual1, qual2, ins1, ins2)

        self.assertEqual(expected_mseq, mseq)

    def testTwoInsertions(self):
        seq1          = 'AGT'  # @IgnorePep8
        seq2          = '---GCA'  # @IgnorePep8
        qual1         = 'JJJ'  # @IgnorePep8
        qual2         = '!!!JJJ'  # @IgnorePep8
        ins1 = {2: ('CCC', 'JJJ')}
        ins2 = {5: ('TTT', 'JJJ')}
        expected_mseq = 'AGCCCTGCTTTA'

        mseq = merge_pairs(seq1, seq2, qual1, qual2, ins1, ins2)

        self.assertEqual(expected_mseq, mseq)

    def testConflictingInsertions(self):
        seq1          = 'AGTGCA'  # @IgnorePep8
        seq2          = 'AGTGCA'  # @IgnorePep8
        qual1         = 'JJJJJJ'  # @IgnorePep8
        qual2         = 'JJJJJJ'  # @IgnorePep8
        ins1 = {2: ('CCC', 'JJJ')}
        ins2 = {2: ('CTC', 'JAJ')}
        expected_mseq = 'AGCCCTGCA'

        mseq = merge_pairs(seq1, seq2, qual1, qual2, ins1, ins2)

        self.assertEqual(expected_mseq, mseq)
