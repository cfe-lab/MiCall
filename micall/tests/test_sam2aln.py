import unittest
from io import StringIO

from micall.core.sam2aln import sam2aln, apply_cigar, merge_pairs, merge_inserts


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
        sam2aln(remap_file, actual_aligned_csv)

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
        sam2aln(remap_file, actual_aligned_csv)

        self.assertMultiLineEqual(expected_aligned_csv,
                                  actual_aligned_csv.getvalue())

    def test_low_mapq(self):
        """ We no longer fail reads because of low mapq.

        When we use more than one reference, reads can receive low mapq if they
        are in a conserved region that matches more than one reference.
        """
        remap_file = StringIO("""\
qname,flag,rname,pos,mapq,cigar,rnext,pnext,tlen,seq,qual
Example_read_1,99,V3LOOP,1,44,32M,=,1,-32,TGTACAAGACCCAACAACAATACAAGAAAAAG,AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
Example_read_1,147,V3LOOP,1,44,32M,=,1,-32,TGTACAAGACCCAACAACAATACAAGAAAAAG,AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
Example_read_2,99,INT,1,8,32M,=,1,-32,TGTACAAGACCCAACAACAATACAAGAAAAAG,AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
Example_read_2,147,INT,1,44,32M,=,1,-32,TGTACAAGACCCAACAACAATACAAGAAAAAG,AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
""")
        expected_aligned_csv = """\
refname,qcut,rank,count,offset,seq
INT,15,0,1,0,TGTACAAGACCCAACAACAATACAAGAAAAAG
V3LOOP,15,0,1,0,TGTACAAGACCCAACAACAATACAAGAAAAAG
"""
        expected_failed_csv = """\
qname,cause
"""
        actual_aligned_csv = StringIO()
        actual_failed_csv = StringIO()
        sam2aln(remap_file,
                actual_aligned_csv,
                failed_csv=actual_failed_csv)

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
        sam2aln(remap_file,
                actual_aligned_csv,
                failed_csv=actual_failed_csv)

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
        sam2aln(remap_file,
                actual_aligned_csv,
                failed_csv=actual_failed_csv)

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
        sam2aln(remap_file,
                actual_aligned_csv,
                failed_csv=actual_failed_csv)

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
        sam2aln(remap_file,
                actual_aligned_csv,
                failed_csv=actual_failed_csv)

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
        sam2aln(remap_file,
                actual_aligned_csv,
                actual_insert_csv)

        self.assertMultiLineEqual(expected_aligned_csv,
                                  actual_aligned_csv.getvalue())
        self.assertMultiLineEqual(expected_insert_csv,
                                  actual_insert_csv.getvalue())

    def test_soft_clipping(self):
        remap_file = StringIO("""\
qname,flag,rname,pos,mapq,cigar,rnext,pnext,tlen,seq,qual
Example_read_1,99,V3LOOP,18,44,7S10M9S,=,1,-32,TGTACAAGACCCAATACAAGAAAAAG,AAAAAAAAAAAAAAAAAAAAAAAAAA
Example_read_1,147,V3LOOP,18,44,3S10M13S,=,1,-32,CAAGACCCAATACAAGAAAAAGCAAC,AAAAAAAAAAAAAAAAAAAAAAAAAA
""")
        expected_aligned_csv = """\
refname,qcut,rank,count,offset,seq
V3LOOP,15,0,1,17,GACCCAATAC
"""
        expected_clipping_csv = """\
refname,pos,count
V3LOOP,11,1
V3LOOP,12,1
V3LOOP,13,1
V3LOOP,14,1
V3LOOP,15,1
V3LOOP,16,1
V3LOOP,17,1
V3LOOP,28,1
V3LOOP,29,1
V3LOOP,30,1
V3LOOP,31,1
V3LOOP,32,1
V3LOOP,33,1
V3LOOP,34,1
V3LOOP,35,1
V3LOOP,36,1
V3LOOP,37,1
V3LOOP,38,1
V3LOOP,39,1
V3LOOP,40,1
"""
        actual_aligned_csv = StringIO()
        actual_clipping_csv = StringIO()
        sam2aln(remap_file,
                actual_aligned_csv,
                clipping_csv=actual_clipping_csv)

        self.assertMultiLineEqual(expected_aligned_csv,
                                  actual_aligned_csv.getvalue())
        self.assertMultiLineEqual(expected_clipping_csv,
                                  actual_clipping_csv.getvalue())


class CigarTest(unittest.TestCase):
    def setUp(self):
        super(CigarTest, self).setUp()
        self.addTypeEqualityFunc(str, self.assertMultiLineEqual)

    def testTrivial(self):
        cigar = '9M'
        inp_sequence = 'AAACAACCA'
        inp__quality = 'BBBBBBBBB'
        exp_sequence = inp_sequence
        exp__quality = inp__quality

        clipped_seq, clipped_quality, inserts = apply_cigar(cigar, inp_sequence, inp__quality)

        self.assertEqual(exp_sequence, clipped_seq)
        self.assertEqual(exp__quality, clipped_quality)
        self.assertEqual({}, inserts)

    def testDeletion(self):
        cigar = '6M3D3M'
        inp_sequence = 'AAACAACCA'
        inp__quality = 'BBBDDDEEE'
        exp_sequence = 'AAACAA---CCA'
        exp__quality = 'BBBDDD   EEE'

        clipped_seq, clipped_quality, inserts = apply_cigar(cigar, inp_sequence, inp__quality)

        self.assertEqual(exp_sequence, clipped_seq)
        self.assertEqual(exp__quality, clipped_quality)
        self.assertEqual({}, inserts)

    def testSoftClip(self):
        cigar = '3S6M'
        inp_sequence = 'AAACAACCA'
        inp__quality = 'BBBDDDEEE'
        expect_sequence = 'CAACCA'
        expect__quality = 'DDDEEE'

        clipped_seq, clipped_quality, inserts = apply_cigar(cigar, inp_sequence, inp__quality)

        self.assertEqual(expect_sequence, clipped_seq)
        self.assertEqual(expect__quality, clipped_quality)
        self.assertEqual({}, inserts)

    def testInsertion(self):
        cigar = '3M3I6M'
        inp_sequence = 'AAACAACCACCC'
        inp__quality = 'BBBDDDEEEFFF'
        exp_sequence = 'AAACCACCC'
        exp__quality = 'BBBEEEFFF'

        clipped_seq, clipped_quality, inserts = apply_cigar(cigar, inp_sequence, inp__quality)

        self.assertEqual(exp_sequence, clipped_seq)
        self.assertEqual(exp__quality, clipped_quality)
        self.assertEqual({3: ('CAA', 'DDD')}, inserts)

    def testInsertionLowQuality(self):
        cigar = '3M3I6M'
        inp_sequence = 'AAACAACCACCC'
        inp__quality = 'BBBD*DEEEFFF'
        exp_sequence = 'AAACCACCC'
        exp__quality = 'BBBEEEFFF'

        clipped_seq, clipped_quality, inserts = apply_cigar(cigar, inp_sequence, inp__quality)

        self.assertEqual(exp_sequence, clipped_seq)
        self.assertEqual(exp__quality, clipped_quality)
        self.assertEqual({3: ('CAA', 'D*D')}, inserts)

    def testLargeToken(self):
        cigar = '12M'
        inp_sequence = 'AAACAACCACCC'
        inp__quality = 'BBBBBBBBBBBB'
        exp_sequence = inp_sequence
        exp__quality = inp__quality

        clipped_seq, clipped_quality, inserts = apply_cigar(cigar, inp_sequence, inp__quality)

        self.assertEqual(exp_sequence, clipped_seq)
        self.assertEqual(exp__quality, clipped_quality)
        self.assertEqual({}, inserts)

    def testPadding(self):
        cigar = '12M'
        inp_sequence = 'AAACAACCACCC'
        inp__quality = 'BBBDDDEEEFFF'
        pos = 3
        exp_sequence = '---AAACAACCACCC'
        exp__quality = '!!!BBBDDDEEEFFF'

        clipped_seq, clipped_quality, inserts = apply_cigar(cigar, inp_sequence, inp__quality, pos)

        self.assertEqual(exp_sequence, clipped_seq)
        self.assertEqual(exp__quality, clipped_quality)
        self.assertEqual({}, inserts)

    def testClipping(self):
        cigar = '12M'
        inp_sequence = 'AAACAACCACCC'
        inp__quality = 'BBBDDDEEEFFF'
        pos = 0
        clip_from = 3
        clip_to = 8
        exp_sequence = 'CAACCA'
        exp__quality = 'DDDEEE'

        clipped_seq, clipped_quality, inserts = apply_cigar(
            cigar,
            inp_sequence,
            inp__quality,
            pos,
            clip_from,
            clip_to)

        self.assertEqual(exp_sequence, clipped_seq)
        self.assertEqual(exp__quality, clipped_quality)
        self.assertEqual({}, inserts)

    def testClipInsertion(self):
        cigar = '6M3I6M'
        inp_sequence = 'AAACAAGGGCCACCC'
        inp__quality = 'BBBDDDHHHEEEFFF'
        pos = 0
        clip_from = 3
        clip_to = 8
        exp_sequence = 'CAACCA'
        exp__quality = 'DDDEEE'

        clipped_seq, clipped_quality, inserts = apply_cigar(
            cigar,
            inp_sequence,
            inp__quality,
            pos,
            clip_from,
            clip_to)

        self.assertEqual(exp_sequence, clipped_seq)
        self.assertEqual(exp__quality, clipped_quality)
        self.assertEqual({3: ('GGG', 'HHH')}, inserts)

    def testClipInsertionLowQuality(self):
        cigar = '6M3I6M'
        inp_sequence = 'AAACAAGGGCCACCC'
        inp__quality = 'BBBDDDHH*EEEFFF'
        pos = 0
        clip_from = 3
        clip_to = 8
        exp_sequence = 'CAACCA'
        exp__quality = 'DDDEEE'

        clipped_seq, clipped_quality, inserts = apply_cigar(
            cigar,
            inp_sequence,
            inp__quality,
            pos,
            clip_from,
            clip_to)

        self.assertEqual(exp_sequence, clipped_seq)
        self.assertEqual(exp__quality, clipped_quality)
        self.assertEqual({3: ('GGG', 'HH*')}, inserts)

    def testInsertionBeforeClip(self):
        cigar = '3M3I9M'
        inp_sequence = 'AAAGGGCAACCACCC'
        inp__quality = 'BBBHHHDDDEEEFFF'
        pos = 0
        clip_from = 3
        clip_to = 8
        exp_sequence = 'CAACCA'
        exp__quality = 'DDDEEE'

        clipped_seq, clipped_quality, inserts = apply_cigar(
            cigar,
            inp_sequence,
            inp__quality,
            pos,
            clip_from,
            clip_to)

        self.assertEqual(exp_sequence, clipped_seq)
        self.assertEqual(exp__quality, clipped_quality)
        self.assertEqual({0: ('GGG', 'HHH')}, inserts)

    def testInsertionInsideClipRegionWithOffset(self):
        cigar = '2M1I2M'
        inp_sequence = 'TAGCT'
        inp__quality = 'AABCC'
        pos = 3
        clip_from = 4
        clip_to = 20
        exp_sequence = 'ACT'
        exp__quality = 'ACC'

        clipped_seq, clipped_quality, inserts = apply_cigar(
            cigar,
            inp_sequence,
            inp__quality,
            pos,
            clip_from,
            clip_to)

        self.assertEqual(exp_sequence, clipped_seq)
        self.assertEqual(exp__quality, clipped_quality)
        self.assertEqual({1: ('G', 'B')}, inserts)

    def testInsertionAfterClipRegionWithOffset(self):
        cigar = '5M1I2M'
        inp_sequence = 'TAGCTCAG'
        inp__quality = 'AAAAABCC'
        pos = 10
        clip_from = 10
        clip_to = 13
        exp_sequence = 'TAGC'
        exp__quality = 'AAAA'
        expected_inserts = {}

        clipped_seq, clipped_quality, inserts = apply_cigar(
            cigar,
            inp_sequence,
            inp__quality,
            pos,
            clip_from,
            clip_to)

        self.assertEqual(exp_sequence, clipped_seq)
        self.assertEqual(exp__quality, clipped_quality)
        self.assertEqual(expected_inserts, inserts)

    def testClippingEverything(self):
        cigar = '12M'
        inp_sequence = 'AAACAACCACCC'
        inp__quality = 'BBBDDDEEEFFF'
        pos = 0
        clip_from = 100
        clip_to = 108
        exp_sequence = ''
        exp__quality = ''

        clipped_seq, clipped_quality, inserts = apply_cigar(
            cigar,
            inp_sequence,
            inp__quality,
            pos,
            clip_from,
            clip_to)

        self.assertEqual(exp_sequence, clipped_seq)
        self.assertEqual(exp__quality, clipped_quality)
        self.assertEqual({}, inserts)

    def testSoftClipPositions(self):
        cigar = '3S6M'
        pos = 4
        inp__sequence = 'AAACAACCA'
        inp___quality = 'BBBDDDEEE'
        exp_sequence = '----CAACCA'
        exp__quality = '!!!!DDDEEE'
        mapped = set()
        soft_clipped = set()
        expected_mapped = {4, 5, 6, 7, 8, 9}
        expected_soft_clipped = {1, 2, 3}

        clipped_seq, clipped_quality, inserts = apply_cigar(
            cigar,
            inp__sequence,
            inp___quality,
            pos=pos,
            mapped=mapped,
            soft_clipped=soft_clipped)

        self.assertEqual(exp_sequence, clipped_seq)
        self.assertEqual(exp__quality, clipped_quality)
        self.assertEqual({}, inserts)
        self.assertEqual(expected_mapped, mapped)
        self.assertEqual(expected_soft_clipped, soft_clipped)

    def testInvalidCigar(self):
        cigar = '3M...6M'
        inp_sequence = 'AAACAACCACCC'
        inp__quality = 'BBBDDDEEEFFF'

        with self.assertRaises(RuntimeError) as result:
            apply_cigar(cigar, inp_sequence, inp__quality)

        self.assertEqual(
            "Invalid CIGAR string: '3M...6M'.",
            result.exception.args[0])

    def testUnsupportedCigarToken(self):
        cigar = '3M3X6M'
        inp_sequence = 'AAACAACCACCC'
        inp__quality = 'BBBDDDEEEFFF'

        with self.assertRaises(RuntimeError) as result:
            apply_cigar(cigar, inp_sequence, inp__quality)

        self.assertEqual(
            "Unsupported CIGAR token: '3X'.",
            result.exception.args[0])

    def testShortCigar(self):
        cigar = '8M'
        inp_sequence = 'AAACAACCA'
        inp__quality = 'BBBDDDEEE'

        with self.assertRaises(RuntimeError) as result:
            apply_cigar(cigar, inp_sequence, inp__quality)

        self.assertEqual(
            "CIGAR string '8M' is too short for sequence 'AAACAACCA'.",
            result.exception.args[0])

    def testLongCigar(self):
        cigar = '10M'
        inp_sequence = 'AAACAACCA'
        inp__quality = 'BBBDDDEEE'

        with self.assertRaises(RuntimeError) as result:
            apply_cigar(cigar, inp_sequence, inp__quality)

        self.assertEqual(
            "CIGAR string '10M' is too long for sequence 'AAACAACCA'.",
            result.exception.args[0])

    def testInsertionAfterClipping(self):
        cigar = '3M3I3M'
        inp_sequence = "ACTTAGAAA"
        inp__quality = 'AAABBBDDD'
        pos = 0
        clip_from = 0
        clip_to = 2
        exp_sequence = 'ACT'
        exp__quality = 'AAA'

        clipped_seq, clipped_quality, inserts = apply_cigar(
            cigar,
            inp_sequence,
            inp__quality,
            pos,
            clip_from,
            clip_to)

        self.assertEqual(exp_sequence, clipped_seq)
        self.assertEqual(exp__quality, clipped_quality)
        self.assertEqual({}, inserts)

    def testInsertionAtEndOfClipping(self):
        cigar = '3M3I3M'
        inp_sequence = "ACTTAGAAA"
        inp__quality = 'AAABBBDDD'
        pos = 0
        clip_from = 0
        clip_to = 3
        exp_sequence = 'ACTA'
        exp__quality = 'AAAD'

        clipped_seq, clipped_quality, inserts = apply_cigar(
            cigar,
            inp_sequence,
            inp__quality,
            pos,
            clip_from,
            clip_to)

        self.assertEqual(exp_sequence, clipped_seq)
        self.assertEqual(exp__quality, clipped_quality)
        self.assertEqual({3: ('TAG', 'BBB')}, inserts)

    def testInsertionAfterInsertion(self):
        cigar = '3M3I3M3I3M'
        (seq,
         quality,
         expected_seq,
         expected_quality) = ('TTTGGGCCCAAATTT',
                              '111222333444555',
                              'TTTCCCTTT',
                              '111333555')
        expected_inserts = {3: ('GGG', '222'), 6: ('AAA', '444')}

        seq, quality, inserts = apply_cigar(
          cigar,
          seq,
          quality)

        self.assertEqual(expected_seq, seq)
        self.assertEqual(expected_quality, quality)
        self.assertEqual(expected_inserts, inserts)

    def testInsertionAfterDeletion(self):
        cigar = '3M3D3M3I3M'
        (seq,
         quality,
         expected_seq,
         expected_quality) = ('TTTCCCAAATTT',
                              '111222333444',
                              'TTT---CCCTTT',
                              '111   222444')
        expected_inserts = {9: ('AAA', '333')}

        seq, quality, inserts = apply_cigar(
          cigar,
          seq,
          quality)

        self.assertEqual(expected_seq, seq)
        self.assertEqual(expected_quality, quality)
        self.assertEqual(expected_inserts, inserts)


class MergePairsTest(unittest.TestCase):
    def setUp(self):
        self.addTypeEqualityFunc(str, self.assertMultiLineEqual)

    def testSimple(self):
        sequence1 = 'ACTGCA'
        sequence2 = 'ACTGCA'
        quality_1 = 'JJJJJJ'
        quality_2 = 'JJJJJJ'
        exp_m_seq = 'ACTGCA'

        mseq = merge_pairs(sequence1, sequence2, quality_1, quality_2)

        self.assertEqual(exp_m_seq, mseq)

    def testDifferentLength(self):
        sequence1 = 'ACTGCATCT'
        sequence2 = 'ACTGCA'
        quality_1 = 'JJJJJJJJJ'
        quality_2 = 'JJJJJJ'
        exp_m_seq = 'ACTGCATCT'

        mseq = merge_pairs(sequence1, sequence2, quality_1, quality_2)

        self.assertEqual(exp_m_seq, mseq)

    def testOffset(self):
        sequence1 = '-CTGCA'
        sequence2 = '---GCATCT'
        quality_1 = '!JJJJJ'
        quality_2 = '!!!JJJJJJ'
        exp_m_seq = '-CTGCATCT'

        mseq = merge_pairs(sequence1, sequence2, quality_1, quality_2)

        self.assertEqual(exp_m_seq, mseq)

    def testForwardDeletion(self):
        sequence1 = 'C-GCA'
        sequence2 = '---CATCT'
        quality_1 = 'J!JJJ'
        quality_2 = '!!!JJJJJ'
        exp_m_seq = 'C-GCATCT'

        mseq = merge_pairs(sequence1, sequence2, quality_1, quality_2)

        self.assertEqual(exp_m_seq, mseq)

    def testReverseDeletion(self):
        sequence1 = 'CTGCA'
        sequence2 = '--GCAT-T'
        quality_1 = 'JJJJJ'
        quality_2 = '!!JJJJ!J'
        exp_m_seq = 'CTGCAT-T'

        mseq = merge_pairs(sequence1, sequence2, quality_1, quality_2)

        self.assertEqual(exp_m_seq, mseq)

    def testDisagreementWithDifferentQualityFirstHigher(self):
        sequence1 = 'AGTGCA'
        sequence2 = 'ACTGCA'
        quality_1 = 'JJJJJJ'
        quality_2 = 'JEJJJJ'
        exp_m_seq = 'AGTGCA'

        mseq = merge_pairs(sequence1, sequence2, quality_1, quality_2)

        self.assertEqual(exp_m_seq, mseq)

    def testDisagreementWithCloseQualityFirstHigher(self):
        sequence1 = 'AGTGCA'
        sequence2 = 'ACTGCA'
        quality_1 = 'JJJJJJ'
        quality_2 = 'JFJJJJ'
        exp_m_seq = 'ANTGCA'

        mseq = merge_pairs(sequence1, sequence2, quality_1, quality_2)

        self.assertEqual(exp_m_seq, mseq)

    def testDisagreementWithDifferentQualitySecondHigher(self):
        sequence1 = 'AGTGCA'
        sequence2 = 'ACTGCA'
        quality_1 = 'JEJJJJ'
        quality_2 = 'JJJJJJ'
        exp_m_seq = 'ACTGCA'

        mseq = merge_pairs(sequence1, sequence2, quality_1, quality_2)

        self.assertEqual(exp_m_seq, mseq)

    def testDisagreementWithCloseQualitySecondHigher(self):
        sequence1 = 'AGTGCA'
        sequence2 = 'ACTGCA'
        quality_1 = 'JFJJJJ'
        quality_2 = 'JJJJJJ'
        exp_m_seq = 'ANTGCA'

        mseq = merge_pairs(sequence1, sequence2, quality_1, quality_2)

        self.assertEqual(exp_m_seq, mseq)

    def testDisagreementWithLowQuality(self):
        sequence1 = 'AGTGCA'
        sequence2 = 'ACTGCA'
        quality_1 = 'J!JJJJ'
        quality_2 = 'J*JJJJ'
        exp_m_seq = 'ANTGCA'

        mseq = merge_pairs(sequence1, sequence2, quality_1, quality_2)

        self.assertEqual(exp_m_seq, mseq)

    def testGap(self):
        sequence1 = 'AGT'
        sequence2 = '------GCA'
        quality_1 = 'JJJ'
        quality_2 = '!!!!!!JJJ'
        exp_m_seq = 'AGTnnnGCA'

        mseq = merge_pairs(sequence1, sequence2, quality_1, quality_2)

        self.assertEqual(exp_m_seq, mseq)

    def testLowQualityInSecondRead(self):
        sequence1 = 'AGT'
        sequence2 = '---GCA'
        quality_1 = 'JJJ'
        quality_2 = '!!!J*J'
        exp_m_seq = 'AGTGNA'

        mseq = merge_pairs(sequence1, sequence2, quality_1, quality_2)

        self.assertEqual(exp_m_seq, mseq)

    def testOneInsertion(self):
        sequence1 = 'AGT'
        sequence2 = '---GCA'
        quality_1 = 'JJJ'
        quality_2 = '!!!JJJ'
        ins1 = {2: ('CCC', 'JJJ')}
        ins2 = {}
        exp_m_seq = 'AGCCCTGCA'

        mseq = merge_pairs(sequence1, sequence2, quality_1, quality_2, ins1, ins2)

        self.assertEqual(exp_m_seq, mseq)

    def testTwoInsertions(self):
        sequence1 = 'AGT'
        sequence2 = '---GCA'
        quality_1 = 'JJJ'
        quality_2 = '!!!JJJ'
        ins1 = {2: ('CCC', 'JJJ')}
        ins2 = {5: ('TTT', 'JJJ')}
        exp_m_seq = 'AGCCCTGCTTTA'

        mseq = merge_pairs(sequence1, sequence2, quality_1, quality_2, ins1, ins2)

        self.assertEqual(exp_m_seq, mseq)

    def testConflictingInsertions(self):
        sequence1 = 'AGTGCA'
        sequence2 = 'AGTGCA'
        quality_1 = 'JJJJJJ'
        quality_2 = 'JJJJJJ'
        ins1 = {2: ('CCC', 'JJJ')}
        ins2 = {2: ('CTC', 'JAJ')}
        exp_m_seq = 'AGCCCTGCA'

        mseq = merge_pairs(sequence1, sequence2, quality_1, quality_2, ins1, ins2)

        self.assertEqual(exp_m_seq, mseq)


class MergeInsertionsTest(unittest.TestCase):
    def setUp(self):
        self.addTypeEqualityFunc(str, self.assertMultiLineEqual)

    def testSeparateInsertions(self):
        ins1 = {2: ('CCC', 'JJJ')}
        ins2 = {10: ('GGG', 'JJJ')}
        expected_merged = {2: 'CCC', 10: 'GGG'}

        merged = merge_inserts(ins1, ins2)

        self.assertEqual(expected_merged, merged)

    def testIdenticalInsertions(self):
        ins1 = {2: ('CCC', 'JJJ')}
        ins2 = {2: ('CCC', 'JJJ')}
        expected_merged = {2: 'CCC'}

        merged = merge_inserts(ins1, ins2)

        self.assertEqual(expected_merged, merged)

    def testConflictingInsertions(self):
        ins1 = {2: ('CCC', 'JJJ')}
        ins2 = {2: ('CTC', 'JAJ')}
        expected_merged = {2: 'CCC'}

        merged = merge_inserts(ins1, ins2)

        self.assertEqual(expected_merged, merged)

    def testConflictingInsertionsTooCloseInQuality(self):
        ins1 = {2: ('CCC', 'JJJ')}
        ins2 = {2: ('CTC', 'JAJ')}
        expected_merged = {2: 'CNC'}

        merged = merge_inserts(ins1, ins2, minimum_q_delta=20)

        self.assertEqual(expected_merged, merged)

    def testInsertionQualityTooLowForward(self):
        ins1 = {2: ('CCC', 'JAJ')}
        ins2 = {10: ('GGG', 'JJJ')}
        expected_merged = {10: 'GGG'}

        merged = merge_inserts(ins1, ins2, q_cutoff=32)

        self.assertEqual(expected_merged, merged)

    def testInsertionQualityTooLowReverse(self):
        ins1 = {2: ('CCC', 'JJJ')}
        ins2 = {10: ('GGG', 'JAJ')}
        expected_merged = {2: 'CCC'}

        merged = merge_inserts(ins1, ins2, q_cutoff=32)

        self.assertEqual(expected_merged, merged)

    def testInsertionQualityJustEnough(self):
        ins1 = {2: ('CCC', 'JJJ')}
        ins2 = {10: ('GGG', 'JBJ')}
        expected_merged = {2: 'CCC', 10: 'GGG'}

        merged = merge_inserts(ins1, ins2, q_cutoff=32)

        self.assertEqual(expected_merged, merged)

    def testNone(self):
        ins1 = None
        ins2 = None
        expected_merged = {}

        merged = merge_inserts(ins1, ins2)

        self.assertEqual(expected_merged, merged)
