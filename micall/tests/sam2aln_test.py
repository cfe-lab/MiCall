import unittest
from StringIO import StringIO

from micall.core import sam2aln


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
Example_read_1,99,V3LOOP,1,44,32M,=,1,-32,TGT,"A,A"
Example_read_1,147,V3LOOP,1,44,32M,=,1,-32,TGT,"A""A"
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
qname,qcut,seq1,qual1,seq2,qual2,prop_N,mseq,mapq1,mapq2
Example_read_2,,TGTACAAGACCCAACAACAATACAAGAAAAAG,AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA,TGTACAAGACCCAACAACAATACAAGAAAAAG,AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA,,,8,44
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
qname,qcut,seq1,qual1,seq2,qual2,prop_N,mseq,mapq1,mapq2
Example_read_1,15,TGTACAAGACCCAACAACAATACAAGAAAAAG,000000000000000000AAAAAAAAAAAAAA,TGTACAAGACCCAACAACAATACAAGAAAAAG,000000000000000000AAAAAAAAAAAAAA,0.5625,NNNNNNNNNNNNNNNNNNAATACAAGAAAAAG,,
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
