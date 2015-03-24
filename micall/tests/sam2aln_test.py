import unittest
import StringIO

import sam2aln


class RemapReaderTest(unittest.TestCase):
    def test_grouping(self):
        remap_file = StringIO.StringIO("""\
qname,flag,rname,pos,mapq,cigar,rnext,pnext,tlen,seq,qual
Example_read_1,99,V3LOOP,1,44,32M,=,1,-32,TGTACAAGACCCAACAACAATACAAGAAAAAG,AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
Example_read_1,147,V3LOOP,1,44,32M,=,1,-32,TGTACAAGACCCAACAACAATACAAGAAAAAG,AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
Example_read_2,99,INT,1,44,32M,=,1,-32,TGTACAAGACCCAACAACAATACAAGAAAAAG,AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
Example_read_2,147,INT,1,44,32M,=,1,-32,TGTACAAGACCCAACAACAATACAAGAAAAAG,AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
""")
        
        reader = sam2aln.RemapReader(remap_file)
        regions = []
        for region, _group in reader.read_groups():
            regions.append(region)

        self.assertEqual(regions, ['V3LOOP', 'INT'])

    def test_escaping(self):
        remap_file = StringIO.StringIO("""\
qname,flag,rname,pos,mapq,cigar,rnext,pnext,tlen,seq,qual
Example_read_1,99,V3LOOP,1,44,32M,=,1,-32,TGT,"A,A"
Example_read_1,147,V3LOOP,1,44,32M,=,1,-32,TGT,"A""A"
""")
        expected_qualities = ['A,A',
                              'A"A']
        
        reader = sam2aln.RemapReader(remap_file)
        qualities = []
        for _region, group in reader.read_groups():
            for row in group:
                qualities.append(row['qual'])
        
        self.assertSequenceEqual(expected_qualities, qualities)
