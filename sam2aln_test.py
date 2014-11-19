import unittest
import StringIO

import sam2aln


class RemapReaderTest(unittest.TestCase):
    def test_grouping(self):
        remap_file = StringIO.StringIO("""\
sample_name,qname,flag,rname,pos,mapq,cigar,rnext,pnext,tlen,seq,qual
1234A,Example_read_1,99,V3LOOP,1,44,32M,=,1,-32,TGTACAAGACCCAACAACAATACAAGAAAAAG,AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
1234A,Example_read_1,147,V3LOOP,1,44,32M,=,1,-32,TGTACAAGACCCAACAACAATACAAGAAAAAG,AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
1234A,Example_read_2,99,INT,1,44,32M,=,1,-32,TGTACAAGACCCAACAACAATACAAGAAAAAG,AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
1234A,Example_read_2,147,INT,1,44,32M,=,1,-32,TGTACAAGACCCAACAACAATACAAGAAAAAG,AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
""")
        
        reader = sam2aln.RemapReader(remap_file)
        sample_names = []
        regions = []
        for sample_name, region, _group in reader.read_groups():
            sample_names.append(sample_name)
            regions.append(region)

        self.assertEqual(regions, ['V3LOOP', 'INT'])
        self.assertEqual(sample_names, ['1234A', '1234A'])
