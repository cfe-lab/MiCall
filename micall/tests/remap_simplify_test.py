import unittest

from micall.utils.remap_simplify import SamBase, SamHeader
from unittest.util import safe_repr

class SamBaseTest(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super(SamBaseTest, self).__init__(*args, **kwargs)
        self.addTypeEqualityFunc(str, self.assertMultiLineEqual)
        self.addTypeEqualityFunc(SamBase, self.assertSamBaseEqual)
        self.addTypeEqualityFunc(SamHeader, self.assertSamHeaderEqual)
        

    def assertAttributesEqual(self, base1, base2, msg, attr_names):
        longMessage = self.longMessage
        self.longMessage = True
        try:
            for attr_name in attr_names:
                attr1 = getattr(base1, attr_name)
                attr2 = getattr(base2, attr_name)
                context_msg = self._formatMessage(msg, attr_name)
                self.assertEqual(attr1, attr2, context_msg)
        finally:
            
            self.longMessage = longMessage

    def assertSamBaseEqual(self, base1, base2, msg=None):
        attr_names = 'nuc quality pos token_type qname flag rname mapq'.split()
        self.assertAttributesEqual(base1, base2, msg, attr_names)

    def assertSamHeaderEqual(self, base1, base2, msg=None):
        attr_names = ['line']
        self.assertAttributesEqual(base1, base2, msg, attr_names)
    
    def assertSplits(self, sam_lines, expected_sam_bases, msg=None):
        sam_bases = SamBase.split_bases(sam_lines)
        
        if len(expected_sam_bases) != len(sam_bases):
            self.assertEqual(expected_sam_bases, sam_bases, msg)
        else:
            longMessage = self.longMessage
            self.longMessage = True
            try:
                for i, (expected_base, base) in enumerate(zip(expected_sam_bases,
                                                              sam_bases)):
                    if expected_base != base:
                        context_msg = self._formatMessage(
                            msg,
                            'base index {} ({} and {})'.format(
                                i,
                                safe_repr(expected_base),
                                safe_repr(base)))
                    self.assertEqual(expected_base, base, context_msg)
            finally:
                self.longMessage = longMessage

    def assertJoins(self, sam_bases, expected_sam_lines, msg=None):
        sam_lines = SamBase.join(sam_bases)
        
        expected_text = ''.join(expected_sam_lines)
        actual_text = ''.join(sam_lines)
        
        # Tabs mess up the failure display, so replace them with spaces
        expected_text = expected_text.replace('\t', ' ')
        actual_text = actual_text.replace('\t', ' ')
        
        self.assertEqual(expected_text, actual_text)
                
    def testSimple(self):
        qname = 'test1'
        flag = '99'
        rname = 'test'
        mapq = '44'
        sam_bases = [SamBase('A', 'J', 1, 'M', qname, flag, rname, mapq),
                     SamBase('C', 'K', 2, 'M', qname, flag, rname, mapq)]
        #SAM:qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual
        sam_lines = ["test1\t99\ttest\t1\t44\t2M\t=\t1\t2\tAC\tJK\n"]

        self.assertJoins(sam_bases, sam_lines)
        self.assertSplits(sam_lines, sam_bases)
    
    def testTwoLines(self):
        qname = 'test1'
        flag_forward = '99'
        flag_reverse = '147'
        rname = 'test'
        mapq = '44'
        sam_bases = [SamBase('A', 'J', 1, 'M', qname, flag_forward, rname, mapq),
                     SamBase('C', 'K', 2, 'M', qname, flag_forward, rname, mapq),
                     SamBase('G', 'M', 1, 'M', qname, flag_reverse, rname, mapq),
                     SamBase('T', 'L', 2, 'M', qname, flag_reverse, rname, mapq)]
        #SAM:qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual
        sam_lines = ["test1\t99\ttest\t1\t44\t2M\t=\t1\t2\tAC\tJK\n",
                     "test1\t147\ttest\t1\t44\t2M\t=\t1\t-2\tGT\tML\n"]
         
        self.assertJoins(sam_bases, sam_lines)
        self.assertSplits(sam_lines, sam_bases)
     
    def testOffset(self):
        qname = 'test1'
        flag = '99'
        rname = 'test'
        mapq = '44'
        sam_bases = [SamBase('A', 'J', 10, 'M', qname, flag, rname, mapq),
                     SamBase('C', 'K', 11, 'M', qname, flag, rname, mapq)]
        #SAM:qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual
        sam_lines = ["test1\t99\ttest\t10\t44\t2M\t=\t10\t2\tAC\tJK\n"]
         
        self.assertJoins(sam_bases, sam_lines)
        self.assertSplits(sam_lines, sam_bases)
     
    def testInsert(self):
        qname = 'test1'
        flag = '99'
        rname = 'test'
        mapq = '44'
        sam_bases = [SamBase('A', 'J', 1, 'M', qname, flag, rname, mapq),
                     SamBase('G', 'F', 2, 'M', qname, flag, rname, mapq),
                     SamBase('T', 'L', None, 'I', qname, flag, rname, mapq),
                     SamBase('C', 'K', 3, 'M', qname, flag, rname, mapq)]
        #SAM:qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual
        sam_lines = ["test1\t99\ttest\t1\t44\t2M1I1M\t=\t1\t4\tAGTC\tJFLK\n"]
         
        self.assertJoins(sam_bases, sam_lines)
        self.assertSplits(sam_lines, sam_bases)
     
    def testDelete(self):
        qname = 'test1'
        flag = '99'
        rname = 'test'
        mapq = '44'
        sam_bases = [SamBase('A', 'J', 1, 'M', qname, flag, rname, mapq),
                     SamBase(None, None, 2, 'D', qname, flag, rname, mapq),
                     SamBase('C', 'K', 3, 'M', qname, flag, rname, mapq)]
        #SAM:qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual
        sam_lines = ["test1\t99\ttest\t1\t44\t1M1D1M\t=\t1\t2\tAC\tJK\n"]
        
        self.assertJoins(sam_bases, sam_lines)
        self.assertSplits(sam_lines, sam_bases)
     
    def testDeleteAtEndOfLine(self):
        qname = 'test1'
        qname2 = 'test2'
        flag = '99'
        rname = 'test'
        mapq = '44'
        sam_bases = [SamBase('A', 'J', 1, 'M', qname, flag, rname, mapq),
                     SamBase(None, None, 2, 'D', qname, flag, rname, mapq),
                     SamBase('C', 'K', 3, 'M', qname2, flag, rname, mapq)]
        #SAM:qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual
        sam_lines = ["test1\t99\ttest\t1\t44\t1M\t=\t1\t1\tA\tJ\n",
                     "test2\t99\ttest\t3\t44\t1M\t=\t3\t1\tC\tK\n"]
        
        self.assertJoins(sam_bases, sam_lines)
     
    def testDeleteAtStartOfLine(self):
        qname = 'test1'
        qname2 = 'test2'
        flag = '99'
        rname = 'test'
        mapq = '44'
        sam_bases = [SamBase('A', 'J', 1, 'M', qname, flag, rname, mapq),
                     SamBase(None, None, 2, 'D', qname2, flag, rname, mapq),
                     SamBase('C', 'K', 3, 'M', qname2, flag, rname, mapq)]
        #SAM:qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual
        sam_lines = ["test1\t99\ttest\t1\t44\t1M\t=\t1\t1\tA\tJ\n",
                     "test2\t99\ttest\t3\t44\t1M\t=\t3\t1\tC\tK\n"]
        
        self.assertJoins(sam_bases, sam_lines)
    
    def testDeleteAfterSoftClip(self):
        qname = 'test1'
        flag = '99'
        rname = 'test'
        mapq = '44'
        sam_bases = [SamBase('A', 'J', None, 'S', qname, flag, rname, mapq),
                     SamBase(None, None, 2, 'D', qname, flag, rname, mapq),
                     SamBase('C', 'K', 3, 'M', qname, flag, rname, mapq)]
        #SAM:qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual
        sam_lines = ["test1\t99\ttest\t3\t44\t1S1M\t=\t3\t2\tAC\tJK\n"]
        
        self.assertJoins(sam_bases, sam_lines)
    
    def testDeleteAfterInsert(self):
        qname = 'test1'
        flag = '99'
        rname = 'test'
        mapq = '44'
        sam_bases = [SamBase('A', 'J', 1, 'M', qname, flag, rname, mapq),
                     SamBase('T', 'L', None, 'I', qname, flag, rname, mapq),
                     SamBase(None, None, 3, 'D', qname, flag, rname, mapq),
                     SamBase('C', 'K', 4, 'M', qname, flag, rname, mapq)]
        #SAM:qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual
        sam_lines = ["test1\t99\ttest\t1\t44\t1M1I1M1D1M\t=\t1\t4\tATAC\tJL#K\n"]
        
        self.assertJoins(sam_bases, sam_lines)
    
    def testInsertAfterDelete(self):
        qname = 'test1'
        flag = '99'
        rname = 'test'
        mapq = '44'
        sam_bases = [SamBase('A', 'J', 1, 'M', qname, flag, rname, mapq),
                     SamBase(None, None, 2, 'D', qname, flag, rname, mapq),
                     SamBase('T', 'L', None, 'I', qname, flag, rname, mapq),
                     SamBase('C', 'K', 4, 'M', qname, flag, rname, mapq)]
        #SAM:qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual
        sam_lines = ["test1\t99\ttest\t1\t44\t1M1D2M\t=\t1\t3\tAAC\tJ#K\n"]
        
        self.assertJoins(sam_bases, sam_lines)
    
    def testSoftClip(self):
        qname = 'test1'
        flag = '99'
        rname = 'test'
        mapq = '44'
        sam_bases = [SamBase('A', 'J', None, 'S', qname, flag, rname, mapq),
                     SamBase('C', 'K', 1, 'M', qname, flag, rname, mapq)]
        #SAM:qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual
        sam_lines = ["test1\t99\ttest\t1\t44\t1S1M\t=\t1\t2\tAC\tJK\n"]
         
        self.assertJoins(sam_bases, sam_lines)
        self.assertSplits(sam_lines, sam_bases)
                
    def testJoinWithGap(self):
        qname = 'test1'
        flag = '99'
        rname = 'test'
        mapq = '44'
        sam_bases = [SamBase('T', 'J', 1, 'M', qname, flag, rname, mapq),
                     SamBase('C', 'K', 3, 'M', qname, flag, rname, mapq)]
        #SAM:qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual
        sam_lines = ["test1\t99\ttest\t1\t44\t3M\t=\t1\t3\tTAC\tJ#K\n"]

        self.assertJoins(sam_bases, sam_lines)
                
    def testHeader(self):
        qname = 'test1'
        flag = '99'
        rname = 'test'
        mapq = '44'
        sam_bases = [SamHeader('@SQ\tSN:test\tLN:1568'),
                     SamBase('C', 'K', 3, 'M', qname, flag, rname, mapq)]
        #SAM:qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual
        sam_lines = ["@SQ\tSN:test\tLN:1568\n",
                     "test1\t99\ttest\t3\t44\t1M\t=\t3\t1\tC\tK\n"]

        self.assertJoins(sam_bases, sam_lines)
        self.assertSplits(sam_lines, sam_bases)
