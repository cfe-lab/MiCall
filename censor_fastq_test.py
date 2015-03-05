import StringIO
import unittest

import censor_fastq

#TODO: remaining tests:
# * multiple directions
# * multiple reads
class CensorTest(unittest.TestCase):
    def setUp(self):
        self.original_text = """\
@M01841:45:000000000-A5FEG:1:1101:5296:13227 1:N:0:9
ACGT
+
AAAA
"""
        self.original_file = StringIO.StringIO(self.original_text)
        self.bad_cycles = []
        self.censored_file = StringIO.StringIO()
     
    def testNoBadCycles(self):
        expected_text = self.original_text
        
        censor_fastq.censor(self.original_file,
                            self.bad_cycles,
                            self.censored_file)
         
        self.assertMultiLineEqual(expected_text, self.censored_file.getvalue())
     
    def testBadCycle(self):
        self.bad_cycles = [{'tile': '1101', 'cycle': '3'}]
        expected_text = """\
@M01841:45:000000000-A5FEG:1:1101:5296:13227 1:N:0:9
ACNT
+
AA#A
"""

        censor_fastq.censor(self.original_file,
                            self.bad_cycles,
                            self.censored_file)
         
        self.assertMultiLineEqual(expected_text, self.censored_file.getvalue())

    def testDifferentTile(self):
        self.bad_cycles = [{'tile': '1102', 'cycle': '3'}]
        expected_text = self.original_text
        
        censor_fastq.censor(self.original_file,
                            self.bad_cycles,
                            self.censored_file)
         
        self.assertMultiLineEqual(expected_text, self.censored_file.getvalue())

    def testDifferentDirection(self):
        self.original_text = """\
@M01841:45:000000000-A5FEG:1:1101:5296:13227 2:N:0:9
ACGT
+
AAAA
"""
        self.original_file = StringIO.StringIO(self.original_text)
        self.bad_cycles = [{'tile': '1101', 'cycle': '3'}]
        expected_text = self.original_text
        
        censor_fastq.censor(self.original_file,
                            self.bad_cycles,
                            self.censored_file)
         
        self.assertMultiLineEqual(expected_text, self.censored_file.getvalue())

    def testReverseDirection(self):
        self.original_text = """\
@M01841:45:000000000-A5FEG:1:1101:5296:13227 2:N:0:9
ACGT
+
AAAA
"""
        self.original_file = StringIO.StringIO(self.original_text)
        self.bad_cycles = [{'tile': '1101', 'cycle': '-3'}]
        expected_text = """\
@M01841:45:000000000-A5FEG:1:1101:5296:13227 2:N:0:9
ACNT
+
AA#A
"""
        
        censor_fastq.censor(self.original_file,
                            self.bad_cycles,
                            self.censored_file)
         
        self.assertMultiLineEqual(expected_text, self.censored_file.getvalue())

    def testTwoReads(self):
        self.original_text = """\
@M01841:45:000000000-A5FEG:1:1101:5296:13227 1:N:0:9
ACGT
+
AAAA
@M01841:45:000000000-A5FEG:1:1102:1234:12345 1:N:0:9
TGCA
+
BBBB
"""
        self.original_file = StringIO.StringIO(self.original_text)
        self.bad_cycles = [{'tile': '1101', 'cycle': '2'},
                           {'tile': '1102', 'cycle': '3'}]
        expected_text = """\
@M01841:45:000000000-A5FEG:1:1101:5296:13227 1:N:0:9
ANGT
+
A#AA
@M01841:45:000000000-A5FEG:1:1102:1234:12345 1:N:0:9
TGNA
+
BB#B
"""
        
        censor_fastq.censor(self.original_file,
                            self.bad_cycles,
                            self.censored_file)
         
        self.assertMultiLineEqual(expected_text, self.censored_file.getvalue())
