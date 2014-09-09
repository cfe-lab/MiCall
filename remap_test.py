import unittest

import remap
from remap import is_first_read

class RemapTest(unittest.TestCase):
    def testSampleName(self):
        filename = '/some/path/2020A-V3LOOP_S4_L001_R1_001.fastq'
        expected_sample_name = '2020A-V3LOOP_S4'
        
        sample_name = remap.calculate_sample_name(filename)
        
        self.assertEqual(sample_name, expected_sample_name)

class IsFirstReadTest(unittest.TestCase):
    def testFirstRead(self):
        flag = '99'
        isFirstExpected = True
        
        isFirst = is_first_read(flag)
        
        self.assertEqual(isFirstExpected, isFirst)

    def testSecondRead(self):
        flag = '147'
        isFirstExpected = False
        
        isFirst = is_first_read(flag)
        
        self.assertEqual(isFirstExpected, isFirst)        
        
    def testSmallFlag(self):
        flag = '3'
        isFirstExpected = False
        
        isFirst = is_first_read(flag)
        
        self.assertEqual(isFirstExpected, isFirst)
