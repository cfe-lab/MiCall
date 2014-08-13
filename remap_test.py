import unittest

import remap

class RemapTest(unittest.TestCase):
    def testSampleName(self):
        filename = '/some/path/2020A-V3LOOP_S4_L001_R1_001.fastq'
        expected_sample_name = '2020A-V3LOOP_S4'
        
        sample_name = remap.calculate_sample_name(filename)
        
        self.assertEqual(sample_name, expected_sample_name)
