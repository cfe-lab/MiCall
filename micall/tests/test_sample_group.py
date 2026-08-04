from unittest import TestCase

from micall.drivers.sample_group import SampleGroup
from micall.drivers.sample import Sample


class SampleGroupTest(TestCase):
    def test_pair(self):
        expected_paths = ['1_R1_001.fastq', '2_R1_001.fastq']
        group = SampleGroup(Sample(fastq1='1_R1_001.fastq'),
                            Sample(fastq1='2_R1_001.fastq'))

        paths = [sample.fastq1 for sample in group]

        self.assertEqual(expected_paths, paths)

    def test_iter(self):
        expected_paths = ['1_R1_001.fastq']
        group = SampleGroup(Sample(fastq1='1_R1_001.fastq'))

        paths = [sample.fastq1 for sample in group]

        self.assertEqual(expected_paths, paths)

    def test_repr(self):
        expected_repr = ("SampleGroup(Sample(fastq1='x_R1_001.fastq'), "
                         "Sample(fastq1='b_R1_001.fastq'))")
        group = SampleGroup(Sample(fastq1='x_R1_001.fastq'),
                            Sample(fastq1='b_R1_001.fastq'))

        self.assertEqual(expected_repr, repr(group))
