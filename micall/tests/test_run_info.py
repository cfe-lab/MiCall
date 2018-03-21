from unittest import TestCase

from micall.drivers.run_info import RunInfo
from micall.drivers.sample import Sample
from micall.drivers.sample_group import SampleGroup


class RunInfoTest(TestCase):
    def test_get_all_samples(self):
        expected_fastq_paths = ['1a_R1_001.fastq',
                                '1b_R1_001.fastq',
                                '2_R1_001.fastq']

        run_info = RunInfo(
            sample_groups=[SampleGroup(Sample(fastq1='1a_R1_001.fastq'),
                                       Sample(fastq1='1b_R1_001.fastq')),
                           SampleGroup(Sample(fastq1='2_R1_001.fastq'))])
        fastq_paths = [sample.fastq1 for sample in run_info.get_all_samples()]

        self.assertEqual(expected_fastq_paths, fastq_paths)
