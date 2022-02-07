import pickle
from unittest import TestCase

import pytest

from micall.drivers.sample import Sample, exclude_extra_seeds


class SampleTest(TestCase):
    def test_explicit_paths(self):
        expected_amino_csv = '/path/to/amino.csv'
        expected_fastq1_path = '/path/to/1234A_S1_R1_001.fastq.gz'

        sample = Sample(fastq1=expected_fastq1_path,
                        amino_csv=expected_amino_csv)
        fastq1_path = sample.fastq1
        amino_path = sample.amino_csv

        self.assertEqual(expected_fastq1_path, fastq1_path)
        self.assertEqual(expected_amino_csv, amino_path)

    def test_default_csv(self):
        scratch_path = '/path/to/scratch'
        expected_amino_path = '/path/to/scratch/amino.csv'

        sample = Sample(fastq1='/path/to/1234A_S1_R1_001.fastq.gz',
                        scratch_path=scratch_path,
                        nuc_csv='/path/to/nuc.csv')
        amino_path = sample.amino_csv

        self.assertEqual(expected_amino_path, amino_path)

    def test_default_fastq(self):
        scratch_path = '/path/to/scratch'
        expected_trimmed1_path = '/path/to/scratch/trimmed1.fastq'

        sample = Sample(fastq1='/path/to/1234A_S1_R1_001.fastq.gz',
                        scratch_path=scratch_path)
        trimmed1_path = sample.trimmed1_fastq

        self.assertEqual(expected_trimmed1_path, trimmed1_path)

    def test_extension_in_middle(self):
        scratch_path = '/path/to/scratch'
        expected_amino_path = '/path/to/scratch/amino_csvision.csv'

        sample = Sample(fastq1='/path/to/1234A_S1_R1_001.fastq.gz',
                        scratch_path=scratch_path,
                        nuc_csv='/path/to/nuc.csv')
        amino_path = sample.amino_csvision_csv

        self.assertEqual(expected_amino_path, amino_path)

    def test_unknown_extension(self):
        scratch_path = '/path/to/scratch'
        expected_path = '/path/to/scratch/example_stuff'

        sample = Sample(fastq1='/path/to/1234A_S1_R1_001.fastq.gz',
                        scratch_path=scratch_path,
                        nuc_csv='/path/to/nuc.csv')
        path = sample.example_stuff

        self.assertEqual(expected_path, path)

    def test_set_path(self):
        scratch_path = '/path/to/scratch'
        expected_path = '/path/to/elsewhere/example.csv'

        sample = Sample(fastq1='/path/to/1234A_S1_R1_001.fastq.gz',
                        scratch_path=scratch_path,
                        nuc_csv='/path/to/nuc.csv')
        sample.example_csv = expected_path
        path = sample.example_csv

        self.assertEqual(expected_path, path)

    def test_explicit_fastq2(self):
        expected_fastq1 = '/path/to/1234Ax_S1_R1_001.fastq.gz'
        expected_fastq2 = '/path/to/1234Ay_S1_R2_001.fastq.gz'

        sample = Sample(fastq1=expected_fastq1,
                        fastq2=expected_fastq2)

        self.assertEqual(expected_fastq1, sample.fastq1)
        self.assertEqual(expected_fastq2, sample.fastq2)

    def test_generate_fastq2(self):
        expected_fastq1 = '/path/to/1234A_S1_R1_001.fastq.gz'
        expected_fastq2 = '/path/to/1234A_S1_R2_001.fastq.gz'

        sample = Sample(fastq1=expected_fastq1)

        self.assertEqual(expected_fastq1, sample.fastq1)
        self.assertEqual(expected_fastq2, sample.fastq2)

    def test_cannot_generate_fastq2(self):
        with self.assertRaisesRegex(
                ValueError,
                r"fastq2 not given, and fastq1 does not contain '_R1_'."):
            Sample(fastq1='/path/to/1234A_S1.fastq.gz')

    def test_no_fastq(self):
        scratch_path = '/path/to/scratch'
        expected_fastq1 = '/path/to/scratch/fastq1'

        sample = Sample(amino_csv='/other/path/amino.csv',
                        scratch_path=scratch_path)

        self.assertEqual(expected_fastq1, sample.fastq1)

    def test_repr(self):
        expected_repr = "Sample(fastq1='1_R1_001.fastq')"
        sample = Sample(fastq1='1_R1_001.fastq',
                        amino_csv='/path/to_amino.csv')

        self.assertEqual(expected_repr, repr(sample))

    def test(self):
        expected_repr = "Sample()"
        sample = Sample(amino_csv='/path/to_amino.csv')

        self.assertEqual(expected_repr, repr(sample))

    def test_str(self):
        expected_str = "Sample 1234A-V3LOOP_S1"
        sample = Sample(fastq1='/path/to/1234A-V3LOOP_S1_L001_R1_001.fastq.gz')

        self.assertEqual(expected_str, str(sample))

    def test_str_without_fastq1(self):
        expected_str = "Sample"
        sample = Sample()

        self.assertEqual(expected_str, str(sample))

    def test_str_with_rank(self):
        expected_str = "Sample 1234A-V3LOOP_S1 (2 of 10)"
        sample = Sample(fastq1='/path/to/1234A-V3LOOP_S1_L001_R1_001.fastq.gz',
                        rank='2 of 10')

        self.assertEqual(expected_str, str(sample))

    def test_name(self):
        expected_sample_name = "1234A-V3LOOP_S1"

        sample = Sample(fastq1='/path/to/1234A-V3LOOP_S1_L001_R1_001.fastq.gz',
                        rank='2 of 10')
        sample_name = sample.name

        self.assertEqual(expected_sample_name, sample_name)

    def test_pickle(self):
        sample = Sample(fastq1='1A_R1_001.fastq',
                        scratch_path='/scratch')

        data = pickle.dumps(sample)

        self.assertNotEqual(b'', data)


@pytest.mark.parametrize(
    'project_code,excluded,expected',
    [('HIVGHA', (), ()),
     (None, (), ["HIV1-CRF06_CPX-GH-AB286851-seed",
                 "HIV1-CRF30_0206-GH-AB286854-seed"]),
     ('HIV', (), ["HIV1-CRF06_CPX-GH-AB286851-seed",
                  "HIV1-CRF30_0206-GH-AB286854-seed"]),
     ('HIV', ["HLA-B-seed"], ["HIV1-CRF06_CPX-GH-AB286851-seed",
                              "HIV1-CRF30_0206-GH-AB286854-seed",
                              "HLA-B-seed"])])
def test_exclude_extra_seeds(project_code, excluded, expected):
    all_excluded = exclude_extra_seeds(excluded, project_code)

    assert all_excluded == expected
