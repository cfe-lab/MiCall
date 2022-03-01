import io
import pickle
from unittest import TestCase

import pytest
import os
from unittest.mock import patch, DEFAULT
import logging
from io import StringIO

from micall.drivers.sample import Sample, exclude_extra_seeds, open_files


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
     (None, (), ["HIV1-CRF02_AG-GH-AB286855-seed",
                 "HIV1-CRF06_CPX-GH-AB286851-seed"]),
     ('HIV', (), ["HIV1-CRF02_AG-GH-AB286855-seed",
                  "HIV1-CRF06_CPX-GH-AB286851-seed"]),
     ('HIV', ["HLA-B-seed"], ["HIV1-CRF02_AG-GH-AB286855-seed",
                              "HIV1-CRF06_CPX-GH-AB286851-seed",
                              "HLA-B-seed"])])
def test_exclude_extra_seeds(project_code, excluded, expected):
    all_excluded = exclude_extra_seeds(excluded, project_code)

    assert all_excluded == expected


def test_context_manager(tmp_path):
    file_path = os.path.join(tmp_path, 'testfile.csv')
    file_info = (file_path, 'w')
    with open_files(testfile=file_info) as opened_files:
        file = opened_files['testfile']
        assert isinstance(file, io.TextIOBase)

    assert file.closed


def test_context_manager_write(tmp_path):
    file_path = os.path.join(tmp_path, 'testfile.csv')
    file_info = (file_path, 'w')
    with open_files(testfile=file_info) as opened_files:
        file = opened_files['testfile']
        file.write("Testing...")

    with open(file_path, 'r') as file_read:
        assert file_read.read() == "Testing..."


def mock_side_effect(_, mode):
    if mode == 'r':
        raise FileNotFoundError
    else:
        return DEFAULT


@patch('builtins.open')
def test_context_manager_fail_to_open(mock_open, tmp_path):
    # Check that the already opened files are closed if one fails to open.
    file_path1 = os.path.join(tmp_path, 'testfile.csv')
    file_path2 = os.path.join(tmp_path, 'testfile2.csv')
    file_info1 = (file_path1, 'w')
    file_info2 = (file_path2, 'r')
    mock_open.side_effect = mock_side_effect
    with pytest.raises(FileNotFoundError):
        with open_files(testfile1=file_info1, testfile2=file_info2) as _:
            pass
    assert mock_open.return_value.close.call_count == 1


@patch('builtins.open')
def test_context_manager_fail_to_close(mock_open, tmp_path, caplog):
    # Check that the file manager tries to close the other files if one fails to close.
    file_path1 = os.path.join(tmp_path, 'testfile.csv')
    file_path2 = os.path.join(tmp_path, 'testfile2.csv')
    file_info1 = (file_path1, 'w')
    file_info2 = (file_path2, 'w')
    expected_log = [
        ('micall.drivers.sample',
         logging.ERROR,
         "The following files could not be closed: ['testfile1', 'testfile2']")]
    mock_open.return_value.close.side_effect = IOError
    with pytest.raises(IOError):
        with open_files(testfile1=file_info1, testfile2=file_info2) as _:
            pass
    assert mock_open.return_value.close.call_count == 2
    assert caplog.record_tuples == expected_log
