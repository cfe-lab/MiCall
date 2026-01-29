"""Tests for the list_fastq_files utility module."""

import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from micall.utils.list_fastq_files import (
    _get_base_calls_path,
    find_fastq_source_folder,
    list_fastq_files,
    list_fastq_file_names
)


class TestGetBaseCallsPath(unittest.TestCase):
    def test_get_base_calls_path(self):
        """Test that _get_base_calls_path returns the standard BaseCalls path."""
        run_path = Path('/path/to/run')
        expected = Path('/path/to/run/Data/Intensities/BaseCalls')
        result = _get_base_calls_path(run_path)
        self.assertEqual(expected, result)

    def test_get_base_calls_path_with_string(self):
        """Test that _get_base_calls_path works with string input."""
        run_path = '/path/to/run'
        expected = Path('/path/to/run/Data/Intensities/BaseCalls')
        result = _get_base_calls_path(run_path)
        self.assertEqual(expected, result)


class TestFindFastqSourceFolder(unittest.TestCase):
    def setUp(self):
        """Set up a temporary directory structure for testing."""
        self.temp_dir = TemporaryDirectory()
        self.run_path = Path(self.temp_dir.name) / 'test_run'
        self.run_path.mkdir()

    def tearDown(self):
        """Clean up the temporary directory."""
        self.temp_dir.cleanup()

    def test_find_in_base_calls(self):
        """Test finding FASTQ files in the standard BaseCalls folder."""
        base_calls = self.run_path / 'Data' / 'Intensities' / 'BaseCalls'
        base_calls.mkdir(parents=True)
        fastq_file = base_calls / 'Sample1_S1_L001_R1_001.fastq.gz'
        fastq_file.touch()

        result = find_fastq_source_folder(self.run_path, 'Sample1*_R1_*')
        self.assertEqual(base_calls, result)

    def test_find_in_alignment_folder(self):
        """Test finding FASTQ files in an Alignment folder."""
        alignment = self.run_path / 'Alignment_1' / 'L001' / 'Fastq'
        alignment.mkdir(parents=True)
        fastq_file = alignment / 'Sample1_S1_L001_R1_001.fastq.gz'
        fastq_file.touch()

        result = find_fastq_source_folder(self.run_path, 'Sample1*_R1_*')
        self.assertEqual(alignment, result)

    def test_prioritize_base_calls_over_alignment(self):
        """Test that BaseCalls is preferred over Alignment folders."""
        # Create files in both locations
        base_calls = self.run_path / 'Data' / 'Intensities' / 'BaseCalls'
        base_calls.mkdir(parents=True)
        (base_calls / 'Sample1_S1_L001_R1_001.fastq.gz').touch()

        alignment = self.run_path / 'Alignment_1' / 'L001' / 'Fastq'
        alignment.mkdir(parents=True)
        (alignment / 'Sample1_S1_L001_R1_001.fastq.gz').touch()

        result = find_fastq_source_folder(self.run_path, 'Sample1*_R1_*')
        self.assertEqual(base_calls, result)

    def test_prioritize_higher_alignment_number(self):
        """Test that higher alignment numbers are preferred."""
        # Create Alignment_5
        alignment5 = self.run_path / 'Alignment_5' / 'L001' / 'Fastq'
        alignment5.mkdir(parents=True)
        (alignment5 / 'Sample1_S1_L001_R1_001.fastq.gz').touch()

        # Create Alignment_20
        alignment20 = self.run_path / 'Alignment_20' / 'L001' / 'Fastq'
        alignment20.mkdir(parents=True)
        (alignment20 / 'Sample1_S1_L001_R1_001.fastq.gz').touch()

        result = find_fastq_source_folder(self.run_path, 'Sample1*_R1_*')
        self.assertEqual(alignment20, result)

    def test_prioritize_lexicographically_largest_subdir(self):
        """Test that lexicographically largest subdirectory is preferred within same alignment."""
        # Create Alignment_10/L001/Fastq
        alignment_l001 = self.run_path / 'Alignment_10' / 'L001' / 'Fastq'
        alignment_l001.mkdir(parents=True)
        (alignment_l001 / 'Sample1_S1_L001_R1_001.fastq.gz').touch()

        # Create Alignment_10/L002/Fastq
        alignment_l002 = self.run_path / 'Alignment_10' / 'L002' / 'Fastq'
        alignment_l002.mkdir(parents=True)
        (alignment_l002 / 'Sample1_S1_L001_R1_001.fastq.gz').touch()

        result = find_fastq_source_folder(self.run_path, 'Sample1*_R1_*')
        self.assertEqual(alignment_l002, result)

    def test_no_matching_files(self):
        """Test that None is returned when no matching files are found."""
        base_calls = self.run_path / 'Data' / 'Intensities' / 'BaseCalls'
        base_calls.mkdir(parents=True)
        # Create a file that doesn't match the pattern
        (base_calls / 'OtherSample_S1_L001_R1_001.fastq.gz').touch()

        result = find_fastq_source_folder(self.run_path, 'Sample1*_R1_*')
        self.assertIsNone(result)

    def test_fallback_to_run_path(self):
        """Test fallback to run_path when files are in the root."""
        fastq_file = self.run_path / 'Sample1_S1_L001_R1_001.fastq.gz'
        fastq_file.touch()

        result = find_fastq_source_folder(self.run_path, 'Sample1*_R1_*')
        self.assertEqual(self.run_path, result)

    def test_no_files_anywhere(self):
        """Test that None is returned when no files exist anywhere."""
        # Don't create any FASTQ files
        result = find_fastq_source_folder(self.run_path, 'Sample1*_R1_*')
        self.assertIsNone(result)


class TestListFastqFiles(unittest.TestCase):
    def setUp(self):
        """Set up a temporary directory structure for testing."""
        self.temp_dir = TemporaryDirectory()
        self.run_path = Path(self.temp_dir.name) / 'test_run'
        self.run_path.mkdir()

    def tearDown(self):
        """Clean up the temporary directory."""
        self.temp_dir.cleanup()

    def test_list_files_in_base_calls(self):
        """Test listing FASTQ files in BaseCalls folder."""
        base_calls = self.run_path / 'Data' / 'Intensities' / 'BaseCalls'
        base_calls.mkdir(parents=True)
        fastq1 = base_calls / 'Sample1_S1_L001_R1_001.fastq.gz'
        fastq2 = base_calls / 'Sample1_S1_L001_R2_001.fastq.gz'
        fastq1.touch()
        fastq2.touch()

        result = list_fastq_files(self.run_path, 'Sample1*_R1_*')
        self.assertEqual(1, len(result))
        self.assertEqual(fastq1, result[0])

    def test_list_files_in_alignment_folder(self):
        """Test listing FASTQ files in Alignment folder."""
        alignment = self.run_path / 'Alignment_10' / 'L001' / 'Fastq'
        alignment.mkdir(parents=True)
        fastq1 = alignment / 'Sample1_S1_L001_R1_001.fastq.gz'
        fastq2 = alignment / 'Sample1_S1_L001_R2_001.fastq.gz'
        fastq1.touch()
        fastq2.touch()

        result = list_fastq_files(self.run_path, 'Sample1*_R1_*')
        self.assertEqual(1, len(result))
        self.assertEqual(fastq1, result[0])

    def test_list_multiple_files(self):
        """Test listing multiple matching FASTQ files."""
        base_calls = self.run_path / 'Data' / 'Intensities' / 'BaseCalls'
        base_calls.mkdir(parents=True)
        fastq1 = base_calls / 'Sample1_S1_L001_R1_001.fastq.gz'
        fastq2 = base_calls / 'Sample1_S1_L002_R1_001.fastq.gz'
        fastq3 = base_calls / 'Sample2_S2_L001_R1_001.fastq.gz'
        fastq1.touch()
        fastq2.touch()
        fastq3.touch()

        result = list_fastq_files(self.run_path, 'Sample1*_R1_*')
        self.assertEqual(2, len(result))
        self.assertIn(fastq1, result)
        self.assertIn(fastq2, result)

    def test_empty_list_when_no_files(self):
        """Test that an empty list is returned when no matching files are found."""
        base_calls = self.run_path / 'Data' / 'Intensities' / 'BaseCalls'
        base_calls.mkdir(parents=True)

        result = list_fastq_files(self.run_path, 'Sample1*_R1_*', 
                                 fallback_to_run_path=False)
        self.assertEqual([], result)

    def test_fallback_to_run_path(self):
        """Test fallback to run_path when enabled."""
        fastq1 = self.run_path / 'Sample1_S1_L001_R1_001.fastq.gz'
        fastq1.touch()

        result = list_fastq_files(self.run_path, 'Sample1*_R1_*')
        self.assertEqual(1, len(result))
        self.assertEqual(fastq1, result[0])


class TestListFastqFileNames(unittest.TestCase):
    def setUp(self):
        """Set up a temporary directory structure for testing."""
        self.temp_dir = TemporaryDirectory()
        self.run_path = Path(self.temp_dir.name) / 'test_run'
        self.run_path.mkdir()

    def tearDown(self):
        """Clean up the temporary directory."""
        self.temp_dir.cleanup()

    def test_list_file_names(self):
        """Test that list_fastq_file_names returns only the filenames."""
        base_calls = self.run_path / 'Data' / 'Intensities' / 'BaseCalls'
        base_calls.mkdir(parents=True)
        fastq1 = base_calls / 'Sample1_S1_L001_R1_001.fastq.gz'
        fastq2 = base_calls / 'Sample1_S1_L002_R1_001.fastq.gz'
        fastq1.touch()
        fastq2.touch()

        result = list_fastq_file_names(self.run_path, 'Sample1*_R1_*')
        self.assertEqual(2, len(result))
        self.assertIn('Sample1_S1_L001_R1_001.fastq.gz', result)
        self.assertIn('Sample1_S1_L002_R1_001.fastq.gz', result)

    def test_file_names_from_alignment_folder(self):
        """Test getting file names from Alignment folder."""
        alignment = self.run_path / 'Alignment_15' / 'L001' / 'Fastq'
        alignment.mkdir(parents=True)
        fastq1 = alignment / 'Sample1_S1_L001_R1_001.fastq.gz'
        fastq1.touch()

        result = list_fastq_file_names(self.run_path, 'Sample1*_R1_*')
        self.assertEqual(['Sample1_S1_L001_R1_001.fastq.gz'], result)
