"""Tests for release_test_setup.py functionality."""

import unittest
from argparse import Namespace
from pathlib import Path
from tempfile import TemporaryDirectory

from micall.utils.release_test_setup import Sample


class TestSampleDirectoryReplication(unittest.TestCase):
    """Test that Sample class correctly replicates source directory structure."""

    def setUp(self):
        """Set up temporary directories for testing."""
        self.temp_dir = TemporaryDirectory()
        self.source_folder = Path(self.temp_dir.name) / 'source'
        self.test_folder = Path(self.temp_dir.name) / 'test'
        self.source_folder.mkdir()
        self.test_folder.mkdir()
        
        # Create a mock config object
        self.config = Namespace(
            test_folder=self.test_folder,
            pipeline_version='0-test',
            no_links=True  # Use copies for easier testing
        )

    def tearDown(self):
        """Clean up the temporary directories."""
        self.temp_dir.cleanup()

    def _create_source_run_with_base_calls(self, run_name):
        """Create a source run with files in the standard BaseCalls structure."""
        run_path = self.source_folder / 'MiSeq' / 'runs' / run_name
        base_calls = run_path / 'Data' / 'Intensities' / 'BaseCalls'
        base_calls.mkdir(parents=True)
        
        # Create FASTQ files
        (base_calls / 'Sample1_S1_L001_R1_001.fastq.gz').write_text('R1 data')
        (base_calls / 'Sample1_S1_L001_R2_001.fastq.gz').write_text('R2 data')
        
        # Create InterOp and metadata files
        interop = run_path / 'InterOp'
        interop.mkdir()
        (interop / 'dummy.bin').write_text('interop data')
        (run_path / 'RunInfo.xml').write_text('<RunInfo/>')
        (run_path / 'SampleSheet.csv').write_text('[Header]\n')
        (run_path / 'needsprocessing').write_text('')
        
        return run_path

    def _create_source_run_with_alignment(self, run_name, alignment_num=10):
        """Create a source run with files in an Alignment_X/*/Fastq structure."""
        run_path = self.source_folder / 'MiSeq' / 'runs' / run_name
        alignment = run_path / f'Alignment_{alignment_num}' / 'L001' / 'Fastq'
        alignment.mkdir(parents=True)
        
        # Create FASTQ files
        (alignment / 'Sample1_S1_L001_R1_001.fastq.gz').write_text('R1 data')
        (alignment / 'Sample1_S1_L001_R2_001.fastq.gz').write_text('R2 data')
        
        # Create InterOp and metadata files
        interop = run_path / 'InterOp'
        interop.mkdir()
        (interop / 'dummy.bin').write_text('interop data')
        (run_path / 'RunInfo.xml').write_text('<RunInfo/>')
        (run_path / 'SampleSheet.csv').write_text('[Header]\n')
        (run_path / 'needsprocessing').write_text('')
        
        return run_path, alignment

    def test_base_calls_structure_replicated(self):
        """Test that BaseCalls structure is replicated in target."""
        run_name = '230101_M01234_0001_000000000-ABCDE'
        self._create_source_run_with_base_calls(run_name)
        
        sample = Sample(run_name, 'Sample1', self.config)
        sample.find(str(self.source_folder))
        
        # Verify the source folder was found
        expected_source = Path(sample.run_name) / 'Data' / 'Intensities' / 'BaseCalls'
        self.assertEqual(expected_source, sample.source_fastq_folder)
        
        # Setup the run and samples
        sample.setup_run()
        sample_paths = sample.setup_samples()
        
        # Verify target structure matches source structure
        target_base_calls = (self.test_folder / 'MiSeq' / 'runs' / run_name / 
                            'Data' / 'Intensities' / 'BaseCalls')
        self.assertTrue(target_base_calls.exists(), 
                       f"Target BaseCalls folder should exist at {target_base_calls}")
        
        # Verify FASTQ files are in the correct location
        target_r1 = target_base_calls / 'Sample1_S1_L001_R1_001.fastq.gz'
        target_r2 = target_base_calls / 'Sample1_S1_L001_R2_001.fastq.gz'
        self.assertTrue(target_r1.exists(), f"R1 should exist at {target_r1}")
        self.assertTrue(target_r2.exists(), f"R2 should exist at {target_r2}")
        
        # Verify returned paths are correct
        self.assertEqual(2, len(sample_paths))
        self.assertIn(str(target_r1), sample_paths)
        self.assertIn(str(target_r2), sample_paths)

    def test_alignment_structure_replicated(self):
        """Test that Alignment_X/*/Fastq structure is replicated in target."""
        run_name = '230101_M01234_0002_000000000-ABCDE'
        _, alignment_path = self._create_source_run_with_alignment(run_name, alignment_num=15)
        
        sample = Sample(run_name, 'Sample1', self.config)
        sample.find(str(self.source_folder))
        
        # Verify the source folder was found
        expected_source = Path(sample.run_name) / 'Alignment_15' / 'L001' / 'Fastq'
        self.assertEqual(expected_source, sample.source_fastq_folder)
        
        # Setup the run and samples
        sample.setup_run()
        sample_paths = sample.setup_samples()
        
        # Verify target structure matches source structure (Alignment, not BaseCalls!)
        target_alignment = (self.test_folder / 'MiSeq' / 'runs' / run_name / 
                           'Alignment_15' / 'L001' / 'Fastq')
        self.assertTrue(target_alignment.exists(), 
                       f"Target Alignment folder should exist at {target_alignment}")
        
        # Verify FASTQ files are in the correct location
        target_r1 = target_alignment / 'Sample1_S1_L001_R1_001.fastq.gz'
        target_r2 = target_alignment / 'Sample1_S1_L001_R2_001.fastq.gz'
        self.assertTrue(target_r1.exists(), f"R1 should exist at {target_r1}")
        self.assertTrue(target_r2.exists(), f"R2 should exist at {target_r2}")
        
        # Verify returned paths are correct
        self.assertEqual(2, len(sample_paths))
        self.assertIn(str(target_r1), sample_paths)
        self.assertIn(str(target_r2), sample_paths)

    def test_higher_alignment_number_preferred(self):
        """Test that higher alignment numbers are preferred and replicated."""
        run_name = '230101_M01234_0003_000000000-ABCDE'
        run_path = self.source_folder / 'MiSeq' / 'runs' / run_name
        
        # Create Alignment_5
        alignment5 = run_path / 'Alignment_5' / 'L001' / 'Fastq'
        alignment5.mkdir(parents=True)
        (alignment5 / 'Sample1_S1_L001_R1_001.fastq.gz').write_text('R1 data v5')
        (alignment5 / 'Sample1_S1_L001_R2_001.fastq.gz').write_text('R2 data v5')
        
        # Create Alignment_20
        alignment20 = run_path / 'Alignment_20' / 'L001' / 'Fastq'
        alignment20.mkdir(parents=True)
        (alignment20 / 'Sample1_S1_L001_R1_001.fastq.gz').write_text('R1 data v20')
        (alignment20 / 'Sample1_S1_L001_R2_001.fastq.gz').write_text('R2 data v20')
        
        # Create InterOp and metadata files
        interop = run_path / 'InterOp'
        interop.mkdir()
        (interop / 'dummy.bin').write_text('interop data')
        (run_path / 'RunInfo.xml').write_text('<RunInfo/>')
        (run_path / 'SampleSheet.csv').write_text('[Header]\n')
        (run_path / 'needsprocessing').write_text('')
        
        sample = Sample(run_name, 'Sample1', self.config)
        sample.find(str(self.source_folder))
        
        # Verify Alignment_20 was chosen
        expected_source = Path(sample.run_name) / 'Alignment_20' / 'L001' / 'Fastq'
        self.assertEqual(expected_source, sample.source_fastq_folder)
        
        # Setup and verify the correct structure was replicated
        sample.setup_run()
        sample.setup_samples()
        
        target_alignment20 = (self.test_folder / 'MiSeq' / 'runs' / run_name / 
                             'Alignment_20' / 'L001' / 'Fastq')
        self.assertTrue(target_alignment20.exists())
        
        # Verify the content is from Alignment_20
        target_r1 = target_alignment20 / 'Sample1_S1_L001_R1_001.fastq.gz'
        self.assertEqual('R1 data v20', target_r1.read_text())

    def test_lexicographically_largest_subdir_preferred(self):
        """Test that lexicographically largest subdirectory is preferred within same alignment."""
        run_name = '230101_M01234_0004_000000000-ABCDE'
        run_path = self.source_folder / 'MiSeq' / 'runs' / run_name
        
        # Create Alignment_10/L001/Fastq
        alignment_l001 = run_path / 'Alignment_10' / 'L001' / 'Fastq'
        alignment_l001.mkdir(parents=True)
        (alignment_l001 / 'Sample1_S1_L001_R1_001.fastq.gz').write_text('R1 L001')
        (alignment_l001 / 'Sample1_S1_L001_R2_001.fastq.gz').write_text('R2 L001')
        
        # Create Alignment_10/L002/Fastq
        alignment_l002 = run_path / 'Alignment_10' / 'L002' / 'Fastq'
        alignment_l002.mkdir(parents=True)
        (alignment_l002 / 'Sample1_S1_L001_R1_001.fastq.gz').write_text('R1 L002')
        (alignment_l002 / 'Sample1_S1_L001_R2_001.fastq.gz').write_text('R2 L002')
        
        # Create InterOp and metadata files
        interop = run_path / 'InterOp'
        interop.mkdir()
        (interop / 'dummy.bin').write_text('interop data')
        (run_path / 'RunInfo.xml').write_text('<RunInfo/>')
        (run_path / 'SampleSheet.csv').write_text('[Header]\n')
        (run_path / 'needsprocessing').write_text('')
        
        sample = Sample(run_name, 'Sample1', self.config)
        sample.find(str(self.source_folder))
        
        # Verify L002 was chosen (lexicographically larger)
        expected_source = Path(sample.run_name) / 'Alignment_10' / 'L002' / 'Fastq'
        self.assertEqual(expected_source, sample.source_fastq_folder)
        
        # Setup and verify the correct structure was replicated
        sample.setup_run()
        sample.setup_samples()
        
        target_alignment_l002 = (self.test_folder / 'MiSeq' / 'runs' / run_name / 
                                'Alignment_10' / 'L002' / 'Fastq')
        self.assertTrue(target_alignment_l002.exists())
        
        # Verify the content is from L002
        target_r1 = target_alignment_l002 / 'Sample1_S1_L001_R1_001.fastq.gz'
        self.assertEqual('R1 L002', target_r1.read_text())
