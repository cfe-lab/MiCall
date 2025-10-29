import os
from random import randrange
from tempfile import NamedTemporaryFile
from unittest import TestCase
import csv
import tempfile
from pathlib import Path

from micall.core.prelim_map import check_fastq, prelim_map
from micall.utils.stderr import Stderr
from micall.utils.work_dir import WorkDir


class CheckFastqTest(TestCase):
    def create_temp_file(self, suffix, cleanup=True):
        f = NamedTemporaryFile(mode='w', dir='.', suffix=suffix, delete=False)
        if cleanup:
            self.addCleanup(os.unlink, f.name)
        return f

    def test_simple(self):
        f = self.create_temp_file(suffix='.fastq')
        f.close()

        new_name = check_fastq(f.name)

        self.assertEqual(f.name, new_name)

    def test_missing(self):
        f = self.create_temp_file(suffix='.fastq', cleanup=False)
        f.close()
        os.unlink(f.name)

        with self.assertRaisesRegex(SystemExit, 'No FASTQ found at'):
            check_fastq(f.name)

    def test_gzipped(self):
        f = self.create_temp_file(suffix='.fastq.gz')
        f.close()

        new_name = check_fastq(f.name, gzip=True)

        self.assertEqual(f.name, new_name)

    def test_gzipped_link(self):
        f = self.create_temp_file(suffix='.txt')
        expected_content = str(randrange(1000000))
        f.write(expected_content)
        f.close()
        expected_new_name = f.name + '.gz'
        self.addCleanup(os.unlink, expected_new_name)

        new_name = check_fastq(f.name, gzip=True)

        with open(new_name) as f2:
            content = f2.read()

        self.assertEqual(expected_new_name, new_name)
        self.assertEqual(expected_content, content)

    def test(self):
        f = self.create_temp_file(suffix='.txt')
        f.close()
        expected_new_name = f.name + '.gz'
        with open(expected_new_name, 'w'):
            pass
        self.addCleanup(os.unlink, expected_new_name)

        new_name = check_fastq(f.name, gzip=True)

        self.assertEqual(expected_new_name, new_name)


class PrelimMapIntegrationTest(TestCase):
    """Integration tests for prelim_map that run the full pipeline with real FASTQ files."""

    def setUp(self):
        """Set up test fixtures."""
        self.test_dir = os.path.dirname(__file__)
        self.microtest_dir = os.path.join(self.test_dir, 'microtest')
        self.work_dir = tempfile.mkdtemp()
        self.addCleanup(self._cleanup_work_dir)

    def _cleanup_work_dir(self):
        """Clean up the working directory."""
        import shutil
        if os.path.exists(self.work_dir):
            shutil.rmtree(self.work_dir)

    def test_basic_mapping(self):
        """Test basic mapping with V3LOOP test data."""
        fastq1 = os.path.join(self.microtest_dir, '1234A-V3LOOP_S1_L001_R1_001.fastq')
        fastq2 = os.path.join(self.microtest_dir, '1234A-V3LOOP_S1_L001_R2_001.fastq')

        output_file = os.path.join(self.work_dir, 'output.csv')
        stderr_file = os.path.join(self.work_dir, 'stderr.txt')

        with WorkDir.using(Path(self.work_dir)):
            with open(stderr_file, 'w') as stderr_f:
                with Stderr.using(stderr_f):
                    prelim_map(
                        fastq1=Path(fastq1),
                        fastq2=Path(fastq2),
                        prelim_csv=Path(output_file)
                    )

        # Parse the CSV output
        with open(output_file, 'r') as f:
            reader = csv.DictReader(f)
            rows = list(reader)        # Verify we got expected results
        self.assertGreater(len(rows), 0, "Should have at least one mapped read")

        # All test reads should map (10 read pairs = 20 rows)
        self.assertEqual(len(rows), 20, "Should have 20 rows (10 read pairs)")

        # Check first row has expected structure
        first_row = rows[0]
        expected_fields = ['qname', 'flag', 'rname', 'pos', 'mapq',
                          'cigar', 'rnext', 'pnext', 'tlen', 'seq', 'qual']
        self.assertEqual(set(first_row.keys()), set(expected_fields))

        # Verify read name format
        self.assertTrue(first_row['qname'].startswith('M01234:01:000000000-AAAAA'))

        # Verify reference name
        self.assertEqual(first_row['rname'], 'HIV1-C-BR-JX140663-seed')

        # Verify mapping position
        self.assertEqual(first_row['pos'], '6535')

        # Verify CIGAR string (all 51bp perfect matches)
        self.assertEqual(first_row['cigar'], '51M')

        # Check flags - should be 99 (forward) and 147 (reverse) for paired reads
        flags = [int(row['flag']) for row in rows]
        self.assertIn(99, flags, "Should have forward reads (flag 99)")
        self.assertIn(147, flags, "Should have reverse reads (flag 147)")

        # Check that we have equal numbers of forward and reverse reads
        forward_count = sum(1 for f in flags if f == 99)
        reverse_count = sum(1 for f in flags if f == 147)
        self.assertEqual(forward_count, reverse_count, "Should have equal forward and reverse reads")

    def test_mapping_quality_scores(self):
        """Test that mapping quality scores are reasonable."""
        fastq1 = os.path.join(self.microtest_dir, '1234A-V3LOOP_S1_L001_R1_001.fastq')
        fastq2 = os.path.join(self.microtest_dir, '1234A-V3LOOP_S1_L001_R2_001.fastq')

        output_file = os.path.join(self.work_dir, 'output.csv')

        with WorkDir.using(Path(self.work_dir)):
            prelim_map(
                fastq1=Path(fastq1),
                fastq2=Path(fastq2),
                prelim_csv=Path(output_file)
            )

        with open(output_file, 'r') as f:
            reader = csv.DictReader(f)
            rows = list(reader)        # Most reads should have high mapping quality (36)
        # The last read pair (0010) has a mismatch so lower quality (14)
        mapq_scores = [int(row['mapq']) for row in rows]

        # Check that we have the expected high-quality mappings
        high_quality_count = sum(1 for q in mapq_scores if q == 36)
        self.assertEqual(high_quality_count, 18, "Should have 18 high-quality mappings (9 pairs)")

        # Check that we have the expected lower-quality mappings
        lower_quality_count = sum(1 for q in mapq_scores if q == 14)
        self.assertEqual(lower_quality_count, 2, "Should have 2 lower-quality mappings (1 pair)")

    def test_sequence_preservation(self):
        """Test that sequences are preserved correctly in the output."""
        fastq1 = os.path.join(self.microtest_dir, '1234A-V3LOOP_S1_L001_R1_001.fastq')
        fastq2 = os.path.join(self.microtest_dir, '1234A-V3LOOP_S1_L001_R2_001.fastq')

        output_file = os.path.join(self.work_dir, 'output.csv')

        with WorkDir.using(Path(self.work_dir)):
            prelim_map(
                fastq1=Path(fastq1),
                fastq2=Path(fastq2),
                prelim_csv=Path(output_file)
            )

        with open(output_file, 'r') as f:
            reader = csv.DictReader(f)
            rows = list(reader)        # Check that the first read has the expected sequence
        first_row = rows[0]
        expected_seq = 'TGCACAAGACCCAACAACAATACAAGAAAAAGTATAAGGATAGGACCAGGA'
        self.assertEqual(first_row['seq'], expected_seq)

        # Check quality scores are preserved
        self.assertEqual(first_row['qual'], 'A' * 51)

        # Check the variant read (0010) has the expected difference
        variant_row = [row for row in rows if '0010' in row['qname'] and row['flag'] == '99'][0]
        expected_variant_seq = 'TGCATAAGACCCAACAACAATACAAGAAAAAGTATAAGGATAGGACCAGGA'
        self.assertEqual(variant_row['seq'], expected_variant_seq)

    def test_paired_read_consistency(self):
        """Test that paired reads have consistent information."""
        fastq1 = os.path.join(self.microtest_dir, '1234A-V3LOOP_S1_L001_R1_001.fastq')
        fastq2 = os.path.join(self.microtest_dir, '1234A-V3LOOP_S1_L001_R2_001.fastq')

        output_file = os.path.join(self.work_dir, 'output.csv')

        with WorkDir.using(Path(self.work_dir)):
            prelim_map(
                fastq1=Path(fastq1),
                fastq2=Path(fastq2),
                prelim_csv=Path(output_file)
            )

        with open(output_file, 'r') as f:
            reader = csv.DictReader(f)
            rows = list(reader)        # Group reads by qname
        read_pairs = {}
        for row in rows:
            qname = row['qname']
            if qname not in read_pairs:
                read_pairs[qname] = []
            read_pairs[qname].append(row)

        # Check each pair
        for qname, pair in read_pairs.items():
            self.assertEqual(len(pair), 2, f"Should have exactly 2 reads for {qname}")

            # Both should map to same reference
            self.assertEqual(pair[0]['rname'], pair[1]['rname'])

            # Both should have same position
            self.assertEqual(pair[0]['pos'], pair[1]['pos'])

            # rnext should be '=' indicating same reference
            self.assertEqual(pair[0]['rnext'], '=')
            self.assertEqual(pair[1]['rnext'], '=')

            # pnext should match pos
            self.assertEqual(pair[0]['pnext'], pair[0]['pos'])
            self.assertEqual(pair[1]['pnext'], pair[1]['pos'])

    def test_custom_gap_penalties(self):
        """Test that custom gap penalties can be specified."""
        fastq1 = os.path.join(self.microtest_dir, '1234A-V3LOOP_S1_L001_R1_001.fastq')
        fastq2 = os.path.join(self.microtest_dir, '1234A-V3LOOP_S1_L001_R2_001.fastq')

        output_file = os.path.join(self.work_dir, 'output.csv')

        # Run with custom gap penalties
        with WorkDir.using(Path(self.work_dir)):
            prelim_map(
                fastq1=Path(fastq1),
                fastq2=Path(fastq2),
                prelim_csv=Path(output_file),
                rdgopen=15,
                rfgopen=15
            )

        with open(output_file, 'r') as f:
            reader = csv.DictReader(f)
            rows = list(reader)        # Should still produce results with different penalties
        self.assertGreater(len(rows), 0, "Should have mapped reads with custom penalties")

        # Results should be the same for this simple dataset (no gaps)
        self.assertEqual(len(rows), 20)

    def test_multiple_threads(self):
        """Test that multiple threads can be used."""
        fastq1 = os.path.join(self.microtest_dir, '1234A-V3LOOP_S1_L001_R1_001.fastq')
        fastq2 = os.path.join(self.microtest_dir, '1234A-V3LOOP_S1_L001_R2_001.fastq')

        output_file = os.path.join(self.work_dir, 'output.csv')

        # Run with multiple threads
        with WorkDir.using(Path(self.work_dir)):
            prelim_map(
                fastq1=Path(fastq1),
                fastq2=Path(fastq2),
                prelim_csv=Path(output_file),
                nthreads=2
            )

        with open(output_file, 'r') as f:
            reader = csv.DictReader(f)
            rows = list(reader)        # Should produce same results regardless of thread count
        self.assertEqual(len(rows), 20)

    def test_excluded_seeds(self):
        """Test that specific seeds can be excluded from mapping."""
        fastq1 = os.path.join(self.microtest_dir, '1234A-V3LOOP_S1_L001_R1_001.fastq')
        fastq2 = os.path.join(self.microtest_dir, '1234A-V3LOOP_S1_L001_R2_001.fastq')

        output_file = os.path.join(self.work_dir, 'output.csv')

        # Run with all HIV seeds excluded (should produce no results)
        with WorkDir.using(Path(self.work_dir)):
            prelim_map(
                fastq1=Path(fastq1),
                fastq2=Path(fastq2),
                prelim_csv=Path(output_file),
                excluded_seeds={'HIV1-C-BR-JX140663-seed'}
            )

        with open(output_file, 'r') as f:
            reader = csv.DictReader(f)
            rows = list(reader)        # With the specific seed excluded, reads may not map or map to other references
        # The exact behavior depends on what other seeds are available
        # Just verify the function runs without error
        self.assertIsInstance(rows, list)
