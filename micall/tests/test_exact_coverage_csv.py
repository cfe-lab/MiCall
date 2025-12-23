"""
Tests for exact_coverage CSV input functionality.
"""
import csv
import tempfile
import unittest
from io import StringIO
from pathlib import Path

from micall.utils.exact_coverage import (
    calculate_exact_coverage_from_csv,
    read_aligned_csv,
    write_coverage_csv,
)


class TestReadAlignedCSV(unittest.TestCase):
    def test_read_aligned_csv_basic(self):
        """Test reading basic aligned CSV"""
        csv_data = StringIO("""\
refname,seq
1-HIV1-seed,ACGTACGT
1-HIV1-seed,GGGGCCCC
""")

        reads = list(read_aligned_csv(csv_data))

        self.assertEqual(len(reads), 2)
        self.assertEqual(reads[0], ('1-HIV1-seed', 'ACGTACGT'))
        self.assertEqual(reads[1], ('1-HIV1-seed', 'GGGGCCCC'))

    def test_read_aligned_csv_empty(self):
        """Test reading empty CSV"""
        csv_data = StringIO("refname,seq\n")

        reads = list(read_aligned_csv(csv_data))

        self.assertEqual(len(reads), 0)

    def test_read_aligned_csv_skip_empty_rows(self):
        """Test that rows with empty refname or seq are skipped"""
        csv_data = StringIO("""\
refname,seq
1-HIV1-seed,ACGTACGT
,GGGGCCCC
1-HIV1-seed,
1-HIV1-seed,TTTTAAAA
""")

        reads = list(read_aligned_csv(csv_data))

        self.assertEqual(len(reads), 2)
        self.assertEqual(reads[0], ('1-HIV1-seed', 'ACGTACGT'))
        self.assertEqual(reads[1], ('1-HIV1-seed', 'TTTTAAAA'))


class TestCalculateExactCoverageFromCSV(unittest.TestCase):
    def test_exact_coverage_from_csv_simple(self):
        """Test calculating exact coverage from CSV input"""
        aligned_csv = StringIO("""\
refname,seq
contig1,ACGTACGTACGT
contig1,TACGTACGTACG
""")

        contigs_csv = StringIO("""\
region,sequence
contig1,ACGTACGTACGTACGTACGTACGT
""")

        coverage, contigs = calculate_exact_coverage_from_csv(
            aligned_csv, contigs_csv, overlap_size=2
        )

        self.assertIn('contig1', coverage)
        self.assertEqual(len(coverage['contig1']), 24)
        # Read ACGTACGTACGT (12 bases) matches at position 0
        # With overlap_size=2, inner portion is positions 2-10
        for i in range(2, 10):
            self.assertGreater(coverage['contig1'][i], 0)

    def test_exact_coverage_from_csv_no_matches(self):
        """Test coverage when reads don't match contig"""
        aligned_csv = StringIO("""\
refname,seq
contig1,TTTTTTTTTTTT
""")

        contigs_csv = StringIO("""\
region,sequence
contig1,ACGTACGTACGT
""")

        coverage, contigs = calculate_exact_coverage_from_csv(
            aligned_csv, contigs_csv, overlap_size=2
        )

        self.assertIn('contig1', coverage)
        # No matches, all coverage should be 0
        for cov in coverage['contig1']:
            self.assertEqual(cov, 0)

    def test_exact_coverage_from_csv_reverse_complement(self):
        """Test that reverse complement matches are found"""
        aligned_csv = StringIO("""\
refname,seq
contig1,ACGTACGTACGT
""")

        # Contig is reverse complement of read
        contigs_csv = StringIO("""\
region,sequence
contig1,ACGTACGTACGT
""")

        coverage, contigs = calculate_exact_coverage_from_csv(
            aligned_csv, contigs_csv, overlap_size=2
        )

        self.assertIn('contig1', coverage)
        # Should find exact match
        for i in range(2, 10):
            self.assertGreater(coverage['contig1'][i], 0)

    def test_exact_coverage_from_csv_multiple_contigs(self):
        """Test coverage across multiple contigs"""
        aligned_csv = StringIO("""\
refname,seq
contig1,AAAAAAAA
contig2,GGGGGGGG
""")

        contigs_csv = StringIO("""\
region,sequence
contig1,AAAAAAAAAAAAAAAA
contig2,GGGGGGGGGGGGGGGG
""")

        coverage, contigs = calculate_exact_coverage_from_csv(
            aligned_csv, contigs_csv, overlap_size=1
        )

        self.assertIn('contig1', coverage)
        self.assertIn('contig2', coverage)

        # Both contigs should have some coverage
        self.assertGreater(sum(coverage['contig1']), 0)
        self.assertGreater(sum(coverage['contig2']), 0)


class TestIntegrationCSV(unittest.TestCase):
    def test_full_pipeline_csv_input(self):
        """Test full pipeline with CSV input"""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create test CSV files
            aligned_csv_path = Path(tmpdir) / "aligned.csv"
            contigs_csv_path = Path(tmpdir) / "contigs.csv"
            output_csv_path = Path(tmpdir) / "output.csv"

            # Write aligned CSV
            with open(aligned_csv_path, 'w') as f:
                f.write("refname,seq\n")
                f.write("1-HIV1-seed,ACGTACGTACGTACGTACGT\n")
                f.write("1-HIV1-seed,CGTACGTACGTACGTACGTA\n")

            # Write contigs CSV
            with open(contigs_csv_path, 'w') as f:
                f.write("region,sequence\n")
                f.write("1-HIV1-seed,ACGTACGTACGTACGTACGTACGTACGT\n")

            # Calculate coverage
            with open(aligned_csv_path, 'r') as aligned_f, \
                 open(contigs_csv_path, 'r') as contigs_f, \
                 open(output_csv_path, 'w') as output_f:

                coverage, contigs = calculate_exact_coverage_from_csv(
                    aligned_f, contigs_f, overlap_size=2
                )
                write_coverage_csv(coverage, contigs, output_f)

            # Verify output
            with open(output_csv_path, 'r') as f:
                reader = csv.DictReader(f)
                rows = list(reader)

            self.assertGreater(len(rows), 0)
            self.assertEqual(rows[0]['contig'], '1-HIV1-seed')

            # Check that some positions have coverage
            coverages = [int(row['exact_coverage']) for row in rows]
            self.assertGreater(sum(coverages), 0)


class TestCSVValidation(unittest.TestCase):
    def test_missing_refname_column(self):
        """Test that missing refname column raises ValueError"""
        csv_data = StringIO("""\
sequence,other
ACGTACGT,data
""")

        with self.assertRaises(ValueError) as ctx:
            list(read_aligned_csv(csv_data))

        self.assertIn("missing required columns", str(ctx.exception).lower())
        self.assertIn("refname", str(ctx.exception))

    def test_missing_seq_column(self):
        """Test that missing seq column raises ValueError"""
        csv_data = StringIO("""\
refname,other
contig1,data
""")

        with self.assertRaises(ValueError) as ctx:
            list(read_aligned_csv(csv_data))

        self.assertIn("missing required columns", str(ctx.exception).lower())
        self.assertIn("seq", str(ctx.exception))

    def test_missing_both_columns(self):
        """Test that missing both columns raises ValueError"""
        csv_data = StringIO("""\
other1,other2
data1,data2
""")

        with self.assertRaises(ValueError) as ctx:
            list(read_aligned_csv(csv_data))

        error_msg = str(ctx.exception).lower()
        self.assertIn("missing required columns", error_msg)
        self.assertIn("refname", str(ctx.exception))
        self.assertIn("seq", str(ctx.exception))

#    def test_no_header_row(self):
#        """Test that CSV without header raises ValueError"""
#        csv_data = StringIO("")
#
#        with self.assertRaises(ValueError) as ctx:
#            list(read_aligned_csv(csv_data))
#
#        self.assertIn("no header", str(ctx.exception).lower())

    def test_empty_refname_skipped(self):
        """Test that rows with empty refname are skipped"""
        csv_data = StringIO("""\
refname,seq
,ACGTACGT
contig2,GGGGCCCC
""")

        reads = list(read_aligned_csv(csv_data))

        self.assertEqual(len(reads), 1)
        self.assertEqual(reads[0], ('contig2', 'GGGGCCCC'))

    def test_empty_seq_skipped(self):
        """Test that rows with empty seq are skipped"""
        csv_data = StringIO("""\
refname,seq
contig1,
contig2,GGGGCCCC
""")

        reads = list(read_aligned_csv(csv_data))

        self.assertEqual(len(reads), 1)
        self.assertEqual(reads[0], ('contig2', 'GGGGCCCC'))

    def test_whitespace_trimmed(self):
        """Test that whitespace is trimmed from refname and seq"""
        csv_data = StringIO("""\
refname,seq
  contig1  ,  ACGTACGT
""")

        reads = list(read_aligned_csv(csv_data))

        self.assertEqual(len(reads), 1)
        self.assertEqual(reads[0], ('contig1', 'ACGTACGT'))

    def test_negative_overlap_size(self):
        """Test that negative overlap_size raises ValueError"""
        aligned_csv = StringIO("refname,seq\ncontig1,ACGT\n")
        contigs_csv = StringIO("region,sequence\ncontig1,ACGTACGT\n")

        with self.assertRaises(ValueError) as ctx:
            calculate_exact_coverage_from_csv(aligned_csv, contigs_csv, overlap_size=-1)

        self.assertIn("non-negative", str(ctx.exception))

    def test_empty_contigs_file(self):
        """Test that empty contigs file raises ValueError"""
        aligned_csv = StringIO("refname,seq\ncontig1,ACGT\n")
        contigs_csv = StringIO("region,sequence\n")

        with self.assertRaises(ValueError) as ctx:
            calculate_exact_coverage_from_csv(aligned_csv, contigs_csv, overlap_size=2)

        self.assertIn("no contigs", str(ctx.exception).lower())
