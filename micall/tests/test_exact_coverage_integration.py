"""
Integration test for exact_coverage tool.
Creates sample data and tests the complete workflow.
"""

import tempfile
import csv
from pathlib import Path


def create_test_data(tmp_dir):
    """Create test FASTQ and FASTA files."""

    # Create a simple contig
    contig_fasta = tmp_dir / "contigs.fasta"
    contig_fasta.write_text(""">contig1
ACGTACGTACGTACGTACGT
>contig2
GGGGCCCCGGGGCCCC
""")

    # Create FASTQ files with reads that match exactly
    # Read pairs that match contig1
    fastq1 = tmp_dir / "reads_R1.fastq"
    fastq1.write_text("""\
@read1
ACGTACGT
+
IIIIIIII
@read2
GTACGTAC
+
IIIIIIII
@read3
GGGGCCCC
+
IIIIIIII
""")

    fastq2 = tmp_dir / "reads_R2.fastq"
    fastq2.write_text("""\
@read1
ACGTACGT
+
IIIIIIII
@read2
GTACGTAC
+
IIIIIIII
@read3
CCCCGGGG
+
IIIIIIII
""")

    return contig_fasta, fastq1, fastq2


def test_exact_coverage_integration():
    """Test the complete exact_coverage workflow."""

    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_path = Path(tmp_dir)

        # Create test data
        contig_fasta, fastq1, fastq2 = create_test_data(tmp_path)
        output_csv = tmp_path / "output.csv"

        # Run the tool
        from micall.utils.exact_coverage import (
            calculate_exact_coverage,
            write_coverage_csv,
        )

        with open(contig_fasta, "r") as fc, open(output_csv, "w") as fo:
            coverage0, coverage, contigs = calculate_exact_coverage(
                fastq1, fastq2, fc, overlap_size=2
            )
            write_coverage_csv(coverage0, coverage, contigs, fo)

        # Verify results
        with open(output_csv, "r") as f:
            reader = csv.DictReader(f)
            rows = list(reader)

        # Check structure
        assert len(rows) > 0, "Output CSV should have rows"
        assert rows[0]["contig"] in ["contig1", "contig2"], "Should have contig names"
        assert "position" in rows[0], "Should have position column"
        assert "exact_coverage" in rows[0], "Should have exact_coverage column"
        # base column should NOT be present
        assert "base" not in rows[0], "Should NOT have base column"

        # Check that we got coverage for both contigs
        contigs_in_output = set(row["contig"] for row in rows)
        assert "contig1" in contigs_in_output, "Should have contig1 in output"
        assert "contig2" in contigs_in_output, "Should have contig2 in output"

        # Check that some positions have non-zero coverage
        coverages = [int(row["exact_coverage"]) for row in rows]
        assert any(c > 0 for c in coverages), "Should have some non-zero coverage"

        print("✓ Integration test passed!")
        print(f"✓ Total positions: {len(rows)}")
        print(f"✓ Positions with coverage > 0: {sum(1 for c in coverages if c > 0)}")
        print(f"✓ Max coverage: {max(coverages)}")
        print(f"✓ Average coverage: {sum(coverages) / len(coverages):.2f}")

        # Print first few rows for inspection
        print("\nFirst 10 rows of output:")
        print("contig,position,exact_coverage")
        for row in rows[:10]:
            print(f"{row['contig']},{row['position']},{row['exact_coverage']}")


def test_exact_coverage_with_csv_contigs():
    """Integration test with CSV contigs file (contigs.csv format)"""
    with tempfile.TemporaryDirectory() as tmpdir:
        # Create test files
        fastq1, fastq2, contig_csv, output_csv = create_test_data_contigs_csv(tmpdir)

        # Run exact coverage
        from micall.utils.exact_coverage import (
            calculate_exact_coverage,
            write_coverage_csv,
        )

        with open(contig_csv, "r") as fc, open(output_csv, "w") as fo:
            coverage0, coverage, contigs = calculate_exact_coverage(
                Path(fastq1), Path(fastq2), fc, overlap_size=2
            )
            write_coverage_csv(coverage0, coverage, contigs, fo)

        # Verify results
        with open(output_csv, "r") as f:
            reader = csv.DictReader(f)
            rows = list(reader)

        # Check structure
        assert len(rows) > 0, "Output CSV should have rows"
        # Should use 'ref' column names since CSV has 'ref' column
        # (priority: region > ref > sample)
        assert any(row["contig"] == "ref1" for row in rows), (
            "Should have ref1 contig name from 'ref' column"
        )
        assert any(row["contig"] == "ref2" for row in rows), (
            "Should have ref2 contig name from 'ref' column"
        )
        assert "position" in rows[0], "Should have position column"
        assert "exact_coverage" in rows[0], "Should have exact_coverage column"
        assert "base" not in rows[0], "Should NOT have base column"

        # Check that some positions have non-zero coverage
        coverages = [int(row["exact_coverage"]) for row in rows]
        assert any(c > 0 for c in coverages), "Should have some non-zero coverage"

        print("✓ CSV contigs integration test passed!")
        print(f"✓ Total positions: {len(rows)}")
        print(f"✓ Positions with coverage > 0: {sum(1 for c in coverages if c > 0)}")


def test_exact_coverage_with_conseq_csv():
    """Integration test with CSV conseq file (conseq.csv format)"""
    with tempfile.TemporaryDirectory() as tmpdir:
        # Create test files
        fastq1, fastq2, conseq_csv, output_csv = create_test_data_conseq_csv(tmpdir)

        # Run exact coverage
        from micall.utils.exact_coverage import (
            calculate_exact_coverage,
            write_coverage_csv,
        )

        with open(conseq_csv, "r") as fc, open(output_csv, "w") as fo:
            coverage0, coverage, contigs = calculate_exact_coverage(
                Path(fastq1), Path(fastq2), fc, overlap_size=2
            )
            write_coverage_csv(coverage0, coverage, contigs, fo)

        # Verify results
        with open(output_csv, "r") as f:
            reader = csv.DictReader(f)
            rows = list(reader)

        # Check structure
        assert len(rows) > 0, "Output CSV should have rows"
        # Should use 'region' or 'sample' column as contig name
        assert any(
            row["contig"] in ["region1", "region2", "sample1"] for row in rows
        ), "Should have contig names from CSV"
        assert "position" in rows[0], "Should have position column"
        assert "exact_coverage" in rows[0], "Should have exact_coverage column"

        # Check that some positions have non-zero coverage
        coverages = [int(row["exact_coverage"]) for row in rows]
        assert any(c > 0 for c in coverages), "Should have some non-zero coverage"

        print("✓ conseq.csv format integration test passed!")
        print(f"✓ Total positions: {len(rows)}")
        print(f"✓ Positions with coverage > 0: {sum(1 for c in coverages if c > 0)}")


def create_test_data_contigs_csv(tmpdir):
    """
    Create test data files for exact coverage with contigs.csv format.
    Returns paths to fastq1, fastq2, contigs_csv, output_csv.
    """
    tmpdir = Path(tmpdir)

    # Create FASTQ files
    fastq1 = tmpdir / "test_R1.fastq"
    fastq2 = tmpdir / "test_R2.fastq"

    # Simple test data: two contigs with overlapping reads
    fastq1.write_text("""\
@read1
ACGTACGTACGT
+
IIIIIIIIIIII
@read2
GGGGCCCCTTTT
+
IIIIIIIIIIII
""")

    fastq2.write_text("""\
@read1
TTTTTTTTTTTT
+
IIIIIIIIIIII
@read2
AAAAAAAAAAA
+
IIIIIIIIIII
""")

    # Create CSV contigs file (contigs.csv format with 'contig' column)
    contig_csv = tmpdir / "contigs.csv"
    contig_csv.write_text("""\
ref,match,group_ref,contig
ref1,1.0,group1,ACGTACGTACGTACGTACGTACGTACGTACGT
ref2,1.0,group2,GGGGCCCCTTTTAAAACCCCGGGGTTTTAAAA
""")

    output_csv = tmpdir / "exact_coverage_csv.csv"

    return str(fastq1), str(fastq2), str(contig_csv), str(output_csv)


def create_test_data_conseq_csv(tmpdir):
    """
    Create test data files for exact coverage with conseq.csv format.
    Returns paths to fastq1, fastq2, conseq_csv, output_csv.
    """
    tmpdir = Path(tmpdir)

    # Create FASTQ files
    fastq1 = tmpdir / "test_R1.fastq"
    fastq2 = tmpdir / "test_R2.fastq"

    # Simple test data
    fastq1.write_text("""\
@read1
ACGTACGTACGT
+
IIIIIIIIIIII
@read2
GGGGCCCCTTTT
+
IIIIIIIIIIII
""")

    fastq2.write_text("""\
@read1
TTTTTTTTTTTT
+
IIIIIIIIIIII
@read2
AAAAAAAAAAA
+
IIIIIIIIIII
""")

    # Create CSV conseq file (conseq.csv format with 'sequence' column)
    conseq_csv = tmpdir / "conseq.csv"
    conseq_csv.write_text("""\
sample,region,q-cutoff,consensus-percent-cutoff,offset,sequence
sample1,region1,15,MAX,0,ACGTACGTACGTACGTACGTACGTACGTACGT
sample1,region2,15,MAX,0,GGGGCCCCTTTTAAAACCCCGGGGTTTTAAAA
""")

    output_csv = tmpdir / "exact_coverage_conseq.csv"

    return str(fastq1), str(fastq2), str(conseq_csv), str(output_csv)
