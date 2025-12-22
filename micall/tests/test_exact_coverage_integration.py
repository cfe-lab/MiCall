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

        with (
            open(fastq1, "r") as f1,
            open(fastq2, "r") as f2,
            open(contig_fasta, "r") as fc,
            open(output_csv, "w") as fo,
        ):
            coverage, contigs = calculate_exact_coverage(
                f1, f2, fc, kmer_size=4, min_overlap=4
            )
            write_coverage_csv(coverage, contigs, fo)

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

        return True
