"""
Tests for exact_coverage integration in aln2counts.
These tests verify that the exact_coverage column is properly populated.
"""

import csv
from io import StringIO
import pytest

from micall.core.aln2counts import aln2counts

# Import fixture
from micall.tests.test_aln2counts_report import default_sequence_report  # noqa: F401


def test_exact_coverage_with_remap_conseq():
    """Test that exact_coverage column is populated when remap_conseq_csv is provided."""
    # Use a seed name that exists in the default project config
    seed_name = "HIV1-B-FR-K03455-seed"

    aligned_csv = StringIO(f"""\
refname,qcut,rank,count,offset,seq
{seed_name},15,0,5,0,AAATTTCCC
{seed_name},15,0,5,0,AAATTTCCC
{seed_name},15,0,5,0,AAATTTCCC
""")

    remap_conseq_csv = StringIO(f"""\
region,sequence
{seed_name},AAATTTCCC
""")

    nuc_csv = StringIO()
    amino_csv = StringIO()
    insertions_csv = StringIO()
    conseq_csv = StringIO()
    failed_align_csv = StringIO()
    coverage_summary_csv = StringIO()

    aln2counts(aligned_csv=aligned_csv,
               nuc_csv=nuc_csv,
               amino_csv=amino_csv,
               insertions_csv=insertions_csv,
               conseq_csv=conseq_csv,
               failed_align_csv=failed_align_csv,
               coverage_summary_csv=coverage_summary_csv,
               remap_conseq_csv=remap_conseq_csv)

    nuc_csv.seek(0)
    reader = csv.DictReader(nuc_csv)
    rows = list(reader)

    # Should have rows with exact_coverage values
    assert len(rows) > 0, "Should have nuc rows"

    # Check that exact_coverage column exists
    assert 'exact_coverage' in rows[0], "Should have exact_coverage column"

    # Check that at least some rows have non-empty exact_coverage
    exact_coverages = [row['exact_coverage'] for row in rows]
    non_empty = [ec for ec in exact_coverages if ec and ec.strip()]

    assert len(non_empty) > 0, f"Should have some non-empty exact_coverage values, got: {exact_coverages}"

    # Check that values are numeric
    for ec in non_empty:
        assert ec.isdigit(), f"exact_coverage should be numeric, got: {ec}"
        assert int(ec) > 0, f"exact_coverage should be positive, got: {ec}"


def test_exact_coverage_without_remap_conseq():
    """Test that exact_coverage column is empty when remap_conseq_csv is NOT provided."""
    # Use a known seed from projects
    aligned_csv = StringIO("""\
refname,qcut,rank,count,offset,seq
HIV1-B-FR-K03455-seed,15,0,5,0,AAATTT
""")

    nuc_csv = StringIO()
    amino_csv = StringIO()
    insertions_csv = StringIO()
    conseq_csv = StringIO()
    failed_align_csv = StringIO()
    coverage_summary_csv = StringIO()

    aln2counts(aligned_csv=aligned_csv,
               nuc_csv=nuc_csv,
               amino_csv=amino_csv,
               insertions_csv=insertions_csv,
               conseq_csv=conseq_csv,
               failed_align_csv=failed_align_csv,
               coverage_summary_csv=coverage_summary_csv,
               remap_conseq_csv=None)  # No remap_conseq_csv

    nuc_csv.seek(0)
    reader = csv.DictReader(nuc_csv)
    rows = list(reader)

    # Should have rows
    assert len(rows) > 0, "Should have nuc rows"

    # Check that exact_coverage column exists but is empty
    assert 'exact_coverage' in rows[0], "Should have exact_coverage column"

    # All exact_coverage values should be empty
    exact_coverages = [row['exact_coverage'] for row in rows]
    assert all(not ec or not ec.strip() for ec in exact_coverages), \
        f"exact_coverage should be empty without remap_conseq_csv, got: {exact_coverages}"


def test_exact_coverage_multiple_contigs():
    """Test exact_coverage with multiple contigs."""
    # Use two different HIV seeds
    seed1 = "HIV1-B-FR-K03455-seed"
    seed2 = "HIV1-CRF02_AG-GH-AB286855-seed"

    aligned_csv = StringIO(f"""\
refname,qcut,rank,count,offset,seq
{seed1},15,0,3,0,AAATTTCCC
{seed1},15,0,3,0,AAATTTCCC
{seed2},15,0,2,0,GGGCCCAAA
{seed2},15,0,2,0,GGGCCCAAA
""")

    remap_conseq_csv = StringIO(f"""\
region,sequence
{seed1},AAATTTCCC
{seed2},GGGCCCAAA
""")

    nuc_csv = StringIO()
    amino_csv = StringIO()
    insertions_csv = StringIO()
    conseq_csv = StringIO()
    failed_align_csv = StringIO()
    coverage_summary_csv = StringIO()

    aln2counts(aligned_csv=aligned_csv,
               nuc_csv=nuc_csv,
               amino_csv=amino_csv,
               insertions_csv=insertions_csv,
               conseq_csv=conseq_csv,
               failed_align_csv=failed_align_csv,
               coverage_summary_csv=coverage_summary_csv,
               remap_conseq_csv=remap_conseq_csv)

    nuc_csv.seek(0)
    reader = csv.DictReader(nuc_csv)
    rows = list(reader)

    # Group by seed (contig)
    by_seed = {}
    for row in rows:
        seed = row['seed']
        if seed not in by_seed:
            by_seed[seed] = []
        by_seed[seed].append(row)

    # Should have both contigs
    assert seed1 in by_seed, f"Should have {seed1}"
    assert seed2 in by_seed, f"Should have {seed2}"

    # Each contig should have some non-empty exact_coverage
    for seed in [seed1, seed2]:
        exact_coverages = [row['exact_coverage'] for row in by_seed[seed]]
        non_empty = [ec for ec in exact_coverages if ec and ec.strip()]
        assert len(non_empty) > 0, f"Contig {seed} should have non-empty exact_coverage"
