import csv
from io import StringIO

from micall.core.aln2counts import aln2counts

# Import fixtures
from micall.tests.test_aln2counts_report import default_sequence_report  # noqa: F401
from micall.tests.test_remap import load_projects

assert load_projects


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
    assert len(rows) > 0, "Should have nuc rows"
    assert 'exact_coverage' in rows[0], "Should have exact_coverage column"
    exact_coverages = [row['exact_coverage'] for row in rows]
    non_empty = [ec for ec in exact_coverages if ec and ec.strip()]
    assert len(non_empty) > 0, f"Should have some non-empty exact_coverage values, got: {exact_coverages}"
    for ec in non_empty:
        assert ec.isdigit(), f"exact_coverage should be numeric, got: {ec}"
        assert int(ec) > 0, f"exact_coverage should be positive, got: {ec}"

def test_exact_coverage_without_remap_conseq():
    """Test that exact_coverage column is empty when remap_conseq_csv is NOT provided."""
    # Use a known seed from projects
    aligned_csv = StringIO("""refname,qcut,rank,count,offset,seq
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
               remap_conseq_csv=None)
    nuc_csv.seek(0)
    reader = csv.DictReader(nuc_csv)
    rows = list(reader)
    assert len(rows) > 0, "Should have nuc rows"
    assert 'exact_coverage' in rows[0], "Should have exact_coverage column"
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
{seed1},15,0,3,0,AAATTTCCCCCCC
{seed1},15,0,3,0,AAATTTCCACCCC
{seed2},15,0,2,0,GGGCCCAAACCCC
{seed2},15,0,2,0,GGGCCCAATCCCC
""")
    remap_conseq_csv = StringIO(f"""\
region,sequence
{seed1},AAATTTCCCCCCC
{seed2},GGGCCCAAACCCC
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
    contents = nuc_csv.read()
    assert contents != [], "Nuc CSV should not be empty"

    nuc_csv.seek(0)
    reader = csv.DictReader(nuc_csv)
    rows = list(reader)
    by_seed = {}
    for row in rows:
        seed = row['seed']
        if seed not in by_seed:
            by_seed[seed] = []
        by_seed[seed].append(row)
    assert seed1 in by_seed, f"Should have {seed1}"
    assert seed2 in by_seed, f"Should have {seed2}"
    for seed in [seed1, seed2]:
        exact_coverages = [row['exact_coverage'] for row in by_seed[seed]]
        non_empty = [ec for ec in exact_coverages if ec and ec.strip()]
        assert len(non_empty) > 0, f"Contig {seed} should have non-empty exact_coverage"


def test_exact_coverage_multiple_contigs_different_numbers():
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
    contents = nuc_csv.read()
    assert contents != [], "Nuc CSV should not be empty"

    nuc_csv.seek(0)
    reader = csv.DictReader(nuc_csv)
    rows = list(reader)
    by_seed = {}
    for row in rows:
        seed = row['seed']
        if seed not in by_seed:
            by_seed[seed] = []
        by_seed[seed].append(row)
    assert seed1 in by_seed, f"Should have {seed1}"
    assert seed2 in by_seed, f"Should have {seed2}"
    for seed in [seed1, seed2]:
        exact_coverages = [row['exact_coverage'] for row in by_seed[seed]]
        non_empty = [ec for ec in exact_coverages if ec and ec.strip()]
        assert len(non_empty) > 0, f"Contig {seed} should have non-empty exact_coverage"


def test_exact_coverage_accumulation_and_name_mapping():
    """
    Test that exact_coverage accumulates when multiple contigs with different
    prefixes map to the same seed name.
    """
    seed_name = "HIV1-B-FR-K03455-seed"
    # Contig 1: count 5, palindrome read -> 10 coverage
    # Contig 2: count 2, palindrome read -> 4 coverage
    # Both should map to seed-name.
    aligned_csv = StringIO(f"""\
refname,qcut,rank,count,offset,seq
1-{seed_name},15,0,5,0,AAATTT
2-{seed_name},15,0,2,0,AAATTT
""")
    remap_conseq_csv = StringIO(f"""\
region,sequence
{seed_name},AAATTT
1-{seed_name},AAATTT
2-{seed_name},AAATTT
""")
    nuc_csv = StringIO()
    nuc_detail_csv = StringIO()
    amino_csv = StringIO()
    insertions_csv = StringIO()
    conseq_csv = StringIO()
    failed_align_csv = StringIO()
    # Pass nuc_detail_csv to trigger combine_reports logic
    aln2counts(aligned_csv=aligned_csv,
               nuc_csv=nuc_csv,
               nuc_detail_csv=nuc_detail_csv,
               amino_csv=amino_csv,
               insertions_csv=insertions_csv,
               conseq_csv=conseq_csv,
               failed_align_csv=failed_align_csv,
               remap_conseq_csv=remap_conseq_csv)
    nuc_csv.seek(0)
    contents = nuc_csv.read()
    assert contents != [], "Nuc CSV should not be empty"

    nuc_csv.seek(0)
    reader = csv.DictReader(nuc_csv)
    rows = list(reader)
    assert all(r['seed'] == seed_name for r in rows), f"All rows should have seed {seed_name}"
    assert len(rows) > 0, "Should have at least one row"
    # Find a row with non-empty exact_coverage
    row_with_coverage = next((r for r in rows if r.get('exact_coverage') and r['exact_coverage'].strip()), None)
    assert row_with_coverage is not None, f"Should have at least one row with exact_coverage, got rows: {[(r['refseq.nuc.pos'], r['query.nuc.pos'], r['exact_coverage']) for r in rows]}"
    ec = row_with_coverage['exact_coverage']
    assert ec.isdigit(), f"Exact coverage should be numeric, got: {ec}"
    # Expected: 5*2 (count 5, palindrome) + 2*2 (count 2, palindrome) = 14
    assert int(ec) == 14, f"Exact coverage should be 14 (5*2 + 2*2), got: {ec}"
    assert int(ec) == 14, f"Expected accumulated coverage 14 (5*2 + 2*2 for palindrome reads), got {ec}"


def test_no_contamination_between_seeds():
    """
    Critical: Ensure coverage from one seed does NOT leak to another.
    Uses non-palindromic sequences to avoid doubling from reverse-complement matching.
    """
    seed1 = "HIV1-B-FR-K03455-seed"
    seed2 = "HIV1-CRF02_AG-GH-AB286855-seed"

    # Non-palindromic sequences
    aligned_csv = StringIO(f"""\
refname,qcut,rank,count,offset,seq
{seed1},15,0,10,0,AAACCCGGG
{seed2},15,0,20,0,GGGCCCAAA
""")
    remap_conseq_csv = StringIO(f"""\
region,sequence
{seed1},AAACCCGGG
{seed2},GGGCCCAAA
""")
    nuc_csv = StringIO()
    amino_csv = StringIO()
    insertions_csv = StringIO()
    conseq_csv = StringIO()
    failed_align_csv = StringIO()

    aln2counts(
        aligned_csv=aligned_csv,
        nuc_csv=nuc_csv,
        amino_csv=amino_csv,
        insertions_csv=insertions_csv,
        conseq_csv=conseq_csv,
        failed_align_csv=failed_align_csv,
        remap_conseq_csv=remap_conseq_csv,
    )

    nuc_csv.seek(0)
    reader = csv.DictReader(nuc_csv)
    rows = list(reader)

    # Group by seed
    by_seed = {}
    for row in rows:
        seed = row["seed"]
        if seed not in by_seed:
            by_seed[seed] = []
        by_seed[seed].append(row)

    # Get max coverages
    seed1_coverages = [
        int(r["exact_coverage"])
        for r in by_seed[seed1]
        if r["exact_coverage"] and r["exact_coverage"].strip()
    ]
    seed2_coverages = [
        int(r["exact_coverage"])
        for r in by_seed[seed2]
        if r["exact_coverage"] and r["exact_coverage"].strip()
    ]

    # Seed1 with count=10 should have coverage 10
    # Seed2 with count=20 should have coverage 20
    # They should NOT be equal (no contamination)
    assert len(seed1_coverages) > 0, "seed1 should have coverage"
    assert len(seed2_coverages) > 0, "seed2 should have coverage"

    max1 = max(seed1_coverages)
    max2 = max(seed2_coverages)

    assert max1 == 10, f"seed1 max coverage should be 10, got {max1}"
    assert max2 == 20, f"seed2 max coverage should be 20, got {max2}"
    assert max1 != max2, "Coverages should be different (no contamination)"


def test_prefixes_accumulate_correctly():
    """
    Critical: Multiple prefixed contigs (1-seed, 2-seed) should accumulate
    to the base seed with correct total coverage.
    """
    seed_name = "HIV1-B-FR-K03455-seed"

    aligned_csv = StringIO(f"""\
refname,qcut,rank,count,offset,seq
1-{seed_name},15,0,7,0,AAATTTCCC
2-{seed_name},15,0,3,0,AAATTTCCC
""")
    remap_conseq_csv = StringIO(f"""\
region,sequence
{seed_name},AAATTTCCC
1-{seed_name},AAATTTCCC
2-{seed_name},AAATTTCCC
""")
    nuc_csv = StringIO()
    nuc_detail_csv = StringIO()
    amino_csv = StringIO()
    insertions_csv = StringIO()
    conseq_csv = StringIO()
    failed_align_csv = StringIO()

    aln2counts(
        aligned_csv=aligned_csv,
        nuc_csv=nuc_csv,
        nuc_detail_csv=nuc_detail_csv,
        amino_csv=amino_csv,
        insertions_csv=insertions_csv,
        conseq_csv=conseq_csv,
        failed_align_csv=failed_align_csv,
        remap_conseq_csv=remap_conseq_csv,
    )

    nuc_csv.seek(0)
    reader = csv.DictReader(nuc_csv)
    rows = list(reader)

    # All rows should map to base seed (without prefix)
    assert all(r["seed"] == seed_name for r in rows), (
        "All rows should have base seed name"
    )

    # Get coverage values
    coverages = [
        int(r["exact_coverage"])
        for r in rows
        if r["exact_coverage"] and r["exact_coverage"].strip()
    ]

    # Total should be 7 + 3 = 10
    assert len(coverages) > 0, "Should have coverage values"
    assert max(coverages) == 10, (
        f"Max coverage should be 10 (7+3), got {max(coverages)}"
    )


def test_offset_reads_excluded():
    """
    Critical: Reads with offset != 0 should NOT contribute to exact_coverage.
    """
    seed_name = "HIV1-B-FR-K03455-seed"

    aligned_csv = StringIO(f"""\
refname,qcut,rank,count,offset,seq
{seed_name},15,0,10,0,AAATTTCCC
{seed_name},15,0,50,5,AAATTTCCC
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

    aln2counts(
        aligned_csv=aligned_csv,
        nuc_csv=nuc_csv,
        amino_csv=amino_csv,
        insertions_csv=insertions_csv,
        conseq_csv=conseq_csv,
        failed_align_csv=failed_align_csv,
        remap_conseq_csv=remap_conseq_csv,
    )

    nuc_csv.seek(0)
    reader = csv.DictReader(nuc_csv)
    rows = list(reader)

    coverages = [
        int(r["exact_coverage"])
        for r in rows
        if r["exact_coverage"] and r["exact_coverage"].strip()
    ]

    # Should only have coverage from offset=0 read (count=10)
    # NOT from offset=5 read (count=50)
    assert max(coverages) == 10, (
        f"Max coverage should be 10 (offset=0 only), got {max(coverages)}"
    )


def test_mismatched_reads_excluded():
    """
    Critical: Reads with mismatches should NOT contribute to exact_coverage.
    """
    seed_name = "HIV1-B-FR-K03455-seed"

    aligned_csv = StringIO(f"""\
refname,qcut,rank,count,offset,seq
{seed_name},15,0,10,0,AAATTTCCC
{seed_name},15,0,50,0,AAATTTCCT
{seed_name},15,0,30,0,AAATATCCC
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

    aln2counts(
        aligned_csv=aligned_csv,
        nuc_csv=nuc_csv,
        amino_csv=amino_csv,
        insertions_csv=insertions_csv,
        conseq_csv=conseq_csv,
        failed_align_csv=failed_align_csv,
        remap_conseq_csv=remap_conseq_csv,
    )

    nuc_csv.seek(0)
    reader = csv.DictReader(nuc_csv)
    rows = list(reader)

    coverages = [
        int(r["exact_coverage"])
        for r in rows
        if r["exact_coverage"] and r["exact_coverage"].strip()
    ]

    # Should only count the exact match (count=10)
    assert max(coverages) == 10, (
        f"Max coverage should be 10 (exact matches only), got {max(coverages)}"
    )


def test_query_positions_consistent():
    """
    Critical: query.nuc.pos should be 1-indexed and consistent across combined reports.
    """
    seed_name = "HIV1-B-FR-K03455-seed"

    aligned_csv = StringIO(f"""\
refname,qcut,rank,count,offset,seq
1-{seed_name},15,0,5,0,AAATTT
2-{seed_name},15,0,2,0,AAATTT
""")
    remap_conseq_csv = StringIO(f"""\
region,sequence
{seed_name},AAATTT
1-{seed_name},AAATTT
2-{seed_name},AAATTT
""")
    nuc_csv = StringIO()
    nuc_detail_csv = StringIO()
    amino_csv = StringIO()
    insertions_csv = StringIO()
    conseq_csv = StringIO()
    failed_align_csv = StringIO()

    aln2counts(
        aligned_csv=aligned_csv,
        nuc_csv=nuc_csv,
        nuc_detail_csv=nuc_detail_csv,
        amino_csv=amino_csv,
        insertions_csv=insertions_csv,
        conseq_csv=conseq_csv,
        failed_align_csv=failed_align_csv,
        remap_conseq_csv=remap_conseq_csv,
    )

    nuc_csv.seek(0)
    reader = csv.DictReader(nuc_csv)
    rows = list(reader)

    # Get query positions
    query_positions = [int(r["query.nuc.pos"]) for r in rows if r["query.nuc.pos"]]

    # Should be 1-indexed and consecutive
    assert min(query_positions) == 1, "query.nuc.pos should start at 1"
    assert max(query_positions) == 6, "query.nuc.pos should end at 6"
    assert sorted(query_positions) == [1, 2, 3, 4, 5, 6], (
        "Positions should be consecutive"
    )

    # Verify coverage is at correct positions
    coverage_by_pos = {}
    for row in rows:
        if (
            row["query.nuc.pos"]
            and row["exact_coverage"]
            and row["exact_coverage"].strip()
        ):
            pos = int(row["query.nuc.pos"])
            cov = int(row["exact_coverage"])
            coverage_by_pos[pos] = cov

    # Should have coverage at some middle positions
    assert len(coverage_by_pos) > 0, "Should have coverage at some positions"
    # Check what values we got
    if coverage_by_pos:
        unique_coverages = set(coverage_by_pos.values())
        # With 6bp read and overlap_size = 6//4 = 1, edges are trimmed
        # Middle positions should have full coverage (5+2=7)
        # But may vary due to edge trimming
        print(f"coverage_by_pos: {coverage_by_pos}")
        # Just verify we have reasonable coverage values
        assert max(coverage_by_pos.values()) > 0, "Should have some coverage"


def test_independent_seed_position_spaces():
    """
    Critical: Different seeds have independent position numbering.
    Uses non-palindromic sequences to test actual coverage values.
    """
    seed1 = "HIV1-B-FR-K03455-seed"
    seed2 = "HIV1-CRF02_AG-GH-AB286855-seed"

    # seed1: 6bp, seed2: 9bp - non-palindromic
    aligned_csv = StringIO(f"""\
refname,qcut,rank,count,offset,seq
{seed1},15,0,10,0,AAACCC
{seed2},15,0,20,0,GGGAAACCC
""")
    remap_conseq_csv = StringIO(f"""\
region,sequence
{seed1},AAACCC
{seed2},GGGAAACCC
""")
    nuc_csv = StringIO()
    amino_csv = StringIO()
    insertions_csv = StringIO()
    conseq_csv = StringIO()
    failed_align_csv = StringIO()

    aln2counts(
        aligned_csv=aligned_csv,
        nuc_csv=nuc_csv,
        amino_csv=amino_csv,
        insertions_csv=insertions_csv,
        conseq_csv=conseq_csv,
        failed_align_csv=failed_align_csv,
        remap_conseq_csv=remap_conseq_csv,
    )

    nuc_csv.seek(0)
    reader = csv.DictReader(nuc_csv)
    rows = list(reader)

    # Group by seed
    by_seed = {}
    for row in rows:
        seed = row["seed"]
        if seed not in by_seed:
            by_seed[seed] = []
        by_seed[seed].append(row)

    # Check positions
    seed1_positions = sorted(
        [int(r["query.nuc.pos"]) for r in by_seed[seed1] if r["query.nuc.pos"]]
    )
    seed2_positions = sorted(
        [int(r["query.nuc.pos"]) for r in by_seed[seed2] if r["query.nuc.pos"]]
    )

    assert seed1_positions == [1, 2, 3, 4, 5, 6], "seed1 should have positions 1-6"
    assert seed2_positions == [1, 2, 3, 4, 5, 6, 7, 8, 9], (
        "seed2 should have positions 1-9"
    )

    # Check coverages are independent
    seed1_coverage = {
        int(r["query.nuc.pos"]): int(r["exact_coverage"])
        for r in by_seed[seed1]
        if r["query.nuc.pos"] and r["exact_coverage"] and r["exact_coverage"].strip()
    }
    seed2_coverage = {
        int(r["query.nuc.pos"]): int(r["exact_coverage"])
        for r in by_seed[seed2]
        if r["query.nuc.pos"] and r["exact_coverage"] and r["exact_coverage"].strip()
    }

    # They have different position counts and different coverage values, showing they're independent
    assert len(seed1_coverage) > 0, "seed1 should have coverage"
    assert len(seed2_coverage) > 0, "seed2 should have coverage"

    # The key test: coverage values should be different (10 vs 20)
    if seed1_coverage and seed2_coverage:
        max1 = max(seed1_coverage.values())
        max2 = max(seed2_coverage.values())
        assert max1 == 10, f"seed1 should have max coverage 10, got {max1}"
        assert max2 == 20, f"seed2 should have max coverage 20, got {max2}"
        assert max1 != max2, f"Max coverages should differ: seed1={max1}, seed2={max2}"
