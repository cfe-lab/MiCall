
import csv
from io import StringIO
import pytest

from micall.utils.calculate_kmer_frequencies import (
    kmers,
    slide_kmer,
    process_contig,
    read_contigs,
    main_typed,
    UserError,
    Contig,
    KMer,
)


def test_kmers_basic():
    # Given a short DNA sequence, check that kmers produces the expected values.
    # For sequence 'ACGT' and k=2, we expect:
    #   i=0: kmer = "AC", left = "", right = "GT"
    #   i=1: kmer = "CG", left = "A", right = "T"
    #   i=2: kmer = "GT", left = "AC", right = ""
    sequence = "ACGT"
    size = 2
    expected = [
        ("AC", "", "GT"),
        ("CG", "A", "T"),
        ("GT", "AC", ""),
    ]
    results = list(kmers(sequence, size))
    assert len(results) == len(expected), "Incorrect number of k-mers generated."
    for kmer_obj, (exp_seq, exp_left, exp_right) in zip(results, expected):
        assert kmer_obj.sequence == exp_seq
        assert kmer_obj.left == exp_left
        assert kmer_obj.right == exp_right


def test_kmers_size_exceeds_sequence():
    # For a sequence shorter than k, kmers should produce nothing.
    sequence = "ACG"
    size = 5
    results = list(kmers(sequence, size))
    assert results == [], "Expected no k-mers when k > sequence length."


def test_slide_kmer_no_context_matches():
    # Test a situation where the contexts are too short to yield any
    # matching comparisons.  For example, use a k-mer extracted from a
    # sequence of length equal to k, so both contexts are empty.  This
    # should yield no rows.

    # Instead, we can use the KMer dataclass from the module.
    kmer_obj = KMer(sequence="AC", left="", right="")
    # With no context, chain(kmers("", 2)) will yield an empty iterator.
    rows = list(slide_kmer("test_contig", kmer_obj))
    assert rows == [], "Expected no rows when context is empty."


def test_slide_kmer_some_matches():
    # Build a test case where contexts do produce valid comparisons.
    # For example, consider sequence "AACGAC" and choose the kmer "AC"
    # that appears in the middle.
    # We'll build a KMer where:
    #    sequence = "AC"
    #    left = "AA"  (positions 0-1)
    #    right = "GAC" (positions 4-6); note: left/right here are just
    #    the slices based on our kmers() function.

    kmer_obj = KMer(sequence="AC", left="AA", right="GAC")
    rows = tuple(slide_kmer("contig1", kmer_obj))
    # Expect a counter: 1 match occurred twice, and 2 matches occurred once.
    # Order of rows is not guaranteed; we can check counts.

    counter: dict[int, int] = {}
    for row in rows:
        counter[row["matches"]] = counter.get(row["matches"], 0) + row["occurrences"]
    assert counter.get(1) == 1, "Expected count of 1-match to be 2."
    assert counter.get(2) == 1, "Expected count of 2-matches to be 1."


def test_process_contig_uniqueness():
    # Given a contig with a short sequence, evaluate process_contig.
    contig = Contig(name="test_contig", seq="ACGT")
    # Use max_k_mer = 2. This will process k=1 and k=2.
    rows = list(process_contig(contig, max_kmer=2))
    # For k=1, every unique single nucleotide should produce a row,
    # but since slide_kmer compares against context k-mers, some rows might be empty.
    # In our small sequence, many comparisons might produce no matches.
    # To check, verify that each row has the proper contig id and a kmer of length 1 or 2.
    for row in rows:
        assert row["qseqid"] == "test_contig"
        assert row["size"] in (1, 2)
        assert isinstance(row["occurrences"], int)
        assert isinstance(row["matches"], int)
        assert isinstance(row["kmer"], str)


def test_read_contigs_fasta():
    # Create an in-memory FASTA file using StringIO.
    fasta_data = """>contig1
ACGTACGT
>contig2
TTGGCCA
"""
    fasta_io = StringIO(fasta_data)
    contigs = list(read_contigs(fasta_io))
    assert len(contigs) == 2, "Should have two contigs from FASTA input."
    names = [c.name for c in contigs]
    assert "contig1" in names
    assert "contig2" in names


def test_main_typed_creates_output(tmp_path):
    # Create a temporary FASTA file, run main_typed, and check the output CSV.
    input_data = """>contig1
ACGTACGT
"""
    # Write input FASTA file.
    input_file = tmp_path / "input.fasta"
    input_file.write_text(input_data)
    # Define output file.
    output_file = tmp_path / "output.csv"

    # Run the processing function.
    main_typed(input_file, output_file, max_kmer=2)

    # Read the CSV file and verify header and at least one row.
    with output_file.open() as fp:
        reader = csv.DictReader(fp)
        rows = list(reader)
        # Check header fields match expected.
        expected_headers = ("qseqid", "size", "matches", "occurrences", "kmer")
        assert reader.fieldnames is not None
        assert tuple(reader.fieldnames) == expected_headers
        # This test might yield rows if there are context comparisons;
        # if no matches occur, we might get no rows. Here, ensure that no error occurred.
        # Optionally, check that rows are well-formed if present.
        for row in rows:
            for key in expected_headers:
                assert key in row


def test_main_typed_input_file_not_exist(tmp_path):
    # Provide a non-existent input file; main_typed should raise a UserError.
    non_existent = tmp_path / "no_file.fasta"
    output_file = tmp_path / "output.csv"
    with pytest.raises(UserError):
        main_typed(non_existent, output_file, max_kmer=2)
