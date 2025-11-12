
import csv
from io import StringIO
import pytest


# Import functions, classes, and types from your module.
# Adjust the module name as necessary.
from micall.utils.calculate_kmer_frequencies import (
    kmers,
    slide_kmer,
    process_contig,
    process_all_contigs,
    counters_to_rows,
    read_contigs,
    main_typed,
    parse_arguments,
    configure_logging,
    UserError,
    Contig,
    KMer,
    KMerWithCounter,
    logger,
)


def test_kmers_basic():
    # For sequence 'ACGT' with k=2 we expect three k-mers:
    #   i=0: kmer "AC", left="", right="GT"
    #   i=1: kmer "CG", left="A", right="T"
    #   i=2: kmer "GT", left="AC", right=""
    sequence = "ACGT"
    size = 2
    expected = [
        ("AC", "", "GT"),
        ("CG", "A", "T"),
        ("GT", "AC", ""),
    ]
    results = list(kmers(sequence, size))
    assert len(results) == len(expected)
    for result, (exp_seq, exp_left, exp_right) in zip(results, expected):
        assert result.sequence == exp_seq
        assert result.left == exp_left
        assert result.right == exp_right


def test_kmers_size_exceeds_sequence():
    # When the k-mer size is greater than the sequence length,
    # kmers() should yield nothing.
    sequence = "ACG"
    size = 5
    results = list(kmers(sequence, size))
    assert results == []


def test_slide_kmer_counts():
    # Create a contrived KMer and assign left and right contexts that will
    # produce predictable matching counts.
    # We want a KMer whose sequence is "AC".
    # We set left context to "ACA" and right context to "ACC".
    #
    # left kmers (size 2 from "ACA"):
    #   i=0: "AC" → compare with "AC" gives 2 matches.
    #   i=1: "CA" → compare "AC" vs "CA": 0 matches.
    #
    # right kmers (size 2 from "ACC"):
    #   i=0: "AC" → compare with "AC" gives 2 matches.
    #   i=1: "CC" → compare "AC" vs "CC": only second position matches (1 match).
    #
    # Expected counter: {2: 2, 1: 1}
    kmer_obj = KMer(sequence="AC", left="ACA", right="ACC")
    with_counter = KMerWithCounter(kmer=kmer_obj, counter=dict())
    # Replace the counter with a defaultdict(int) so slide_kmer can add counts.
    from collections import defaultdict
    with_counter = KMerWithCounter(kmer=kmer_obj, counter=defaultdict(int))
    slide_kmer(with_counter)
    counts = dict(with_counter.counter)
    assert counts.get(2) == 2, f"Expected two occurrences of 2 matches, got {counts.get(2)}"
    assert counts.get(1) == 1, f"Expected one occurrence of 1 match, got {counts.get(1)}"
    # Optionally, ensure that no unexpected keys are present.
    assert set(counts.keys()) == {1, 2}


def test_process_contig_deduplication():
    # Set up a simple contig.
    contig = Contig(name="test1", seq="ACGT", reads_count=None)
    pool = {}
    # Use max_kmer=2 so both 1-mers and 2-mers are processed.
    process_contig(pool, contig, max_kmer=2)
    # Check that pool contains the unique kmer sequences
    # For k=1: expected unique kmers: "A", "C", "G", "T"
    # For k=2: expected unique kmers: "AC", "CG", "GT"
    expected_keys = {"A", "C", "G", "T", "AC", "CG", "GT"}
    assert set(pool.keys()) == expected_keys
    # Check that counters are non-negative integers (they might be zero if no match was found)
    for entry in pool.values():
        for k, count in entry.counter.items():
            assert isinstance(k, int)
            assert isinstance(count, int)
            assert count >= 0


def test_process_all_contigs_multiple():
    contig1 = Contig(name="c1", seq="ACGT", reads_count=None)
    contig2 = Contig(name="c2", seq="CGTA", reads_count=None)
    contigs = [contig1, contig2]
    pool = {}
    process_all_contigs(pool, contigs, max_kmer=1)
    # For k=1, across both contigs, pool keys should be a subset of {A, C, G, T}
    for key in pool.keys():
        assert key in {"A", "C", "G", "T"}


def test_counters_to_rows():
    # Create a fake pool with one kmer and a counter.
    from collections import defaultdict
    km = KMer(sequence="AC", left="A", right="C")
    count_dict = defaultdict(int)
    count_dict[1] = 5
    count_dict[2] = 2
    pool = {"AC": KMerWithCounter(kmer=km, counter=count_dict)}
    rows = list(counters_to_rows(pool))
    # There should be two rows (one for each distinct match count)
    assert len(rows) == 2
    # Verify the structure of each row.
    for row in rows:
        assert "size" in row
        assert "matches" in row
        assert "occurrences" in row
        assert "kmer" in row
        # size should equal len("AC") i.e. 2.
        assert row["size"] == 2
        assert row["kmer"] == "AC"


def test_read_contigs():
    fasta_data = """>contig1
ACGTACGT
>contig2
TTGGCCA
"""
    fasta_io = StringIO(fasta_data)
    contigs = list(read_contigs(fasta_io))
    assert len(contigs) == 2
    names = {contig.name for contig in contigs}
    assert names == {"contig1", "contig2"}


def test_main_typed_creates_csv(tmp_path):
    # Create a temporary FASTA file with multiple contigs.
    fasta_data = """>c1
ACGTAC
>c2
GTACGT
"""
    input_file = tmp_path / "input.fasta"
    input_file.write_text(fasta_data)
    output_file = tmp_path / "output.csv"
    # Run main_typed with, say, max_kmer=2.
    main_typed(input_file, output_file, max_kmer=2)
    # Open and read the CSV and check header and that at least one row is present.
    with output_file.open() as fp:
        reader = csv.DictReader(fp)
        # FIELDNAMES should be: ["size", "matches", "occurrences", "kmer"]
        expected_fields = ("size", "matches", "occurrences", "kmer")
        assert tuple(reader.fieldnames) == expected_fields
        rows = list(reader)
        # We cannot predict exact row values, but we expect rows to be present.
        assert len(rows) > 0
        # Verify each row has the expected keys.
        for row in rows:
            for key in expected_fields:
                assert key in row


def test_parse_arguments_and_logging():
    test_args = ["input.fasta", "output.csv", "--max", "5", "--debug"]
    args = parse_arguments(test_args)
    # Ensure arguments are parsed correctly.
    assert args.input.name == "input.fasta"
    assert args.output.name == "output.csv"
    assert args.max == 5
    assert args.debug is True

    # Set logging level based on args.
    configure_logging(args)
    # When debug flag is set, logger level should be DEBUG.
    assert logger.getEffectiveLevel() == 10  # 10 corresponds to DEBUG


def test_main_typed_missing_file(tmp_path):
    bad_input = tmp_path / "nonexistent.fasta"
    output_file = tmp_path / "output.csv"
    with pytest.raises(UserError):
        main_typed(bad_input, output_file, max_kmer=2)
