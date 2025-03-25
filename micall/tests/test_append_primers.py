
import pytest
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from micall.utils.append_primers import (
    append_primers_to_record,
    entry,
    parse_arguments,
    append_primers,
    UserError,
    DEFAULT_FORWARD_PRIMER,
    DEFAULT_REVERSE_PRIMER,
)


def test_append_primers_to_record():
    original_seq = "ATGCGT"
    record = SeqRecord(Seq(original_seq), id="test", description="dummy")
    fwd_primer = "AAA"
    rev_primer = "TTT"
    new_record = append_primers_to_record(record, fwd_primer, rev_primer)
    expected_seq = fwd_primer + original_seq + rev_primer
    assert str(new_record.seq) == expected_seq
    assert new_record.id == record.id
    assert new_record.description == record.description


def test_parse_arguments_defaults():
    # simulate command-line arguments,
    # note that the first two arguments are input and output files.
    test_args = [
        "input.fa",
        "output.fa",
    ]
    args = parse_arguments(test_args)
    # Check that default primers are set.
    assert args.input_fasta == Path("input.fa")
    assert args.output_fasta == Path("output.fa")
    assert args.forward_primer == DEFAULT_FORWARD_PRIMER
    assert args.reverse_primer == DEFAULT_REVERSE_PRIMER


def test_parse_arguments_custom_primers():
    test_args = [
        "in.fa",
        "out.fa",
        "--forward-primer", "CUSTOMFWD",
        "--reverse-primer", "CUSTOMREV"
    ]
    args = parse_arguments(test_args)
    assert args.forward_primer == "CUSTOMFWD"
    assert args.reverse_primer == "CUSTOMREV"


@pytest.fixture
def temp_fasta_file(tmp_path: Path) -> Path:
    # Create a temporary FASTA file with two records.
    fasta_path = tmp_path / "input.fa"
    records = [
        SeqRecord(Seq("ATGC"), id="rec1", description="first"),
        SeqRecord(Seq("GATTACA"), id="rec2", description="second"),
    ]
    with fasta_path.open("w") as f:
        SeqIO.write(records, f, "fasta")
    return fasta_path


def test_append_primers_default(temp_fasta_file: Path, tmp_path: Path):
    # Use default primers.
    output_path = tmp_path / "output.fa"
    append_primers(
        input=temp_fasta_file,
        output=output_path,
        forward_primer=DEFAULT_FORWARD_PRIMER,
        reverse_primer=DEFAULT_REVERSE_PRIMER,
    )
    # Now read back the output.
    records = list(SeqIO.parse(str(output_path), "fasta"))
    # Check that the primers were appended.
    for rec in records:
        # The record sequence should start with DEFAULT_FORWARD_PRIMER
        # and end with DEFAULT_REVERSE_PRIMER.
        seq_str = str(rec.seq)
        assert seq_str.startswith(DEFAULT_FORWARD_PRIMER)
        assert seq_str.endswith(DEFAULT_REVERSE_PRIMER)
        # The middle part should be the original sequence. For instance,
        # for rec1, the original sequence was "ATGC".
        # We can extract the original: sequence[len(forward): -len(reverse)]
        original_seq = seq_str[len(DEFAULT_FORWARD_PRIMER): -len(DEFAULT_REVERSE_PRIMER)]
        # Check that the original sequence exists in the input file.
        # Without building a dict from the input file, we can simply check that
        # the extracted sequence is one of the known ones.
        assert original_seq in {"ATGC", "GATTACA"}


def test_append_primers_default_main(temp_fasta_file: Path,
                                     tmp_path: Path,
                                     monkeypatch):
    # Use default primers.
    output_path = tmp_path / "output.fa"

    monkeypatch.setattr("sys.argv", ["append_primers.py",
                                     str(temp_fasta_file),
                                     str(output_path),
                                     ])

    with pytest.raises(SystemExit) as exinfo:
        entry()

    ex: SystemExit = exinfo.value
    assert ex.code == 0

    # Now read back the output.
    records = list(SeqIO.parse(str(output_path), "fasta"))
    # Check that the primers were appended.
    for rec in records:
        # The record sequence should start with DEFAULT_FORWARD_PRIMER
        # and end with DEFAULT_REVERSE_PRIMER.
        seq_str = str(rec.seq)
        assert seq_str.startswith(DEFAULT_FORWARD_PRIMER)
        assert seq_str.endswith(DEFAULT_REVERSE_PRIMER)
        # The middle part should be the original sequence. For instance,
        # for rec1, the original sequence was "ATGC".
        # We can extract the original: sequence[len(forward): -len(reverse)]
        original_seq = seq_str[len(DEFAULT_FORWARD_PRIMER): -len(DEFAULT_REVERSE_PRIMER)]
        # Check that the original sequence exists in the input file.
        # Without building a dict from the input file, we can simply check that
        # the extracted sequence is one of the known ones.
        assert original_seq in {"ATGC", "GATTACA"}


def test_append_primers_custom(temp_fasta_file: Path, tmp_path: Path):
    output_path = tmp_path / "output_custom.fa"
    custom_fwd = "CUSTOMFWD"
    custom_rev = "CUSTOMREV"
    append_primers(
        input=temp_fasta_file,
        output=output_path,
        forward_primer=custom_fwd,
        reverse_primer=custom_rev,
    )
    records = list(SeqIO.parse(str(output_path), "fasta"))
    for rec in records:
        seq_str = str(rec.seq)
        assert seq_str.startswith(custom_fwd)
        assert seq_str.endswith(custom_rev)


def test_append_primers_missing_input(tmp_path: Path):
    # Create a non-existent input file path.
    missing_file = tmp_path / "non_existent.fa"
    output_file = tmp_path / "output.fa"
    with pytest.raises(UserError):
        append_primers(
            input=missing_file,
            output=output_file,
            forward_primer=DEFAULT_FORWARD_PRIMER,
            reverse_primer=DEFAULT_REVERSE_PRIMER,
        )
