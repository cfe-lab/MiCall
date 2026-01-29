import random
import pytest
from pathlib import Path
from io import StringIO

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from micall.utils.randomize_fastq import (
    introduce_errors,
    process_records,
    process_fastq,
    main,
    entry,
    NUCLEOTIDES,
)


# --- Tests for introduce_errors ---
def test_introduce_errors_no_error():
    seq = "ACGTACGT"
    qualities = [30] * len(seq)
    # no errors: all error rates set to 0.
    rng = random.Random(42)
    new_seq, new_quals = introduce_errors(seq, qualities,
                                          subst_rate=0.0,
                                          ins_rate=0.0,
                                          del_rate=0.0,
                                          ins_quality=20,
                                          rng=rng)
    assert new_seq == seq
    assert new_quals == qualities


def test_introduce_errors_full_deletion():
    seq = "ACGTACGT"
    qualities = [30] * len(seq)
    rng = random.Random(42)
    new_seq, new_quals = introduce_errors(seq, qualities,
                                          subst_rate=0.0,
                                          ins_rate=0.0,
                                          del_rate=1.0,
                                          ins_quality=20,
                                          rng=rng)
    # With full deletion, expect empty sequence and empty quality list.
    assert new_seq == ""
    assert new_quals == []


def test_introduce_errors_full_substitution():
    seq = "ACGT"
    qualities = [30, 31, 32, 33]
    # Set subst_rate=1 so every base is substituted.
    rng = random.Random(42)
    new_seq, new_quals = introduce_errors(seq, qualities,
                                          subst_rate=1.0,
                                          ins_rate=0.0,
                                          del_rate=0.0,
                                          ins_quality=20,
                                          rng=rng)
    assert len(new_seq) == len(seq)
    # for each base, new base should be different than original.
    for original, new in zip(seq, new_seq):
        assert new.upper() in NUCLEOTIDES
        assert new.upper() != original.upper()
    assert new_quals == qualities


def test_introduce_errors_full_insertion():
    seq = "ACGT"
    qualities = [30] * len(seq)
    rng = random.Random(42)
    new_seq, new_quals = introduce_errors(seq, qualities,
                                          subst_rate=0.0,
                                          ins_rate=0.99,
                                          del_rate=0.0,
                                          ins_quality=15,
                                          rng=rng)
    # Every base should be followed by an inserted base: length 2*len(seq)
    assert len(new_seq) > 2 * len(seq)


# --- Test for process_records using an in-memory FASTQ string ---
def test_process_records_no_error():
    fastq_str = (
        "@read1\n"
        "ACGT\n"
        "+\n"
        "!!!!\n"
    )

    # in FASTQ, qualities are usually encoded by ascii, but BioPython converts
    # them to a list of phred quality scores.
    handle = StringIO(fastq_str)
    # All error rates set to 0, so no changes.
    records = list(process_records(handle,
                                   subst_rate=0.0,
                                   ins_rate=0.0,
                                   del_rate=0.0,
                                   ins_quality=20,
                                   rng=random.Random(42)))
    assert len(records) == 1
    rec = records[0]
    assert str(rec.seq) == "ACGT"

    # Check phred qualities:
    # '!' converts to 0 in Phred.
    assert rec.letter_annotations["phred_quality"] == [0]*4


@pytest.fixture
def tmp_fastq_file(tmp_path: Path) -> Path:
    fastq_path = tmp_path / "input.fastq"
    records = [
        SeqRecord(Seq("ACGTACGT"),
                  id="read1",
                  description="",
                  letter_annotations={"phred_quality": [30]*8}),
    ]
    with fastq_path.open("w") as fh:
        SeqIO.write(records, fh, "fastq")
    return fastq_path


def test_process_fastq_no_error(tmp_fastq_file: Path, tmp_path: Path):
    output_path = tmp_path / "output.fastq"
    # Run process_fastq with error rates set to 0.
    process_fastq(in_fastq=tmp_fastq_file,
                  out_fastq=output_path,
                  subst_rate=0.0,
                  ins_rate=0.0,
                  del_rate=0.0,
                  ins_quality=20,
                  rng=random.Random(42))
    # Read back the output file.
    recs = list(SeqIO.parse(str(output_path), "fastq"))
    assert len(recs) == 1
    # Without errors, the sequence should be identical.
    assert str(recs[0].seq) == "ACGTACGT"


def test_main_integration(tmp_fastq_file: Path, tmp_path: Path, monkeypatch):
    # Build a fake command line to run main.
    out_path = tmp_path / "output2.fastq"
    args = [
        str(tmp_fastq_file),
        str(out_path),
        "--subst_rate", "0.0",
        "--ins_rate", "0.0",
        "--del_rate", "0.0",
        "--ins_quality", "20",
        "--seed", "42"
    ]
    # Call main directly.
    exit_code = main(args)
    assert exit_code == 0
    recs = list(SeqIO.parse(str(out_path), "fastq"))
    assert len(recs) == 1
    assert str(recs[0].seq) == "ACGTACGT"


def test_entry_integration(tmp_fastq_file: Path, tmp_path: Path, monkeypatch):
    out_path = tmp_path / "output3.fastq"
    monkeypatch.setattr("sys.argv", [
        "error_introducer.py",
        str(tmp_fastq_file),
        str(out_path),
        "--subst_rate", "0.0",
        "--ins_rate", "0.0",
        "--del_rate", "0.0",
        "--ins_quality", "20",
        "--seed", "42"
    ])

    # entry() calls sys.exit; we can catch the SystemExit exception.
    with pytest.raises(SystemExit) as e:
        entry()

    assert e.value.code == 0

    recs = list(SeqIO.parse(str(out_path), "fastq"))
    assert len(recs) == 1
    assert str(recs[0].seq) == "ACGTACGT"
