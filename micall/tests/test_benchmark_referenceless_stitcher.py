import csv
import gzip
from pathlib import Path
from unittest import mock

import pytest

from micall.utils.benchmark_referenceless_stitcher import (
    resolve_megacompare_run,
    prepare_benchmark_inputs,
    run_benchmark_worker,
    _fingerprint_contigs,
)


def _make_fake_megacompare_run(tmp_path: Path, run_id: str, sample: str, project: str, run_name: str):
    run_dir = tmp_path / "megacompare" / "runs" / run_id
    run_dir.mkdir(parents=True)
    # info.csv
    info_path = run_dir / f"{sample}_info.csv"
    with open(info_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["sample", "project", "run_name"])
        writer.writeheader()
        writer.writerow({"sample": sample, "project": project, "run_name": run_name})
    # unstitched contigs
    unstitched = run_dir / f"unstitched_contigs_{run_id}.csv"
    with open(unstitched, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["ref", "contig"])
        writer.writeheader()
        writer.writerow({"ref": "r", "contig": "AAAACCCC"})
        writer.writerow({"ref": "r", "contig": "CCCCGGGG"})
    return run_dir, unstitched, info_path


def _make_fake_raw_data(tmp_path: Path, run_name: str, sample: str):
    # raw_data_root is .../MiSeq/runs
    raw_root = tmp_path / "raw_data" / "MiSeq" / "runs"
    full_run = raw_root / f"{run_name}_000000000-TEST"
    full_run.mkdir(parents=True)
    basecalls = full_run / "Data" / "Intensities" / "BaseCalls"
    basecalls.mkdir(parents=True)
    r1 = basecalls / f"{sample}_L001_R1_001.fastq.gz"
    r2 = basecalls / f"{sample}_L001_R2_001.fastq.gz"
    # Small synthetic FASTQs
    for p, seq in [(r1, "ACGTACGTACGT"), (r2, "TGCATGCATGCA")]:
        with gzip.open(p, "wt") as f:
            f.write(f"@r1\n{seq}\n+\n{'I'*len(seq)}\n")
    return raw_root, r1, r2


def test_resolve_run_directory(tmp_path):
    run_dir, _, _ = _make_fake_megacompare_run(tmp_path, "12345", "SAMPLE1", "HCV", "160101_M01234")
    raw_root, _, _ = _make_fake_raw_data(tmp_path, "160101_M01234", "SAMPLE1")
    resolved = resolve_megacompare_run(run_dir, raw_root)
    assert resolved["run_id"] == "12345"
    assert resolved["sample"] == "SAMPLE1"
    assert resolved["run_name"] == "160101_M01234"


def test_resolve_nested_file(tmp_path):
    run_dir, _, _ = _make_fake_megacompare_run(tmp_path, "12345", "SAMPLE1", "HCV", "160101_M01234")
    _make_fake_raw_data(tmp_path, "160101_M01234", "SAMPLE1")
    # Create a nested file inside run_dir
    nested = run_dir / "subdir" / "file.txt"
    nested.parent.mkdir(parents=True, exist_ok=True)
    nested.write_text("hi")
    raw_root = tmp_path / "raw_data" / "MiSeq" / "runs"
    resolved = resolve_megacompare_run(nested, raw_root)
    assert resolved["run_id"] == "12345"


def test_sample_info_mapping(tmp_path):
    run_dir, _, _ = _make_fake_megacompare_run(tmp_path, "12345", "MYSAMPLE", "MyProject", "170202_M99999")
    raw_root, _, _ = _make_fake_raw_data(tmp_path, "170202_M99999", "MYSAMPLE")
    resolved = resolve_megacompare_run(run_dir, raw_root)
    assert resolved["sample"] == "MYSAMPLE"
    assert resolved["project"] == "MyProject"
    assert resolved["run_name"] == "170202_M99999"


def test_raw_r1_r2_resolution(tmp_path):
    run_dir, _, _ = _make_fake_megacompare_run(tmp_path, "12345", "SAMPLE1", "HCV", "160101_M01234")
    raw_root, r1, r2 = _make_fake_raw_data(tmp_path, "160101_M01234", "SAMPLE1")
    resolved = resolve_megacompare_run(run_dir, raw_root)
    assert resolved["raw1"] == r1
    assert resolved["raw2"] == r2


def test_missing_run_dir(tmp_path):
    with pytest.raises(FileNotFoundError):
        resolve_megacompare_run(tmp_path / "nonexistent", tmp_path / "raw_data" / "MiSeq" / "runs")


def test_ambiguous_info(tmp_path):
    run_dir = tmp_path / "megacompare" / "runs" / "99999"
    run_dir.mkdir(parents=True)
    # Create two info files -> ambiguous
    for sample in ["SAMPLE1", "SAMPLE2"]:
        info = run_dir / f"{sample}_info.csv"
        with open(info, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=["sample", "project", "run_name"])
            writer.writeheader()
            writer.writerow({"sample": sample, "project": "P", "run_name": "RUN"})
    # Also create unstitched
    (run_dir / "unstitched_contigs_99999.csv").write_text("ref,contig\nr,AAA\n")
    raw_root = tmp_path / "raw_data" / "MiSeq" / "runs"
    raw_root.mkdir(parents=True)
    with pytest.raises(FileNotFoundError, match="exactly one.*_info"):
        resolve_megacompare_run(run_dir, raw_root)


def _fake_trim(orig_files, bad_csv, trimmed_files, **kw):
    import gzip

    for src, dst in zip(orig_files, trimmed_files):
        # raw is gzipped, trimmed is expected uncompressed
        with gzip.open(src, "rt") as fin, open(dst, "w") as fout:
            fout.write(fin.read())


def test_conversion_to_fasta(tmp_path):
    run_dir, _, _ = _make_fake_megacompare_run(tmp_path, "12345", "SAMPLE1", "HCV", "160101_M01234")
    raw_root, _, _ = _make_fake_raw_data(tmp_path, "160101_M01234", "SAMPLE1")
    resolved = resolve_megacompare_run(run_dir, raw_root)
    work_dir = tmp_path / "work"
    with mock.patch("micall.utils.benchmark_referenceless_stitcher.trim", side_effect=_fake_trim):
        prep = prepare_benchmark_inputs(resolved, work_dir)
    assert prep["fasta"].exists()
    assert prep["n_contigs"] == 2
    assert prep["input_fp"] == _fingerprint_contigs(["AAAACCCC", "CCCCGGGG"])


def test_preparation_outside_repeated_runs(tmp_path):
    run_dir, _, _ = _make_fake_megacompare_run(tmp_path, "12345", "SAMPLE1", "HCV", "160101_M01234")
    raw_root, _, _ = _make_fake_raw_data(tmp_path, "160101_M01234", "SAMPLE1")
    resolved = resolve_megacompare_run(run_dir, raw_root)
    work_dir = tmp_path / "work"
    with mock.patch("micall.utils.benchmark_referenceless_stitcher.trim", side_effect=_fake_trim):
        prep1 = prepare_benchmark_inputs(resolved, work_dir)
    mtime1 = prep1["trimmed1"].stat().st_mtime
    # Second preparation should reuse
    with mock.patch("micall.utils.benchmark_referenceless_stitcher.trim", side_effect=_fake_trim):
        prep2 = prepare_benchmark_inputs(resolved, work_dir)
    mtime2 = prep2["trimmed1"].stat().st_mtime
    assert mtime1 == mtime2  # not re-trimmed


def test_min_depth_zero_skips_index(tmp_path):
    run_dir, _, _ = _make_fake_megacompare_run(tmp_path, "12345", "SAMPLE1", "HCV", "160101_M01234")
    raw_root, _, _ = _make_fake_raw_data(tmp_path, "160101_M01234", "SAMPLE1")
    resolved = resolve_megacompare_run(run_dir, raw_root)
    work_dir = tmp_path / "work"
    with mock.patch("micall.utils.benchmark_referenceless_stitcher.trim", side_effect=_fake_trim):
        prep = prepare_benchmark_inputs(resolved, work_dir)
    out = work_dir / "out.fasta"
    res = run_benchmark_worker(prep["fasta"], prep["trimmed1"], prep["trimmed2"], 0, 150, out)
    assert res["read_index_wall"] == 0.0
    assert res["read_index_meta"]["buckets"] == 0
    assert out.exists()


def test_fingerprint_stability(tmp_path):
    assert _fingerprint_contigs(["AAA", "CCC"]) == _fingerprint_contigs(["CCC", "AAA"])
    assert _fingerprint_contigs(["AAA"]) != _fingerprint_contigs(["AAA", "CCC"])


def test_json_result_shape(tmp_path):
    run_dir, _, _ = _make_fake_megacompare_run(tmp_path, "12345", "SAMPLE1", "HCV", "160101_M01234")
    raw_root, _, _ = _make_fake_raw_data(tmp_path, "160101_M01234", "SAMPLE1")
    resolved = resolve_megacompare_run(run_dir, raw_root)
    work_dir = tmp_path / "work"
    with mock.patch("micall.utils.benchmark_referenceless_stitcher.trim", side_effect=_fake_trim):
        prep = prepare_benchmark_inputs(resolved, work_dir)
    out = work_dir / "out.fasta"
    res = run_benchmark_worker(prep["fasta"], prep["trimmed1"], prep["trimmed2"], 1, 150, out)
    # Simulate JSON structure as in main
    json_out = {
        "micall_commit": "abc",
        "run_id": resolved["run_id"],
        "sample": resolved["sample"],
        "resolved": {k: str(v) for k, v in resolved.items()},
        "prepared": {"work_dir": str(work_dir), "n_contigs": prep["n_contigs"]},
        "params": {"repeats": 1, "min_depth": 1, "read_length": 150},
        "per_repetition": [res],
        "output_fp": res["output_fp"],
    }
    assert "micall_commit" in json_out
    assert "run_id" in json_out
    assert "resolved" in json_out
    assert "per_repetition" in json_out
    assert "output_fp" in json_out
    assert isinstance(res["output_fp"], str) and len(res["output_fp"]) == 64


def test_prepared_cache_not_reused_after_failed_trim(tmp_path):
    run_dir, _, _ = _make_fake_megacompare_run(tmp_path, "12345", "SAMPLE1", "HCV", "160101_M01234")
    raw_root, _, _ = _make_fake_raw_data(tmp_path, "160101_M01234", "SAMPLE1")
    resolved = resolve_megacompare_run(run_dir, raw_root)
    work_dir = tmp_path / "work"

    def fake_trim_fail(orig_files, bad_csv, trimmed_files, **kw):
        for dst in trimmed_files:
            Path(dst).write_text("partial")
        raise RuntimeError("trim failed")

    with mock.patch("micall.utils.benchmark_referenceless_stitcher.trim", side_effect=fake_trim_fail):
        with pytest.raises(RuntimeError, match="trim failed"):
            prepare_benchmark_inputs(resolved, work_dir)
        # Manifest must not exist after failure
        assert not (work_dir / "prepared.json").exists()
        # Stale trimmed files exist but should be ignored
        assert (work_dir / "trimmed_R1.fastq").exists()

    # Second call must not reuse the failed outputs; it must call trim again
    call_count: list[int] = []

    def fake_trim_ok(orig_files, bad_csv, trimmed_files, **kw):
        call_count.append(1)
        for src, dst in zip(orig_files, trimmed_files):
            with gzip.open(src, "rt") as fin, open(dst, "w") as fout:
                fout.write(fin.read())

    with mock.patch("micall.utils.benchmark_referenceless_stitcher.trim", side_effect=fake_trim_ok):
        prep = prepare_benchmark_inputs(resolved, work_dir)
        assert len(call_count) == 1
        assert (work_dir / "prepared.json").exists()
        assert prep["trimmed1"].exists()
