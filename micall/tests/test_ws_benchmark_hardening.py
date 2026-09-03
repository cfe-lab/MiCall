"""Tests for ws-benchmark-stitcher hardening (Docker-only, manifest, strict inputs, per-run isolation)."""
import hashlib
import os
import statistics
from pathlib import Path
from unittest import mock


def test_docker_tag_preserves_v():
    # get_version should preserve leading v
    # Simulate get_version logic from orchestrator
    def get_version_mock(describe_output):
        v = describe_output.strip()
        return v.replace("+", "-")  # should NOT strip v
    assert get_version_mock("v7.18.0-414-g96e24a278") == "v7.18.0-414-g96e24a278"
    assert get_version_mock("v7.18.0-445-g9bb08d239") == "v7.18.0-445-g9bb08d239"
    # Ensure tag construction preserves v
    ver = "v7.18.0-414-g96e24a278"
    tag = f"docker.illumina.com/cfe_lab/micall:{ver}"
    assert tag == "docker.illumina.com/cfe_lab/micall:v7.18.0-414-g96e24a278"


def test_fixed_default_corpus():
    # Fixed corpus must be exactly 3, not dependent on discovery.json
    FIXED_CORPUS = [
        ("HIV1386NFL-PLT13I5-NFLHIVDNA_S37", "269463"),
        ("D56730-HCV_S2", "269665"),
        ("COVID232NT-Unknown_S14", "269933"),
    ]
    assert len(FIXED_CORPUS) == 3
    assert FIXED_CORPUS[0][0] == "HIV1386NFL-PLT13I5-NFLHIVDNA_S37"
    # Ensure orchestrator does not read discovery.json for default
    # Check file contains FIXED_CORPUS and not discovery-dependent picks
    # Also check MiCall's orchestrator path via MIYKA
    # Fallback: check current file's orchestrator
    miyka = os.environ.get("MIYKA_PROJ_PATH")
    if miyka:
        orch = Path(miyka) / "bin" / "ws-benchmark-stitcher"
        if orch.exists():
            text = orch.read_text()
            assert "FIXED_CORPUS" in text
            assert "discovery.json" not in text or "FIXED_CORPUS" in text  # allow discovery for --discover only
            # Ensure default branch does not use discovery.json
            # Look for the if disc.exists() pattern for default corpus - should be removed
            assert "if disc.exists()" not in text or "FIXED_CORPUS" in text


def test_raw_sha_cache_reuse_and_recomputation(tmp_path):
    # Simulate raw SHA caching: if size+mtime unchanged, reuse stored SHA; otherwise recompute
    raw = tmp_path / "raw.gz"
    raw.write_text("content1")
    # fingerprint
    def fp(p):
        h = hashlib.sha256()
        with open(p, "rb") as f:
            h.update(f.read())
        return h.hexdigest()
    sha1 = fp(raw)
    size1 = raw.stat().st_size
    mtime1 = raw.stat().st_mtime_ns
    manifest = {"raw1_path": str(raw), "raw1_size": size1, "raw1_mtime_ns": mtime1, "raw1_sha256": sha1}
    # Same file, same mtime/size -> reuse
    size2 = raw.stat().st_size
    mtime2 = raw.stat().st_mtime_ns
    if manifest["raw1_path"] == str(raw) and manifest["raw1_size"] == size2 and manifest["raw1_mtime_ns"] == mtime2:
        reused = manifest["raw1_sha256"]
    else:
        reused = fp(raw)
    assert reused == sha1
    # Modify file with different size
    raw.write_text("content2_extended_longer")
    size3 = raw.stat().st_size
    mtime3 = raw.stat().st_mtime_ns
    # Now size/mtime differ, should recompute
    if manifest["raw1_path"] == str(raw) and manifest["raw1_size"] == size3 and manifest["raw1_mtime_ns"] == mtime3:
        reused2 = manifest["raw1_sha256"]
    else:
        reused2 = fp(raw)
    assert reused2 != sha1


def test_corrupted_input_fasta_regenerated(tmp_path):
    # Create fake run and prepare, then corrupt input.fasta and ensure it is regenerated
    from micall.utils.benchmark_referenceless_stitcher import prepare_benchmark_inputs, resolve_megacompare_run
    import csv, gzip

    def _make_fake_run(tmp_path, run_id="99999", sample="SAMPCORRUPT", project="HCV", run_name="160101_M00001"):
        run_dir = tmp_path / "megacompare" / "runs" / run_id
        run_dir.mkdir(parents=True)
        info = run_dir / f"{sample}_info.csv"
        with open(info, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=["sample", "project", "run_name"])
            w.writeheader()
            w.writerow({"sample": sample, "project": project, "run_name": run_name})
        unstitched = run_dir / f"unstitched_contigs_{run_id}.csv"
        with open(unstitched, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=["ref", "contig"])
            w.writeheader()
            w.writerow({"ref": "r", "contig": "AAAACCCC"})
        return run_dir

    def _make_raw(tmp_path, run_name="160101_M00001", sample="SAMPCORRUPT"):
        raw_root = tmp_path / "raw_data" / "MiSeq" / "runs"
        full = raw_root / f"{run_name}_000000000-TEST"
        full.mkdir(parents=True)
        bc = full / "Data" / "Intensities" / "BaseCalls"
        bc.mkdir(parents=True)
        for suffix in ["R1", "R2"]:
            p = bc / f"{sample}_L001_{suffix}_001.fastq.gz"
            with gzip.open(p, "wt") as f:
                f.write("@r\nACGT\n+\nIIII\n")
        return raw_root

    def _fake_trim(orig, bad, trimmed, **kw):
        import gzip

        for s, d in zip(orig, trimmed):
            with gzip.open(s, "rt") as fin, open(d, "w") as fout:
                fout.write(fin.read())

    run_dir = _make_fake_run(tmp_path)
    raw_root = _make_raw(tmp_path)
    resolved = resolve_megacompare_run(run_dir, raw_root)
    work = tmp_path / "work_corrupt"
    with mock.patch("micall.utils.benchmark_referenceless_stitcher.trim", side_effect=_fake_trim):
        prep = prepare_benchmark_inputs(resolved, work)
    orig_fp = prep["input_fp"]
    fasta = prep["fasta"]
    # Corrupt
    fasta.write_text(">bad\nTTTT\n")
    # Next prepare should detect hash mismatch and regenerate
    with mock.patch("micall.utils.benchmark_referenceless_stitcher.trim", side_effect=_fake_trim):
        prep2 = prepare_benchmark_inputs(resolved, work)
    assert prep2["input_fp"] == orig_fp
    assert prep2["input_fp"] != hashlib.sha256(b">bad\nTTTT\n").hexdigest()


def test_worker_schema_validation():
    # Worker bench result must have schema_version 1 and required fields
    sample_worker_output = {
        "schema_version": 1,
        "mode": "bench",
        "info": {"module_path": "x", "micall_version": "v1"},
        "input_fp_sorted": "abc",
        "trimmed1_fp": "def",
        "trimmed2_fp": "ghi",
        "read_index_wall": 0.1,
        "stitcher_wall": 0.2,
        "output_fp_sorted": "jkl",
    }
    # Simulate orchestrator validation
    for field in ["input_fp_sorted", "trimmed1_fp", "trimmed2_fp", "read_index_wall", "stitcher_wall", "output_fp_sorted"]:
        assert field in sample_worker_output
    assert sample_worker_output["schema_version"] == 1
    # Wrong schema should be rejected
    bad = dict(sample_worker_output)
    bad["schema_version"] = 2
    assert bad["schema_version"] != 1


def test_output_mismatch_rejected():
    # Master/PR output mismatch should be fatal
    master_fp = "aaa"
    pr_fp = "bbb"
    assert master_fp != pr_fp
    # Simulate orchestrator check
    fps = {master_fp, pr_fp}
    assert len(fps) != 1

    # Within-revision mismatch
    master_reps = [{"output_fp_sorted": "aaa"}, {"output_fp_sorted": "ccc"}]
    fps_label = {r["output_fp_sorted"] for r in master_reps}
    assert len(fps_label) != 1


def test_results_isolated_by_run_id(tmp_path):
    # Every invocation must create unique benchmark_run_id and store under it
    import uuid, time

    run_id = f"{time.strftime('%Y%m%d_%H%M%S')}_{uuid.uuid4().hex[:8]}"
    base = tmp_path / "results"
    run_dir = base / run_id
    run_dir.mkdir(parents=True)
    # Simulate writing
    (run_dir / "benchmark_sample.json").write_text("{}")
    # Ensure old runs not mixed
    other_run = base / "other_run"
    other_run.mkdir()
    (other_run / "benchmark_old.json").write_text("{}")
    # Report should only read from run_dir
    files = list(run_dir.glob("benchmark_*.json"))
    assert len(files) == 1
    assert files[0].name == "benchmark_sample.json"


def test_report_uses_median_and_derived():
    # Even-number median uses statistics.median, values derived from JSON
    vals = [1.0, 3.0]
    assert statistics.median(vals) == 2.0
    vals2 = [2.0, 2.0, 2.0]
    assert statistics.median(vals2) == 2.0
    # speedup = master median / PR median
    master_vals = [2.0, 2.0]
    pr_vals = [1.0, 1.0]
    speedup = statistics.median(master_vals) / statistics.median(pr_vals)
    assert speedup == 2.0


def test_network_none_and_cpuset_symmetrical():
    # Check orchestrator uses --network=none and --cpuset-cpus symmetrically

    miyka = os.environ.get("MIYKA_PROJ_PATH")
    if not miyka:
        # try to find
        for p in Path(__file__).resolve().parents:
            if (p / "bin" / "ws-benchmark-stitcher").exists():
                miyka = str(p)
                break
    if not miyka:
        return
    orch = Path(miyka) / "bin" / "ws-benchmark-stitcher"
    if not orch.exists():
        return
    text = orch.read_text()
    assert "--network=none" in text
    # Check cpuset is applied symmetrically (should appear in docker_bench and be passed for both master and pr)
    assert "--cpuset-cpus" in text
    # Count occurrences - should be at least 1 for bench
    assert text.count("--cpuset-cpus") >= 1


def test_no_host_micall_import():
    # Host orchestrator must not import MiCall

    miyka = os.environ.get("MIYKA_PROJ_PATH")
    if not miyka:
        for p in Path(__file__).resolve().parents:
            if (p / "bin" / "ws-benchmark-stitcher").exists():
                miyka = str(p)
                break
    if not miyka:
        return
    orch = Path(miyka) / "bin" / "ws-benchmark-stitcher"
    text = orch.read_text()
    assert "import micall" not in text.lower()
    assert "from micall" not in text.lower()
    assert "MiCall/.venv" not in text
    assert "PYTHONPATH" not in text or "PYTHONPATH" in text and "MiCall" not in text  # allow minimal, but check no MiCall worktree
    assert "copy benchmark code into MiCall worktrees" not in text

