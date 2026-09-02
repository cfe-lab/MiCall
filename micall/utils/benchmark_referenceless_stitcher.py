"""Benchmark tool for the referenceless contig stitcher with real megacompare inputs."""

import argparse
import csv
import hashlib
import json
import statistics
import subprocess
import sys
import time
from pathlib import Path
from typing import Dict, List, Tuple
import resource

from micall.core.trim_fastqs import trim
from micall.utils.referenceless_contig_stitcher import build_read_index, referenceless_contig_stitcher


def resolve_megacompare_run(input_path: Path, raw_data_root: Path = Path("/media/raw_data/MiSeq/runs")) -> Dict:
    """Resolve a path anywhere inside megacompare/runs/<run_id>/ to benchmark inputs.

    Walks upward until a directory whose parent is `runs` is found.
    Expects to find exactly one `<sample>_info.csv` and one `unstitched_contigs_<run_id>.csv`.
    Returns dict with run_id, sample, project, run_name, raw_fastq1/2, unstitched_csv.
    """
    p = Path(input_path).resolve()
    if p.is_file():
        p = p.parent
    # Walk upward
    run_dir = None
    for parent in [p] + list(p.parents):
        if parent.name == "runs":
            # Should not happen, we need the child
            continue
        if parent.parent.name == "runs":
            # parent is like runs/<run_id>
            # Check if it looks like a run directory (contains *_info.csv)
            infos = list(parent.glob("*_info.csv"))
            if infos:
                run_dir = parent
                break
            # Also check if the original input_path was inside runs/<run_id> but we are at the run_id itself
            if parent == p and parent.is_dir() and any(parent.glob("*_info.csv")):
                run_dir = parent
                break
        # Also handle case where p itself is the run dir
        if p.parent.name == "runs" and p.is_dir() and any(p.glob("*_info.csv")):
            run_dir = p
            break
    if run_dir is None:
        # Fallback: walk all parents and look for runs/<run_id> pattern
        for parent in [Path(input_path).resolve()] + list(Path(input_path).resolve().parents):
            if parent.is_dir() and parent.parent.name == "runs" and list(parent.glob("unstitched_contigs_*.csv")):
                run_dir = parent
                break
        if run_dir is None:
            raise FileNotFoundError(f"Could not find megacompare run directory for {input_path} (expected runs/<run_id> with *_info.csv)")
    run_id = run_dir.name
    # Find unstitched contigs
    unstitched_candidates = list(run_dir.glob(f"unstitched_contigs_{run_id}.csv"))
    if not unstitched_candidates:
        # Fallback: any unstitched_contigs_*.csv
        unstitched_candidates = list(run_dir.glob("unstitched_contigs_*.csv"))
    if len(unstitched_candidates) != 1:
        raise FileNotFoundError(f"Expected exactly one unstitched_contigs csv in {run_dir}, found {unstitched_candidates}")
    unstitched_csv = unstitched_candidates[0]
    # Find sample info
    info_candidates = list(run_dir.glob("*_info.csv"))
    if len(info_candidates) != 1:
        # Filter out the case where there are multiple samples? Should be one per run, but check
        # If multiple, try to pick the one matching the run's sample? For now fail if ambiguous
        raise FileNotFoundError(f"Expected exactly one *_info.csv in {run_dir}, found {info_candidates}")
    info_csv = info_candidates[0]
    with open(info_csv) as f:
        reader = csv.DictReader(f)
        rows = list(reader)
        if len(rows) != 1:
            raise ValueError(f"Expected exactly one row in {info_csv}, found {len(rows)}")
        row = rows[0]
        sample = row["sample"]
        project = row.get("project", "")
        run_name = row["run_name"]
    if not sample or not run_name:
        raise ValueError(f"Missing sample/run_name in {info_csv}: {row}")
    # Find full MiSeq run directory
    # raw_data_root is /media/raw_data/MiSeq/runs by default, so parent is MiSeq/runs, need to handle correctly
    # Actually raw_data_root is the runs directory itself, so we need to handle both cases
    if raw_data_root.name == "runs":
        # raw_data_root is .../MiSeq/runs
        candidates = list(raw_data_root.glob(f"{run_name}*"))
    else:
        candidates = list(Path(raw_data_root).glob(f"{run_name}*"))
    if len(candidates) != 1:
        raise FileNotFoundError(f"Expected exactly one MiSeq run directory for {run_name} in {raw_data_root}, found {candidates}")
    full_run = candidates[0]
    # Raw FASTQs
    r1_candidates = list(full_run.glob(f"Data/Intensities/BaseCalls/{sample}_*_R1_*.fastq.gz"))
    r2_candidates = list(full_run.glob(f"Data/Intensities/BaseCalls/{sample}_*_R2_*.fastq.gz"))
    # Fallback: search recursively if not in BaseCalls
    if not r1_candidates:
        r1_candidates = list(full_run.rglob(f"{sample}_*R1*.fastq.gz"))
    if not r2_candidates:
        r2_candidates = list(full_run.rglob(f"{sample}_*R2*.fastq.gz"))
    if len(r1_candidates) != 1 or len(r2_candidates) != 1:
        raise FileNotFoundError(f"Expected exactly one R1/R2 for {sample} in {full_run}, found R1={r1_candidates} R2={r2_candidates}")
    raw1 = r1_candidates[0]
    raw2 = r2_candidates[0]
    if not raw1.exists() or raw1.stat().st_size == 0:
        raise FileNotFoundError(f"Raw R1 missing or empty: {raw1}")
    if not raw2.exists() or raw2.stat().st_size == 0:
        raise FileNotFoundError(f"Raw R2 missing or empty: {raw2}")
    if not unstitched_csv.exists() or unstitched_csv.stat().st_size == 0:
        # Allow empty? But we require at least header
        pass
    return {
        "run_dir": run_dir,
        "run_id": run_id,
        "sample": sample,
        "project": project,
        "run_name": run_name,
        "full_run": full_run,
        "raw1": raw1,
        "raw2": raw2,
        "unstitched_csv": unstitched_csv,
    }


def _fingerprint_contigs(contigs) -> str:
    h = hashlib.sha256()
    for seq in sorted(contigs):
        h.update(seq.encode())
        h.update(b"\0")
    return h.hexdigest()


def _fingerprint_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()


def _convert_unstitched_csv_to_fasta(csv_path: Path, fasta_path: Path) -> int:
    count = 0
    with open(csv_path) as csvfile, open(fasta_path, "w") as out:
        reader = csv.DictReader(csvfile)
        # Support both 'contig' and 'sequence' columns
        for i, row in enumerate(reader):
            seq = row.get("contig") or row.get("sequence") or ""
            seq = seq.strip()
            if not seq:
                continue
            out.write(f">contig_{i}\n{seq}\n")
            count += 1
    return count


def _get_input_contig_stats(fasta_path: Path) -> Tuple[int, int, int]:
    total = 0
    mx = 0
    n = 0
    with open(fasta_path) as f:
        seq = ""
        for line in f:
            if line.startswith(">"):
                if seq:
                    total += len(seq)
                    mx = max(mx, len(seq))
                    n += 1
                    seq = ""
            else:
                seq += line.strip()
        if seq:
            total += len(seq)
            mx = max(mx, len(seq))
            n += 1
    return n, total, mx


def _manifest_path(work_dir: Path) -> Path:
    return work_dir / "prepared.json"


def _raw_file_fingerprint(path: Path) -> str:
    # Use size + mtime for quick check, but also fingerprint for exactness if needed
    # For manifest, store both fingerprint and mtime/size for robustness
    try:
        stat = path.stat()
        return f"{stat.st_size}:{stat.st_mtime_ns}:{_fingerprint_file(path)}"
    except OSError:
        return ""


def prepare_benchmark_inputs(resolved: Dict, work_dir: Path, verbose: bool = False) -> Dict:
    """Prepare work dir outside timed section: FASTA + trimmed FASTQs.

    Returns dict with prepared paths and metadata, and records preparation wall time.
    Keeps trimmed FASTQs so repeated runs do not retrim, but only when a
    `prepared.json` manifest exists and matches the current sources/configuration.
    The manifest is written only after all preparation succeeds, so interrupted
    or failed trimming never leaves a valid cache.
    """
    start = time.perf_counter()
    work_dir = Path(work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = _manifest_path(work_dir)
    fasta_path = work_dir / "input_contigs.fasta"
    trimmed1 = work_dir / "trimmed_R1.fastq"
    trimmed2 = work_dir / "trimmed_R2.fastq"
    # Try to reuse via manifest
    if manifest_path.exists():
        try:
            manifest = json.loads(manifest_path.read_text())
            # Check all relevant source/configuration
            # Use raw file fingerprints (size+mtime+fingerprint) to detect changes
            # For unstitched CSV, also check its fingerprint
            expected_bad = ""
            for cand in [
                Path(resolved["full_run"]) / "Data" / "Intensities" / "BaseCalls" / "bad_cycles.csv",
                Path(resolved["full_run"]) / "bad_cycles.csv",
                work_dir / "bad_cycles.csv",
            ]:
                if cand.exists():
                    expected_bad = str(cand)
                    break
            if not expected_bad:
                expected_bad = str(work_dir / "bad_cycles_missing.csv")
            # Compare manifest
            if (
                manifest.get("raw1") == str(resolved["raw1"])
                and manifest.get("raw2") == str(resolved["raw2"])
                and manifest.get("raw1_fp") == _raw_file_fingerprint(Path(resolved["raw1"]))
                and manifest.get("raw2_fp") == _raw_file_fingerprint(Path(resolved["raw2"]))
                and manifest.get("unstitched_csv") == str(resolved["unstitched_csv"])
                and manifest.get("unstitched_fp") == _fingerprint_file(Path(resolved["unstitched_csv"]))
                and manifest.get("project") == resolved["project"]
                and manifest.get("bad_csv") == expected_bad
                and manifest.get("bad_fp") == (_fingerprint_file(Path(expected_bad)) if Path(expected_bad).exists() else "")
                and Path(manifest.get("trimmed1", "")).exists()
                and Path(manifest.get("trimmed2", "")).exists()
                and _fingerprint_file(Path(manifest["trimmed1"])) == manifest.get("trimmed1_fp")
                and _fingerprint_file(Path(manifest["trimmed2"])) == manifest.get("trimmed2_fp")
                and Path(manifest.get("fasta", "")).exists()
                and _fingerprint_file(Path(manifest["fasta"])) == manifest.get("fasta_fp")
            ):
                if verbose:
                    print(f"Reusing prepared inputs via {manifest_path}", file=sys.stderr)
                # Recompute stats for return (cheap, not trimmed)
                n_contigs, total_len, max_len = _get_input_contig_stats(fasta_path)
                # Also need trimmed_count and fingerprints (already in manifest, but recompute for return)
                trimmed_count = manifest.get("trimmed_count")
                if trimmed_count is None:
                    try:
                        with open(trimmed1) as tf:
                            trimmed_count = sum(1 for _ in tf) // 4
                    except OSError:
                        trimmed_count = None
                prep_wall = time.perf_counter() - start
                return {
                    "fasta": fasta_path,
                    "trimmed1": trimmed1,
                    "trimmed2": trimmed2,
                    "n_contigs": n_contigs,
                    "total_len": total_len,
                    "max_len": max_len,
                    "input_fp": manifest.get("input_fp"),
                    "prep_wall": prep_wall,
                    "trimmed_count": trimmed_count,
                    "trimmed1_fp": manifest.get("trimmed1_fp"),
                    "trimmed2_fp": manifest.get("trimmed2_fp"),
                }
        except Exception:  # noqa: BLE001, S110
            # Any manifest read/parse error -> treat as stale and re-prepare
            pass
    # Need to (re)prepare: convert and trim
    # 1. Convert unstitched CSV to FASTA (always, but will be overwritten)
    n = _convert_unstitched_csv_to_fasta(resolved["unstitched_csv"], fasta_path)
    if verbose:
        print(f"Converted {n} contigs to {fasta_path}", file=sys.stderr)
    n_contigs, total_len, max_len = _get_input_contig_stats(fasta_path)
    # 2. Trim raw FASTQs (nearest-recoverable: bad_cycles=[] if missing, skip=(), project_code)
    # Find bad_cycles file
    bad_candidates = [
        Path(resolved["full_run"]) / "Data" / "Intensities" / "BaseCalls" / "bad_cycles.csv",
        Path(resolved["full_run"]) / "bad_cycles.csv",
        work_dir / "bad_cycles.csv",
    ]
    bad_csv = ""
    for cand in bad_candidates:
        if cand.exists():
            bad_csv = str(cand)
            break
    if not bad_csv:
        bad_csv = str(work_dir / "bad_cycles_missing.csv")
    trim(
        (str(resolved["raw1"]), str(resolved["raw2"])),
        bad_csv,
        (str(trimmed1), str(trimmed2)),
        use_gzip=True,
        skip=(),
        project_code=resolved["project"],
    )
    if verbose:
        print(f"Trimmed {resolved['raw1']} -> {trimmed1}", file=sys.stderr)
    prep_wall = time.perf_counter() - start
    # Input fingerprints
    with open(fasta_path) as f:
        contigs = []
        seq = ""
        for line in f:
            if line.startswith(">"):
                if seq:
                    contigs.append(seq)
                    seq = ""
            else:
                seq += line.strip()
        if seq:
            contigs.append(seq)
    input_fp = _fingerprint_contigs(contigs)
    # Trimmed read count (cheap, but outside timed section)
    trimmed_count = None
    try:
        with open(trimmed1) as tf:
            trimmed_count = sum(1 for _ in tf) // 4
    except OSError:
        pass
    # Fingerprints for trimmed FASTQs (outside timed section, hash once)
    trimmed1_fp = _fingerprint_file(trimmed1) if trimmed1.exists() else None
    trimmed2_fp = _fingerprint_file(trimmed2) if trimmed2.exists() else None
    # Write manifest only after all preparation succeeds (atomic)
    # Recompute expected bad_csv for manifest
    expected_bad = ""
    for cand in [
        Path(resolved["full_run"]) / "Data" / "Intensities" / "BaseCalls" / "bad_cycles.csv",
        Path(resolved["full_run"]) / "bad_cycles.csv",
        work_dir / "bad_cycles.csv",
    ]:
        if cand.exists():
            expected_bad = str(cand)
            break
    if not expected_bad:
        expected_bad = str(work_dir / "bad_cycles_missing.csv")
    manifest = {
        "raw1": str(resolved["raw1"]),
        "raw2": str(resolved["raw2"]),
        "raw1_fp": _raw_file_fingerprint(Path(resolved["raw1"])),
        "raw2_fp": _raw_file_fingerprint(Path(resolved["raw2"])),
        "unstitched_csv": str(resolved["unstitched_csv"]),
        "unstitched_fp": _fingerprint_file(Path(resolved["unstitched_csv"])),
        "project": resolved["project"],
        "bad_csv": expected_bad,
        "bad_fp": _fingerprint_file(Path(expected_bad)) if Path(expected_bad).exists() else "",
        "fasta": str(fasta_path),
        "fasta_fp": _fingerprint_file(fasta_path),
        "trimmed1": str(trimmed1),
        "trimmed2": str(trimmed2),
        "trimmed1_fp": trimmed1_fp,
        "trimmed2_fp": trimmed2_fp,
        "input_fp": input_fp,
        "trimmed_count": trimmed_count,
        "n_contigs": n_contigs,
        "total_len": total_len,
        "max_len": max_len,
    }
    # Write atomically via temp file + rename
    tmp_manifest = manifest_path.with_suffix(".tmp")
    tmp_manifest.write_text(json.dumps(manifest, indent=2))
    tmp_manifest.rename(manifest_path)
    return {
        "fasta": fasta_path,
        "trimmed1": trimmed1,
        "trimmed2": trimmed2,
        "n_contigs": n_contigs,
        "total_len": total_len,
        "max_len": max_len,
        "input_fp": input_fp,
        "prep_wall": prep_wall,
        "trimmed_count": trimmed_count,
        "trimmed1_fp": trimmed1_fp,
        "trimmed2_fp": trimmed2_fp,
    }


def run_benchmark_worker(
    fasta_path: Path,
    trimmed1: Path,
    trimmed2: Path,
    min_depth: int,
    read_length: int,
    output_fasta: Path,
) -> Dict:
    """Timed worker: build_read_index + stitcher, returns measurements and output stats."""
    # Wall and CPU for read index
    t0_wall = time.perf_counter()
    t0_cpu = time.process_time()
    if min_depth == 0:
        read_index = None
        idx_wall = time.perf_counter() - t0_wall
        idx_cpu = time.process_time() - t0_cpu
        # Still need to count distinct reads? For metadata, read_index is None
        read_index_meta = {"buckets": 0, "min_L": None, "max_L": None, "distinct": 0, "total_mult": 0}
    else:
        read_index = build_read_index(trimmed1, trimmed2)
        idx_wall = time.perf_counter() - t0_wall
        idx_cpu = time.process_time() - t0_cpu
        # Metadata
        buckets = len(read_index)
        if buckets:
            min_L = min(read_index.keys())
            max_L = max(read_index.keys())
            distinct = sum(len(c) for c in read_index.values())
            total_mult = sum(sum(c.values()) for c in read_index.values())
        else:
            min_L = None
            max_L = None
            distinct = 0
            total_mult = 0
        read_index_meta = {"buckets": buckets, "min_L": min_L, "max_L": max_L, "distinct": distinct, "total_mult": total_mult}
    # Stitcher
    t1_wall = time.perf_counter()
    t1_cpu = time.process_time()
    with open(fasta_path) as fin, open(output_fasta, "w") as fout:
        referenceless_contig_stitcher(
            fin,
            fout,
            read_index=read_index if min_depth != 0 else None,
            minimum_read_depth=min_depth,
            read_length=read_length,
        )
    stitch_wall = time.perf_counter() - t1_wall
    stitch_cpu = time.process_time() - t1_cpu
    # Peak RSS (resident set size) for this process (includes both phases, but we report per repetition)
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss  # KB on Linux
    # Output stats
    with open(output_fasta) as f:
        out_contigs = []
        seq = ""
        for line in f:
            if line.startswith(">"):
                if seq:
                    out_contigs.append(seq)
                    seq = ""
            else:
                seq += line.strip()
        if seq:
            out_contigs.append(seq)
    out_fp = _fingerprint_contigs(out_contigs)
    out_n = len(out_contigs)
    out_total = sum(len(s) for s in out_contigs)
    out_max = max((len(s) for s in out_contigs), default=0)
    return {
        "read_index_wall": idx_wall if min_depth != 0 else 0.0,
        "read_index_cpu": idx_cpu if min_depth != 0 else 0.0,
        "stitcher_wall": stitch_wall,
        "stitcher_cpu": stitch_cpu,
        "combined_wall": (idx_wall if min_depth != 0 else 0.0) + stitch_wall,
        "peak_rss_kb": rss,
        "read_index_meta": read_index_meta,
        "output_n": out_n,
        "output_total": out_total,
        "output_max": out_max,
        "output_fp": out_fp,
    }


def _run_worker_subprocess(
    fasta: Path, trimmed1: Path, trimmed2: Path, min_depth: int, read_length: int, out_fasta: Path
) -> Dict:
    """Run one repetition in a fresh subprocess."""
    # Use python -c to avoid argparse path consumption issues with --worker
    code = (
        "from pathlib import Path; "
        "from micall.utils.benchmark_referenceless_stitcher import run_benchmark_worker; "
        "import json, sys; "
        f"res=run_benchmark_worker(Path({str(fasta)!r}), Path({str(trimmed1)!r}), Path({str(trimmed2)!r}), {min_depth}, {read_length}, Path({str(out_fasta)!r})); "
        "json.dump(res, sys.stdout)"
    )
    cmd = [sys.executable, "-c", code]
    proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if proc.returncode != 0:
        raise RuntimeError(f"Worker failed: {proc.stderr}\n{proc.stdout}")
    return json.loads(proc.stdout)


def _worker_main(args: List[str]):
    # args: fasta trimmed1 trimmed2 out_fasta min_depth read_length
    fasta = Path(args[0])
    t1 = Path(args[1])
    t2 = Path(args[2])
    out = Path(args[3])
    min_depth = int(args[4])
    read_length = int(args[5])
    res = run_benchmark_worker(fasta, t1, t2, min_depth, read_length, out)
    # Convert Path etc. not needed
    json.dump(res, sys.stdout)


def _get_git_commit() -> str:
    try:
        out = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=Path(__file__).parent.parent.parent, text=True)
        return out.strip()
    except Exception:  # noqa: BLE001
        return "unknown"


def main():
    parser = argparse.ArgumentParser(description="Benchmark referenceless stitcher with megacompare inputs")
    parser.add_argument("path", nargs="?", default=None, help="Path anywhere inside megacompare/runs/<run_id>/ (file or directory)")
    parser.add_argument("--repeats", type=int, default=1, help="Number of timed repetitions (each in fresh subprocess)")
    parser.add_argument("--work-dir", type=str, default=None, help="Work directory for prepared inputs (default: <run_dir>/benchmark_work)")
    parser.add_argument("--raw-data-root", type=str, default="/media/raw_data/MiSeq/runs", help="Raw data root for R1/R2 lookup")
    parser.add_argument("--minimum-read-depth", type=int, default=1, dest="min_depth", help="Minimum read depth (0 disables)")
    parser.add_argument("--read-length", type=int, default=150, help="Read length for window")
    parser.add_argument("--json", action="store_true", help="Print machine-readable JSON to stdout")
    parser.add_argument("--json-out", type=str, default=None, help="Write JSON to file")
    parser.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("worker_args", nargs=argparse.REMAINDER, help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.worker:
        _worker_main(args.worker_args)
        return
    if not args.path:
        parser.error("path is required unless --worker is used")
    if args.repeats < 1:
        parser.error("--repeats must be >=1")
    if args.min_depth < 0:
        parser.error("--minimum-read-depth must be >=0")
    if args.read_length < 1:
        parser.error("--read-length must be >=1")
    # Resolve
    try:
        resolved = resolve_megacompare_run(Path(args.path), Path(args.raw_data_root))
    except Exception as e:  # noqa: BLE001
        print(f"Error resolving {args.path}: {e}", file=sys.stderr)
        sys.exit(1)
    # Work dir
    work_dir = Path(args.work_dir) if args.work_dir else Path(resolved["run_dir"]) / "benchmark_work"
    # Preparation (outside timed)
    prep_start = time.perf_counter()
    try:
        prep = prepare_benchmark_inputs(resolved, work_dir)
    except Exception as e:  # noqa: BLE001
        print(f"Preparation failed: {e}", file=sys.stderr)
        sys.exit(1)
    prep_wall = time.perf_counter() - prep_start
    # Also include prep_wall from prepare_benchmark_inputs (which already measured conversion+trim)
    # Use the larger? We'll report both
    # Input stats
    n_contigs, total_len, max_len = prep["n_contigs"], prep["total_len"], prep["max_len"]
    input_fp = prep["input_fp"]
    # Historical contigs for comparison
    hist_contigs_path = resolved["run_dir"] / f"contigs_{resolved['run_id']}.csv"
    hist_fp = None
    hist_match = None
    if hist_contigs_path.exists():
        try:
            hist_seqs = []
            with open(hist_contigs_path) as cf:
                reader = csv.DictReader(cf)
                for row in reader:
                    s = row.get("contig") or row.get("sequence") or ""
                    s = s.strip()
                    if s:
                        hist_seqs.append(s)
            hist_fp = _fingerprint_contigs(hist_seqs)
        except Exception:  # noqa: BLE001, S110
            pass
    # Benchmark repetitions
    results: List[Dict] = []
    for i in range(args.repeats):
        out_fasta = work_dir / f"output_{i}.fasta"
        res = _run_worker_subprocess(prep["fasta"], prep["trimmed1"], prep["trimmed2"], args.min_depth, args.read_length, out_fasta)
        results.append(res)
    # Aggregate if N>1
    def _agg(key):
        vals = [r[key] for r in results]
        return {"min": min(vals), "median": statistics.median(vals), "mean": sum(vals)/len(vals), "max": max(vals)} if vals else None
    agg = {}
    if len(results) > 1:
        for k in ["read_index_wall", "read_index_cpu", "stitcher_wall", "stitcher_cpu", "combined_wall", "peak_rss_kb"]:
            agg[k] = _agg(k)
    # Historical match
    if hist_fp is not None and results:
        hist_match = (results[-1]["output_fp"] == hist_fp)
    # Git commit
    git_commit = _get_git_commit()
    # Build JSON
    json_out = {
        "micall_commit": git_commit,
        "run_id": resolved["run_id"],
        "sample": resolved["sample"],
        "project": resolved["project"],
        "run_name": resolved["run_name"],
        "resolved": {k: str(v) for k, v in resolved.items()},
        "prepared": {
            "work_dir": str(work_dir),
            "fasta": str(prep["fasta"]),
            "trimmed1": str(prep["trimmed1"]),
            "trimmed2": str(prep["trimmed2"]),
            "n_contigs": n_contigs,
            "total_len": total_len,
            "max_len": max_len,
            "input_fp": input_fp,
            "prep_wall": prep_wall,
            "prep_detail_wall": prep["prep_wall"],
            "trimmed_count": prep.get("trimmed_count"),
            "trimmed1_fp": prep.get("trimmed1_fp"),
            "trimmed2_fp": prep.get("trimmed2_fp"),
        },
        "params": {"repeats": args.repeats, "min_depth": args.min_depth, "read_length": args.read_length},
        "per_repetition": results,
        "aggregate": agg if agg else None,
        "output_fp": results[-1]["output_fp"] if results else None,
        "historical_fp": hist_fp,
        "historical_match": hist_match,
    }
    # Human-readable summary
    if not args.json or args.json_out is None:
        # Print human if not pure json stdout
        pass
    # Always print human to stderr if json to stdout, or to stdout if not json
    out_stream = sys.stdout if not args.json else sys.stderr
    print(f"Sample: {resolved['sample']}  Run: {resolved['run_id']}  Project: {resolved['project']}", file=out_stream)
    print(f"Inputs: {prep['n_contigs']} contigs, total {prep['total_len']} max {prep['max_len']}, fingerprint {input_fp[:12]}", file=out_stream)
    print(f"Prepared: {work_dir} (prep wall {prep_wall:.2f}s, detail {prep['prep_wall']:.2f}s)", file=out_stream)
    print(f"Params: min_depth={args.min_depth} read_length={args.read_length} repeats={args.repeats}", file=out_stream)
    if results and results[0].get("read_index_meta"):
        meta = results[0]["read_index_meta"]
        print(f"Read index: buckets {meta['buckets']} min {meta['min_L']} max {meta['max_L']} distinct {meta['distinct']} total {meta['total_mult']} trimmed_count {prep.get('trimmed_count')}", file=out_stream)
    for idx, r in enumerate(results):
        print(f"Rep {idx}: idx wall {r['read_index_wall']:.3f}s cpu {r['read_index_cpu']:.3f}s | stitcher wall {r['stitcher_wall']:.3f}s cpu {r['stitcher_cpu']:.3f}s | combined {r['combined_wall']:.3f}s RSS {r['peak_rss_kb']/1024:.1f} MB | out {r['output_n']} contigs total {r['output_total']} max {r['output_max']} fp {r['output_fp'][:12]}", file=out_stream)
    if agg:
        print(f"Aggregate combined wall: min {agg['combined_wall']['min']:.3f} median {agg['combined_wall']['median']:.3f} mean {agg['combined_wall']['mean']:.3f} max {agg['combined_wall']['max']:.3f}", file=out_stream)
    if hist_fp is not None:
        print(f"Historical {hist_fp[:12]} match: {hist_match} (informational, trimmed FASTQs not historically bit-identical)", file=out_stream)
    print(f"Git commit: {git_commit}", file=out_stream)
    # JSON output
    if args.json:
        json.dump(json_out, sys.stdout, indent=2)
    if args.json_out:
        with open(args.json_out, "w") as jf:
            json.dump(json_out, jf, indent=2)
        print(f"JSON written to {args.json_out}", file=out_stream)


if __name__ == "__main__":
    main()
