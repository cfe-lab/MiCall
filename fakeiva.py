#!/usr/bin/env python3

"""
fakeiva.py - pretend to be IVA

Expected invocation pattern (from your pipeline):
    IVA --fr <joined.fastq> -t <threads> [--contigs <seeds.fasta> --make_new_seeds] <iva_out_dir>

Behavior:
  * Computes MD5 of <joined.fastq> and <seeds.fasta>
  * Looks up: /data/assembly-ios/<md5(joined)>__<md5(seeds)>/contigs.fasta
  * Copies that file into <iva_out_dir>/contigs.fasta

Notes:
  * Only a minimal subset of IVA flags is parsed. Unknown flags are ignored.
  * Requires that the directory structure was created previously by your copy script.
"""

from __future__ import annotations
import argparse
import os
import shutil
import subprocess as sp
import sys
from pathlib import Path

ASSEMBLY_IOS_ROOT = Path("/data/assembly-ios")

def run_capture(cmd: list[str]) -> str:
    proc = sp.run(cmd, check=True, text=True, capture_output=True)
    return proc.stdout

def md5sum(path: Path | None) -> str:
    """
    Prefer the system md5sum for exact compatibility; fallback to hashlib if unavailable.
    """
    if path is None:
        path = Path("/dev/null")
    out = run_capture(["md5sum", str(path)])
    return out.strip().split()[0]

def parse_args(argv: list[str]) -> tuple[Path, Path | None, Path]:
    # Parse a minimal subset and tolerate unknown options.
    p = argparse.ArgumentParser(add_help=False)
    p.add_argument("--fr", dest="joined", required=True)
    p.add_argument("--contigs", dest="seeds", required=False)
    p.add_argument("-t", dest="threads", required=False)  # ignored
    p.add_argument("--make_new_seeds", action="store_true")  # ignored
    # Collect everything else, with the expectation that the LAST token is the output dir
    args, rest = p.parse_known_args(argv)

    if not rest:
        print("[iva-stub] missing output directory argument", file=sys.stderr)
        sys.exit(2)

    iva_out_dir = Path(rest[-1])

    joined_path = Path(args.joined)
    if not joined_path.is_file():
        print(f"[iva-stub] --fr file not found: {joined_path}", file=sys.stderr)
        sys.exit(2)

    seeds_path: Path | None = None
    if args.seeds is not None:
        seeds_path = Path(args.seeds)
        if not seeds_path.is_file():
            print(f"[iva-stub] --contigs file not found: {seeds_path}", file=sys.stderr)
            sys.exit(2)

    return joined_path, seeds_path, iva_out_dir

def main(argv: list[str]) -> None:
    joined_path, seeds_path, out_dir = parse_args(argv)

    # Hash inputs
    md5_joined = md5sum(joined_path)
    md5_seeds = md5sum(seeds_path)

    src_dir = ASSEMBLY_IOS_ROOT / f"{md5_joined}__{md5_seeds}"
    src_contigs = src_dir / "contigs.fasta"

    if not src_contigs.is_file():
        print(f"[iva-stub] expected contigs at: {src_contigs}", file=sys.stderr)
        sys.exit(1)

    # Ensure output directory exists
    try:
        out_dir.mkdir(parents=True, exist_ok=True)
    except Exception as e:
        print(f"[iva-stub] cannot create output dir {out_dir}: {e}", file=sys.stderr)
        sys.exit(1)

    dst_contigs = out_dir / "contigs.fasta"

    try:
        shutil.copy2(src_contigs, dst_contigs)
    except Exception as e:
        print(f"[iva-stub] copy failed: {src_contigs} -> {dst_contigs}: {e}", file=sys.stderr)
        sys.exit(1)

    # Mimic a chatty tool on stdout (your pipeline captures it anyway)
    print(f"[iva-stub] copied {src_contigs} -> {dst_contigs}")

if __name__ == "__main__":
    # argv[0] is the program name; IVA gets called like:
    #   iva_stub.py --fr <joined> -t 2 --contigs <seeds> --make_new_seeds <iva_out_dir>
    main(sys.argv[1:])
