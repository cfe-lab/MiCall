#!/usr/bin/env python3
"""Verify the cached-IVA docker build (adopted from 15a2b99).

Needs docker and raw MiSeq data, e.g.::

    build/verify_fake_iva.py --image micall:tag \\
        --run /media/raw_data/MiSeq/runs/240216_M04401_0292_000000000-L5WBN \\
        --sample NEG-R763-V3-3-V3LOOP_S132 \\
        --entry d41d8cd98f00b204e9800998ecf8427e__d41d8cd98f00b204e9800998ecf8427e

Checks: /bin/iva is fakeiva.py, /data/assembly-ios holds 512 entries,
and denovo on the sample yields byte-identical precomputed contigs.
"""

from __future__ import annotations

import argparse
import hashlib
import subprocess
import sys
import tempfile
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent


def sh(image: str, script: str, mounts: list[str]) -> subprocess.CompletedProcess[str]:
    cmd = ["docker", "run", "--rm", "--entrypoint", "sh"]
    for mount in mounts:
        cmd += ["-v", mount]
    return subprocess.run(cmd + [image, "-c", script],
                          capture_output=True, text=True)


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--image", required=True)
    parser.add_argument("--run", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--entry", required=True,
                        help="<md5(joined)>__<md5(seeds)> entry name")
    return parser.parse_args(argv)


def main(argv: list[str]) -> int:
    args = parse_args(argv)
    expected_iva = hashlib.md5((REPO_ROOT / "fakeiva.py").read_bytes()).hexdigest()
    expected = (REPO_ROOT / "build" / "assembly-ios" / args.entry / "contigs.fasta").read_bytes()

    proc = sh(args.image, "md5sum /bin/iva && ls /data/assembly-ios | wc -l", [])
    assert proc.returncode == 0, proc.stderr
    tokens = proc.stdout.split()
    iva_md5, count = tokens[0], tokens[-1]
    assert iva_md5 == expected_iva, f"/bin/iva mismatch: {iva_md5}"
    assert count == "512", f"expected 512 entries, got {count}"
    print(f"PASS image: /bin/iva is fakeiva.py, {count} cached entries")

    basecalls = f"{args.run}/Data/Intensities/BaseCalls"
    with tempfile.TemporaryDirectory(prefix="verify-fake-iva-") as work:
        mounts = [f"{basecalls}:{basecalls}:ro", f"{work}:{work}"]
        cmd = ("/opt/venv/bin/python /opt/micall/micall/core/denovo.py "
               f"{basecalls}/{args.sample}_L001_R1_001.fastq.gz "
               f"{basecalls}/{args.sample}_L001_R2_001.fastq.gz {work}/out.fasta")
        proc = sh(args.image, cmd, mounts)
        assert proc.returncode == 0, proc.stderr[-2000:]
        out = Path(work, "out.fasta").read_bytes()

    assert out == expected, "denovo output differs from precomputed contigs"
    print(f"PASS denovo: {args.sample} -> {len(expected)} cached bytes, byte-identical")
    print("ALL VERIFICATION CHECKS PASSED")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
