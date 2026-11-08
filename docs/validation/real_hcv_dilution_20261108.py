"""Generate private, lower-fraction real-read source-origin HCV controls.

Uses only original, already sequenced paired-end reads. It changes neither
nucleotides nor qualities, and shuffles paired records in lockstep.
The read-origin label is known; biological genotype is not independently
confirmed. Outputs go to a protected, Git-ignored local directory.
"""
import argparse
import gzip
import json
from pathlib import Path
import random


def pairs(folder):
    first, = folder.glob('*_L001_R1_001.fastq.gz')
    second, = folder.glob('*_L001_R2_001.fastq.gz')
    with gzip.open(first, 'rt') as a, gzip.open(second, 'rt') as b:
        while True:
            left = tuple(a.readline() for _ in range(4))
            right = tuple(b.readline() for _ in range(4))
            if not any(left) and not any(right):
                return
            if not all(left) or not all(right):
                raise ValueError('Incomplete FASTQ record')
            if left[0].split()[0] != right[0].split()[0]:
                raise ValueError('Mates not synchronized')
            yield left, right


def create(source, target, source_minor_count, seed,
           expected_major=50000, expected_minor=500):
    parent = Path(source)
    selected_major = []
    minor = []
    for left, right in pairs(parent):
        if left[0].startswith('@HCV1a|'):
            selected_major.append((left, right))
        elif left[0].startswith('@HCV2b|'):
            minor.append((left, right))
        else:
            raise ValueError('Unexpected source-origin label')
    if len(selected_major) != expected_major or len(minor) != expected_minor:
        raise ValueError('Source pair counts differ from the expected archive')
    if source_minor_count > len(minor) or source_minor_count < 0:
        raise ValueError('Invalid minor count')
    rng = random.Random(seed)
    reads = selected_major + rng.sample(minor, source_minor_count)
    rng.shuffle(reads)
    out = Path(target)
    out.mkdir(parents=True, exist_ok=False)
    names = (
        out / f'{out.name}_L001_R1_001.fastq.gz',
        out / f'{out.name}_L001_R2_001.fastq.gz',
    )
    try:
        with gzip.open(names[0], 'wt', compresslevel=3) as first, \
             gzip.open(names[1], 'wt', compresslevel=3) as second:
            for left, right in reads:
                first.writelines(left)
                second.writelines(right)
        info = {
            'source_major_pairs': len(selected_major),
            'source_minor_pairs': source_minor_count,
            'total_pairs': len(reads),
            'actual_minor_source_fraction': source_minor_count / len(reads),
            'shuffle_seed': seed,
            'note': ('Source identity by original read provenance only. '
                     'Neither independent HCV genotype nor minority '
                     'biological molecule fractions are certified.'),
        }
        (out / 'source_manifest.json').write_text(
            json.dumps(info, indent=2) + '\n')
        return info
    except BaseException:
        for file in names:
            file.unlink(missing_ok=True)
        out.rmdir()
        raise


def main():
    cli = argparse.ArgumentParser()
    cli.add_argument('--source', type=Path, required=True)
    cli.add_argument('--output', type=Path, required=True)
    cli.add_argument('--minor', type=int, required=True)
    cli.add_argument('--seed', type=int, required=True)
    args = cli.parse_args()
    print(json.dumps(create(args.source, args.output,
                            args.minor, args.seed), indent=2))


if __name__ == '__main__':
    main()
