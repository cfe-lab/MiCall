"""Audit the nine real HIV clone reads against 17 real biological references.

Nine Sanger-verified HIV9 sequences, five independently sequenced HIV5
isolates and three archived real HCV contigs. No sequence is invented.
"""
import argparse
import csv
from itertools import groupby
import json
from pathlib import Path
import subprocess
import time

from micall.core.denovo_multimap import select_best_pairs

def build(study, panels, bowtie):
    fasta = study / '17-real-references.fasta'
    records = []
    for kind, path in panels:
        with path.open() as f:
            for i, row in enumerate(csv.DictReader(f), 1):
                records.append((f'{kind}-{i}', row['contig']))
    assert len(records) == 17
    with fasta.open('w') as dest:
        for name, sequence in records:
            dest.write(f'>{name}\n{sequence}\n')
    subprocess.run([str(bowtie / 'bowtie2-build-s'), '-q',
                    str(fasta), str(fasta)], check=True,
                   capture_output=True, timeout=90)
    return fasta


def compare(label, flags, study, bowtie, fasta, fastq1, fastq2):
    sam = study / f'{label}.sam'
    args = [str(bowtie / 'bowtie2-align-s'), '--wrapper', 'micall-0',
            '--quiet', '-x', str(fasta), '--rdg', '10,3',
            '--rfg', '10,3', '-1', str(fastq1),
            '-2', str(fastq2), '--no-hd',
            '--local', '-X', '1200', '-p', '4', '--reorder', *flags]
    start = time.perf_counter()
    with sam.open('w') as dest:
        subprocess.run(args, stdout=dest, check=True,
                       stderr=subprocess.PIPE, timeout=240)
    selected = study / f'{label}.best.sam'
    diagnostics = {}
    with sam.open() as source, selected.open('w') as dest:
        select_best_pairs(source, dest,
                          limit=None if label == 'all' else 8,
                          diagnostics=diagnostics)
    refs = {}
    with selected.open() as source:
        for q, lines in groupby((line.strip().split('\t') for line in source
                                 if not line.startswith('@')),
                                lambda row: row[0]):
            refs[q] = sorted({row[2] for row in lines
                              if row[2] != '*' and not int(row[1]) & 4})
    return {'seconds': round(time.perf_counter() - start, 3),
            'diagnostics': diagnostics, 'refs': refs}


def count_fastq_pairs(path):
    with path.open() as f:
        lines = sum(1 for _ in f)
    if lines % 4:
        raise ValueError('Input FASTQ length is not a multiple of four')
    return lines // 4


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--study', type=Path, required=True)
    parser.add_argument('--hiv9-contigs', type=Path, required=True)
    parser.add_argument('--hiv5-contigs', type=Path, required=True)
    parser.add_argument('--hcv-contigs', type=Path, required=True)
    parser.add_argument('--read1', type=Path, required=True)
    parser.add_argument('--read2', type=Path, required=True)
    parser.add_argument('--bowtie-bin', type=Path, required=True)
    args = parser.parse_args()
    args.study.mkdir(parents=True, exist_ok=True)
    panels = (('HIV9', args.hiv9_contigs),
              ('HIV5', args.hiv5_contigs),
              ('HCV', args.hcv_contigs))
    fasta = build(args.study, panels, args.bowtie_bin)
    bounded = compare('k8', ['-k', '8'], args.study, args.bowtie_bin,
                      fasta, args.read1, args.read2)
    deeper = compare('all', ['-a'], args.study, args.bowtie_bin,
                     fasta, args.read1, args.read2)
    if (bounded['diagnostics'].get('oversized_pair_groups', 0)
            or deeper['diagnostics'].get('oversized_pair_groups', 0)):
        raise ValueError('Incomplete large-reference audit: oversized SAM group')
    failures = [
        len(set(refs) - set(bounded['refs'].get(q, [])))
        for q, refs in deeper['refs'].items()
    ]
    report = {
        'references': 17,
        'original_real_hiv9_pairs': count_fastq_pairs(args.read1),
        'bounded_seconds': bounded['seconds'],
        'all_seconds': deeper['seconds'],
        'bounded_diagnostics': bounded['diagnostics'],
        'all_diagnostics': deeper['diagnostics'],
        'pairs_losing_references_at_k8': sum(n > 0 for n in failures),
        'missing_reference_alternatives': sum(failures),
        'warning': 'All-hit Bowtie2 is still heuristic and not a truth oracle.',
    }
    (args.study / 'aggregate.json').write_text(
        json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
