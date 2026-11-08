"""Prepare reproducible public biological-truth challenges for MiCall.

Truth panels come from independently resequenced HIV strains or
Sanger-confirmed HIV plasmids. Original SRA reads and derived FASTQ
remain git-ignored; never produce patient-specific reports.
"""
import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path
import random

from Bio import SeqIO


BENCHMARKS = {
    'HIV5': {
        'accession': 'SRR961514',
        'read_count': 714994,
        'target_pairs': 60000,
        'truth_fasta': '5-virus-mix/data/REF.fasta',
        'raw_dir': 'hiv5-real/reads',
    },
    'HIV9': {
        'accession': 'SRR12042290',
        'read_count': 2783356,
        'target_pairs': 240000,
        'truth_fasta': 'CliqueSNV-validation/relevant_haplotypes/HIV9exp.fasta',
        'raw_dir': 'hiv9-sanger/reads',
    },
}


def variant_base(original):
    return {'A': 'C', 'C': 'G', 'G': 'T', 'T': 'A'}[original]


def build_panels(root, label):
    config = BENCHMARKS[label]
    truth_source = root / config['truth_fasta']
    refs = [(r.id, str(r.seq).upper()) for r in SeqIO.parse(truth_source, 'fasta')]
    assert len(refs) in (5, 9)
    assert all(set(seq) <= set('ACGT') for _, seq in refs)
    folder = root / ('study-' + label)
    folder.mkdir(parents=True, exist_ok=True)
    log = []
    for kind in ('truth', 'perturbed'):
        with (folder / (kind + '-contigs.csv')).open('w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=['ref', 'match', 'group_ref',
                                                   'contig'], lineterminator='\n')
            writer.writeheader()
            for i, (name, seq) in enumerate(refs, 1):
                bases = list(seq)
                if kind == 'perturbed':
                    # Hold back a sequence-verified substitution every few
                    # hundred bases, away from primer/amplicon termini.
                    n = 10 if label == 'HIV5' else 5
                    for j in range(n):
                        position = (len(bases) * (j + 1)) // (n + 1)
                        position += (i % 3) - 1
                        original = bases[position]
                        mutant = variant_base(original)
                        bases[position] = mutant
                        log.append({'refname': f"{i}-HIV1-B-FR-K03455-seed",
                                    'strain': name, 'position': position + 1,
                                    'original': original, 'mutated': mutant})
                writer.writerow({'ref': 'HIV1-B-FR-K03455-seed', 'match': '1.0',
                                 'group_ref': 'HIV1-B-FR-K03455-seed',
                                 'contig': ''.join(bases)})
    (folder / 'verified_reference.fasta').write_text(
        ''.join(f">{i}-HIV1-B-FR-K03455-seed {name}\n{seq}\n"
                for i, (name, seq) in enumerate(refs, 1)))
    (folder / 'perturbations.json').write_text(json.dumps(log, indent=2) + '\n')
    return folder


def names_consistent(a, b):
    return a.split()[0].removeprefix('@').removesuffix('/1') == (
           b.split()[0].removeprefix('@').removesuffix('/2'))


def sample_reads(root, label, dest):
    conf = BENCHMARKS[label]
    accession = conf['accession']
    raw_dir = root / conf['raw_dir']
    r1 = raw_dir / (accession + '_1.fastq.gz')
    r2 = raw_dir / (accession + '_2.fastq.gz')
    assert r1.is_file() and r2.is_file()
    outdir = dest / label
    outdir.mkdir(parents=True, exist_ok=True)
    sampled1 = outdir / (label + '_L001_R1_001.fastq.gz')
    sampled2 = outdir / (label + '_L001_R2_001.fastq.gz')
    if sampled1.exists() or sampled2.exists():
        raise FileExistsError('Refusing to overwrite already sampled control')
    rng = random.Random(0x42A319 + len(label))
    threshold = conf['target_pairs'] / conf['read_count']
    chosen = total = 0
    with gzip.open(r1, 'rt') as in1, gzip.open(r2, 'rt') as in2, \
         gzip.open(sampled1, 'wt', compresslevel=3) as out1, \
         gzip.open(sampled2, 'wt', compresslevel=3) as out2:
        while True:
            a = [in1.readline() for _ in range(4)]
            b = [in2.readline() for _ in range(4)]
            if not any(a) and not any(b):
                break
            if not all(a) or not all(b) or not names_consistent(a[0], b[0]):
                raise RuntimeError(f'Paired FASTQ loss/desynchronization at {total}')
            if not a[0].startswith('@') or a[2][:1] != '+' or b[2][:1] != '+':
                raise RuntimeError(f'Invalid original FASTQ layout at {total}')
            total += 1
            if rng.random() < threshold:
                out1.writelines(a)
                out2.writelines(b)
                chosen += 1
    if total != conf['read_count']:
        raise RuntimeError(f'ENA provided {total} pairs, expected {conf["read_count"]}')
    digest=lambda f:hashlib.sha256(f.read_bytes()).hexdigest()
    manifest = {'benchmark': label, 'accession': accession,
                'selected_pairs': chosen, 'original_pairs': total,
                'truth_source': conf['truth_fasta'],
                'truth_sha256': digest(root / conf['truth_fasta']),
                'r1_sha256': digest(sampled1), 'r2_sha256': digest(sampled2)}
    (outdir / 'selection.json').write_text(json.dumps(manifest,indent=2)+'\n')
    return manifest


def run():
    parser=argparse.ArgumentParser()
    parser.add_argument('--benchmark', choices=tuple(BENCHMARKS), required=True)
    parser.add_argument('--panel-only',action='store_true')
    args=parser.parse_args()
    cwd=Path.cwd()
    root=cwd/'.validation-work/independent-truth'
    panel=build_panels(root,args.benchmark)
    print('PANEL',panel,args.benchmark,flush=True)
    if not args.panel_only:
        result=sample_reads(root,args.benchmark,root/'input-samples')
        print('SAMPLED',json.dumps(result),flush=True)


if __name__=='__main__':
    run()
