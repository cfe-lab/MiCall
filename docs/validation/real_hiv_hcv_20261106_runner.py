"""Controlled A/B remapping of archived real sequencing samples.

This script must run in a new interpreter with PYTHONPATH pointed at the
unmodified MiCall checkout (upstream arm) or candidate checkout.
No raw reads or patient-identifying metadata are written into Git.
"""
import argparse
import csv
import gzip
import os
import shutil
from contextlib import ExitStack
import hashlib
import json
from pathlib import Path
import time
import traceback

from micall.core.remap import map_to_contigs
from micall.core.sam2aln import sam2aln
from micall.core.aln2counts import aln2counts


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for part in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(part)
    return h.hexdigest()


def run(args):
    source = Path(args.source_root) / args.sample
    gz1, = source.glob('*_L001_R1_001.fastq.gz')
    gz2, = source.glob('*_L001_R2_001.fastq.gz')
    # The direct Bowtie2 2.2.8 binary in this environment requires plain
    # FASTQ input. Stage decompressed originals once, outside Git tracking.
    staged = Path(args.input_cache) / args.sample
    staged.mkdir(parents=True, exist_ok=True)
    r1, r2 = staged / 'R1.fastq', staged / 'R2.fastq'
    gz_hashes = (digest(gz1), digest(gz2))
    for compressed, plain, fingerprint in (
            (gz1, r1, gz_hashes[0]), (gz2, r2, gz_hashes[1])):
        receipt = staged / (plain.name + '.compressed-sha256')
        if (plain.is_file() and plain.stat().st_size
                and receipt.is_file()
                and receipt.read_text().strip() == fingerprint):
            continue
        # Atomic writes make concurrent upstream/candidate staging safe.
        temporary = plain.with_name(plain.name + '.' + str(os.getpid()))
        with gzip.open(compressed, 'rb') as upstream, temporary.open('wb') as dest:
            shutil.copyfileobj(upstream, dest)
        os.replace(temporary, plain)
        temp_receipt = receipt.with_name(receipt.name + '.' + str(os.getpid()))
        temp_receipt.write_text(fingerprint + '\n')
        os.replace(temp_receipt, receipt)
    contigs = (Path(args.contigs) if args.contigs
               else next(source.glob('unstitched_contigs_*.csv')))
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)
    work = out / 'work'
    work.mkdir(exist_ok=True)
    provenance = {'alias': args.alias, 'arm': args.arm,
                  'sha256_r1': gz_hashes[0], 'sha256_r2': gz_hashes[1],
                  'sha256_contigs': digest(contigs),
                  'baseline_commit': 'cba1898a',
                  'candidate_commit': args.candidate_commit}
    with contigs.open() as f:
        provenance['contigs_input_count'] = sum(1 for _ in csv.DictReader(f))
    (out / 'provenance.json').write_text(json.dumps(provenance, indent=2))
    times = {}
    try:
        start = time.perf_counter()
        with ExitStack() as stack:
            inp = stack.enter_context(contigs.open())
            remap = stack.enter_context((out/'remap.csv').open('w'))
            counts = stack.enter_context((out/'remap_counts.csv').open('w'))
            refs = stack.enter_context((out/'remap_conseq.csv').open('w'))
            um1 = stack.enter_context((out/'unmapped_R1.fastq').open('w'))
            um2 = stack.enter_context((out/'unmapped_R2.fastq').open('w'))
            log = stack.enter_context((out/'bowtie2.log').open('w'))
            map_to_contigs(str(r1), str(r2), inp, remap, counts, refs,
                           um1, um2, str(work), stderr=log, gzip=False)
        times['remap_s'] = time.perf_counter() - start
        start = time.perf_counter()
        with ExitStack() as stack:
            remap = stack.enter_context((out/'remap.csv').open())
            aligned = stack.enter_context((out/'aligned.csv').open('w'))
            ins = stack.enter_context((out/'conseq_ins.csv').open('w'))
            failed = stack.enter_context((out/'failed.csv').open('w'))
            clipping = stack.enter_context((out/'clipping.csv').open('w'))
            extra = {}
            if args.arm == 'nonstealing':
                extra = {'depth_csv': stack.enter_context(
                    (out/'multimap_depth.csv').open('w')),
                    'depth_reference_csv': stack.enter_context(
                        (out/'remap_conseq.csv').open())}
            sam2aln(remap, aligned, ins, failed, clipping, **extra)
        times['sam2aln_s'] = time.perf_counter() - start
        start = time.perf_counter()
        with ExitStack() as stack:
            outputs = {name: stack.enter_context(
                (out / (name+'.csv')).open('w')) for name in (
                    'nuc', 'amino', 'insertions', 'conseq', 'failed_align',
                    'conseq_all', 'nuc_detail', 'genome_coverage')}
            inputs = {name: stack.enter_context(path.open())
                      for name, path in (
                          ('aligned', out/'aligned.csv'),
                          ('clipping', out/'clipping.csv'),
                          ('conseq_ins', out/'conseq_ins.csv'),
                          ('remap_conseq', out/'remap_conseq.csv'),
                          ('contigs', contigs))}
            aln2counts(inputs['aligned'], outputs['nuc'], outputs['amino'],
                       outputs['insertions'], outputs['conseq'],
                       outputs['failed_align'], clipping_csv=inputs['clipping'],
                       conseq_ins_csv=inputs['conseq_ins'],
                       remap_conseq_csv=inputs['remap_conseq'],
                       contigs_csv=inputs['contigs'],
                       conseq_all_csv=outputs['conseq_all'],
                       nuc_detail_csv=outputs['nuc_detail'],
                       genome_coverage_csv=outputs['genome_coverage'])
        times['aln2counts_s'] = time.perf_counter() - start
        times['total_s'] = sum(times.values())
        status = {'ok': True, 'timings': times, 'provenance': provenance}
    except Exception as e:  # noqa: BLE001 - preserve evidence of all runner failures
        status = {'ok': False, 'timings': times, 'provenance': provenance,
                      'error': str(e), 'traceback': traceback.format_exc()}
    (out/'status.json').write_text(json.dumps(status, indent=2)+'\n')
    print(json.dumps({'alias':args.alias,'arm':args.arm,
                      'ok':status['ok'],**times,
                      'error':status.get('error','')}),flush=True)
    if not status['ok']:
        raise RuntimeError(status['error'])


if __name__ == '__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--source-root', required=True)
    parser.add_argument('--sample',required=True)
    parser.add_argument('--alias',required=True)
    parser.add_argument('--arm',choices=('upstream','nonstealing'),required=True)
    parser.add_argument('--output',required=True)
    parser.add_argument('--candidate-commit',required=True)
    parser.add_argument('--input-cache', default=str(Path('.validation-work/input-cache').resolve()))
    parser.add_argument('--contigs')
    run(parser.parse_args())
