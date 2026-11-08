"""Audit real-source HCV dilution experiments against original paired-read IDs.

Known input source-library membership is independent of the remapper.
These are NOT lab-certified genotypes or calibrated molecular abundances.
"""
from collections import Counter
import csv
import gzip
import json
from pathlib import Path

from hcv_origin_eval import measure


ROOT = Path(__file__).resolve().parents[2]
STUDY = ROOT / '.validation-work/independent-truth/hcv-origin'
SOURCES = STUDY / 'input-samples'
RESULTS = STUDY / 'results'
PUBLIC = ROOT / 'docs/validation'
CASES = (
    ('HCV-origin-zero', 0),
    ('HCV-origin-0p2-A', 100),
    ('HCV-origin-0p2-B', 100),
    ('HCV-origin-0p4-A', 200),
    ('HCV-origin-0p4-B', 200),
    ('HCV-origin-1pct', 500),
)


def count_input_pairs(case):
    folder = SOURCES / case
    first, = folder.glob('*_L001_R1_001.fastq.gz')
    second, = folder.glob('*_L001_R2_001.fastq.gz')
    counts = Counter()
    with gzip.open(first, 'rt') as r1, gzip.open(second, 'rt') as r2:
        while True:
            one = tuple(r1.readline() for _ in range(4))
            two = tuple(r2.readline() for _ in range(4))
            if not any(one) and not any(two):
                break
            if not all(one) or not all(two):
                raise ValueError(f'{case}: truncated FASTQ stream')
            if one[0].split()[0] != two[0].split()[0]:
                raise ValueError(f'{case}: mates out of synchronization')
            source = one[0].split('|')[0].removeprefix('@')
            if source not in ('HCV1a', 'HCV2b'):
                raise ValueError(f'{case}: unexpected source library')
            counts[source] += 1
    return dict(counts)


def fractional_decoy_evidence(folder):
    by_ref = {}
    with (folder / 'conseq.csv').open(newline='') as stream:
        for row in csv.DictReader(stream):
            if (row['q-cutoff'], row['consensus-percent-cutoff']) != (
                    '15', 'MAX'):
                continue
            if not row['region'].endswith('HCV-1b'):
                continue
            offset = int(row['offset'])
            by_ref[row['region']] = {
                offset + i + 1 for i, base in enumerate(row['sequence'])
                if base in 'ACGT'}
    low_fraction = high_fraction = 0
    sidecar = folder / 'multimap_depth.csv'
    if not sidecar.is_file():
        raise ValueError('Candidate fractional-depth output missing')
    with sidecar.open(newline='') as stream:
        for row in csv.DictReader(stream):
            if int(row['position']) not in by_ref.get(row['refname'], ()):
                continue
            if float(row['avg_depth']) >= 100:
                high_fraction += 1
            else:
                low_fraction += 1
    if high_fraction + low_fraction != sum(map(len, by_ref.values())):
        raise ValueError('Missing decoy depth for a callable position')
    return {'fractional_depth_at_least_100': high_fraction,
            'fractional_depth_below_100': low_fraction}


def analyze():
    complete = {}
    for sample, expected_minor in CASES:
        source = count_input_pairs(sample)
        if (source.get('HCV1a') != 50000 or
                source.get('HCV2b', 0) != expected_minor):
            raise ValueError(f'{sample}: unexpected actual source fraction')
        observed = {}
        for arm in ('upstream',
                    'nonstealing-all-p4' if sample in (
                        'HCV-origin-zero', 'HCV-origin-1pct')
                    else 'nonstealing'):
            result = measure(sample, arm)
            if not result.get('complete'):
                raise ValueError(f'{sample}/{arm}: not a completed experiment')
            original = json.loads(
                (RESULTS / sample / arm / 'status.json').read_text())
            observed[arm] = result
            observed[arm]['sha256_compressed_inputs'] = [
                original['provenance'][f'sha256_r{i}'] for i in (1, 2)]
            observed[arm]['sha256_contigs'] = (
                original['provenance']['sha256_contigs'])
        a, b = observed.values()
        if (a['sha256_compressed_inputs'] != b['sha256_compressed_inputs']
                or a['sha256_contigs'] != b['sha256_contigs']):
            raise ValueError(f'{sample}: experimental inputs differ')
        for result in observed.values():
            del result['sha256_compressed_inputs']
            del result['sha256_contigs']
        candidate = 'nonstealing-all-p4' if sample in (
            'HCV-origin-zero', 'HCV-origin-1pct') else 'nonstealing'
        observed[candidate]['decoy_fractional_evidence'] = (
            fractional_decoy_evidence(RESULTS / sample / candidate))
        complete[sample] = {
            'input_source_pairs': source,
            'source_minor_fraction': expected_minor/(50000+expected_minor),
            'arms': observed,
        }
    return complete


def plot_callable_positions(data):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    names = list(data)
    labels = [name.removeprefix('HCV-origin-').replace('0p', '0.')
              for name in names]
    minor = []
    decoy = []
    for entry in data.values():
        candidate = entry['arms'].get(
            'nonstealing', entry['arms'].get('nonstealing-all-p4'))
        found = candidate['per_reference']
        minor.append(sum(x['callable_q15_max'] for ref, x in found.items()
                         if ref.endswith('HCV-2b')))
        decoy.append(sum(x['callable_q15_max'] for ref, x in found.items()
                         if ref.endswith('HCV-1b')))
    places = range(len(names))
    fig, ax = plt.subplots(figsize=(9.2, 4.8))
    ax.bar([i-0.19 for i in places], minor, 0.38,
           label='Expected minority-origin contig')
    ax.bar([i+0.19 for i in places], decoy, 0.38,
           label='Absent-source genotype decoy')
    ax.set_xticks(list(places), labels)
    ax.set_xlabel('Actual source-library mixture (0, 0.2, 0.4, ~1%)')
    ax.set_ylabel('Callable Q15 MAX contig positions')
    ax.set_title('Real-read minority dilution: consensus vs absent-source decoy')
    ax.legend()
    fig.tight_layout()
    dest = PUBLIC / 'figures/real_hcv_dilution_20261108.svg'
    fig.savefig(dest)
    plt.close(fig)
    dest.write_text('\n'.join(line.rstrip() for line in
                               dest.read_text().splitlines()) + '\n')
    return dest


def make_report(data):
    header = [
        '# Real HCV read-source minority dilution: 0%, 0.2%, 0.4%, 0.99%',
        '',
        ('**Interpretation: known source-library identity, not independently'
        ' certified HCV genotype or full minority-haplotype truth.**'),
        '',
        ('This experiment retains original unaltered reads from archived HCV'
        ' libraries and resamples source-tagged, synchronized R1/R2 pairs.'
        ' The same 50,000 majority-library pairs and fixed three-contig'
        ' mapping panel are used throughout. Two distinct subsamples of 100'
        ' minority-library pairs and two of 200 minority-library pairs'
        ' test empirical variability. The 500-pair control is inherited'
        ' from the earlier experiment.'),
        '',
        '## A/B results',
        '',
        ('| Source-minority fraction | Arm | Known-origin pairs mapped to'
        ' minority contig | Minority Q15 MAX callable positions |'
        ' Absent-library decoy callable | Decoy fractional ≥100 | Runtime |'),
        '|:--|:--|--:|--:|--:|--:|--:|',
    ]
    for sample, info in data.items():
        num = info['input_source_pairs'].get('HCV2b', 0)
        for arm, values in info['arms'].items():
            true_pair = values['counts'].get('HCV2b', {}).get(
                'paired_correct', 0)
            calls = values['per_reference']
            minor = sum(x['callable_q15_max'] for ref, x in calls.items()
                        if ref.endswith('HCV-2b'))
            decoy = sum(x['callable_q15_max'] for ref, x in calls.items()
                        if ref.endswith('HCV-1b'))
            fractional = values.get('decoy_fractional_evidence')
            header.append(
                f"| {sample.replace('HCV-origin-', '')}"
                f" ({num}/{50000+num}) | {arm} | {true_pair}/{num} |"
                f" {minor} | {decoy} | "
                f"{fractional['fractional_depth_at_least_100'] if fractional else 'N/A'}"
                f" | {values['timings']['total_s']:.1f}s |")
    header.extend([
        '',
        '![Real-read source dilution and absent-source decoy](figures/real_hcv_dilution_20261108.svg)',
        '',
        '## Interpretation and release implications',
        '',
        ('- At 0.2% measured read-source fraction (100 input minority pairs),'
        ' both 0.2% subsamples mapped 88 minority-origin pairs correctly'
        ' but produced **zero callable positions** in the minority contig.'
        ' Under a 100× per-base minimum, this is expected when fewer than'
        ' 100 usable original pairs support the entire contig.'),
        ('- At 0.4%, two subsamples mapped 176/200 and 178/200 minority'
        ' pairs correctly but yielded just **1 versus 22 callable'
        ' positions**. This illustrates sequencing and coverage variance.'),
        ('- At 0.99%, 428/500 minority-library source pairs mapped to the'
        ' expected contig and 173 bases were callable in both arms.'),
        ('- The absent-source HCV1b decoy was callable at 26 sites in the'
        ' zero, 0.2% and 0.4% cases, and 32 sites in the 0.99% case.'
        ' All decoy callable sites in the candidate met fractional 100×'
        ' depth; average depth is **not** an adequate genotype-specificity'
        ' filter in this control.'),
        ('- The candidate preserved the same minority-contig callable'
        ' counts and correctly paired minority-source counts as upstream'
        ' across all cases. This shows no new *source-origin* regression'
        ' here but does not prove clinical genotype specificity.'),
        ('- The repeated subsamples share the same underlying archived'
        ' libraries and are statistically correlated; do not pool their'
        ' success counts as if they were independent biological specimens.'),
        ('- No independent LiPA confirmation, molecular barcode abundance'
        ' truth, full-length haplotype phasing, and complete clinical MiCall'
        ' reporting was available. A genotype-presence decision must be'
        ' tested on independently adjudicated negative controls.'),
        '',
        '## Reproducibility and provenance',
        '',
        ('The majority and minority libraries were selected from the'
         ' protected 198-sample real FASTQ archive. The source pool contained'
         ' 50,000 original HCV1a-labelled library pairs and 500'
         ' HCV2b-labelled library pairs with quality scores retained.'
         ' Four conditions were created by synchronized whole-pair'
         ' resampling and shuffling with deterministic seeds:'),
        '',
        '- 0.2% replicate A: 100 minority pairs, seed 1001.',
        '- 0.2% replicate B: 100 minority pairs, seed 2002.',
        '- 0.4% replicate A: 200 minority pairs, seed 3003.',
        '- 0.4% replicate B: 200 minority pairs, seed 4004.',
        '',
        ('The script real_hcv_dilution_20261108.py produces the private'
         ' resampled FASTQs. Every candidate/upstream A/B arm has matching'
         ' compressed R1/R2 and fixed-contig SHA-256 fingerprints,'
         ' verified by real_hcv_dilution_eval_20261108.py. Private'
         ' sequence reads, original specimen identities and per-read'
         ' mapping outputs are not tracked in Git.'),
        '',
        ('Both MiCall branches ran map_to_contigs, sam2aln and'
         ' aln2counts with Bowtie2 2.2.8; upstream source cba1898a,'
         ' current candidate algorithm source ccb56a87.'
         ' Candidate uses four alignment threads versus one upstream.'
         ' The reported times were collected concurrently and cannot'
         ' establish controlled performance equivalence.'),
        '',
        ('The independent evaluator recounts both original mate streams,'
         ' matches the source-library labels and distinguishes fully'
         ' mapped R1/R2 pairs from names with only one mapped end.'
         ' Five unit regressions test this distinction and the'
         ' conservation of original reads and qualities.'),
        '',
        ('**Production verdict: NO-GO.** These controlled real reads establish'
        ' useful behavior but do not demonstrate sensitivity and specificity'
        ' at 0.2% true biological haplotype abundance.'),
        '',
    ])
    return '\n'.join(header)


def main():
    data = analyze()
    plot_callable_positions(data)
    (PUBLIC / 'real_hcv_dilution_20261108_metrics.json').write_text(
        json.dumps(data, indent=2)+'\n')
    (PUBLIC / 'real_hcv_dilution_20261108.md').write_text(make_report(data))
    for case, stats in data.items():
        result = {}
        for arm, values in stats['arms'].items():
            origin = values['counts'].get('HCV2b', {})
            result[arm] = {
                'true_paired': origin.get('paired_correct', 0),
                'callable': sum(v['callable_q15_max'] for ref, v in
                                values['per_reference'].items()
                                if ref.endswith('HCV-2b'))}
        print(case, result)


if __name__ == '__main__':
    main()
