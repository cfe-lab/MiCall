"""Prepare reproducible aggregate reports for exact-source HIV9 clone mixtures."""
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from known_hiv9_mixture_eval_20261108 import evaluate

ROOT = Path(__file__).resolve().parents[2]
PRIVATE = ROOT / '.validation-work/independent-truth'
REPORTS = ROOT / 'docs/validation'
FRACTIONS = {'0': '0%', '0p2': '0.2%', '1': '1%'}


def analyze():
    ref = PRIVATE / 'study-HIV9/verified_reference.fasta'
    output = {}
    for fraction in FRACTIONS:
        sample = 'HIV9-synthetic-' + fraction
        inputs = PRIVATE / 'input-samples' / sample
        arms = {}
        checksums = []
        for arm in ('upstream', 'nonstealing'):
            folder = PRIVATE / 'results' / sample / arm
            result = evaluate(ref, inputs, folder, fraction)
            arms[arm] = result
            status = json.loads((folder/'status.json').read_text())
            checksums.append({key: status['provenance'][key]
                              for key in ('sha256_r1', 'sha256_r2',
                                          'sha256_contigs')})
        if checksums[0] != checksums[1]:
            raise ValueError('Compared arms used different inputs')
        manifest = json.loads((inputs/'synthetic_manifest.json').read_text())
        if manifest['input_pairs'] != arms['upstream']['input_pairs']:
            raise ValueError('Synthetic source manifest does not match input')
        phase = json.loads((inputs/'phase_private.json').read_text())
        if phase['original_pairs'] != arms['upstream']['input_pairs']:
            raise ValueError('Incomplete alignment-independent paired FASTQ phase audit')
        private = phase['per_clone']['9-HIV1-B-FR-K03455-seed']
        output[fraction] = {
            'source_pairs': arms['upstream']['source_labels'],
            'generator_seed': manifest['rng_seed'],
            'read_length': manifest['read_length'],
            'insert_length': manifest['insert_length'],
            'independent_local_phase': {
                'one_or_more': private['at_least_one'],
                'linked_two_or_more': private['at_least_two'],
                'linked_three': private['at_least_three'],
                'marker_positions': private['private_positions'],
            },
            'arms': arms,
        }
    return output


def plot(data):
    x = list(range(len(FRACTIONS)))
    labels = list(FRACTIONS)
    fig, ax = plt.subplots(figsize=(8.4, 4.6))
    for arm, shift in (('upstream', -.18), ('nonstealing', .18)):
        n = [data[k]['arms'][arm]['truth_calls']['clone9']['callable']
             for k in labels]
        ax.bar([i + shift for i in x], n, .36, label=arm)
    for i, label in enumerate(labels):
        clone_pairs = data[label]['source_pairs']['clone9']
        ax.text(i, max(data[label]['arms'][a]['truth_calls']['clone9']['callable']
                       for a in ('upstream', 'nonstealing'))+8,
                str(clone_pairs)+' originating rare-clone pairs', ha='center')
    ax.set_xticks(x, [FRACTIONS[k] for k in labels])
    ax.set_ylim(0, 525)
    ax.set_ylabel('Callable rare-clone Q15 MAX bases')
    ax.set_xlabel('Exact rare-clone source fraction')
    ax.set_title('Independent sequence truth, controlled synthetic source')
    ax.legend()
    fig.tight_layout()
    path = REPORTS/'figures/known_hiv9_mixture_20261108.svg'
    fig.savefig(path)
    plt.close(fig)
    path.write_text('\n'.join(line.rstrip() for line in
                              path.read_text().splitlines())+'\n')


def markdown(data):
    out = [
        '# Independently sequenced HIV9 clones: exact-source minority control',
        '',
        ('**Verdict: NO-GO for clinical production. This is an error-free'
        ' synthetic control, not an independently measured biological mixture.**'),
        '',
        '## Study design',
        '',
        ('Nine independently Sanger-verified HIV clone sequences are used'
        ' as genomic truth. 20,000 error-free 250 nt paired fragments from'
        ' 450 nt inserts are generated using fixed random seeds. Clone9'
        ' has exactly 0, 40, or 200 originating read pairs for 0%,'
        ' 0.2%, or 1% tests; clone1 is adjusted to maintain total depth.'
        ' All other clone counts are held constant. Both MiCall arms'
        ' use identical original FASTQ and fixed Sanger reference panels,'
        ' verified by compressed-input and contig SHA-256.'
        ' The generator RNG seeds are 4500 (0%), 4501 (0.2%) and'
        ' 4502 (1%); the exact input mixture and seed are tracked'
        ' in an ignored source manifest and in aggregate metrics.'),
        '',
        ('The upstream source revision is cba1898a; the revised'
         ' non-stealing mapping code is equivalent to local commit'
         ' fad064abb503774103b9ac16f810c762285de97f.'
         ' The original source FASTQ and contig checksums matched'
         ' between each A/B pair; the challenge used Sanger-truth'
         ' reference sequences as mapping seeds.'),
        '',
        '## Rare-clone source attribution versus consensus',
        '',
        ('| Minority source fraction | Original pairs | Mapper |'
        ' Pairs correctly assigned | Exclusive correct pairs |'
        ' Callable Q15 MAX bases | Wrong Q15 MAX bases | Wall seconds |'),
        '|:--|--:|:--|--:|--:|--:|--:|--:|',
    ]
    for fraction,item in data.items():
        for arm,result in item['arms'].items():
            stats=result['source_assignment']['clone9']
            truth=result['truth_calls']['clone9']
            out.append(
                f'| {FRACTIONS[fraction]} | {item["source_pairs"]["clone9"]}'
                f' | {arm} | {stats.get("paired_to_correct",0)} |'
                f' {stats.get("paired_exclusively_correct",0)} |'
                f' {truth["callable"]} | {truth["wrong"]} |'
                f' {result["run_seconds"]:.1f} |')
    out.extend(['',
                '## Alignment-independent local marker linkage',
                '',
                ('Source-labeled clone9 private 21-mers were counted directly'
                 ' in the original unaligned read pairs, including their'
                 ' reverse complements. Multiple private markers on one'
                 ' sequenced pair demonstrate local linkage without relying'
                 ' on MiCall remapping. This cannot demonstrate full'
                 ' clone-length phase.'),
                '',
                ('| Rare source fraction | One or more private markers |'
                 ' At least two linked | Three linked |'),
                '|:--|--:|--:|--:|'])
    for fraction,item in data.items():
        evidence=item['independent_local_phase']
        out.append(
            f'| {FRACTIONS[fraction]} | {evidence["one_or_more"]} |'
            f' {evidence["linked_two_or_more"]} |'
            f' {evidence["linked_three"]} |')
    out.extend(['',
                ('At a controlled 0.2% fraction, 16 original pairs link'
                 ' at least two of the three clone9-private markers and'
                 ' eight link all three, despite zero callable rare-clone'
                 ' consensus bases. At 0% no clone9-private pair'
                 ' markers were observed; at 1% 91 pairs link two or more.'
                 ' These are alignment-independent synthetic signals, not'
                 ' evidence of a validated biological detection limit.'),
                '',
                ('![Minority clone read origin and callable bases]'
                '(figures/known_hiv9_mixture_20261108.svg)'),
                '',
                '## Full-panel truth concordance and correlated depth',
                '',
                ('| Rare clone fraction | Mapper | Correct Q15 MAX bases |'
                ' Wrong Q15 MAX bases | True clone8 pairs |'
                ' Callable clone8 bases |'),
                '|:--|:--|--:|--:|--:|--:|'])
    for fraction,item in data.items():
        for arm,result in item['arms'].items():
            calls=result['truth_calls']
            out.append(
                f'| {FRACTIONS[fraction]} | {arm} |'
                f' {sum(v["correct"] for v in calls.values()):,} |'
                f' {sum(v["wrong"] for v in calls.values()):,} |'
                f' {item["source_pairs"]["clone8"]} |'
                f' {calls["clone8"]["callable"]} |')
    out.extend([
        '',
        '## Conclusions and limitations',
        '',
        ('- No clone9 reads in the zero-minority control gave zero'
        ' callable clone9 positions for either mapper. This is a'
        ' favorable negative control only in the error-free simulation.'),
        ('- With a true 0.2% fraction (40 read pairs), both implementations'
        ' mapped all 40 pairs to clone9 but produced zero callable'
        ' clone9 consensus bases under the 100x requirement.'),
        ('- With 1% rare-clone origin (200 read pairs), the candidate'
        ' recovered all 200 source pairs and 453 callable clone9'
        ' bases; upstream recovered 196 pairs and 450 callable bases.'),
        ('- Across all three arms neither mapper produced any nucleotide'
        ' disagreement with the matching Sanger references. Because'
        ' those very references are used as contig seeds, this is not'
        ' independent clinical variant adjudication.'),
        ('- Clone8 has just 80 original source read pairs in every'
        ' mixture yet more than 550 Q15 MAX callable bases. Compatible'
        ' reads from other clones may support shared positions,'
        ' so callable maximum-depth bases do not measure'
        ' independent clone-specific molecules or full haplotype phase.'),
        ('- This single constructed replicate per fraction cannot'
        ' establish a clinical detection limit. There are no sequencing'
        ' errors, PCR artifacts, primer bias, or independent measurements'
        ' of a biological mixture. A 450 nt insert does not demonstrate'
        ' full-length linkage of a 1074 nt clone.'),
        ('- The real HCV read-origin dilution study separately shows'
        ' source-pair mapping at 0.2% without any callable'
        ' minority-contig position and apparent support on an absent'
        ' reference. Large-panel real HIV reference loss remains.'),
        '',
        ('Reproduction: the scripts known_hiv9_mixture_20261108.py,'
        ' known_hiv9_mixture_eval_20261108.py and'
        ' known_hiv9_mixture_report_20261108.py regenerate synthetic'
        ' original pairs, independently verify their source labels,'
        ' and summarize MiCall A/B outputs. Private FASTQ and'
        ' per-read alignment files remain Git-ignored.'),
        '',
        '**Clinical production decision: NO-GO.**',
        '',
    ])
    return '\n'.join(out)


def main():
    results = analyze()
    plot(results)
    (REPORTS/'known_hiv9_mixture_20261108_metrics.json').write_text(
        json.dumps(results,indent=2)+'\n')
    (REPORTS/'known_hiv9_mixture_20261108.md').write_text(
        markdown(results))
    print('Completed',len(results),'known-source HIV9 mixture A/B conditions')


if __name__ == '__main__':
    main()
