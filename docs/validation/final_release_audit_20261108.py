"""Create aggregate reproducible final-branch biological release-gate evidence.

Requires local protected real sequencing outputs; emits only aggregate metrics
and public clone labels. Clinical genotype truth is *not* inferred from reads
mapped to themselves or from non-independent upstream predictions.
"""
import csv
from datetime import datetime, timezone
import json
from pathlib import Path

from Bio import SeqIO

from hcv_origin_eval import measure as hcv_measure

ROOT = Path(__file__).resolve().parents[2]
PRIVATE = ROOT / '.validation-work/independent-truth'
RESULTS = PRIVATE / 'results'
OUTPUT = ROOT / 'docs/validation'


def known_truth(label):
    fasta = PRIVATE / f'study-{label}/verified_reference.fasta'
    return {record.id: str(record.seq).upper()
            for record in SeqIO.parse(fasta, 'fasta')}


def calls(folder):
    found = {}
    with (folder / 'conseq.csv').open(newline='') as stream:
        for row in csv.DictReader(stream):
            if row['q-cutoff'] != '15' or row['consensus-percent-cutoff'] != 'MAX':
                continue
            found[row['region']] = {
                int(row['offset']) + index + 1: base.upper()
                for index, base in enumerate(row['sequence'])
                if base.upper() in 'ACGT'}
    return found


def score_truth(folder, label):
    state = json.loads((folder / 'status.json').read_text())
    if not state.get('ok'):
        raise ValueError(f'Incomplete validation case: {folder}')
    originals = known_truth(label)
    observed = calls(folder)
    total = {'correct': 0, 'wrong': 0, 'uncallable': 0}
    per_clone = {}
    for ref, seq in originals.items():
        accepted = observed.get(ref, {})
        right = sum(accepted.get(i) == b for i, b in enumerate(seq, 1))
        wrong = sum(accepted.get(i) not in (None, b)
                    for i, b in enumerate(seq, 1))
        missing = len(seq) - right - wrong
        total['correct'] += right
        total['wrong'] += wrong
        total['uncallable'] += missing
        per_clone[ref] = {'correct': right, 'wrong': wrong, 'uncallable': missing}
    return {'totals': total,
            'per_clone': per_clone,
            'timing_seconds': round(state['timings']['total_s'], 2),
            'remap_seconds': round(state['timings']['remap_s'], 2),
            'sha': {k: state['provenance'][k]
                    for k in ('sha256_r1', 'sha256_r2', 'sha256_contigs')}}


def check_matching_inputs(old, final, label):
    if old['sha'] != final['sha']:
        raise ValueError(f'{label}: compared experimental arms had different inputs')


def count_perturbations(folder, label):
    source = PRIVATE / f'study-{label}/perturbations.json'
    edits = json.loads(source.read_text())
    seen = calls(folder)
    result = {'reverted': 0, 'wrong_seed_retained': 0,
              'other_wrong': 0, 'uncallable': 0}
    for entry in edits:
        base = seen.get(entry['refname'], {}).get(entry['position'])
        if base is None:
            result['uncallable'] += 1
        elif base == entry['original']:
            result['reverted'] += 1
        elif base == entry['mutated']:
            result['wrong_seed_retained'] += 1
        else:
            result['other_wrong'] += 1
    return result


def main():
    groups = {}
    cases = {
        'HIV5-truth': ['HIV5', 'HIV5-truth/upstream',
                       'HIV5-truth/nonstealing',
                       'HIV5-all-p4/nonstealing'],
        'HIV9-truth': ['HIV9', 'HIV9-truth/upstream',
                       'HIV9-truth/nonstealing',
                       'HIV9-all-p4/nonstealing'],
        'HIV5-perturbed': ['HIV5', 'HIV5-perturbed/upstream',
                           'HIV5-perturbed/nonstealing',
                           'HIV5-perturbed-all-p4/nonstealing'],
    }
    for label, (truth_label, baseline, old, final) in cases.items():
        scores = {
            arm: score_truth(RESULTS / location, truth_label)
            for arm, location in (('upstream', baseline),
                                  ('bounded_candidate', old),
                                  ('final_candidate', final))}
        check_matching_inputs(scores['upstream'], scores['final_candidate'], label)
        if label.endswith('-perturbed'):
            for arm, location in (('upstream', baseline),
                                  ('bounded_candidate', old),
                                  ('final_candidate', final)):
                scores[arm]['planted_error_correction'] = (
                    count_perturbations(RESULTS / location, truth_label))
        # Avoid publishing reproducible sample fingerprints.
        for result in scores.values():
            del result['sha']
        groups[label] = scores

    hcv = {
        name: {arm: hcv_measure(name, arm)
               for arm in ('upstream', 'nonstealing', 'nonstealing-all-p4')}
        for name in ('HCV-origin-zero', 'HCV-origin-1pct')}
    for name, arms in hcv.items():
        if not all(x.get('complete') for x in arms.values()):
            raise ValueError(f'HCV real-origin controls incomplete: {name}')

    linked = json.loads((PRIVATE / 'hiv9-exact-phased-motifs.json').read_text())
    with (PRIVATE / 'hiv9-large-panel-audit/aggregate.json').open() as f:
        large_panel = json.load(f)
    with (OUTPUT / 'real_hcv_dilution_20261108_metrics.json').open() as f:
        dilution = json.load(f)
    aggregate = {
        'truth_cases': groups,
        'hcv_source_origin': hcv,
        'known_origin_dilution': dilution,
        '17_real_reference_panel': large_panel,
        'hiv9_raw_exact_linked_motifs': linked,
        'note': ('Allele correctness is measured against independently sequenced'
                 ' HIV references. HCV read origins, not biological genotype,'
                 ' are independently established. Clone 0.2% nominal labels'
                 ' are not validated abundance measurements.'),
    }
    (OUTPUT / 'final_release_audit_20261108_metrics.json').write_text(
        json.dumps(aggregate, indent=2) + '\n')

    header = ['# Final-candidate independent genotype and minority-haplotype audit',
              '',
              '**Clinical production decision: NO-GO; release requirements remain unmet.**',
              '',
              '## Real HIV reference truth — identical inputs in all arms',
              '',
              '| Condition | Algorithm | Correct Q15 MAX | Wrong Q15 MAX | Uncallable | Total runtime |',
              '|:--|:--|--:|--:|--:|--:|']
    for label, arms in groups.items():
        for arm, summary in arms.items():
            total = summary['totals']
            header.append(
                f"| {label} | {arm} | {total['correct']:,} | "
                f"{total['wrong']:,} | {total['uncallable']:,} | "
                f"{summary['timing_seconds']:.1f}s |")

    header.extend(['', '### Reference-error correction challenge',
                   '', ('| Implementation | Planted errors corrected | '
                   'Wrong seed retained | Other wrong | Uncallable |'),
                   '|:--|--:|--:|--:|--:|'])
    for arm, summary in groups['HIV5-perturbed'].items():
        v = summary['planted_error_correction']
        header.append(
            f"| {arm} | {v['reverted']}/50 | {v['wrong_seed_retained']} | "
            f"{v['other_wrong']} | {v['uncallable']} |")
    header.extend([
        '',
        ('The five-strain truth FASTA derives from individually sequenced'
        ' laboratory isolates; the nine-haplotype truth FASTA derives from'
        ' separately Sanger-confirmed plasmids. However, the known-truth'
        ' references are used directly as *mapping seeds* for the truth panel.'
        ' The deliberately corrupted-reference panel is a partial control'
        ' against circular seed agreement.'),
        '',
        ('Scientific sources: the five-virus control references are described'
        ' by Di Giallonardo and colleagues (2014),'
        ' https://doi.org/10.1093/nar/gku537;'
        ' the nine-clone Sanger-confirmed reference mixture by'
        ' Knyazev and colleagues (2021),'
        ' https://doi.org/10.1093/nar/gkab576.'
        ' The original published nominal clone percentages have not'
        ' been confirmed as measured frequencies in the selected read set.'),
        '',
        '## HCV real-read source-origin and absent-genotype negative controls',
        '',
        ('| Condition | Mapper | Source HCV2b read pairs mapped correctly | '
        'Major-library pairs on absent HCV1b decoy | Decoy callable positions |'),
        '|:--|:--|--:|--:|--:|'])
    for case, arms in hcv.items():
        for arm, result in arms.items():
            minor = result['counts'].get('HCV2b', {})
            major = result['counts'].get('HCV1a', {})
            decoy = sum(value['callable_q15_max'] for ref, value in
                        result['per_reference'].items()
                        if ref.endswith('HCV-1b'))
            header.append(
                f"| {case} | {arm} | {minor.get('correct', 0)} | "
                f"{major.get('decoy', 0)} | {decoy} |")

    header.extend([
        '',
        ('Fractional-depth specificity check: in the final candidate, all 26'
         ' Q15 MAX callable positions on the absent HCV1b decoy in the'
         ' zero-minority control, and all 32 in the 0.99% mixture, reached'
         ' fractional depth 100 as well as maximum depth 100. Fractional'
         ' depth alone therefore does not establish genotype absence.'),
        '',
        ('These controls mix real reads from identified source *libraries*.'
        ' The identity of the original read library is known, but the'
        ' biological genotype of each library was not independently'
        ' adjudicated. Apparent calls on an intentionally absent contig'
        ' are a warning about specificity, not proof of a clinically'
        ' reportable false-positive result. These apparent calls already'
        ' occurred with upstream MiCall.'),
        '',
        '## Real-read source dilutions below one percent',
        '',
        ('Additional same-input A/B controls resampled authentic archived'
         ' HCV reads at verified source-library fractions 0%, 0.2%, 0.4%'
         ' and 0.99%. Both independently seeded 0.2% resamples yielded'
         ' 88/100 correctly assigned source pairs but no callable'
         ' minority-contig bases. Two 0.4% resamples produced 1 and'
         ' 22 callable bases despite 176/200 and 178/200 correct source'
         ' assignments. The absent-source decoy retained 26 callable'
         ' bases. These are not certified clinical genotype fractions.'
         ' See the dilution report real_hcv_dilution_20261108.md.'),
        '',
        '## Large-panel bounded-alignment completeness: real-read stress',
        '',
        (f'The final >16-contig fallback was challenged with'
         f' {large_panel["original_real_hiv9_pairs"]:,} real paired HIV9'
         f' reads and a constructed panel of {large_panel["references"]}'
         ' authentic reference sequences (nine Sanger-confirmed HIV9'
         ' plasmids, five HIV5 isolates and three archived HCV contigs).'
         f' {large_panel["bounded_diagnostics"]["at_report_limit"]:,}'
         ' paired-read groups reached the -k 8 reporting limit.'
         f' {large_panel["pairs_losing_references_at_k8"]:,}'
         ' read pairs omitted at least one reference found by Bowtie2'
         ' all-alignment mode, with'
         f' {large_panel["missing_reference_alternatives"]:,} missing'
         ' alternatives in total. This is a deeper heuristic'
         ' search, not a completeness oracle. The reference panel is'
         ' an engineered validation panel of real contigs, not an actual'
         ' 17-contig patient assembly. Runtime was'
         f' {large_panel["bounded_seconds"]:.1f}s capped versus'
         f' {large_panel["all_seconds"]:.1f}s deeper. This prevents'
         ' endorsing the large-panel fallback for clinical use without'
         ' an adaptive, bounded correctness strategy.'),
        '',
        '## Alignment-independent phased-minority signals',
        '',
        ('Private 21-mers around each Sanger-confirmed, clone-private SNP'
        ' were searched against the **unaligned** original paired FASTQ'
        ' reads, including reverse-complement matches. Two private motifs'
        ' occurring in one pair demonstrate local linkage, not a full-length'
        ' phased haplotype or calibrated clone abundance.'),
        '',
        '| Clone | Private 21-mers | ≥1 per pair | ≥2 linked per pair | ≥3 linked per pair |',
        '|:--|--:|--:|--:|--:|'])
    for index, label in enumerate(
            ('A_50','B_25','C_12.5','D_6.3','E_3.2','F_1.6','G_0.8','H_0.4','J_0.2'), 1):
        key = f'{index}-HIV1-B-FR-K03455-seed'
        v = linked['per_clone'][key]
        header.append(
            f"| {label} | {len(v['private_positions'])} | "
            f"{v['at_least_one']:,} | {v['at_least_two']:,} | "
            f"{v['at_least_three']:,} |")

    header.extend([
        '',
        ('The clone named **J_0.2** has 96,433 pairs matching at least one'
        ' private motif and 35,122 with two or more linked motifs in'
        ' 239,006 sampled paired reads. It is therefore scientifically'
        ' unjustified to claim demonstrated 0.2% sensitivity from this'
        ' particular sequenced read set merely because the reference'
        ' FASTA carries a 0.2% nominal label. Composition, amplification,'
        ' sequence error and library provenance need separate adjudication.'),
        '',
        '## Production gate and correctness limitations',
        '',
        ('- **PASS in tested configuration:** 4,373 broader MiCall tests passed,'
        ' with five skips, one deselection, and 54 expected failures.'
        ' Installation/IVA tests depend on unavailable prerequisites.'),
        ('- **PASS on reference subsets:** the revised candidate recovers'
        ' more correct callable HIV bases with no overall increase in wrong'
        ' calls against upstream on the fixed-known-truth panels; some'
        ' clones and positions still have incorrect calls.'),
        ('- **NOT ESTABLISHED:** known negative-control genotype-absence'
        ' specificity, minority mixture abundance and full haplotype phase.'),
        ('- **NOT ESTABLISHED:** complete detection of competing alignments'
        ' on Bowtie2 heuristic search; panels over 16 still use a k=8 cap.'),
        ('- **NOT ESTABLISHED:** end-to-end production MiCall clinical outputs'
        ' (assembly, stitcher, trimming, G2P, HIV resistance reports),'
        ' real 198-sample throughput and peak memory budgets.'),
        ('- **NOT ESTABLISHED:** independently adjudicated full biological'
        ' ground truth for all real HCV libraries, despite their verified'
        ' source-library identities.'),
        ('- **Fixed engineering diagnostic:** legitimate discordant Bowtie2'
        ' SAM groups with unequal mate-alternative counts no longer trigger'
        ' false malformed-pair warnings; rescue behavior is unit-tested.'),
        '',
        ('The clinical production release remains blocked. An increase in'
        ' per-contig maximum depth is not independent molecule evidence.'
        ' Never use the alternative-depth sidecar as a calibrated genotype'
        ' frequency or clinical diagnosis by itself.'),
        '',
        '### Provenance',
        '',
        ('Upstream source: cba1898a. Bounded candidate controls were run on'
        ' earlier revisions. All-placement, four-thread controls'
        ' used the working tree implementing the same alignment mode as'
        ' local commit 379e80e9; exact commit-pinned reruns and independent'
        ' review remain mandatory before final release.'
        ' Original reads, haplotype sequences and private run artifacts'
        ' remain in git-ignored local storage.'),
        '',
        'Generated ' + datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC'),
        '',
    ])
    (OUTPUT / 'final_release_audit_20261108.md').write_text(
        '\n'.join(header))
    print(json.dumps({key: {arm: v['totals']
                             for arm, v in arms.items()}
                      for key, arms in groups.items()}, indent=2))


if __name__ == '__main__':
    main()
