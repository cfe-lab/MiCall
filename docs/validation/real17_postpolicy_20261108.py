"""Summarize the reference-resource policy on genuine HIV paired reads.

17 combined independently sequenced viral reference contigs are an engineered
competition stress panel, not one patient's authentic de novo assembly.
"""
import csv
import json
from pathlib import Path

from independent_truth_eval import truth_data, calls_by_region


ROOT = Path(__file__).resolve().parents[2]
PRIVATE = ROOT / '.validation-work/independent-truth'
PUBLIC = ROOT / 'docs/validation'


def load(arm):
    path = PRIVATE / 'results' / ('HIV9-real17-'+arm) / 'nonstealing'
    status = json.loads((path/'status.json').read_text())
    if not status.get('ok'):
        raise ValueError('Incomplete real17 mapping arm')
    calls = calls_by_region(path / 'conseq.csv')
    metrics = {'timings': status['timings']}
    metrics['input_fingerprint'] = {
        item: status['provenance'][item] for item in (
            'sha256_r1', 'sha256_r2', 'sha256_contigs')}
    metrics['reported_final_placements'] = 0
    metrics['capped_original_pair_groups'] = 0
    with (path/'remap_counts.csv').open() as f:
        for row in csv.DictReader(f):
            if row['type'].startswith('remap-final '):
                metrics['reported_final_placements'] += int(row['count'])
            if row['type'] == 'multimap-capped-pairs':
                metrics['capped_original_pair_groups'] += int(row['count'])
    truth, names, _ = truth_data('HIV9')
    metrics['clone_calls'] = {}
    for ref,seq in truth.items():
        assigned = calls.get(ref,{})
        correct = sum(assigned.get(i) == b for i,b in enumerate(seq,1))
        wrong = sum(assigned.get(i) not in (None,b)
                    for i,b in enumerate(seq,1))
        metrics['clone_calls'][names[ref]] = {
            'correct': correct, 'wrong': wrong,
            'uncallable': len(seq)-correct-wrong,
        }
    return metrics


def main():
    old = load('before')
    new = load('after')
    if old['input_fingerprint'] != new['input_fingerprint']:
        raise ValueError('Input FASTQ or fixed contig panel changed')
    del old['input_fingerprint']
    del new['input_fingerprint']
    data = {'references': 17, 'authentic_hiv_pair_inputs': 1200,
            'status': ('Controlled fixed-contig mapping experiment;'
                       ' no clinical genotype-presence truth'),
            'previous_candidate': old, 'resource_aware_candidate': new}
    (PUBLIC/'real17_postpolicy_20261108_metrics.json').write_text(
        json.dumps(data,indent=2)+'\n')
    heading = [
        '# Real HIV9 fixed-contig 17-reference remapping comparison',
        '',
        ('**This is a constructed, 17-authentic-reference competition panel,'
        ' not a 17-contig clinical specimen or a full pipeline validation.**'),
        '',
        ('The same 1,200 real paired HIV9 reads and 17 fixed authentic'
        ' references (nine Sanger-verified HIV clones, five independently'
        ' sequenced HIV strains and three archived HCV contigs) were mapped'
        ' with two versions of the candidate. The former version capped'
        ' large panels at -k 8; the new version permits all Bowtie2'
        ' reported hits for 17–32 references whose total length is at most'
        ' 65,536 nt. FASTQ and contig SHA-256 checksums matched.'),
        '',
        '## Algorithm revision provenance',
        '',
        ('The previous branch was at ec50f74c and used the fixed'
         ' eight-hit search for 17 contigs. The revised mapping source'
         ' was subsequently committed unchanged as fad064abb503774103b9ac16f810c762285de97f.'
         ' The private runner labels the two execution environments'
         ' search-policy-before and search-policy-after. The exact'
         ' original FASTQ and reference CSV SHA-256 values were compared'
         ' locally but are not published to avoid linking private inputs.'
         ' These runs were concurrent and not a controlled time trial.'),
        '',
        '## Independently sequenced HIV9 clone truth',
        '',
        ('| Clone | Before correct | After correct | Delta |'
        ' Before wrong | After wrong |'),
        '|:--|--:|--:|--:|--:|--:|',
    ]
    for name in old['clone_calls']:
        b=old['clone_calls'][name]; a=new['clone_calls'][name]
        heading.append(
            f'| {name} | {b["correct"]} | {a["correct"]} |'
            f' {a["correct"]-b["correct"]:+} | {b["wrong"]} |'
            f' {a["wrong"]} |')
    for title,field in (
        ('True HIV9 correct callable bases','correct'),
        ('Incorrect callable bases','wrong')):
        a=sum(x[field] for x in old['clone_calls'].values())
        b=sum(x[field] for x in new['clone_calls'].values())
        heading.append(f'- **{title}:** {a:,} before versus {b:,} after.')
    heading.extend([
        '',
        '## Alignment accounting and observed runtime',
        '',
        '| Metric | Previous candidate | Resource-aware candidate |',
        '|:--|--:|--:|',
        (f'| Final mapped-mate placement rows |'
        f' {old["reported_final_placements"]:,} |'
        f' {new["reported_final_placements"]:,} |'),
        (f'| Saturated original read-pair groups |'
        f' {old["capped_original_pair_groups"]:,} |'
        f' {new["capped_original_pair_groups"]:,} |'),
        (f'| Total wall seconds | {old["timings"]["total_s"]:.2f} |'
        f' {new["timings"]["total_s"]:.2f} |'),
        '',
        ('**Nuance:** total correct bases improved but not monotonically on'
        ' each clone. The best-sequence set changes when additional'
        ' alignment alternatives are searched; the 1% and 0.2%'
        ' sample frequencies in the published HIV9 benchmark are'
        ' not independently established for this real-read subsample.'),
        '',
        ('No variant correctness claim beyond the Sanger-verified'
        ' contig-level truth is supported. Bowtie2 all-hit search is'
        ' heuristic, and the larger-than-32 or >65,536 nt fallback'
        ' retains a finite -k 8 limit. The 30-minute processing budget'
        ' has not been established on broad 198-sample production panels.'),
        '',
        '**Clinical production decision: NO-GO.**',
        '',
    ])
    (PUBLIC/'real17_postpolicy_20261108.md').write_text(
        '\n'.join(heading))
    print('17-ref final correct',
          sum(v['correct'] for v in new['clone_calls'].values()),
          'wrong',sum(v['wrong'] for v in new['clone_calls'].values()))


if __name__ == '__main__':
    main()
