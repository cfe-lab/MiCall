"""Genotype read-origin challenge using two real HCV libraries and absent decoy.

The exact source *library* of each pair is ground truth by construction.
Source library genotype names are only MiCall-inferred, not external
LiPA-confirmed genotypes, so this is not a clinical genotype validation.
"""
from collections import Counter, defaultdict
import csv
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
STUDY=ROOT/'.validation-work/independent-truth/hcv-origin'
OUT=STUDY/'results'


def read_csv(file):
    with file.open(newline='') as stream:
        yield from csv.DictReader(stream)


def measure(sample,arm):
    folder=OUT/sample/arm
    status=folder/'status.json'
    if not status.exists():
        return {'complete':False}
    obj=json.loads(status.read_text())
    if not obj['ok']:
        return {'complete':False,'error':obj.get('error')}
    assignments=defaultdict(set)
    for row in read_csv(folder/'remap.csv'):
        qname=row['qname']
        if row['rname'] and row['rname']!='*' and not int(row['flag'])&0x4:
            assignments[qname].add(row['rname'])
    counts=defaultdict(Counter)
    for name,refs in assignments.items():
        source=name.split('|',1)[0]
        if source not in ('HCV1a','HCV2b'):
            raise AssertionError('Read lost its source-of-origin label')
        expected={'HCV1a':'HCV-1a','HCV2b':'HCV-2b'}[source]
        correctly=any(ref.endswith(expected) for ref in refs)
        decoy=any(ref.endswith('HCV-1b') for ref in refs)
        wrong=any(not ref.endswith(expected) for ref in refs)
        counts[source]['mapped']+=1
        counts[source]['correct']+=correctly
        counts[source]['decoy']+=decoy
        counts[source]['other_reference']+=wrong
        counts[source]['correct_only']+=correctly and not wrong
        counts[source]['wrong_only']+=not correctly and wrong
        counts[source]['multiple_contigs']+=len(refs)>1
    count_regions=Counter()
    callable_bases={}
    for row in read_csv(folder/'conseq.csv'):
        if row['q-cutoff']=='15' and row['consensus-percent-cutoff']=='MAX':
            ref=row['region']
            count_regions[ref]=len(row['sequence'])
            callable_bases[ref]=sum(b in 'ACGT' for b in row['sequence'])
    candidates={}
    for ref,callable_count in callable_bases.items():
        candidates[ref]={'callable_q15_max':callable_count,
                         'span_bases':count_regions[ref]}
    return {'complete': True, 'counts': {k:dict(v) for k,v in counts.items()},
                'per_reference': candidates,'timings': obj.get('timings',{})}


def main():
    metrics={}
    for sample in ('HCV-origin-zero','HCV-origin-1pct'):
        metrics[sample]={arm:measure(sample,arm)
                         for arm in ('upstream','nonstealing')}
    (STUDY/'hcv_origin_metrics.json').write_text(
        json.dumps(metrics,indent=2)+'\n')
    for s,arms in metrics.items():
        print(s)
        for arm,result in arms.items():
            print(' ',arm,'complete',result['complete'],
                  'by_source',result.get('counts'),
                  'callable',result.get('per_reference'))


if __name__=='__main__':
    main()
