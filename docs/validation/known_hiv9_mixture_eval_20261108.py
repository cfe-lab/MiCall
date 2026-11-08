"""Score exact-source synthetic HIV9 mixture controls without truth leakage.

Original per-pair clone names are independent of MiCall alignment selection.
Biological clone references are independently verified, but simulated reads
are error-free and represent a favorable, nonclinical performance test.
"""
import argparse
from collections import Counter, defaultdict
import csv
import json
from pathlib import Path
import re

from Bio import SeqIO

READ_ORIGIN = re.compile(r'clone([1-9])_pair[0-9]+\Z')


def consensus(out):
    found = {}
    with (out/'conseq.csv').open(newline='') as f:
        for row in csv.DictReader(f):
            if row['q-cutoff'] != '15' or row['consensus-percent-cutoff'] != 'MAX':
                continue
            found[row['region']] = {
                int(row['offset'])+i+1:base.upper()
                for i,base in enumerate(row['sequence'])
                if base.upper() in 'ACGT'}
    return found


def score_calls(truth, out):
    observed=consensus(out)
    per_ref={}
    for ref,seq in truth.items():
        idx=int(ref.split('-', 1)[0])
        if not 1<=idx<=9:
            raise ValueError('Unexpected Sanger-truth clone reference')
        accepted=observed.get(ref,{})
        correct=sum(accepted.get(i)==b for i,b in enumerate(seq,1))
        wrong=sum(accepted.get(i) not in (None,b) for i,b in enumerate(seq,1))
        per_ref[f'clone{idx}']={
            'correct':correct,'wrong':wrong,
            'uncallable':len(seq)-correct-wrong,
            'callable':correct+wrong,
        }
    return per_ref


def pair_origin(qname):
    m=READ_ORIGIN.fullmatch(qname)
    if not m:
        raise ValueError('Missing independently encoded clone-source label')
    return int(m.group(1))


def pair_support(out, input_counts):
    refs=defaultdict(lambda:defaultdict(int))
    with (out/'remap.csv').open(newline='') as f:
        for row in csv.DictReader(f):
            flag=int(row['flag'])
            name=row['qname']
            pair_origin(name)
            if flag&0x4 or row['rname']=='*' or not row['rname']:
                continue
            if flag&0x40:refs[name][row['rname']]|=1
            if flag&0x80:refs[name][row['rname']]|=2

    stats=defaultdict(Counter)
    for qname,candidates in refs.items():
        source=pair_origin(qname)
        paired={ref for ref,mask in candidates.items() if mask==3}
        expected=f'{source}-HIV1-B-FR-K03455-seed'
        summary=stats[f'clone{source}']
        summary['mapped_with_any_end']+=1
        summary['full_paired']+=bool(paired)
        summary['paired_to_correct']+=expected in paired
        summary['paired_exclusively_correct']+=paired=={expected}
        summary['paired_to_other']+=any(ref!=expected for ref in paired)
        summary['paired_exclusively_wrong']+=bool(paired) and expected not in paired
        summary['unpaired_only']+=not paired
        summary['multi_reference']+=len(paired)>1
    for clone,count in input_counts.items():
        seen=stats[clone]
        seen['source_pairs_input']=count
        seen['not_fully_mapped']=count-seen['full_paired']
        if seen['not_fully_mapped']<0:
            raise ValueError('More mapped pairs than synthetic input pairs')
    return {name:dict(values) for name,values in stats.items()}


def evaluate(truth, source, out, label):
    seqs={r.id:str(r.seq).upper() for r in SeqIO.parse(truth,'fasta')}
    with (source/'synthetic_manifest.json').open() as f:
        manifest=json.load(f)
    expected=manifest['clone_pairs']
    if set(expected)!={f'clone{i}' for i in range(1,10)}:
        raise ValueError('Malformed synthetic source manifest')
    if sum(expected.values())!=manifest['input_pairs']:
        raise ValueError('Input mixture counts do not add up')
    status=json.loads((out/'status.json').read_text())
    if not status.get('ok'):
        raise ValueError('Incomplete MiCall A/B condition')
    return {
        'mixture':label,
        'input_pairs':manifest['input_pairs'],
        'source_labels':expected,
        'truth_calls':score_calls(seqs,out),
        'source_assignment':pair_support(out,expected),
        'run_seconds':round(status['timings']['total_s'],3),
        'note':('Exact-source and Q15 MAX scoring; sequences simulated'
                ' error-free from Sanger-verified clone references.'),
    }


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--truth',type=Path,required=True)
    parser.add_argument('--input-dir',type=Path,required=True)
    parser.add_argument('--output-dir',type=Path,required=True)
    parser.add_argument('--mixture',required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    result=evaluate(args.truth,args.input_dir,args.output_dir,args.mixture)
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print('mixture',result['mixture'],'run_s',result['run_seconds'])
    for clone,stats in result['truth_calls'].items():
        print(clone,'pairs',result['source_labels'][clone],
              'correct',stats['correct'],'wrong',stats['wrong'],
              'assigned_correct',result['source_assignment'][clone].get(
                  'paired_to_correct',0))


if __name__=='__main__':
    main()
