"""Score MiCall A/B against independently established laboratory HIV truth.

The five-virus references were individually resequenced before mixing.
The nine-plasmid haplotypes were individually Sanger-confirmed.
Mapping the same truth sequence as reference is not enough to claim
accuracy: perturbed-reference challenges test actual correction.
"""
import csv
import json
from pathlib import Path

from Bio import SeqIO

ROOT = Path(__file__).resolve().parents[2]
PRIVATE = ROOT / '.validation-work/independent-truth'
RESULTS = PRIVATE / 'results'
REPORT = ROOT / 'docs/validation'


def csv_rows(path):
    with path.open(newline='') as f:
        yield from csv.DictReader(f)


def truth_data(control):
    folder = PRIVATE / ('study-' + control)
    sequences = {record.id: str(record.seq).upper() for record
                 in SeqIO.parse(folder / 'verified_reference.fasta', 'fasta')}
    strain_names = {record.id:record.description.split(' ', 1)[1]
                    for record in SeqIO.parse(
                        folder / 'verified_reference.fasta', 'fasta')}
    edits = json.loads((folder / 'perturbations.json').read_text())
    return sequences, strain_names, edits


def calls_by_region(path):
    # Coordinates in conseq.csv are reference-relative (offset + index).
    calls = {}
    for row in csv_rows(path):
        if row['q-cutoff'] == '15' and row['consensus-percent-cutoff'] == 'MAX':
            calls[row['region']] = {
                int(row['offset']) + i + 1: base.upper()
                for i, base in enumerate(row['sequence'])
                if base.upper() in ('A', 'C', 'G', 'T')
            }
    return calls


def marker_positions(truth):
    if len({len(seq) for seq in truth.values()}) != 1:
        return {}
    return {ref: [i + 1 for i, base in enumerate(seq)
                  if all(base != alt[i] for ref2, alt in truth.items()
                         if ref2 != ref)]
            for ref, seq in truth.items()}


def depths_by_region(path, markers):
    if not path.is_file():
        return {}
    marker_positions = {
        ref: set(positions) for ref, positions in markers.items()}
    result = {}
    for row in csv_rows(path):
        ref = row['refname']
        position = int(row['position'])
        if position in marker_positions.get(ref, ()):
            result[(ref,position)] = (
                int(row['max_depth']),float(row['avg_depth']))
    return result


def score_one(control, panel, arm, sequences, edits):
    folder = RESULTS / (control + '-' + panel) / arm
    status_path = folder / 'status.json'
    # The nine-haplotype comparison may run in eight verified contiguous
    # partitions to avoid a prohibitively long single-thread Bowtie2 run.
    # The merged shard output is analyzed only after every partition has
    # completed and all original read bytes have been accounted for.
    if (control == 'HIV9' and panel == 'truth' and arm == 'nonstealing'
            and not status_path.is_file()):
        alternate = RESULTS / 'HIV9-truth-sharded' / 'nonstealing'
        if (alternate / 'status.json').is_file():
            folder = alternate
            status_path = folder / 'status.json'
    if not status_path.is_file():
        return {'complete': False, 'reason': 'missing status'}
    status = json.loads(status_path.read_text())
    if not status['ok']:
        return {'complete': False, 'reason': status.get('error', 'failed')}
    calls = calls_by_region(folder / 'conseq.csv')
    exact = wrong = missing = 0
    per_strain = {}
    for ref, truth in sequences.items():
        observed = calls.get(ref, {})
        good = sum(observed.get(i) == base for i, base
                   in enumerate(truth,1))
        bad = sum(observed.get(i) not in (None,base) for i,base
                  in enumerate(truth,1))
        no_call = len(truth)-good-bad
        exact += good
        wrong += bad
        missing += no_call
        per_strain[ref] = {'correct': good, 'incorrect': bad,
                               'not_callable': no_call,'total': len(truth)}
    mutation_states={'reverted_to_truth':0,'retained_wrong_seed':0,
                     'another_wrong_base':0,'not_callable':0}
    if panel == 'perturbed':
        for row in edits:
            observed = calls.get(row['refname'], {}).get(row['position'])
            if observed == row['original']:
                mutation_states['reverted_to_truth'] += 1
            elif observed == row['mutated']:
                mutation_states['retained_wrong_seed'] += 1
            elif observed is None:
                mutation_states['not_callable'] += 1
            else:
                mutation_states['another_wrong_base'] += 1
    markers=marker_positions(sequences)
    depth=depths_by_region(folder/'multimap_depth.csv',markers)
    marker_result={}
    for ref, positions in markers.items():
        count = len(positions)
        called = sum(calls.get(ref,{}).get(i)==sequences[ref][i-1]
                     for i in positions)
        wrong_call = sum(
            calls.get(ref,{}).get(i) not in
            (None,sequences[ref][i-1]) for i in positions)
        marker_result[ref] = {
            'loci': count,'truth_called': called,'wrong_called': wrong_call,
            'missing': count-called-wrong_call,
            'max_ge100': sum(depth.get((ref,i),(0,0))[0]>=100
                          for i in positions),
            'fractional_ge100': sum(depth.get((ref,i),(0,0))[1]>=100
                                 for i in positions)}
    return {'complete': True,'totals': {'correct': exact,'incorrect': wrong,
                                          'not_callable': missing,
                                          'total': sum(map(len,sequences.values()))},
                'edits': mutation_states,
                'per_strain': per_strain,
                'markers': marker_result,
                'timing': status.get('timings',{}),
                'input_hashes': {k:status['provenance'][k]
                              for k in ('sha256_r1','sha256_r2',
                                        'sha256_contigs')}}


def compare_panel(baseline, candidate,truth):
    if not baseline.get('complete') or not candidate.get('complete'):
        return {}
    # Compare truth correctness at every expected genomic reference base,
    # not only positions called by both arms.
    return {
        'correct_gain': candidate['totals']['correct']-baseline['totals']['correct'],
        'wrong_gain': candidate['totals']['incorrect']-baseline['totals']['incorrect'],
        'callable_gain': sum(candidate['totals'][k]-baseline['totals'][k]
                          for k in ('correct','incorrect')),
        'corrected_seeds_gain': (
            candidate['edits']['reverted_to_truth']-
            baseline['edits']['reverted_to_truth'])}


def run():
    results={}
    for control in ('HIV5','HIV9'):
        sequences,names,edits=truth_data(control)
        results[control]={'truth_sources':names,
                          'total_haplotype_bases':sum(map(len,sequences.values())),
                          'known_seed_errors':len(edits),'panels':{}}
        for panel in ('truth','perturbed'):
            arms={arm:score_one(control,panel,arm,sequences,edits)
                  for arm in ('upstream','nonstealing')}
            if (all(arms[a].get('complete') for a in ('upstream','nonstealing'))
                    and any(arms['upstream']['input_hashes'][k] !=
                            arms['nonstealing']['input_hashes'][k]
                            for k in ('sha256_r1','sha256_r2','sha256_contigs'))):
                raise ValueError('A/B source read or reference hash mismatch')
            arms['diff']=compare_panel(arms['upstream'],arms['nonstealing'],
                                       sequences)
            results[control]['panels'][panel]=arms
    (PRIVATE/'independent_truth_metrics.json').write_text(
        json.dumps(results,indent=2)+'\n')
    for control,data in results.items():
        for panel,arms in data['panels'].items():
            print(control,panel)
            for arm in ('upstream','nonstealing'):
                r=arms[arm]
                print(' ',arm,'complete',r.get('complete'),
                      'totals',r.get('totals'),'edits',r.get('edits'))
            print('  A/B delta',arms['diff'])
            if control=='HIV9' and panel=='truth':
                for arm in ('upstream','nonstealing'):
                    print('  ',arm,'HIV9 marker detection',
                          arms[arm].get('markers'))


if __name__=='__main__':
    run()
