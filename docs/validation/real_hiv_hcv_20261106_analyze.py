"""Analyze controlled real HIV/HCV A/B outputs without exposing FASTQ.

Inputs: .validation-work/ab/<alias>/<arm> and private-manifest.json.
Outputs are aggregate-only (no patient IDs or base sequences).
"""
from array import array
from collections import Counter, defaultdict
import csv
import json
from itertools import accumulate
from pathlib import Path
import re

ROOT=Path(__file__).resolve().parents[2]
WORK=ROOT/'.validation-work'
def read_manifest():
    """Read private specimen lookup only during an actual cohort run."""
    return json.loads((WORK/'private-manifest.json').read_text())

OUTPUT=ROOT/'docs/validation'
ARMS=('upstream','nonstealing')
CONSENSUS_FIELDS=('region','q-cutoff','consensus-percent-cutoff','offset','sequence')


def rows(path):
    with path.open(newline='') as f:
        yield from csv.DictReader(f)


def refs_for(out):
    return {r['region']:r['sequence'] for r in rows(out/'remap_conseq.csv')}


def consensus(out):
    return {r['region']:(int(r['offset']),r['sequence'])
            for r in rows(out/'conseq.csv') if r['q-cutoff']=='15'
            and r['consensus-percent-cutoff']=='MAX'}


def remap_stats(out):
    counts=Counter()
    reference_by_pair=defaultdict(set)
    pair_positions=defaultdict(set)
    flag_by_pair=defaultdict(int)
    total_records=0
    for r in rows(out/'remap.csv'):
        total_records+=1
        flag=int(r['flag'])
        name=r['qname']
        if flag & 0x40 and not flag & 0x4:
            counts[name]+=1
            reference_by_pair[name].add(r['rname'])
            pair_positions[name].add((r['rname'],r['pos']))
            flag_by_pair[name]|=1
        if flag & 0x80 and not flag & 0x4:
            flag_by_pair[name]|=2
    count_rows={r['type']:int(r['count']) for r in rows(out/'remap_counts.csv')}
    valid={name for name,mask in flag_by_pair.items() if mask==3}
    multi=sum(counts[name]>1 for name in valid)
    cross=sum(len(reference_by_pair[name])>1 for name in valid)
    within=sum(len(pair_positions[name])>1 and len(reference_by_pair[name])==1
               for name in valid)
    return {'original_pairs': count_rows.get('raw',0)//2,
                'mapped_pairs': len(valid), 'multi_pairs': multi,
                'cross_contig_pairs': cross, 'within_contig_repeat_pairs': within,
                'aligned_mate_rows': total_records,
                'selected_placements': sum(counts.values()),
                'remap_unique_reported': count_rows.get('remap-unique'),
                'unmapped_read_rows': count_rows.get('unmapped',0),
                'remap_final_mate_rows': sum(v for k,v in count_rows.items()
                                          if k.startswith('remap-final '))}


def depth_arrays(out, refseqs):
    diff={ref:array('q',[0])*(len(seq)+2) for ref,seq in refseqs.items()}
    aligned_fragments=0
    out_of_bounds=0
    for row in rows(out/'aligned.csv'):
        if row['qcut']!='15':
            continue
        ref=row['refname']
        aligned_fragments+=int(row['count'])
        if ref not in diff:
            out_of_bounds+=1
            continue
        offset=int(row['offset'])
        for m in re.finditer('[ACGT]+',row['seq']):
            first=offset+m.start()+1
            end=offset+m.end()+1
            if first<1 or end>len(refseqs[ref])+1:
                out_of_bounds+=1
                continue
            diff[ref][first]+=int(row['count'])
            diff[ref][end]-=int(row['count'])
    depths={}
    for ref,values in diff.items():
        values=values[1:len(refseqs[ref])+1]
        depths[ref]=array('q',accumulate(values))
    return depths,aligned_fragments,out_of_bounds


def calc_consensus_delta(baseline,candidate):
    refs=set(baseline)|set(candidate)
    before,after=0,0
    gained,lost,changed=0,0,0
    for ref in refs:
        off1,seq1=baseline.get(ref,(0,''))
        off2,seq2=candidate.get(ref,(0,''))
        d1={off1+i:b for i,b in enumerate(seq1) if b in 'ACGT'}
        d2={off2+i:b for i,b in enumerate(seq2) if b in 'ACGT'}
        before+=len(d1);after+=len(d2)
        gained+=len(d2.keys()-d1.keys())
        lost+=len(d1.keys()-d2.keys())
        changed+=sum(d1[p]!=d2[p] for p in d1.keys() & d2.keys())
    return {'consensus_regions_upstream': len(baseline),
                'consensus_regions_nonstealing': len(candidate),
                'callable_upstream': before,'callable_nonstealing': after,
                'gained_callable': gained,'lost_callable': lost,'discordant_bases': changed}


def calc_coverage_delta(baseline,candidate):
    increase=decrease=unchanged=0
    cross100_up=cross100_down=0
    before_bases=after_bases=0
    sum_before=sum_after=0
    for ref in baseline:
        b,c=baseline[ref],candidate[ref]
        assert len(b)==len(c),(ref,len(b),len(c))
        for left,right in zip(b,c):
            increase+=right>left
            decrease+=right<left
            unchanged+=right==left
            cross100_up+=(left<100 and right>=100)
            cross100_down+=(left>=100 and right<100)
            before_bases+=(left>0)
            after_bases+=(right>0)
            sum_before+=left;sum_after+=right
    return {'depth_positions_increased': increase,'depth_positions_decreased': decrease,
                'depth_positions_unchanged': unchanged,
                'positions_cross_100_up': cross100_up,
                'positions_cross_100_down': cross100_down,
                'positions_with_depth_upstream': before_bases,
                'positions_with_depth_nonstealing': after_bases,
                'depth_sum_upstream': sum_before,'depth_sum_nonstealing': sum_after}


def verify_sidecar(path,derived):
    matched=0
    violations=0
    low_avg_vs100=0
    sum_avg=0.0
    positive_multi=0
    for r in rows(path):
        ref,pos=r['refname'],int(r['position'])
        maxd=int(r['max_depth']);avg=float(r['avg_depth'])
        wanted=derived.get(ref)
        if wanted is None or not 1<=pos<=len(wanted):
            violations+=1;continue
        matched+=1
        if maxd!=wanted[pos-1] or not -1e-8<=avg<=maxd+1e-8:
            violations+=1
        low_avg_vs100+=(maxd>=100 and avg<100)
        positive_multi+=(maxd>avg+1e-7)
        sum_avg+=avg
    # A sidecar that silently omitted positions must not pass merely
    # because all of its *present* rows match reconstructed coverage.
    expected_rows = sum(len(values) for values in derived.values())
    violations += abs(expected_rows - matched)
    return {'sidecar_rows': matched,'sidecar_violations': violations,
                'max_ge100_avg_lt100_positions': low_avg_vs100,
                'positions_fractionally_ambiguous': positive_multi,
                'fractional_depth_sum': round(sum_avg,4)}


def case_metrics(case):
    alias=case['alias']
    outroot=WORK/'ab'/alias
    status={}
    for arm in ARMS:
        f=outroot/arm/'status.json'
        status[arm]=json.loads(f.read_text()) if f.is_file() else None
    if any(st is None or not st.get('ok') for st in status.values()):
        return {'alias': alias,'group': case['cohort'],'panel': case['panel'],
                    'complete': False,'reason': 'missing or unsuccessful A/B arm'}
    paths={arm:outroot/arm for arm in ARMS}
    provenances=[status[arm]['provenance'] for arm in ARMS]
    if any(provenances[0][f]!=provenances[1][f]
           for f in ('sha256_r1','sha256_r2','sha256_contigs')):
        return {'alias': alias,'group': case['cohort'],'panel': case['panel'],
                    'complete': False,'reason': 'A/B input provenance mismatch'}
    # Reject contaminated partial writes, including the discarded pilot.
    if any(b'\\0' in (paths[a]/'remap.csv').open('rb').read(2_000_000)
           for a in ARMS):
        return {'alias': alias,'group': case['cohort'],'panel': case['panel'],
                    'complete': False,'reason': 'NUL bytes in remap stream'}
    refs={arm:refs_for(paths[arm]) for arm in ARMS}
    if refs['upstream']!=refs['nonstealing']:
        return {'alias': alias,'group': case['cohort'],'panel': case['panel'],
                    'complete': False,'reason': 'A/B reference content differs'}
    stats={a:remap_stats(paths[a]) for a in ARMS}
    cov={}
    failures={}
    for arm in ARMS:
        cov[arm],n,invalid=depth_arrays(paths[arm],refs[arm])
        stats[arm]['aligned_q15_fragments']=n
        failures[arm]=invalid
    delta=calc_coverage_delta(cov['upstream'],cov['nonstealing'])
    delta.update(calc_consensus_delta(
        consensus(paths['upstream']),consensus(paths['nonstealing'])))
    sidecar=verify_sidecar(paths['nonstealing']/'multimap_depth.csv',
                           cov['nonstealing'])
    old,new=stats['upstream'],stats['nonstealing']
    result={'alias': alias,'group': case['cohort'],'panel': case['panel'],
                'complete': True,'reference_count': len(refs['upstream']),
                'reference_bases': sum(len(x) for x in refs['upstream'].values()),
                'input_contigs': provenances[0]['contigs_input_count'],
                'input_hash_short': provenances[0]['sha256_r1'][:12],
                'original_read_pairs': old['original_pairs'],
                'baseline_mapped_pairs': old['mapped_pairs'],
                'new_mapped_pairs': new['mapped_pairs'],
                'baseline_placements': old['selected_placements'],
                'new_placements': new['selected_placements'],
                'multi_mapped_pairs': new['multi_pairs'],
                'cross_contig_pairs': new['cross_contig_pairs'],
                'within_contig_repeat_pairs': new['within_contig_repeat_pairs'],
                'mapped_pairs_delta': new['mapped_pairs']-old['mapped_pairs'],
                'placement_inflation': (round(new['selected_placements']/old['selected_placements'],5)
                                     if old['selected_placements'] else None),
                'remap_unique_reported': new['remap_unique_reported'],
                'reported_unique_matches': (new['remap_unique_reported'] is None or
                                          new['remap_unique_reported']//2==new['mapped_pairs']),
                'upstream_remap_seconds': round(status['upstream']['timings']['remap_s'],3),
                'nonstealing_remap_seconds': round(status['nonstealing']['timings']['remap_s'],3),
                'upstream_total_seconds': round(status['upstream']['timings']['total_s'],3),
                'nonstealing_total_seconds': round(status['nonstealing']['timings']['total_s'],3),
                'baseline_out_of_ref_segments': failures['upstream'],
                'candidate_out_of_ref_segments': failures['nonstealing']}
    result.update(delta)
    result.update(sidecar)
    return result


def main():
    metrics=[case_metrics(c) for c in read_manifest()]
    (WORK/'analysis-private.json').write_text(
        json.dumps(metrics,indent=2)+'\n')
    complete=[m for m in metrics if m['complete']]
    if complete:
        # Keep sample-input fingerprints only in the ignored provenance
        # store, not in the published case-level aggregate CSV.
        names = [name for name in complete[0] if name != 'input_hash_short']
        with (OUTPUT/'real_hiv_hcv_20261106_metrics.csv').open('w') as f:
            writer=csv.DictWriter(f,fieldnames=names,extrasaction='ignore',lineterminator='\n')
            writer.writeheader(); writer.writerows(complete)
    print(json.dumps({'total': len(metrics),'complete': len(complete),
                          'incomplete': [(m['alias'],m.get('reason')) for m in metrics
                                      if not m['complete']],
                          'present': [m['alias'] for m in complete],
                          'discordant': sum(m['discordant_bases'] for m in complete),
                          'gained': sum(m['gained_callable'] for m in complete),
                          'lost': sum(m['lost_callable'] for m in complete),
                          'depth_decreased': sum(m['depth_positions_decreased'] for m in complete),
                          'sidecar_violations': sum(m['sidecar_violations'] for m in complete)},
                     indent=2))

if __name__=='__main__': main()
