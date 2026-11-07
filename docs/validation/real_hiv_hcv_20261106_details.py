"""Private, position-level audit of Q15 MAX consensus changes."""
import json

from real_hiv_hcv_20261106_analyze import (
    WORK, read_manifest, consensus, refs_for, depth_arrays, rows,
)


def as_calls(record):
    offset,sequence=record
    return {offset+i:b for i,b in enumerate(sequence)
            if b in 'ACGT'}


def run(case):
    alias=case['alias']
    output=WORK/'ab'/alias
    if not all((output/arm/'status.json').is_file()
               and json.loads((output/arm/'status.json').read_text()).get('ok')
               for arm in ('upstream','nonstealing')):
        return {'alias': alias,'complete': False}
    baseline=output/'upstream';candidate=output/'nonstealing'
    refs=refs_for(candidate)
    old=consensus(baseline);new=consensus(candidate)
    b_cov,_,_=depth_arrays(baseline,refs)
    n_cov,_,_=depth_arrays(candidate,refs)
    fraction={}
    for row in rows(candidate/'multimap_depth.csv'):
        key=(row['refname'],int(row['position'])-1)
        fraction[key]=float(row['avg_depth'])
    regions=[]
    variants=[]
    for ref in sorted(old.keys()|new.keys()):
        before=as_calls(old.get(ref,(0,'')))
        after=as_calls(new.get(ref,(0,'')))
        added=sorted(after.keys()-before.keys())
        lost=sorted(before.keys()-after.keys())
        changed=sorted(p for p in (before.keys() & after.keys())
                       if before[p]!=after[p])
        fraction_low=sum(
            p>=0 and ref in n_cov and p<len(n_cov[ref])
            and n_cov[ref][p]>=100
            and fraction.get((ref,p),0)<100 for p in added)
        # Distinguish callable MAX recovery backed by fractional >=100.
        fraction_ge=sum(
            p>=0 and ref in n_cov and p<len(n_cov[ref])
            and fraction.get((ref,p),0)>=100 for p in added)
        regs={'refname': ref,
                  'before_callable': len(before),
                  'after_callable': len(after),
                  'gained': len(added),'lost': len(lost),'discordant': len(changed),
                  'gained_max100_avg_below100': fraction_low,
                  'gained_avg_at_least100': fraction_ge}
        regions.append(regs)
        for p in changed:
            variants.append({'ref': ref,'pos': p+1,
                                 'before': before[p],'after': after[p],
                                 'old_depth': (int(b_cov[ref][p])
                                            if p<len(b_cov.get(ref,())) else None),
                                 'new_max': (int(n_cov[ref][p])
                                          if p<len(n_cov.get(ref,())) else None),
                                 'new_avg': fraction.get((ref,p))})
    return {'alias': alias,'complete': True,'regions': regions,
                'discordant_variant_positions': variants,
                'gained_max100_avg_below100': sum(r['gained_max100_avg_below100']
                                               for r in regions),
                'gained_avg_at_least100': sum(r['gained_avg_at_least100']
                                          for r in regions)}


if __name__=='__main__':
    data=[run(case) for case in read_manifest()]
    (WORK/'per_contig_private.json').write_text(json.dumps(data,indent=2)+'\n')
    for item in data:
        if item['complete']:
            print(item['alias'],
                  'gained_low_fraction',item['gained_max100_avg_below100'],
                  'gained_avg>=100',item['gained_avg_at_least100'],
                  'discordant',len(item['discordant_variant_positions']))
            if item['discordant_variant_positions']:
                print('  diagnostics',item['discordant_variant_positions'][:4])
