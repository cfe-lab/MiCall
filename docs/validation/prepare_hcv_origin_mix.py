"""Real HCV read-origin dilution and absent-genotype decoy controls.

Original paired FASTQ molecules are never synthesized or edited (apart
from identifying QNAME prefixes). The source of each added pair is known,
not the independent laboratory genotype of the source specimen.
"""
import argparse
from collections import Counter
import csv
import gzip
import json
from pathlib import Path
import random


SOURCE_ROOT = Path('.')
SOURCES = {}


def configure_private_sources(root, manifest):
    global SOURCE_ROOT, SOURCES
    SOURCE_ROOT = root
    SOURCES = json.loads(manifest.read_text())
    if set(SOURCES) != {'HCV1a', 'HCV2b', 'HCV1b_DECOY'}:
        raise ValueError('Expected three privately mapped source roles')
    for sample in SOURCES.values():
        if not (SOURCE_ROOT / sample).is_dir():
            raise FileNotFoundError('A private source dataset is absent')




def read_pairs(sample):
    folder=SOURCE_ROOT / SOURCES[sample]
    r1,=folder.glob('*_L001_R1_001.fastq.gz')
    r2,=folder.glob('*_L001_R2_001.fastq.gz')
    with gzip.open(r1,'rt') as a, gzip.open(r2,'rt') as b:
        while True:
            one=[a.readline() for _ in range(4)]
            two=[b.readline() for _ in range(4)]
            if not any(one) and not any(two):break
            if not all(one) or not all(two):
                raise ValueError('Unpaired FASTQ')
            if one[0].split()[0]!=two[0].split()[0]:
                raise ValueError('Mate names differ')
            yield one,two


def pick(sample,count,seed):
    # Reservoir sampling is uniform over the original real read pairs.
    rng=random.Random(seed)
    selected=[]
    for i,pair in enumerate(read_pairs(sample)):
        if i<count:
            selected.append(pair)
        else:
            j=rng.randint(0,i)
            if j<count:selected[j]=pair
    if len(selected)!=count:
        raise ValueError(f'Only {len(selected)} real pairs available for {sample}')
    return selected


def prefix_pair(tag,pair):
    renamed=[]
    for record in pair:
        head=record[0].strip()
        if not head.startswith('@'):
            raise ValueError('Malformed FASTQ name')
        new='@'+tag+'|'+head[1:]+'\n'
        renamed.append([new,*record[1:]])
    return renamed


def input_contigs(folder):
    out=folder/'fixed-mapping-contigs.csv'
    with out.open('w',newline='') as stream:
        writer=csv.DictWriter(stream,['ref','match','group_ref','contig'],
                              lineterminator='\n')
        writer.writeheader()
        for source in SOURCES.values():
            sample=SOURCE_ROOT/source
            source_file=next(sample.glob('unstitched_contigs_*.csv'))
            candidates=[r for r in csv.DictReader(source_file.open())
                        if r['ref']!='unknown' and float(r['match'])>=0.25]
            if len(candidates)!=1:
                raise ValueError(f'{source} has {len(candidates)} usable contigs')
            row=candidates[0]
            writer.writerow({key:row[key] for key
                             in ('ref','match','group_ref','contig')})
    return out


def write_mix(root,name,minority):
    folder=root/name
    folder.mkdir(parents=True,exist_ok=True)
    r1=folder/(name+'_L001_R1_001.fastq.gz')
    r2=folder/(name+'_L001_R2_001.fastq.gz')
    if r1.exists() or r2.exists():
        raise FileExistsError('Do not overwrite existing real experiment reads')
    main=pick('HCV1a',50000,2519)
    rare=pick('HCV2b',minority,3872) if minority else []
    # Interleave every minority pair among known HCV1a reads so the assay
    # does not rely on an artifact of FASTQ record order.
    placed=Counter()
    with gzip.open(r1,'wt',compresslevel=3) as out1, \
         gzip.open(r2,'wt',compresslevel=3) as out2:
        for i,body in enumerate(main):
            if i < len(rare):
                for part,output in zip(prefix_pair('HCV2b',rare[i]),
                                       (out1,out2)):
                    output.writelines(part)
                placed['HCV2b']+=1
            for part,output in zip(prefix_pair('HCV1a',body),
                                   (out1,out2)):
                output.writelines(part)
            placed['HCV1a']+=1
    print(name,dict(placed),'minority_rate',
          round(minority/(50000+minority),6),flush=True)
    return folder


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--name',choices=('HCV-origin-zero','HCV-origin-1pct'),
                    required=True)
    ap.add_argument('--source-root', type=Path, required=True)
    ap.add_argument('--private-sources', type=Path, required=True)
    args=ap.parse_args()
    configure_private_sources(args.source_root,args.private_sources)
    root=Path('.validation-work/independent-truth/hcv-origin')
    root.mkdir(parents=True,exist_ok=True)
    input_contigs(root)
    write_mix(root/'input-samples',args.name,
              0 if args.name.endswith('zero') else 500)


if __name__=='__main__':
    main()
