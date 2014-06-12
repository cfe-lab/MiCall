#! /usr/bin/python

"""
Shipyard-style MiSeq pipeline, step 2
Takes preliminary SAM as CSV input.  Iterative re-mapping of reads from
original FASTQ files.
Also report the number of reads mapped before and after processing.
Dependencies:
    bowtie2-build
    bowtie2-align
    samtools
    settings.py
"""

import argparse
import subprocess
import os
import itertools
import sys
import re

indel_re = re.compile('[+-][0-9]+')

from settings import *  # settings.py is a CodeResourceDependency


def pileup_to_conseq (handle, qCutoff):
    conseq = ''
    to_skip = 0
    last_pos = 0
    for line in handle:
        if to_skip > 0:
            to_skip -= 1
            continue

        label, pos, en, depth, astr, qstr = line.strip('\n').split('\t')
        pos = int(pos)  # position in the pileup, 1-index
        if (pos - last_pos) > 1:
            conseq += 'N' * (pos - last_pos - 1)
        last_pos = pos
        alist = []  # alist stores all bases at a given coordinate
        i = 0       # Current index for astr
        j = 0       # Current index for qstr

        while i < len(astr):
            if astr[i] == '^':
                q = ord(qstr[j])-33
                base = astr[i+2] if q >= qCutoff else 'N'
                alist.append(base.upper())
                i += 3
                j += 1
            elif astr[i] in '*':
                alist.append('-')
                i += 1
            elif astr[i] == '$':
                i += 1
            elif i < len(astr)-1 and astr[i+1] in '+-':
                m = indel_re.match(astr[i+1:])
                indel_len = int(m.group().strip('+-'))
                left = i+1 + len(m.group())
                insertion = astr[left:(left+indel_len)]
                q = ord(qstr[j])-33
                base = astr[i].upper() if q >= qCutoff else 'N'
                token = base + m.group() + insertion.upper()
                if astr[i+1] == '+':
                    alist.append(token)
                else:
                    alist.append(base)
                i += len(token)
                j += 1
            else:
                # Operative case: sequence matches reference (And no indel ahead)
                q = ord(qstr[j])-33
                base = astr[i].upper() if q >= qCutoff else 'N'
                alist.append(base)
                i += 1
                j += 1

        atypes = set(alist)
        intermed = []
        for atype in atypes:
            intermed.append((alist.count(atype), atype))
        intermed.sort(reverse=True)
        token = intermed[0][1]
        if '+' in token:
            m = indel_re.findall(token)[0] # \+[0-9]+
            conseq += token[0] + token[1+len(m):]
        elif token == '-':
            pass
        else:
            conseq += token
    handle.close()
    return conseq


parser = argparse.ArgumentParser('Iterative remapping of bowtie2 by reference.')

parser.add_argument('fastq1', help='<input> FASTQ containing forward reads')
parser.add_argument('fastq2', help='<input> FASTQ containing reverse reads')
parser.add_argument('sam_csv', help='<input> SAM output of bowtie2 in CSV format')
parser.add_argument('ref', help='<input> initial set of references in FASTA format')
parser.add_argument('output_csv', help='<output> CSV containing remap output (modified SAM)')
parser.add_argument('stats_csv', help='<output> CSV containing numbers of mapped reads')
parser.add_argument('conseq_csv', help='<output> CSV containing mapping consensus sequences')

args = parser.parse_args()

# check that the inputs exist
if not os.path.exists(args.fastq1):
    print 'No FASTQ found at', args.fastq1
    sys.exit(1)

if not os.path.exists(args.fastq2):
    print 'No FASTQ found at', args.fastq2
    sys.exit(1)

# check that we have access to bowtie2
try:
    p = subprocess.Popen(['bowtie2', '-h'], stdout=subprocess.PIPE)
except OSError:
    print 'bowtie2 not found; check if it is installed and in $PATH\n'
    raise

# check that the output paths are valid
for path in [args.output_csv, args.stats_csv, args.conseq_csv]:
    output_path = os.path.split(path)[0]
    if not os.path.exists(output_path) and output_path != '':
        print 'Output path does not exist:', output_path
        sys.exit(1)


# get the raw read count
p = subprocess.Popen(['wc', '-l', args.fastq1], stdout=subprocess.PIPE)
raw_count = int(p.stdout.readline().split()[0]) / 2  # 4 lines per record in FASTQ, paired

stat_file = open(args.stats_csv, 'w')
stat_file.write('raw,%d\n' % raw_count)


handle = open(args.sam_csv, 'rU')

# group CSV stream by first item
prelim_count = 0
map_counts = {}
refnames = []
for refname, group in itertools.groupby(handle, lambda x: x.split(',')[2]):
    refnames.append(refname)
    # reconstitute region-specific SAM files
    tmpfile = open('%s.sam' % refname, 'w')
    count = 0
    for line in group:
        tmpfile.write('\t'.join(line.split(',')))
        prelim_count += 1
        count += 1
    stat_file.write('prelim %s,%d\n' % (refname, count))
    map_counts.update({refname: count})
    tmpfile.close()

# settings for iterative remapping
ref = args.ref
n_remaps = 0
mapping_efficiency = float(prelim_count) / raw_count
fnull = open(os.devnull, 'w')
frozen = []  # which regions to stop re-mapping
tmpfile = 'temp.sam'  # temporary bowtie2-align output

conseqs = {}

while n_remaps < max_remaps:
    if len(frozen) == len(refnames):
        # every region is frozen
        break

    for refname in refnames:
        if refname in frozen:
            # don't attempt to re-map this region
            continue

        samfile = refname+'.sam'
        bamfile = refname+'.bam'

        # convert SAM to BAM
        with open(bamfile, 'w') as f:
            p = subprocess.Popen(['samtools', 'view', '-b', '-T', ref, samfile], stdout=f, stderr=fnull)
            p.wait()
        p = subprocess.Popen(['samtools', 'sort', bamfile, refname])  # overwrite
        p.wait()

        # BAM to pileup
        with open(bamfile+'.pileup', 'w') as f:
            p = subprocess.Popen(['samtools', 'mpileup', '-A', bamfile], stdout=f, stderr=fnull)
            p.wait()

        # pileup to consensus sequence
        with open(bamfile+'.pileup', 'rU') as f:
            conseqs[refname] = pileup_to_conseq(f, consensus_q_cutoff)

        # consensus to *.bt2
        p = subprocess.Popen(['bowtie2-build', '-c', '-q', conseqs[refname], refname])
        p.wait()

        # map raw data to new reference - overwrite SAM file
        p = subprocess.Popen(['bowtie2-align', '--quiet', '-p', str(bowtie_threads), '--local',
                              '-x', refname, '-1', args.fastq1, '-2', args.fastq2,
                              '--no-unal', '-S', tmpfile], stdout=fnull)
        p.wait()

        # how many reads did we map?
        p = subprocess.Popen(['wc', '-l', tmpfile], stdout=subprocess.PIPE)
        count = int(p.stdout.readline().split()[0])

        if count <= map_counts[refname]:
            # failed to improve the number of mapped reads
            frozen.append(refname)
            continue

        # overwrite previous SAM file
        os.rename(tmpfile, samfile)
        map_counts[refname] = count

    n_remaps += 1
    mapping_efficiency = sum(map_counts.values()) / float(raw_count)
    if mapping_efficiency > min_mapping_efficiency:
        break  # a sufficient fraction of raw data has been mapped

fnull.close()


seqfile = open(args.conseq_csv, 'w')  # record consensus sequences for later use
outfile = open(args.output_csv, 'w')  # combine SAM files into single CSV output

for refname in refnames:
    stat_file.write('remap %s,%d\n' % (refname, map_counts[refname]))
    seqfile.write('%s,%s\n' % (refname, conseqs[refname]))
    handle = open(refname+'.sam', 'rU')
    for line in handle:
        if line.startswith('@'):
            continue  # omit SAM header lines
        items = line.strip('\n').split('\t')[:11]
        items[2] = refname  # replace '0' due to passing conseq to bowtie2-build on cmd line
        outfile.write(','.join(items) + '\n')
    handle.close()

outfile.close()
stat_file.close()
seqfile.close()
