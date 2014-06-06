#! /usr/bin/python

"""
Shipyard-style MiSeq pipeline, step 2
Takes preliminary SAM as CSV input.
Requires: bowtie2-build, bowtie2-align, samtools
"""

import argparse
import subprocess
import os
import itertools
import sys

from settings import *  # settings.py is a CodeResourceDependency

parser = argparse.ArgumentParser('Iterative remapping of bowtie2 by reference.')

parser.add_argument('sam_csv', help='<input> SAM output of bowtie2 in CSV format')
parser.add_argument('init_ref_fasta', help='<input> initial set of references in FASTA or *.fai format')
parser.add_argument('output_csv', help='<output> CSV containing remap output (modified SAM)')

args = parser.parse_args()


# check that we have access to bowtie2
try:
    p = subprocess.Popen(['bowtie2', '-h'], stdout=subprocess.PIPE)
except OSError:
    print 'bowtie2 not found; check if it is installed and in $PATH\n'
    raise

# check that the output path is valid
output_path, output_filename = os.path.split(args.output_csv)
if not os.path.exists(output_path) and output_path != '':
    print 'Output path does not exist:', output_path
    sys.exit(1)


handle = open(args.sam_csv, 'rU')

# group CSV stream by first item
for refname, group in itertools.groupby(handle, lambda x: x.split(',')[1]):
    # reconstitute temporary SAM file
    tmpfile = open('temp.sam', 'w')
    for line in group:
        tmpfile.write('\t'+'\t'.join(line.split(',')))  # use an empty qname
    tmpfile.close()

    attempts = 0
    while attempts < max_remaps or prop_mapped < min_mapping_efficiency:

        attempts += 1
    break
