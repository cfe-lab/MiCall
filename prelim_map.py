#! /usr/bin/python

"""
Shipyard-style bowtie2
Run bowtie2 on paired-end FASTQ data sets with user-supplied *.bt2
bowtie2 SAM format output to <stdout> for redirection via subprocess.Popen
Sort outputs by refname.
Convert to CSV format and write to file.

Dependencies: settings.py (derived from settings_default.py)
              *.bt2 files produced by bowtie2-align - pass as metapackage
"""

import argparse
import os
import shutil
import subprocess
import sys
import tempfile

import settings  # settings.py is a CodeResourceDependency

parser = argparse.ArgumentParser('Map contents of FASTQ R1 and R2 data sets to references using bowtie2.')

parser.add_argument('fastq1', help='<input> FASTQ containing forward reads')
parser.add_argument('fastq2', help='<input> FASTQ containing reverse reads')
parser.add_argument('ref', help='<input> initial set of references in FASTA format')
parser.add_argument('sam_csv', help='<output> CSV containing bowtie2 output (modified SAM)')

# TODO: pass number of threads and --local to bowtie2 as a CodeResourceDependency
args = parser.parse_args()


# check that we have access to bowtie2
try:
    p = subprocess.Popen(['bowtie2', '-h'], stdout=subprocess.PIPE)
except OSError:
    print 'bowtie2 not found; check if it is installed and in $PATH\n'
    raise

# check that the inputs exist
if not os.path.exists(args.fastq1):
    print 'No FASTQ found at', args.fastq1
    sys.exit(1)

if not os.path.exists(args.fastq2):
    print 'No FASTQ found at', args.fastq2
    sys.exit(1)

if not os.path.exists(args.ref):
    print 'No reference sequences found at', args.ref
    sys.exit(1)

# check that the SAM output path is valid
output_path = os.path.split(args.sam_csv)[0]
if not os.path.exists(output_path) and output_path != '':
    print 'SAM output path does not exist:', output_path
    sys.exit(1)

# Create a temp folder to hold the reference sequences.
referencesFolder = tempfile.mkdtemp(prefix='tmpReferences')

try:
    referenceFilenameTemplate = os.path.join(referencesFolder, 'reference')
    p = subprocess.check_call(['bowtie2-build',
                               '--quiet', 
                               '-f', 
                               args.ref,
                               referenceFilenameTemplate])
    
    # do preliminary mapping
    output = {}
    
    # stream output from bowtie2
    bowtie_args = ['bowtie2', 
                   '--quiet', 
                   '-x', referenceFilenameTemplate, 
                   '-1', args.fastq1, 
                   '-2', args.fastq2, 
                   '--no-unal', 
                   '--local', 
                   '-p', str(settings.bowtie_threads)]
    p = subprocess.Popen(bowtie_args, stdout=subprocess.PIPE)
    with p.stdout:
        for line in p.stdout:
            if line.startswith('@'):
                # skip header line
                continue
            refname = line.split('\t')[2]  # read was mapped to this reference
            if not refname in output:
                output.update({refname: []})
            output[refname].append('\t'.join(line.split('\t')[:11]))  # discard optional items
            
    if p.returncode:
        raise subprocess.CalledProcessError(p.returncode, bowtie_args)
finally:
    shutil.rmtree(referencesFolder)

# lines grouped by refname
with open(args.sam_csv, 'w') as outfile:
    for refname, lines in output.iteritems():
        for line in lines:
            outfile.write(line.replace('\t', ',') + '\n')
