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
import logging
import os
import subprocess
import sys

import miseq_logging
import settings  # settings.py is a CodeResourceDependency

logger = miseq_logging.init_logging_console_only(logging.DEBUG)

def main():
    parser = argparse.ArgumentParser('Map contents of FASTQ R1 and R2 data sets to references using bowtie2.')
    
    parser.add_argument('fastq1', help='<input> FASTQ containing forward reads')
    parser.add_argument('fastq2', help='<input> FASTQ containing reverse reads')
    parser.add_argument('sam_csv', help='<output> CSV containing bowtie2 output (modified SAM)')
    
    args = parser.parse_args()


    # check that we have access to bowtie2
    try:
        subprocess.check_output(['bowtie2', '-h'])
    except OSError:
        raise RuntimeError('bowtie2 not found; check if it is installed and in $PATH\n')

    # check that the inputs exist
    if not os.path.exists(args.fastq1):
        logger.error('No FASTQ found at %s', args.fastq1)
        sys.exit(1)

    if not os.path.exists(args.fastq2):
        logger.error('No FASTQ found at %s', args.fastq2)
        sys.exit(1)

    # check that the SAM output path is valid
    output_path = os.path.split(args.sam_csv)[0]
    if not os.path.exists(output_path) and output_path != '':
        logger.error('SAM output path does not exist: %s', output_path)
        sys.exit(1)

    # generate initial reference files
    is_ref_found = False
    possible_refs = ('cfe.fasta', '../reference_sequences/cfe.fasta')
    for ref in possible_refs:
        if not os.path.isfile(ref):
            continue
        is_ref_found = True
        log_call(['samtools', 'faidx', ref])
        break
    if not is_ref_found:
        raise RuntimeError('No reference sequences found in {!r}'.format(
            possible_refs))
    reffile_template = 'reference'
    log_call(['bowtie2-build',
              '--quiet',
              '-f',
              ref,
              reffile_template])

    # do preliminary mapping
    output = {}

    # stream output from bowtie2
    bowtie_args = ['bowtie2',
                   '--quiet',
                   '-x', reffile_template,
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

    # lines grouped by refname
    with open(args.sam_csv, 'w') as outfile:
        for refname, lines in output.iteritems():
            for line in lines:
                outfile.write(line.replace('\t', ',') + '\n')


def log_call(args, format_string='%s'):
    """ Launch a subprocess, and log any output to the debug logger.
    
    Raise an exception if the return code is not zero. This assumes only a
    small amount of output, and holds it all in memory before logging it.
    @param args: A list of arguments to pass to subprocess.Popen().
    @param format_string: A template for the debug message that will have each
    line of output formatted with it.
    """
    output = subprocess.check_output(args, stderr=subprocess.STDOUT)
    for line in output.splitlines():
        logger.debug(format_string, line)

if __name__ == '__main__':
    main()
