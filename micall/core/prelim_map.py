#! /usr/bin/env python

"""
Shipyard-style bowtie2
Run bowtie2 on paired-end FASTQ data sets with user-supplied *.bt2
bowtie2 SAM format output to <stdout> for redirection via subprocess.Popen
Sort outputs by refname.
Convert to CSV format and write to file.
"""

import argparse
import csv
import logging
import os
import sys

import miseq_logging
import project_config
from micall.utils.externals import Bowtie2, Bowtie2Build, LineCounter

BOWTIE_THREADS = 4    # Bowtie performance roughly scales with number of threads
BOWTIE_VERSION = '2.2.8'        # version of bowtie2, used for version control
BOWTIE_PATH = 'bowtie2-align-s'         # path to executable, so you can install more than one version
BOWTIE_BUILD_PATH = 'bowtie2-build-s'
# Read and reference gap open/extension penalties.
READ_GAP_OPEN = 10
READ_GAP_EXTEND = 3
REF_GAP_OPEN = 10
REF_GAP_EXTEND = 3

logger = miseq_logging.init_logging_console_only(logging.DEBUG)
line_counter = LineCounter()


def prelim_map(fastq1,
               fastq2,
               prelim_csv,
               nthreads=BOWTIE_THREADS,
               callback=None,
               rdgopen=READ_GAP_OPEN,
               rfgopen=REF_GAP_OPEN,
               stderr=sys.stderr,
               gzip=False):
    """ Run the preliminary mapping step.

    @param fastq1: the file name for the forward reads in FASTQ format
    @param fastq2: the file name for the reverse reads in FASTQ format
    @param prelim_csv: an open file object for the output file - all the reads
        mapped to references in CSV version of the SAM format
    @param nthreads: the number of threads to use.
    @param callback: a function to report progress with three optional
        parameters - callback(message, progress, max_progress)
    @param rdgopen: a penalty for opening a gap in the read sequence.
    @param rfgopen: a penalty for opening a gap in the reference sequence.
    @param stderr: where to write the standard error output from bowtie2 calls.
    """
    try:
        bowtie2 = Bowtie2(BOWTIE_VERSION, BOWTIE_PATH)
        bowtie2_build = Bowtie2Build(BOWTIE_VERSION,
                                     BOWTIE_BUILD_PATH,
                                     logger)
    except:
        bowtie2 = Bowtie2(BOWTIE_VERSION, BOWTIE_PATH + '-' + BOWTIE_VERSION)
        bowtie2_build = Bowtie2Build(BOWTIE_VERSION,
                                     BOWTIE_BUILD_PATH + '-' + BOWTIE_VERSION,
                                     logger)

    # check that the inputs exist
    if not os.path.exists(fastq1):
        logger.error('No FASTQ found at %s', fastq1)
        sys.exit(1)

    if not os.path.exists(fastq2):
        logger.error('No FASTQ found at %s', fastq2)
        sys.exit(1)

    # append .gz extension if necessary
    if gzip:
        if not fastq1.endswith('.gz'):
            try:
                os.symlink(fastq1, fastq1+'.gz')
            except OSError:
                # symbolic link already exists
                pass
            fastq1 += '.gz'

        if not fastq2.endswith('.gz'):
            try:
                os.symlink(fastq2, fastq2+'.gz')
            except OSError:
                # symbolic link already exists
                pass
            fastq2 += '.gz'

    if callback:
        # four lines per read, two files
        total_reads = line_counter.count(fastq1, gzip=gzip) / 2
        callback(message='... preliminary mapping',
                 progress=0,
                 max_progress=total_reads)

    # generate initial reference files
    projects = project_config.ProjectConfig.loadDefault()
    ref_path = 'micall.fasta'
    with open(ref_path, 'w') as ref:
        projects.writeSeedFasta(ref)
    reffile_template = 'reference'
    bowtie2_build.build(ref_path, reffile_template)

    # do preliminary mapping
    output = {}
    read_gap_open_penalty = rdgopen
    ref_gap_open_penalty = rfgopen

    # stream output from bowtie2
    bowtie_args = ['--wrapper', 'micall-0',
                   '--quiet',
                   '-x', reffile_template,
                   '-1', fastq1,
                   '-2', fastq2,
                   '--rdg', "{},{}".format(read_gap_open_penalty,
                                           READ_GAP_EXTEND),
                   '--rfg', "{},{}".format(ref_gap_open_penalty,
                                           REF_GAP_EXTEND),
                   '--no-unal',  # don't report reads that failed to align
                   '--no-hd',  # no header lines (start with @)
                   '-X', '1200',
                   '-p', str(nthreads)]

    for i, line in enumerate(bowtie2.yield_output(bowtie_args, stderr=stderr)):
        if callback and i % 1000 == 0:
            callback(progress=i)
        refname = line.split('\t')[2]  # read was mapped to this reference
        if refname not in output:
            output.update({refname: []})
        output[refname].append(line.split('\t')[:11])  # discard optional items

    fieldnames = ['qname',
                  'flag',
                  'rname',
                  'pos',
                  'mapq',
                  'cigar',
                  'rnext',
                  'pnext',
                  'tlen',
                  'seq',
                  'qual']
    writer = csv.DictWriter(prelim_csv, fieldnames, lineterminator=os.linesep)
    writer.writeheader()

    # lines grouped by refname
    for refname, lines in output.iteritems():
        for line in lines:
            writer.writerow(dict(zip(fieldnames, line)))

    if callback:
        # Track progress for second half
        callback(progress=total_reads)


def main():
    parser = argparse.ArgumentParser(
        description='Map contents of FASTQ R1 and R2 data sets to references using bowtie2.')

    parser.add_argument('fastq1', help='<input> FASTQ containing forward reads')
    parser.add_argument('fastq2', help='<input> FASTQ containing reverse reads')
    parser.add_argument('prelim_csv',
                        type=argparse.FileType('w'),
                        help='<output> CSV containing preliminary mapping from bowtie2 (modified SAM)')
    parser.add_argument("--rdgopen", default=None, help="<optional> read gap open penalty")
    parser.add_argument("--rfgopen", default=None, help="<optional> reference gap open penalty")
    parser.add_argument("--gzip", action='store_true', help="<optional> FASTQs are compressed")

    args = parser.parse_args()
    prelim_map(fastq1=args.fastq1,
               fastq2=args.fastq2,
               prelim_csv=args.prelim_csv,
               rdgopen=args.rdgopen,
               rfgopen=args.rfgopen,
               gzip=args.gzip)  # defaults to False


if __name__ == '__main__':
    main()
