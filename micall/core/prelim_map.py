"""
Kive-style bowtie2
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

from micall.core import project_config
from micall.utils.externals import Bowtie2, Bowtie2Build, LineCounter

BOWTIE_THREADS = 1    # Bowtie performance roughly scales with number of threads
# Read and reference gap open/extension penalties.
READ_GAP_OPEN = 10
READ_GAP_EXTEND = 3
REF_GAP_OPEN = 10
REF_GAP_EXTEND = 3

logger = logging.getLogger(__name__)
line_counter = LineCounter()


def prelim_map(fastq1,
               fastq2,
               prelim_csv,
               nthreads=BOWTIE_THREADS,
               rdgopen=READ_GAP_OPEN,
               rfgopen=REF_GAP_OPEN,
               stderr=sys.stderr,
               gzip=False,
               work_path='',
               excluded_seeds=None):
    """ Run the preliminary mapping step.

    @param fastq1: the file name for the forward reads in FASTQ format
    @param fastq2: the file name for the reverse reads in FASTQ format
    @param prelim_csv: an open file object for the output file - all the reads
        mapped to references in CSV version of the SAM format
    @param nthreads: the number of threads to use.
    @param rdgopen: a penalty for opening a gap in the read sequence.
    @param rfgopen: a penalty for opening a gap in the reference sequence.
    @param stderr: where to write the standard error output from bowtie2 calls.
    @param gzip: True if FASTQ files are in gzip format
    @param work_path:  optional path to store working files
    @param excluded_seeds: a list of seed names to exclude from mapping
    """

    bowtie2 = Bowtie2()
    bowtie2_build = Bowtie2Build()
    bowtie2_build.set_logger(logger)

    # check that the inputs exist
    fastq1 = check_fastq(fastq1, gzip)
    fastq2 = check_fastq(fastq2, gzip)

    # generate initial reference files
    projects = project_config.ProjectConfig.loadDefault()
    ref_path = os.path.join(work_path, 'micall.fasta')
    all_excluded_seeds = {project_config.G2P_SEED_NAME}
    if excluded_seeds:
        all_excluded_seeds.update(excluded_seeds)
    with open(ref_path, 'w') as ref:
        projects.writeSeedFasta(ref, all_excluded_seeds)
    reffile_template = os.path.join(work_path, 'reference')
    bowtie2_build.build(ref_path, reffile_template)

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
    writer = csv.writer(prelim_csv, lineterminator=os.linesep)
    writer.writerow(fieldnames)

    # do preliminary mapping
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
                   '--no-hd',  # no header lines (start with @)
                   '-X', '1200',
                   '-p', str(nthreads)]

    for i, line in enumerate(bowtie2.yield_output(bowtie_args, stderr=stderr)):
        writer.writerow(line.split('\t')[:11])  # discard optional items


def check_fastq(filename, gzip=False):
    if not os.path.exists(filename):
        sys.exit('No FASTQ found at ' + filename)
    if gzip:
        if not filename.endswith('.gz'):
            new_filename = filename + '.gz'
            try:
                os.symlink(filename, new_filename)
            except FileExistsError:
                pass
            filename = new_filename
    return filename


def main():
    parser = argparse.ArgumentParser(
        description='Map contents of FASTQ R1 and R2 data sets to references using bowtie2.')

    parser.add_argument('fastq1', help='<input> FASTQ containing forward reads')
    parser.add_argument('fastq2', help='<input> FASTQ containing reverse reads')
    parser.add_argument('prelim_csv',
                        type=argparse.FileType('w'),
                        help='<output> CSV containing preliminary mapping from bowtie2 (modified SAM)')
    parser.add_argument("--rdgopen", default=READ_GAP_OPEN, help="<optional> read gap open penalty")
    parser.add_argument("--rfgopen", default=REF_GAP_OPEN, help="<optional> reference gap open penalty")
    parser.add_argument("--gzip", action='store_true', help="<optional> FASTQs are compressed")

    args = parser.parse_args()
    work_path = os.path.dirname(args.prelim_csv.name)
    prelim_map(fastq1=args.fastq1,
               fastq2=args.fastq2,
               prelim_csv=args.prelim_csv,
               rdgopen=args.rdgopen,
               rfgopen=args.rfgopen,
               gzip=args.gzip,  # defaults to False
               work_path=work_path)


if __name__ == '__main__':
    main()
