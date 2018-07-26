
import argparse
import os
import sys
import csv

from micall.core.parse_interop import read_errors, write_phix_csv
from micall.core.filter_quality import report_bad_cycles
from micall.core.censor_fastq import censor
from micall.core.prelim_map import prelim_map

def parseArgs():
    parser = argparse.ArgumentParser(
        description='Use the MiCall pipeline to process gzip-comprssed FASTQ '
                    'file(s) for one sample.'
    )
    parser.add_argument('--fastq1', '-1', type=argparse.FileType('rb'), required=True,
                        help='Path to *.R1.fastq.gz file of paired set, or to an '
                             'unpaired *.fastq.gz file.')
    parser.add_argument('--fastq2', '-2', type=argparse.FileType('rb'),
                        required=False,
                        help='Path to *.R2.fastq.gz file of paired set.  Unused if'
                             ' processing an unpaired sample.')
    parser.add_argument('--outdir', '-d', default=None, required=False,
                        help='<optional> Path to write output files.')
    parser.add_argument('--unzipped', '-u', action='store_true', required=False,
                        help='Set if the FASTQ file is not compressed.')
    parser.add_argument('--interop', '-i', type=argparse.FileType('rb'),
                        required=False,
                        help='<optional> Path to ErrorMetricsOut.bin interop file.')
    parser.add_argument('--readlen', '-l', type=int, default=251,
                        help='<optional> Read length (default: 251nt).')
    parser.add_argument('--index', '-x', type=int, default=8,
                        help='<optional> Index length (default: 8nt).')

    return parser.parse_args()


def get_prefix(args):
    """
    Get filename prefix for sample
    :param args:
    :return:
    """
    fn = os.path.basename(args.fastq1.name)
    prefix = None
    if args.unzipped:
        if fn.endswith('.fastq'):
            prefix = fn.replace('.fastq', '')
    else:
        if fn.endswith('.fastq.gz'):
            prefix = fn.replace('.fastq.gz', '')

    if prefix is None:
        print('Error: {} has non-standard file extension'.format(args.fastq1))
        sys.exit()

    return prefix


def censor_fastqs(args, prefix):
    """
    The Illumina system generates a set of binary-encoded (InterOp) files
    that contain useful information about the run.  One of these files, called
    ErrorMetricsOut.bin, contains information about the tile-cycle error rates
    as estimated from the phiX174 control sample.  We use this report to
    identify tile-cycle combinations with excessive error rates.
    :param args:  return value from argparse.ArgumentParser()
    :param prefix:  filename stem
    :return:  a new <args> object with .fastq1 and (optionally) .fastq2
              replaced by read-only file objects to censored FASTQs
    """
    # parse ErrorMetricsOut.bin
    lengths = [args.readlen, args.index, args.index, args.readlen]
    records = read_errors(args.interop)
    quality_csv = os.path.join(args.outdir, prefix + '.quality.csv')
    with open(quality_csv, 'w') as handle:
        write_phix_csv(out_file=handle, records=records, read_lengths=lengths)

    # find bad tile-cycle combinations
    bad_cycles_csv = os.path.join(args.outdir, prefix + '.bad_cycles.csv')
    with open(quality_csv, 'r') as f1, open(bad_cycles_csv, 'w') as f2:
        report_bad_cycles(f1, f2)

    bad_cycles = csv.DictReader(open(bad_cycles_csv, 'r'))
    cfastq1 = os.path.relpath(args.fastq1.name.replace('.fastq', '.censor.fastq'))
    censor(src=args.fastq1,
           bad_cycles_reader=bad_cycles,
           dest=open(cfastq1, 'wb'),
           use_gzip=not args.unzipped)
    args.fastq1 = open(cfastq1, 'rb')  # replace original file

    if args.fastq2:
        cfastq2 = os.path.relpath(args.fastq2.name.replace('.fastq', '.censor.fastq'))
        censor(args.fastq2, bad_cycles, open(cfastq2, 'wb'), not args.unzipped)
        args.fastq2 = open(cfastq2, 'rb')

    return args




if __name__ == '__main__':
    args = parseArgs()
    if args.outdir is None:
        args.outdir = os.path.dirname(args.fastq1.name)
    prefix = get_prefix(args)

    if args.interop:
        args = censor_fastqs(args, prefix)

    prelim_csv = os.path.join(args.outdir, prefix+'.prelim.csv')
    with open(prelim_csv, 'w') as handle:
        prelim_map(args.fastq1, args.fastq2, handle, gzip=not args.unzipped)
