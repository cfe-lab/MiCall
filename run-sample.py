
import argparse
import os
import sys

from micall.core.parse_interop import read_errors, write_phix_csv
from micall.core.filter_quality import report_bad_cycles

def parseArgs():
    parser = argparse.ArgumentParser(
        description='Use the MiCall pipeline to process gzip-comprssed FASTQ '
                    'file(s) for one sample.'
    )
    parser.add_argument('--fastq1', '-1', type=argparse.FileType('r'),
                        help='Path to *.R1.fastq.gz file of paired set, or to an '
                             'unpaired *.fastq.gz file.')
    parser.add_argument('--fastq2', '-2', type=argparse.FileType('r'),
                        required=False,
                        help='Path to *.R2.fastq.gz file of paired set.  Unused if'
                             ' processing an unpaired sample.')
    parser.add_argument('--outdir', '-d', default=None, required=False,
                        help='Path to write output files.')
    parser.add_argument('--unzipped', '-u', action='store_true', required=False,
                        help='Set if the FASTQ file is not compressed.')
    parser.add_argument('--interop', '-i', type=argparse.FileType('rb'),
                        required=False,
                        help='<optional> Path to ErrorMetricsOut.bin interop file.')
    parser.add_argument('--readlen', '-l', type=int, default=250,
                        help='<optional> Read length (default: 250nt).')
    parser.add_argument('--index', '-x', type=int, default=8,
                        help='<optional> Index length (default: 8nt).')

    args = parser.parse_args()

    if args.outdir is None:
        args.outdir = os.path.dirname(args.fastq1)
    return args


def main():
    args = parseArgs()
    fn = os.path.basename(args.fastq1)

    # get filename prefix for sample
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

    if args.interop:
        records = read_errors(args.interop)
        quality_csv = os.join(args.outdir, prefix+'.quality.csv')
        with open(quality_csv, 'w') as handle:
            write_phix_csv(handle, records, [300, 8, 8, 300], {})

        report_bad_cycles(quality_csv)


if __name__ == '__main__':
    main()
