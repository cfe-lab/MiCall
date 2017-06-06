import os
from argparse import ArgumentParser, FileType
import csv
from datetime import datetime
from gzip import GzipFile
from itertools import groupby
from operator import itemgetter

from micall.utils.translation import reverse_and_complement


def parse_args():
    parser = ArgumentParser(description='Load FASTQ files from 454 data.')
    parser.add_argument('source',
                        type=FileType('r'),
                        help='CSV file with 454 data.')
    parser.add_argument('dest',
                        help='Folder to create FASTQ files in.')
    return parser.parse_args()


def main():
    args = parse_args()
    with args.source as source:
        reader = csv.DictReader(source)
        for (run, sample), rows in groupby(reader, itemgetter('run', 'enum')):
            sample_name = format_sample_name(run, sample)
            filename1 = os.path.join(args.dest, sample_name + '_R1_001.fastq.gz')
            filename2 = os.path.join(args.dest, sample_name + '_R2_001.fastq.gz')
            print(filename1)
            with open(filename1, 'wb') as dest1, open(filename2, 'wb') as dest2:
                dest1_zip = GzipFile(fileobj=dest1)
                dest2_zip = GzipFile(fileobj=dest2)
                for i, row in enumerate(rows):
                    seq = row['string'].replace('-', '')
                    for j in range(3):
                        # Three duplicates so that G2P doesn't ignore it.
                        prefix = '@M454:01:000000000-AAAAA:1:1101:{}:{}'.format(
                            10*i + j,
                            row['count'])
                        dest1_zip.write(prefix + ' 1:N:0:1\n')
                        dest2_zip.write(prefix + ' 2:N:0:1\n')
                        dest1_zip.write(seq + '\n')
                        dest2_zip.write(reverse_and_complement(seq) + '\n')
                        dest1_zip.write('+\n')
                        dest2_zip.write('+\n')
                        quality = 'A' * len(seq)
                        dest1_zip.write(quality + '\n')
                        dest2_zip.write(quality + '\n')
                dest1_zip.close()
                dest2_zip.close()
    print('Done.')


def format_sample_name(run, sample):
    run_fields = run.split('.')
    run_date = datetime.strptime(run_fields[0], '%d-%b-%Y')
    if len(run_fields) == 1:
        run_suffix = ''
    else:
        run_suffix = '.' + '.'.join(run_fields[1:])
    sample_name = '{}{}-{}'.format(run_date.strftime('%Y-%m-%d'),
                                   run_suffix,
                                   sample)
    return sample_name


if __name__ == '__main__':
    main()
