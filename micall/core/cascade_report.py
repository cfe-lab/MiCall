#!/usr/bin/env python3.6

from argparse import ArgumentParser, FileType
from csv import DictWriter, DictReader
from operator import itemgetter
import os


def parse_args():
    parser = ArgumentParser(
        description='Summarize the cascade of data through the pipeline.')

    parser.add_argument('g2p_summary_csv',
                        help='G2P call and read counts',
                        type=FileType('r'))
    parser.add_argument('remap_counts_csv',
                        help='how many reads mapped to each reference',
                        type=FileType('r'))
    parser.add_argument('aligned_csv',
                        help='aligned reads and their counts',
                        type=FileType('r'))
    parser.add_argument('cascade_csv',
                        help='count of reads at each step',
                        type=FileType('w'))
    return parser.parse_args()


class CascadeReport:
    def __init__(self, cascade_csv):
        """ Initialize a report object.

        :param cascade_csv: an open CSV file to write to
        """
        self.g2p_summary_csv = self.remap_counts_csv = self.aligned_csv = None
        self.unique_aligned_csv = None
        self.counts = None
        self.cascade_csv = cascade_csv

    def generate(self):
        field_names = ['demultiplexed',
                       'v3loop',
                       'g2p',
                       'prelim_map',
                       'remap',
                       'aligned']
        self.counts = dict.fromkeys(field_names, 0)

        self.read_g2p()
        self.read_remap()
        self.read_aligned()
        
        writer = DictWriter(self.cascade_csv,
                            lineterminator=os.linesep,
                            fieldnames=field_names)
        writer.writeheader()
        writer.writerow(self.counts)

    def read_g2p(self):
        if self.g2p_summary_csv is None:
            return

        reader = DictReader(self.g2p_summary_csv)
        row = next(reader)
        mapped = int(row['mapped'])
        self.counts['v3loop'] = mapped
        self.counts['g2p'] += int(row['valid'])

    def read_remap(self):
        if self.remap_counts_csv is None:
            return

        reader = DictReader(self.remap_counts_csv)
        fragment_count = None
        for row in reader:
            row_type = row['type']
            count = int(row['count']) // 2
            if row_type == 'raw':
                self.counts['demultiplexed'] += count
            elif row_type == 'prelim *':
                pass
            elif row_type.startswith('prelim '):
                self.counts['prelim_map'] += count
            elif row_type.startswith('remap-final '):
                # Reference-based remaps report read-end counts per seed.
                # De novo remaps also include an explicit fragment total.
                self.counts['remap'] += count
            elif row_type == 'remap-fragments':
                fragment_count = int(row['count'])
        if fragment_count is not None:
            # Prefer the exact fragment count regardless of row order.
            self.counts['remap'] = fragment_count

    def read_aligned(self):
        if self.unique_aligned_csv is not None:
            # Each accepted alternative placement appears separately in
            # aligned.csv. Do not count one original physical read pair
            # multiple times in MiCall's clinical cascade QC.
            reader = DictReader(self.unique_aligned_csv)
            record = next(reader, None)
            if record is None or next(reader, None) is not None:
                raise ValueError(
                    'Unique aligned fragment report must have one data row')
            count = int(record['unique_aligned_fragments'])
            if not 0 <= count <= self.counts['remap']:
                raise ValueError(
                    'Aligned original fragments cannot exceed mapped fragments')
            self.counts['aligned'] = count
            return
        if self.aligned_csv is None:
            return

        reader = DictReader(self.aligned_csv)
        self.counts['aligned'] = sum(map(int, map(itemgetter('count'), reader)))


def main():
    args = parse_args()
    report = CascadeReport(args.cascade_csv)
    report.g2p_summary_csv = args.g2p_summary_csv
    report.remap_counts_csv = args.remap_counts_csv
    report.aligned_csv = args.aligned_csv

    report.generate()


if __name__ == '__main__':
    main()
