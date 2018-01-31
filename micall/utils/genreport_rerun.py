import csv
import os
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from collections import namedtuple
from csv import DictReader, DictWriter
from itertools import groupby
from operator import itemgetter

from micall.hivdb.genreport import gen_report
from micall.hivdb.hivdb import hivdb
from micall.utils.sample_sheet_parser import sample_sheet_parser

SampleGroup = namedtuple('SampleGroup', 'enum names')


def parse_args():
    parser = ArgumentParser(
        description='Rerun resistance interpretations on a run folder.',
        formatter_class=ArgumentDefaultsHelpFormatter)

    parser.add_argument('--source',
                        '-s',
                        help='source results folder')
    parser.add_argument('--working',
                        '-w',
                        help='working folder')
    return parser.parse_args()


def find_groups(working_paths, source_path):
    sample_sheet_path = os.path.join(source_path, '../../SampleSheet.csv')
    with open(sample_sheet_path) as sample_sheet_file:
        run_info = sample_sheet_parser(sample_sheet_file)

    midi_files = {row['sample']: row['filename']
                  for row in run_info['DataSplit']
                  if row['project'] == 'MidHCV'}
    wide_names = {row['filename']: row['sample']
                  for row in run_info['DataSplit']
                  if row['project'] == 'HCV'}
    for path in working_paths:
        wide_file = os.path.basename(path)
        sample_name = wide_names.get(wide_file)
        if sample_name is None:
            # Not an HCV sample.
            continue
        midi_file = midi_files.get(sample_name + 'MIDI')
        yield SampleGroup(sample_name, (wide_file, midi_file))


def rewrite_file(filename):
    backup_csv = filename + '.original.csv'
    os.rename(filename, backup_csv)
    with open(backup_csv) as source, open(filename, 'w') as dest:
        reader = DictReader(source)
        writer = DictWriter(dest, reader.fieldnames)
        writer.writeheader()
        for row in reader:
            yield row  # Caller can modify it.
            if row:  # Skips a row that got cleared.
                writer.writerow(row)


def main():
    args = parse_args()
    working_paths = split_files(args)

    sorted_working_paths = sorted(working_paths)
    groups = list(find_groups(sorted_working_paths, args.source))
    for group in groups:
        working_path = os.path.join(args.working, group.names[0])
        midi_path = working_path if group.names[1] is None else os.path.join(
            args.working,
            group.names[1])
        print(working_path)
        with open(os.path.join(working_path, 'amino.csv')) as amino_csv, \
                open(os.path.join(midi_path, 'amino.csv')) as midi_amino_csv, \
                open(os.path.join(working_path, 'resistance.csv'), 'w') as resistance_csv, \
                open(os.path.join(working_path, 'mutations.csv'), 'w') as mutations_csv, \
                open(os.path.join(working_path, 'resistance_fail.csv'), 'w') as resistance_fail_csv:
            hivdb(amino_csv,
                  midi_amino_csv,
                  resistance_csv,
                  mutations_csv,
                  resistance_fail_csv)
        sample_name = os.path.basename(working_path)
        with open(os.path.join(working_path, 'resistance.csv')) as resistance_csv, \
                open(os.path.join(working_path, 'mutations.csv')) as mutations_csv, \
                open(os.path.join(working_path, 'resistance_report.pdf'), 'wb') as resistance_report_csv:
            gen_report(resistance_csv,
                       mutations_csv,
                       resistance_report_csv,
                       sample_name=sample_name)

    for file_name in ('resistance.csv', 'mutations.csv', 'resistance_fail.csv'):
        with open(os.path.join(args.working, file_name), 'w') as dest:
            dest_writer = csv.writer(dest)
            for i, group in enumerate(groups):
                working_path = os.path.join(args.working, group.names[0])
                sample_name = os.path.basename(working_path)
                with open(os.path.join(working_path, file_name), 'r') as source:
                    source_reader = csv.reader(source)
                    for j, row in enumerate(source_reader):
                        if j != 0:
                            row.insert(0, sample_name)
                        elif i == 0:
                            row.insert(0, 'sample')
                        else:
                            continue
                        dest_writer.writerow(row)


def split_files(args):
    working_paths = set()
    file_name = 'amino.csv'
    file_path = os.path.join(args.source, file_name)
    with open(file_path) as f:
        reader = DictReader(f)
        for sample, rows in groupby(reader, itemgetter('sample')):
            working_path = os.path.join(args.working, sample)
            working_paths.add(working_path)
            if __name__ == '__live_coding__':
                if len(working_paths) > 20:
                    break
                continue
            os.makedirs(working_path, exist_ok=True)
            target_path = os.path.join(working_path, file_name)
            with open(target_path, 'w') as target_csv:
                writer = DictWriter(target_csv,
                                    reader.fieldnames[1:])
                writer.writeheader()
                for row in rows:
                    del row['sample']
                    writer.writerow(row)
    return working_paths


main()
