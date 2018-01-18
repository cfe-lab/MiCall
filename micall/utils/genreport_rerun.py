import csv
import os
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from collections import namedtuple
from csv import DictReader, DictWriter
from itertools import groupby
from operator import itemgetter

from micall.core.aln2counts import AMINO_ALPHABET
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


def combine_files(base_path, groups):
    for group in groups:
        if group.names[1] is not None:
            combine_midi(base_path, group.names[0], group.names[1])
        yield os.path.join(base_path, group.names[0])


def combine_midi(base_path, wide_name, midi_name):
    amino_columns = list(AMINO_ALPHABET) + ['del', 'coverage']
    src_filename = os.path.join(base_path,
                                midi_name,
                                'coverage_scores.csv')
    midi_covered_seeds = set()
    with open(src_filename) as src:
        reader = DictReader(src)
        for row in reader:
            if (row['region'].endswith('-NS5b') and
                    row['project'] == 'MidHCV' and
                    row['on.score'] == '4'):
                midi_covered_seeds.add(row['seed'])
                break
    dest_filename = os.path.join(base_path,
                                 wide_name,
                                 'coverage_scores.csv')
    has_good_coverage = False
    for row in rewrite_file(dest_filename):
        if (row['region'].endswith('-NS5b') and
                row['on.score'] == '4'):
            if row['seed'] in midi_covered_seeds:
                has_good_coverage = True
            else:
                row['on.score'] = '1'
    if has_good_coverage:
        dest_filename = os.path.join(base_path,
                                     wide_name,
                                     'amino.csv')
        src_filename = os.path.join(base_path,
                                    midi_name,
                                    'amino.csv')
        with open(src_filename) as src:
            reader = DictReader(src)
            source_rows = {(row['region'], row['refseq.aa.pos']): row
                           for row in reader
                           if row['region'].endswith('-NS5b')}
        for row in rewrite_file(dest_filename):
            source_row = source_rows.get((row['region'], row['refseq.aa.pos']))
            if source_row is not None:
                pos = int(row['refseq.aa.pos'])
                wide_coverage = int(row['coverage'])
                midi_coverage = int(row['coverage'])
                if pos > 335 or (pos >= 226 and midi_coverage > wide_coverage):
                    for column in amino_columns:
                        row[column] = source_row[column]


def main():
    args = parse_args()
    working_paths = split_files(args)

    sorted_working_paths = sorted(working_paths)
    groups = find_groups(sorted_working_paths, args.source)
    combined_working_paths = list(combine_files(args.working, groups))
    failed_working_paths = set(combined_working_paths)
    for working_path in combined_working_paths:
        print(working_path)
        with open(os.path.join(working_path, 'amino.csv')) as amino_csv, \
                open(os.path.join(working_path, 'coverage_scores.csv')) as coverage_scores_csv, \
                open(os.path.join(working_path, 'resistance.csv'), 'w') as resistance_csv, \
                open(os.path.join(working_path, 'mutations.csv'), 'w') as mutations_csv:
            hivdb(amino_csv,
                  coverage_scores_csv,
                  resistance_csv,
                  mutations_csv)
        sample_name = os.path.basename(working_path)
        with open(os.path.join(working_path, 'resistance.csv')) as resistance_csv, \
                open(os.path.join(working_path, 'mutations.csv')) as mutations_csv, \
                open(os.path.join(working_path, 'resistance_report.pdf'), 'wb') as resistance_report_csv:
            gen_report(resistance_csv,
                       mutations_csv,
                       resistance_report_csv,
                       sample_name=sample_name)

    for file_name in ('resistance.csv', 'mutations.csv'):
        with open(os.path.join(args.working, file_name), 'w') as dest:
            dest_writer = csv.writer(dest)
            for i, working_path in enumerate(combined_working_paths):
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
                        if j == 1:
                            failed_working_paths.discard(working_path)
                        dest_writer.writerow(row)

    with open(os.path.join(args.working, 'failed.csv'), 'w') as dest:
        dest_writer = csv.writer(dest)
        dest_writer.writerow(['sample'])
        for working_path in sorted(failed_working_paths):
            dest_writer.writerow([os.path.basename(working_path)])


def split_files(args):
    working_paths = set()
    for file_name in ('amino.csv', 'coverage_scores.csv'):
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
