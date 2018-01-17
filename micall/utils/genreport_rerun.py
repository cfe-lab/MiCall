import csv
import os
import re
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from collections import namedtuple, defaultdict
from csv import DictReader, DictWriter
from itertools import groupby
from operator import itemgetter, attrgetter

from micall.core.aln2counts import AMINO_ALPHABET
from micall.hivdb.genreport import gen_report
from micall.hivdb.hivdb import hivdb

SampleInfo = namedtuple('SampleInfo', 'enum suffix project snum name')
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


def find_groups(working_paths):
    groups = defaultdict(list)
    for path in working_paths:
        basename = os.path.basename(path)
        sample_name, snum = basename.split('_')
        parts = sample_name.split('-')
        for project in ('HCV', 'MidHCV'):
            try:
                project_index = parts.index(project)
            except ValueError:
                continue
            extraction = parts[project_index-1]
            if extraction.endswith('MIDI'):
                extraction_num = extraction[:-4]
                suffix = 'MIDI'
            else:
                extraction_num = extraction
                suffix = ''
            groups[extraction_num].append(SampleInfo(extraction_num,
                                                     suffix,
                                                     project,
                                                     snum,
                                                     basename))
            break
    for extraction_num, samples in sorted(groups.items()):
        if len(samples) == 2:
            names = tuple(sample.name
                          for sample in sorted(samples,
                                               key=attrgetter('project')))
            yield SampleGroup(extraction_num, names)
        else:
            print("Couldn't group:", samples)


def parse_sample_info(sample_name):
    head, snum = sample_name.split('_')
    head, project = head.split('-')
    match = re.match(r'([A-Z]*\d+)(.*$)', head)
    return SampleInfo(enum=match.group(1),
                      suffix=match.group(2),
                      project=project,
                      snum=snum,
                      name=sample_name)


def combine_files(base_path, groups):
    amino_columns = list(AMINO_ALPHABET) + ['del', 'coverage']
    for group in groups:
        src_filename = os.path.join(base_path,
                                    group.names[1],
                                    'coverage_scores.csv')
        with open(src_filename) as src:
            reader = DictReader(src)
            has_good_coverage = False
            for row in reader:
                if row['region'].endswith('-NS5b') and row['on.score'] == '4':
                    has_good_coverage = True
                    break
        if has_good_coverage:
            dest_filename = os.path.join(base_path,
                                         group.names[0],
                                         'amino.csv')
            src_filename = os.path.join(base_path,
                                        group.names[1],
                                        'amino.csv')
            with open(src_filename) as src:
                reader = DictReader(src)
                source_rows = {(row['region'], row['refseq.aa.pos']): row
                               for row in reader
                               if row['region'].endswith('-NS5b')}
            dest_copyname = dest_filename + '.orig.csv'
            os.rename(dest_filename, dest_copyname)
            with open(dest_copyname) as src, open(dest_filename, 'w') as dest:
                reader = DictReader(src)
                writer = DictWriter(dest, reader.fieldnames)
                writer.writeheader()
                for row in reader:
                    source_row = source_rows.get((row['region'],
                                                  row['refseq.aa.pos']),
                                                 {})
                    for column in amino_columns:
                        dest_count = int(row[column])
                        source_count = int(source_row.get(column, '0'))
                        row[column] = dest_count + source_count
                    writer.writerow(row)
        yield os.path.join(base_path, group.names[0])


def main():
    args = parse_args()
    working_paths = split_files(args)

    sorted_working_paths = sorted(working_paths)
    groups = find_groups(sorted_working_paths)
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
