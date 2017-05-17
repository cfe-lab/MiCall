""" Compare result files in shared folder with previous release. """
from __future__ import print_function
from argparse import ArgumentParser
import csv
from difflib import Differ, SequenceMatcher
from itertools import zip_longest
from glob import glob
from operator import itemgetter
import os


from micall.settings import pipeline_version, DONE_PROCESSING


# set((sample, seed)) that had a coverage score better than 1 in source or target.
scored_samples = set()


def parse_args():
    parser = ArgumentParser(description='Compare sample results for testing a new release.')
    parser.add_argument('source_folder',
                        help='Main RAWDATA folder with results from previous version.')
    parser.add_argument('target_folder',
                        help='Testing RAWDATA folder to compare with.')
    return parser.parse_args()


def select_rows(csv_path,
                all_sample_names=None,
                filter_sample_names=None):
    """ Yield rows from a CSV file reader.

    :param csv_path: path to the CSV file
    :param all_sample_names: a set to receive all values for the 'sample'
        column, or None
    :param filter_sample_names: a set of sample names to include in the
        output, or None to not filter
    """
    with open(csv_path, 'rU') as f:
        reader = csv.DictReader(f)
        for row in reader:
            if all_sample_names is not None:
                all_sample_names.add(row['sample'])
            if filter_sample_names is None or row['sample'] in filter_sample_names:
                yield row


def select_columns(rows, *column_names):
    """ Yield lines from a CSV file with only some of the columns.

    :param rows: rows from a CSV file reader
    :param column_names: columns to include
    """
    for row in rows:
        score = row.get('on.score', None)
        if score == '1':
            continue
        if score is not None:
            scored_samples.add((row['sample'], row['seed']))
        x4_pct = row.get('X4pct', None)
        if x4_pct not in (None, ''):
            row['X4pct'] = str(int(round(float(x4_pct))))
        yield ',   '.join(itemgetter(*column_names)(row)) + '\n'


def compare_files(source_path, target_path, filename, *columns):
    sample_names = set()
    target_rows = select_rows(os.path.join(target_path, filename),
                              all_sample_names=sample_names)
    target_columns = sorted(select_columns(target_rows, *columns))
    source_rows = select_rows(os.path.join(source_path, filename),
                              filter_sample_names=sample_names)
    source_columns = sorted(select_columns(source_rows, *columns))
    compare_columns(source_columns, target_columns, target_path, filename)


def compare_columns(source_columns, target_columns, target_path, filename):
    differ = Differ()
    for i, diff in enumerate(differ.compare(source_columns, target_columns)):
        if i == 0:
            print('{} changes in {}:'.format(filename, target_path))
        print(diff, end='')


def format_key(key):
    return ',  '.join(key)


def compare_consensus(source_path, target_path):
    filename = 'conseq.csv'
    sample_names = set()
    keygetter = itemgetter('sample', 'region', 'consensus-percent-cutoff')
    target_rows = sorted((row
                          for row in select_rows(os.path.join(target_path, filename),
                                                 all_sample_names=sample_names)
                          if (row['sample'], row['region']) in scored_samples),
                         key=keygetter)
    source_rows = sorted((row
                          for row in select_rows(os.path.join(source_path, filename),
                                                 filter_sample_names=sample_names)
                          if (row['sample'], row['region']) in scored_samples),
                         key=keygetter)
    target_keys = list(map(keygetter, target_rows))
    source_keys = list(map(keygetter, source_rows))
    for source_row, target_row in zip(source_rows, target_rows):
        source_offset = int(source_row['offset'])
        target_offset = int(target_row['offset'])
        min_offset = min(source_offset, target_offset)
        source_row['sequence'] = '-' * (source_offset-min_offset) + source_row['sequence']
        target_row['sequence'] = '-' * (target_offset-min_offset) + target_row['sequence']
        source_row['offset'] = target_row['offset'] = str(min_offset)
    source_lines = list(select_columns(source_rows,
                                       'sample',
                                       'region',
                                       'consensus-percent-cutoff',
                                       'offset',
                                       'sequence'))
    target_lines = list(select_columns(target_rows,
                                       'sample',
                                       'region',
                                       'consensus-percent-cutoff',
                                       'offset',
                                       'sequence'))
    print('{} changes in {}:'.format(filename, target_path))
    matcher = SequenceMatcher(a=source_keys, b=target_keys, autojunk=False)
    for tag, alo, ahi, blo, bhi in matcher.get_opcodes():
        if tag == 'replace':
            for source_index, target_index in zip_longest(range(alo, ahi), range(blo, bhi)):
                if source_index is None:
                    print('+ ' + format_key(target_keys[target_index]))
                elif target_index is None:
                    print('- ' + format_key(source_keys[source_index]))
                else:
                    source_key = source_keys[source_index]
                    target_key = target_keys[target_index]
                    if source_key != target_key:
                        source_line = format_key(source_key) + '\n'
                        target_line = format_key(target_key) + '\n'
                    else:
                        source_line = source_lines[source_index]
                        target_line = target_lines[target_index]
                    diff = line_diff(source_line, target_line)
                    print(''.join(diff), end='')
        elif tag == 'delete':
            for key in source_keys[alo:ahi]:
                print('- ' + ',  '.join(key))
        elif tag == 'insert':
            for key in target_keys[blo:bhi]:
                print('+ ' + ',  '.join(key))
        else:
            assert tag == 'equal', tag
            for key in source_keys[alo:ahi]:
                print('  ' + ',  '.join(key) + '==')


def line_diff(source_line, target_line):
    cruncher = SequenceMatcher(a=source_line, b=target_line, autojunk=False)
    lines = ['- ', '? ', '+ ', '? ']
    for tag, alo, ahi, blo, bhi in cruncher.get_opcodes():
        adiff = bdiff = ' '
        if tag == 'replace':
            adiff = bdiff = '^'
        elif tag == 'delete':
            adiff = '-'
        elif tag == 'insert':
            bdiff = '+'
        else:
            assert tag == 'equal', tag
        lines[0] += source_line[alo:ahi]
        lines[2] += target_line[blo:bhi]
        lines[1] += adiff * (ahi - alo)
        lines[3] += bdiff * (bhi - blo)

    lines[1] += '\n'
    lines[3] += '\n'
    return lines


def main():
    args = parse_args()
    run_paths = glob(os.path.join(args.target_folder, 'MiSeq', 'runs', '*'))
    run_paths.sort()
    for run_path in run_paths:
        run_name = os.path.basename(run_path)
        target_path = os.path.join(run_path,
                                   'Results',
                                   'version_' + pipeline_version)
        done_path = os.path.join(target_path, DONE_PROCESSING)
        if not os.path.exists(done_path):
            print('Not done: ' + run_name)
            continue
        source_results_path = os.path.join(args.source_folder,
                                           'MiSeq',
                                           'runs',
                                           run_name,
                                           'Results')
        source_versions = os.listdir(source_results_path)
        source_versions.sort()
        source_path = os.path.join(source_results_path, source_versions[-1])
        compare_files(source_path,
                      target_path,
                      'coverage_scores.csv',
                      'sample',
                      'seed',
                      'region',
                      'project',
                      'on.score')
        compare_files(source_path,
                      target_path,
                      'g2p_summary.csv',
                      'sample',
                      'X4pct',
                      'final')
        compare_consensus(source_path, target_path)
        print('Done: ' + run_name)
    print('Done.')


if __name__ == '__main__':
    main()
