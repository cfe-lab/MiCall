""" Compare result files in shared folder with previous release. """
from argparse import ArgumentParser
import csv
from difflib import SequenceMatcher
from itertools import zip_longest
from glob import glob
from operator import itemgetter
import os

from micall.core.aln2counts import AMINO_ALPHABET
from micall.settings import pipeline_version, DONE_PROCESSING

MICALL_VERSION = '7.7'
# set((sample, seed)) that had a coverage score better than 1 in source or target.
scored_samples = set()
# set((sample, target_seed, region)) that had lowest coverage on a partial deletion that got fixed.
fixed_deletions = set()


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


def select_columns(rows, column_names, extra_names=''):
    """ Build lines and lists from a CSV file with only some of the columns.

    :param rows: rows from a CSV file reader
    :param column_names: columns to include in lines, space delimited string
    :param extra_names: columns to include in list, space delimited string
    :return: [line], [columns], sorted by lines, columns include all columns
    """
    main_names = column_names.split()
    all_names = main_names + extra_names.split()
    main_getter = itemgetter(*main_names)
    all_getter = itemgetter(*all_names)
    combined_results = []
    for row in rows:
        score = row.get('on.score', None)
        if score is not None and score != '1':
            scored_samples.add((row['sample'], row['seed']))
        x4_pct = row.get('X4pct', None)
        if x4_pct not in (None, ''):
            row['X4pct'] = str(int(round(float(x4_pct))))

        if MICALL_VERSION == '7.7':
            # Version 7.7 switched HIV seed
            for field in ('seed', 'region'):
                field_value = row.get(field, '')
                if field_value.startswith('HIV1') and field_value.endswith('-seed'):
                    row[field] = 'HIV1???-seed'
        line = ',   '.join(main_getter(row)) + '\n'
        fields = list(all_getter(row))
        combined_results.append((line, fields))
    combined_results.sort()
    return zip(*combined_results)


# noinspection PyUnusedLocal
def should_never_skip(source_fields, target_fields):
    """ Dummy function that reports all differences. """
    return False


def compare_files(source_path,
                  target_path,
                  filename,
                  column_names,
                  extra_names='',
                  should_skip=should_never_skip):
    """ Compare columns within a file.

    Build text lines using column_names, and pass all columns in column_names
    and extra_names to should_skip() if differences are found.
    :param source_path: source folder
    :param target_path: target folder
    :param filename: the file name to read from both folders
    :param column_names: space delimited string
    :param extra_names: space delimited string, may be blank
    :param should_skip: function to check if a difference is small enough to
        skip.
    """
    sample_names = set()
    target_rows = select_rows(os.path.join(target_path, filename),
                              all_sample_names=sample_names)
    target_columns = select_columns(target_rows, column_names, extra_names)
    source_rows = select_rows(os.path.join(source_path, filename),
                              filter_sample_names=sample_names)
    source_columns = select_columns(source_rows, column_names, extra_names)
    compare_columns(source_columns, target_columns, target_path, filename, should_skip)


def compare_columns(source_columns, target_columns, target_path, filename, should_skip):
    print('{} changes in {}:'.format(filename, target_path))
    source_lines, source_fields = source_columns
    target_lines, target_fields = target_columns
    matcher = SequenceMatcher(a=source_lines, b=target_lines, autojunk=False)
    for tag, alo, ahi, blo, bhi in matcher.get_opcodes():
        if tag == 'replace':
            for source_index, target_index in zip_longest(range(alo, ahi), range(blo, bhi)):
                if source_index is None:
                    if not should_skip(None, target_fields[target_index]):
                        print('+ ' + target_lines[target_index], end='')
                elif target_index is None:
                    if not should_skip(source_fields[source_index], None):
                        print('- ' + source_lines[source_index], end='')
                else:
                    diff = line_diff(source_lines[source_index],
                                     target_lines[target_index])
                    if not should_skip(source_fields[source_index],
                                       target_fields[target_index]):
                        print(''.join(diff), end='')
        elif tag == 'delete':
            for line, fields in zip(source_lines[alo:ahi],
                                    source_fields[alo:ahi]):
                if not should_skip(fields, None):
                    print('- ' + line, end='')
        elif tag == 'insert':
            for line, fields in zip(target_lines[blo:bhi],
                                    target_fields[blo:bhi]):
                if not should_skip(None, fields):
                    print('+ ' + line, end='')
        else:
            assert tag == 'equal', tag


def format_key(key):
    return ',  '.join(key)


def compare_consensus(source_path, target_path):
    filename = 'conseq.csv'
    sample_names = set()
    target_rows = [row
                   for row in select_rows(os.path.join(target_path, filename),
                                          all_sample_names=sample_names)
                   if (row['sample'], row['region']) in scored_samples]
    source_rows = [row
                   for row in select_rows(os.path.join(source_path, filename),
                                          filter_sample_names=sample_names)
                   if (row['sample'], row['region']) in scored_samples]

    if MICALL_VERSION == '7.7':
        # Version 7.7 made such big changes to consensus that they're not worth comparing.
        field_names = 'sample region consensus-percent-cutoff'
    else:
        field_names = 'sample region consensus-percent-cutoff offset sequence'
    source_columns = select_columns(source_rows, field_names)
    target_columns = select_columns(target_rows, field_names)

    compare_columns(source_columns,
                    target_columns,
                    target_path,
                    filename,
                    should_skip_consensus)


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


def should_skip_consensus(source_fields, target_fields):
    if MICALL_VERSION != '7.7':
        return False
    # Version 7.7 started eliminating consensus for low coverage.
    if target_fields is None:
        return True
    if source_fields is None:
        return False
    return source_fields == target_fields


def should_skip_coverage_score(source_fields, target_fields):
    source_score = source_fields and int(source_fields[4])
    target_score = target_fields and int(target_fields[4])
    if source_fields is None:
        return target_score == 1
    if target_fields is None:
        return source_score == 1
    if source_score == 1 and target_score == 1:
        return True

    source_coverage = int(source_fields[5])
    target_coverage = int(target_fields[5])
    coverage_diff = source_coverage - target_coverage
    if source_coverage != 0:
        coverage_change = coverage_diff / float(source_coverage)
    else:
        coverage_change = 1.0
    if coverage_diff < 10 or coverage_change < 0.01 or source_score == target_score:
        source_fields[4:6] = target_fields[4:6]
    if source_fields == target_fields:
        return True
    if MICALL_VERSION != '7.7':
        return False
    # Version 7.7 corrected some partial deletions to full deletions.
    sample = source_fields[0]
    region = source_fields[2]
    target_seed = target_fields[1]
    if (target_score > source_score and
            (sample, target_seed, region) in fixed_deletions):
        source_fields[4] = target_fields[4]
    return source_fields == target_fields


def should_skip_g2p_summary(source_fields, target_fields):
    if MICALL_VERSION == '7.7':
        # Version 7.7 stopped making calls on low coverage.
        target_call = target_fields[2]
        if target_call.strip() == '':
            return True

    # We can ignore small changes in %X4.
    source_percent = int(source_fields[1]) if source_fields[1] else -1000
    target_percent = int(target_fields[1]) if target_fields[1] else -1000
    if abs(source_percent - target_percent) <= 1:
        source_fields[1] = target_fields[1]
    return source_fields == target_fields


def check_partial_deletions(source_path, target_path):
    fixed_deletions.clear()
    with open(os.path.join(source_path, 'coverage_scores.csv')) as old_scores:
        min_coverages = {(row['sample'],
                          'HIV' if row['seed'].startswith('HIV') else row['seed'],
                          row['region'],
                          int(row['which.key.pos'])): int(row['min.coverage'])
                         for row in csv.DictReader(old_scores)}
    coverage_fields = list(AMINO_ALPHABET) + ['del']
    with open(os.path.join(target_path, 'amino.csv')) as new_amino:
        for row in csv.DictReader(new_amino):
            current_pos = int(row['refseq.aa.pos'])
            new_deletions = int(row['del'])
            if new_deletions > 0:
                new_coverage = sum(int(row[field])
                                   for field in coverage_fields)
                seed = 'HIV' if row['seed'].startswith('HIV') else row['seed']
                for pos in range(current_pos-10, current_pos+10):
                    old_coverage = min_coverages.get(
                        (row['sample'], seed, row['region'], pos))
                    if old_coverage is not None:
                        if new_deletions >= (new_coverage - old_coverage) * 0.9:
                            fixed_deletions.add((row['sample'], row['seed'], row['region']))


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
        scored_samples.clear()
        if MICALL_VERSION == '7.7':
            check_partial_deletions(source_path, target_path)
        compare_files(source_path,
                      target_path,
                      'coverage_scores.csv',
                      'sample seed region project on.score',
                      'min.coverage',
                      should_skip=should_skip_coverage_score)
        compare_files(source_path,
                      target_path,
                      'g2p_summary.csv',
                      'sample X4pct final',
                      should_skip=should_skip_g2p_summary)
        compare_consensus(source_path, target_path)
        print('Done: ' + run_name)
    print('Done.')


if __name__ == '__main__':
    main()
