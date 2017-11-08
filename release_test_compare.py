""" Compare result files in shared folder with previous release. """
from argparse import ArgumentParser
import csv
from collections import Counter, namedtuple, defaultdict
from difflib import SequenceMatcher
from itertools import zip_longest, groupby
from glob import glob
from multiprocessing.pool import Pool
from operator import itemgetter
import os

from micall.core.aln2counts import AMINO_ALPHABET
from micall.settings import pipeline_version, DONE_PROCESSING

MICALL_VERSION = '7.8'
# set((sample, seed)) that had a coverage score better than 1 in source or target.
scored_samples = set()
# set((sample, target_seed, region)) that had lowest coverage on a partial deletion that got fixed.
fixed_deletions = set()


MiseqRun = namedtuple('MiseqRun', 'source_path target_path is_done')
MiseqRun.__new__.__defaults__ = (None,) * 3
SampleFiles = namedtuple(
    'SampleFiles',
    'cascade coverage_scores g2p_summary consensus remap_counts')
SampleFiles.__new__.__defaults__ = (None,) * 5
Sample = namedtuple('Sample', 'run name source_files target_files')


def parse_args():
    parser = ArgumentParser(description='Compare sample results for testing a new release.')
    parser.add_argument('source_folder',
                        help='Main RAWDATA folder with results from previous version.')
    parser.add_argument('target_folder',
                        help='Testing RAWDATA folder to compare with.')
    return parser.parse_args()


class ChangeLogger:
    def __init__(self):
        self.file = None
        self.is_changed = False
        self.is_run_changed = False
        self.unchanged_runs = []
        self.scenario_counts = Counter()

    def start(self, target_file):
        old_run = self.file and os.path.dirname(self.file)
        self.file = target_file
        self.is_changed = False
        new_run = self.file and os.path.dirname(self.file)
        if old_run and new_run != old_run:
            if not self.is_run_changed:
                run_name = '_'.join(os.path.basename(os.path.dirname(
                    os.path.dirname(old_run))).split('_')[:2])
                self.unchanged_runs.append(run_name)
            self.is_run_changed = False

    def report(self, diff, end='\n'):
        if not self.is_changed:
            filename = os.path.basename(self.file)
            folder = os.path.dirname(self.file)
            print('{} changes in {}:'.format(filename, folder))
            self.is_changed = True
            self.is_run_changed = True
        print(diff, end=end)

    def report_scenario(self, label, count):
        self.scenario_counts[label] += count

    def summarize(self):
        self.start(None)
        print('Unchanged runs: ' + ', '.join(self.unchanged_runs))
        print('Scenario counts:')
        for label, count in self.scenario_counts.most_common():
            print(count, label)


change_logger = ChangeLogger()


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
    change_logger.start(csv_path)
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
    if combined_results:
        combined_results.sort()
        return zip(*combined_results)
    return [], []


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
    compare_columns(source_columns, target_columns, should_skip)


def compare_columns(source_columns, target_columns, should_skip):
    source_lines, source_fields = source_columns
    target_lines, target_fields = target_columns
    matcher = SequenceMatcher(a=source_lines, b=target_lines, autojunk=False)
    for tag, alo, ahi, blo, bhi in matcher.get_opcodes():
        if tag == 'replace':
            for source_index, target_index in zip_longest(range(alo, ahi), range(blo, bhi)):
                if source_index is None:
                    if not should_skip(None, target_fields[target_index]):
                        change_logger.report('+ ' + target_lines[target_index], end='')
                elif target_index is None:
                    if not should_skip(source_fields[source_index], None):
                        change_logger.report('- ' + source_lines[source_index], end='')
                else:
                    diff = line_diff(source_lines[source_index],
                                     target_lines[target_index])
                    if not should_skip(source_fields[source_index],
                                       target_fields[target_index]):
                        change_logger.report(''.join(diff), end='')
        elif tag == 'delete':
            for line, fields in zip(source_lines[alo:ahi],
                                    source_fields[alo:ahi]):
                if not should_skip(fields, None):
                    change_logger.report('- ' + line, end='')
        elif tag == 'insert':
            for line, fields in zip(target_lines[blo:bhi],
                                    target_fields[blo:bhi]):
                if not should_skip(None, fields):
                    change_logger.report('+ ' + line, end='')
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
    if MICALL_VERSION == '7.8':
        for rows in (source_rows, target_rows):
            for row in rows:
                row['sequence'] = row['sequence'].rstrip('-')
    source_columns = select_columns(source_rows, field_names)
    target_columns = select_columns(target_rows, field_names)

    compare_columns(source_columns, target_columns, should_skip_consensus)


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
    change_logger.report_scenario('consensuses compared', 1)
    if MICALL_VERSION != '7.8':
        return False
    if target_fields is None:
        return False
    if source_fields is None:
        return False
    # 7.8 renamed the HIV seed used for V3LOOP alignment.
    if target_fields[1] == 'HIV1-CON-XX-Consensus-seed':
        target_fields[1] = source_fields[1] = 'HIV-???'
        target_fields[3] = source_fields[3] = '???'

    old_consensus, new_consensus = source_fields[4], target_fields[4]
    dash_count = single_sub_count = big_count = 0
    opcodes = SequenceMatcher(a=old_consensus,
                              b=new_consensus,
                              autojunk=False).get_opcodes()
    for tag, old1, old2, new1, new2 in opcodes:
        if tag == 'equal':
            continue
        old_length = old2 - old1
        new_length = new2 - new1
        max_length = max(old_length, new_length)
        old_context = get_consensus_context(old_consensus, old1, old2)
        new_context = get_consensus_context(new_consensus, new1, new2)
        has_dash = '-' in old_context + new_context
        if has_dash:
            dash_count += 1
        elif max_length == 1:
            single_sub_count += 1
        else:
            big_count += 1

    if old_consensus != new_consensus:
        source_fields[4] = target_fields[4] = 'DIFF_IGNORED'
        source_fields[3] = target_fields[3] = 'DIFF_IGNORED'
    is_match = source_fields == target_fields
    if is_match:
        change_logger.report_scenario('consensus change near dashes', dash_count)
        change_logger.report_scenario('consensus single substitutions',
                                      single_sub_count)
        change_logger.report_scenario('consensus bigger substitutions',
                                      big_count)
    return is_match


def get_consensus_context(seq, start, end):
    margin = 5
    return seq[max(0, start-margin):end+margin]


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
    if MICALL_VERSION != '7.8':
        return False
    # Version 7.8 added a new seed
    if target_fields[1] == 'HIV1-CON-XX-Consensus-seed':
        source_fields[1] = target_fields[1]
        if source_fields == target_fields:
            change_logger.report_scenario(
                'coverage score switched to ' + target_fields[1],
                1)
    return source_fields == target_fields


def should_skip_g2p_summary(source_fields, target_fields):
    if MICALL_VERSION == '7.7':
        # Version 7.7 stopped making calls on low coverage.
        target_call = target_fields[2]
        if target_call.strip() == '':
            return True

    if source_fields is None or target_fields is None:
        return False

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


def find_runs(source_folder, target_folder):
    run_paths = glob(os.path.join(target_folder, 'MiSeq', 'runs', '*'))
    run_paths.sort()
    for run_path in run_paths:
        run_name = os.path.basename(run_path)
        target_path = os.path.join(run_path,
                                   'Results',
                                   'version_' + pipeline_version)
        done_path = os.path.join(target_path, DONE_PROCESSING)
        is_done = os.path.exists(done_path)
        source_results_path = os.path.join(source_folder,
                                           'MiSeq',
                                           'runs',
                                           run_name,
                                           'Results')
        source_versions = os.listdir(source_results_path)
        source_versions.sort()
        source_path = os.path.join(source_results_path, source_versions[-1])
        yield MiseqRun(source_path, target_path, is_done)


def report_source_versions(runs):
    version_runs = defaultdict(list)  # {version: [source_path]}
    for run in runs:
        version = os.path.basename(run.source_path)
        version_runs[version].append(run.source_path)
        yield run
    max_count = max(len(paths) for paths in version_runs.values())
    for version, paths in sorted(version_runs.items()):
        if len(paths) == max_count:
            print(version, max_count)
        else:
            print(version, paths)


def read_samples(runs):
    missing_sources = []
    missing_targets = []
    missing_files = []
    for run in runs:
        if run.source_path is None:
            missing_sources.append(run.target_path)
            continue
        if run.target_path is None:
            missing_targets.append(run.source_path)
            continue
        try:
            source_cascades = group_samples(os.path.join(run.source_path,
                                                         'cascade.csv'))
            target_cascades = group_samples(os.path.join(run.target_path,
                                                         'cascade.csv'))
            source_g2ps = group_samples(os.path.join(run.source_path,
                                                     'g2p_summary.csv'))
            target_g2ps = group_samples(os.path.join(run.target_path,
                                                     'g2p_summary.csv'))
            source_coverages = group_samples(os.path.join(run.source_path,
                                                          'coverage_scores.csv'))
            target_coverages = group_samples(os.path.join(run.target_path,
                                                          'coverage_scores.csv'))
            source_counts = group_samples(os.path.join(run.source_path,
                                                       'remap_counts.csv'))
            target_counts = group_samples(os.path.join(run.target_path,
                                                       'remap_counts.csv'))
        except FileNotFoundError as ex:
            missing_files.append(str(ex))
            continue

        sample_names = sorted(target_cascades.keys())
        for sample_name in sample_names:
            yield Sample(run,
                         sample_name,
                         SampleFiles(source_cascades.get(sample_name),
                                     source_coverages.get(sample_name),
                                     source_g2ps.get(sample_name),
                                     None,
                                     source_counts.get(sample_name)),
                         SampleFiles(target_cascades.get(sample_name),
                                     target_coverages.get(sample_name),
                                     target_g2ps.get(sample_name),
                                     None,
                                     target_counts.get(sample_name)))
    if missing_targets:
        print('Missing targets: ', missing_targets)
    if missing_sources:
        print('Missing sources: ', missing_sources)
    print('\n'.join(missing_files))


def group_samples(output_file):
    with open(output_file) as f:
        reader = csv.DictReader(f)
        return {key: list(rows)
                for key, rows in groupby(reader, itemgetter('sample'))}


def compare_g2p(sample, diffs):
    source_fields = sample.source_files.g2p_summary
    target_fields = sample.target_files.g2p_summary
    if source_fields == target_fields:
        return
    assert len(source_fields) == 1, source_fields
    assert len(target_fields) == 1, target_fields
    run_name = os.path.basename(
        os.path.dirname(os.path.dirname(sample.run.target_path)))
    source_x4_pct = source_fields[0]['X4pct']
    target_x4_pct = target_fields[0]['X4pct']
    source_final = source_fields[0].get('final')
    target_final = target_fields[0].get('final')
    if source_final != target_final:
        diffs.append('{}:{} G2P: {} {} => {} {}'.format(run_name,
                                                        sample.name,
                                                        source_final,
                                                        source_x4_pct,
                                                        target_final,
                                                        target_x4_pct))
        return
    if source_x4_pct == target_x4_pct:
        return
    try:
        x4_pct_diff = abs(float(target_x4_pct) - float(source_x4_pct))
        if x4_pct_diff < 2.0:
            return
    except ValueError:
        pass
    diffs.append('{}:{} G2P: {} => {}'.format(run_name,
                                              sample.name,
                                              source_x4_pct,
                                              target_x4_pct))


def map_coverage(coverage_scores):
    if coverage_scores is None:
        return {}
    return {(score['project'], score['region']): (score['on.score'],
                                                  score.get('seed'),
                                                  score.get('which.key.pos'))
            for score in coverage_scores}


def map_remap_counts(counts):
    if counts is None:
        return {}
    row_parts = [row['type'].split() for row in counts]
    count_parts = [(parts[0].split('-'), parts[1])
                   for parts in row_parts
                   if len(parts) == 2]
    seed_counts = {parts[1]: int(parts[0][1])
                   for parts in count_parts
                   if len(parts[0]) == 2 and parts[0][1] != 'final'}
    return seed_counts


def compare_coverage(sample, diffs, scenarios):
    if sample.source_files.coverage_scores == sample.target_files.coverage_scores:
        return
    source_scores = map_coverage(sample.source_files.coverage_scores)
    target_scores = map_coverage(sample.target_files.coverage_scores)
    if source_scores == target_scores:
        return
    source_counts = map_remap_counts(sample.source_files.remap_counts)
    target_counts = map_remap_counts(sample.target_files.remap_counts)
    run_name = os.path.basename(
        os.path.dirname(os.path.dirname(sample.run.target_path)))
    keys = sorted(set(source_scores.keys()) | target_scores.keys())
    for key in keys:
        (source_score,
         source_seed,
         source_key_pos) = source_scores.get(key, ('-', None, None))
        (target_score,
         target_seed,
         target_key_pos) = target_scores.get(key, ('-', None, None))
        source_compare = '-' if source_score == '1' else source_score
        target_compare = '-' if target_score == '1' else target_score
        if source_compare != target_compare:
            project, region = key
            message = '{}:{} coverage: {} {} {} => {}'.format(
                run_name,
                sample.name,
                project,
                region,
                source_score,
                target_score)
            scenario = '  ' + message + '\n'
            if source_counts != target_counts:
                scenarios['different remap counts'].append(scenario)
            elif (MICALL_VERSION == '7.8' and
                  (region, source_key_pos) == ('RT', '318')):
                scenarios['key pos removed RT 318'].append(scenario)
            elif (MICALL_VERSION == '7.8' and
                  (region, target_key_pos) == ('INT', '263')):
                scenarios['key pos added INT 263'].append(scenario)
            elif (MICALL_VERSION == '7.8' and
                  (region, target_key_pos) == ('INT', '163')):
                scenarios['key pos added INT 163'].append(scenario)
            else:
                diffs.append(message)


def compare_sample(sample):
    scenarios = defaultdict(list)
    diffs = []
    compare_g2p(sample, diffs)
    compare_coverage(sample, diffs, scenarios)
    diffs.append('')
    return '\n'.join(diffs), scenarios


def main():
    print('Starting.')
    args = parse_args()
    pool = Pool()
    runs = find_runs(args.source_folder, args.target_folder)
    runs = report_source_versions(runs)
    samples = read_samples(runs)
    results = pool.imap(compare_sample,
                        samples,
                        chunksize=50)
    scenario_summaries = defaultdict(list)
    i = None
    for i, (report, scenarios) in enumerate(results):
        print(report, end='')
        for key, messages in scenarios.items():
            scenario_summaries[key] += scenarios[key]
    for key, messages in sorted(scenario_summaries.items()):
        if messages:
            sample_names = {message.split()[0] for message in messages}
            print(key, len(messages), 'changes in', len(sample_names), 'samples')
            print(''.join(messages), end='')
    print('Finished {} samples.'.format(i))
    # run_paths = glob(os.path.join(args.target_folder, 'MiSeq', 'runs', '*'))
    # run_paths.sort()
    # for run_path in run_paths:
    #     run_name = os.path.basename(run_path)
    #     target_path = os.path.join(run_path,
    #                                'Results',
    #                                'version_' + pipeline_version)
    #     done_path = os.path.join(target_path, DONE_PROCESSING)
    #     if not os.path.exists(done_path):
    #         print('Not done: ' + run_name)
    #         continue
    #     source_results_path = os.path.join(args.source_folder,
    #                                        'MiSeq',
    #                                        'runs',
    #                                        run_name,
    #                                        'Results')
    #     source_versions = os.listdir(source_results_path)
    #     source_versions.sort()
    #     source_path = os.path.join(source_results_path, source_versions[-1])
    #     scored_samples.clear()
    #     if MICALL_VERSION == '7.7':
    #         check_partial_deletions(source_path, target_path)
    #     compare_files(source_path,
    #                   target_path,
    #                   'coverage_scores.csv',
    #                   'sample seed region project on.score',
    #                   'min.coverage',
    #                   should_skip=should_skip_coverage_score)
    #     compare_files(source_path,
    #                   target_path,
    #                   'g2p_summary.csv',
    #                   'sample X4pct final',
    #                   should_skip=should_skip_g2p_summary)
    #     compare_consensus(source_path, target_path)
    # change_logger.summarize()


if __name__ == '__main__':
    main()
