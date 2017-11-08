""" Compare result files in shared folder with previous release. """
from argparse import ArgumentParser
import csv
from collections import namedtuple, defaultdict
from difflib import Differ
from itertools import groupby
from glob import glob
from multiprocessing.pool import Pool
from operator import itemgetter
import os

from micall.settings import pipeline_version, DONE_PROCESSING

MICALL_VERSION = '7.8'

MiseqRun = namedtuple('MiseqRun', 'source_path target_path is_done')
MiseqRun.__new__.__defaults__ = (None,) * 3
SampleFiles = namedtuple(
    'SampleFiles',
    'cascade coverage_scores g2p_summary consensus remap_counts')
SampleFiles.__new__.__defaults__ = (None,) * 5
Sample = namedtuple('Sample', 'run name source_files target_files')

differ = Differ()


def parse_args():
    parser = ArgumentParser(description='Compare sample results for testing a new release.')
    parser.add_argument('source_folder',
                        help='Main RAWDATA folder with results from previous version.')
    parser.add_argument('target_folder',
                        help='Testing RAWDATA folder to compare with.')
    return parser.parse_args()


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
            source_consensus = group_samples(os.path.join(run.source_path,
                                                          'conseq.csv'))
            target_consensus = group_samples(os.path.join(run.target_path,
                                                          'conseq.csv'))
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
                                     source_consensus.get(sample_name),
                                     source_counts.get(sample_name)),
                         SampleFiles(target_cascades.get(sample_name),
                                     target_coverages.get(sample_name),
                                     target_g2ps.get(sample_name),
                                     target_consensus.get(sample_name),
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


def get_run_name(sample):
    return os.path.basename(
        os.path.dirname(os.path.dirname(sample.run.target_path)))


def compare_g2p(sample, diffs):
    source_fields = sample.source_files.g2p_summary
    target_fields = sample.target_files.g2p_summary
    if source_fields == target_fields:
        return
    assert len(source_fields) == 1, source_fields
    assert len(target_fields) == 1, target_fields
    run_name = get_run_name(sample)
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
    run_name = get_run_name(sample)
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


def adjust_region(region):
    if MICALL_VERSION == '7.8' and region.startswith('HIV1-'):
        return 'HIV1-???'
    return region


def adjust_offset(offset, region):
    if MICALL_VERSION == '7.8' and region.startswith('HIV1-'):
        return '???'
    return offset


def is_consensus_interesting(row):
    if row['region'].startswith('HLA-'):
        return row['consensus-percent-cutoff'] == '0.250'
    return row['consensus-percent-cutoff'] == 'MAX'


def map_consensus_sequences(rows):
    return {(adjust_region(row['region']),
             row['consensus-percent-cutoff']): (adjust_offset(row['offset'],
                                                              row['region']),
                                                row['sequence'].rstrip('-'))
            for row in rows
            if is_consensus_interesting(row)}


def display_consensus(fields):
    if fields is None:
        return []
    offset, seq = fields
    return [offset + ' ' + seq + '\n']


def compare_consensus(sample, diffs):
    if sample.source_files.consensus == sample.target_files.consensus:
        return
    run_name = get_run_name(sample)
    source_seqs = map_consensus_sequences(sample.source_files.consensus)
    target_seqs = map_consensus_sequences(sample.target_files.consensus)
    keys = sorted(set(source_seqs.keys()) | target_seqs.keys())
    for key in keys:
        region, cutoff = key
        source_fields = source_seqs.get(key)
        target_fields = target_seqs.get(key)
        if source_fields == target_fields:
            continue
        diffs.append('{}:{} consensus: {} {}'.format(run_name,
                                                     sample.name,
                                                     region,
                                                     cutoff))
        diff = list(differ.compare(display_consensus(source_fields),
                                   display_consensus(target_fields)))
        diffs.extend(line.rstrip() for line in diff)


def compare_sample(sample):
    scenarios = defaultdict(list)
    diffs = []
    compare_g2p(sample, diffs)
    compare_coverage(sample, diffs, scenarios)
    compare_consensus(sample, diffs)
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
    report_count = 0
    for i, (report, scenarios) in enumerate(results):
        if report:
            report_count += 1
            if report_count > 100:
                break
        print(report, end='')
        for key, messages in scenarios.items():
            scenario_summaries[key] += scenarios[key]
    for key, messages in sorted(scenario_summaries.items()):
        if messages:
            sample_names = {message.split()[0] for message in messages}
            print(key, len(messages), 'changes in', len(sample_names), 'samples')
            print(''.join(messages), end='')
    print('Finished {} samples.'.format(i))


if __name__ == '__main__':
    main()
