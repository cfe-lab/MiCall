""" Compare result files in shared folder with previous release. """
from argparse import ArgumentParser
import csv
from collections import namedtuple, defaultdict, Counter
from difflib import Differ
from enum import IntEnum
from functools import partial
from itertools import groupby
from glob import glob
from multiprocessing.pool import Pool
from operator import itemgetter
import os

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import Levenshtein

from micall.settings import pipeline_version, DONE_PROCESSING

MICALL_VERSION = '7.8'

MiseqRun = namedtuple('MiseqRun', 'source_path target_path is_done')
MiseqRun.__new__.__defaults__ = (None,) * 3
SampleFiles = namedtuple(
    'SampleFiles',
    'cascade coverage_scores g2p_summary consensus remap_counts')
SampleFiles.__new__.__defaults__ = (None,) * 5
Sample = namedtuple('Sample', 'run name source_files target_files')
ConsensusDistance = namedtuple('ConsensusDistance', 'target_seed cutoff distance')
SampleComparison = namedtuple('SampleComparison',
                              ['diffs',  # multi-line string
                               'scenarios',  # {Scenarios: [description]}
                               'consensus_distances'])  # [ConsensusDistance]


class Scenarios(IntEnum):
    NONE = 0
    REMAP_COUNTS_CHANGED = 1
    MAIN_CONSENSUS_CHANGED = 2
    OTHER_CONSENSUS_CHANGED = 4
    V78_KEY_POS_REMOVED_RT318 = 8
    V78_KEY_POS_ADDED_INT263 = 16
    V78_KEY_POS_ADDED_INT163 = 32


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
    bad_files = []
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
        except Exception as ex:
            bad_files.append((ex, run.target_path))
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
    if bad_files:
        print(bad_files)


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


def compare_coverage(sample, diffs, scenarios_reported, scenarios):
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
            if (source_counts != target_counts and
                    scenarios_reported & Scenarios.REMAP_COUNTS_CHANGED):
                scenarios[Scenarios.REMAP_COUNTS_CHANGED].append(scenario)
            elif (MICALL_VERSION == '7.8' and
                  (region, source_key_pos) == ('RT', '318') and
                  scenarios_reported & Scenarios.V78_KEY_POS_REMOVED_RT318):
                scenarios[Scenarios.V78_KEY_POS_REMOVED_RT318].append(scenario)
            elif (MICALL_VERSION == '7.8' and
                  (region, target_key_pos) == ('INT', '263') and
                  scenarios_reported & Scenarios.V78_KEY_POS_ADDED_INT263):
                scenarios[Scenarios.V78_KEY_POS_ADDED_INT263].append(scenario)
            elif (MICALL_VERSION == '7.8' and
                  (region, target_key_pos) == ('INT', '163') and
                  scenarios_reported & Scenarios.V78_KEY_POS_ADDED_INT163):
                scenarios[Scenarios.V78_KEY_POS_ADDED_INT163].append(scenario)
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


def is_consensus_interesting(region, cutoff):
    if region.startswith('HLA-'):
        return cutoff == '0.250'
    return cutoff == 'MAX'


def map_consensus_sequences(rows):
    return {(adjust_region(row['region']),
             row['consensus-percent-cutoff']): (adjust_offset(row['offset'],
                                                              row['region']),
                                                row['sequence'].rstrip('-'))
            for row in rows}


def display_consensus(fields):
    if fields is None:
        return []
    offset, seq = fields
    return [offset + ' ' + seq + '\n']


def find_duplicates(sample):
    rows = sample.target_files.consensus
    if rows is None:
        return
    counts = Counter((row['region'], row['consensus-percent-cutoff'])
                     for row in rows)
    duplicate_keys = sorted(key
                            for key, count in counts.items()
                            if count > 1)
    return ['{}:{} duplicate consensus: {} {}'.format(get_run_name(sample),
                                                      sample.name,
                                                      region,
                                                      cutoff)
            for region, cutoff in duplicate_keys]


def calculate_distance(region, cutoff, source_fields, target_fields):
    if source_fields is None or target_fields is None:
        return
    offset1, sequence1 = source_fields
    offset2, sequence2 = target_fields
    try:
        offset1 = int(offset1)
        offset2 = int(offset2)
    except ValueError:
        offset1 = offset2 = 0
    if offset1 > offset2:
        sequence1 = '-'*(offset1-offset2) + sequence1
    else:
        sequence2 = '-'*(offset2-offset1) + sequence2
    if len(sequence1) > len(sequence2):
        sequence2 += '-' * (len(sequence1) - len(sequence2))
    elif len(sequence2) > len(sequence1):
        sequence1 += '-' * (len(sequence2) - len(sequence1))
    distance = Levenshtein.distance(sequence1, sequence2)
    if False and distance == 0:
        return None
    return ConsensusDistance(target_seed=region[:3],
                             cutoff=cutoff,
                             distance=distance)


def compare_consensus(sample,
                      diffs,
                      scenarios_reported,
                      scenarios):
    consensus_distances = []
    duplicates = find_duplicates(sample)
    if duplicates:
        diffs.extend(duplicates)
        return consensus_distances
    if sample.source_files.consensus == sample.target_files.consensus:
        return consensus_distances
    run_name = get_run_name(sample)
    source_seqs = map_consensus_sequences(sample.source_files.consensus)
    target_seqs = map_consensus_sequences(sample.target_files.consensus)
    keys = sorted(set(source_seqs.keys()) | target_seqs.keys())
    for key in keys:
        region, cutoff = key
        source_fields = source_seqs.get(key)
        target_fields = target_seqs.get(key)
        consensus_distance = calculate_distance(region, cutoff, source_fields, target_fields)
        if consensus_distance is not None:
            consensus_distances.append(consensus_distance)
        if source_fields == target_fields:
            continue
        is_main = is_consensus_interesting(region, cutoff)
        if is_main and scenarios_reported & Scenarios.MAIN_CONSENSUS_CHANGED:
            scenarios[Scenarios.MAIN_CONSENSUS_CHANGED].append('.')
        elif (not is_main and
                scenarios_reported & Scenarios.OTHER_CONSENSUS_CHANGED):
            scenarios[Scenarios.OTHER_CONSENSUS_CHANGED].append('.')
        else:
            diffs.append('{}:{} consensus: {} {}'.format(run_name,
                                                         sample.name,
                                                         region,
                                                         cutoff))
            diff = list(differ.compare(display_consensus(source_fields),
                                       display_consensus(target_fields)))
            diffs.extend(line.rstrip() for line in diff)
    return consensus_distances


def compare_sample(sample,
                   scenarios_reported=Scenarios.NONE):
    scenarios = defaultdict(list)
    diffs = []
    compare_g2p(sample, diffs)
    compare_coverage(sample, diffs, scenarios_reported, scenarios)
    consensus_distances = compare_consensus(sample,
                                            diffs,
                                            scenarios_reported,
                                            scenarios)
    diffs.append('')
    return SampleComparison(diffs='\n'.join(diffs),
                            scenarios=scenarios,
                            consensus_distances=consensus_distances)


def format_cutoff(row):
    cutoff = row['cutoff']
    count = row['count']
    try:
        cutoff = float(cutoff)
        cutoff = int(100 * cutoff)
        cutoff = str(cutoff)
    except ValueError:
        pass
    return cutoff + '_' + str(count)


def plot_distances(distance_data, filename, title):
    seeds = sorted(set(distance_data['target_seed']))
    distance_data = distance_data.sort_values(['target_seed', 'cutoff'])
    sns.set()
    num_plots = len(seeds)
    figure, axes_sets = plt.subplots(nrows=num_plots, ncols=1)
    for ax, seed in zip(axes_sets, seeds):
        seed_data = distance_data[distance_data['target_seed'] == seed]
        seed_data = seed_data.assign(
            count=lambda df: df['cutoff'].map(
                df.groupby(by=['cutoff'])['distance'].count()))
        seed_data['cutoff_n'] = seed_data.apply(format_cutoff, 'columns')

        sns.violinplot(x='cutoff_n',
                       y='distance',
                       data=seed_data,
                       cut=0,
                       alpha=0.7,
                       ax=ax)
        plt.setp(ax.lines, zorder=100)
        plt.setp(ax.collections, zorder=100)
        sns.swarmplot(x='cutoff_n',
                      y='distance',
                      data=seed_data,
                      color='k',
                      ax=ax)
        ax.set_ylabel(seed + ' distance')
        ax.get_yaxis().set_label_coords(-0.07, 0.5)
    axes_sets[0].set_title(title)
    plt.savefig(filename)


def main():
    print('Starting.')
    args = parse_args()
    pool = Pool()
    runs = find_runs(args.source_folder, args.target_folder)
    runs = report_source_versions(runs)
    samples = read_samples(runs)
    # noinspection PyTypeChecker
    results = pool.imap(partial(compare_sample,
                                scenarios_reported=sum(Scenarios)),
                        samples,
                        chunksize=50)
    scenario_summaries = defaultdict(list)
    i = None
    all_consensus_distances = []
    report_count = 0
    for i, (report, scenarios, consensus_distances) in enumerate(results):
        if report:
            report_count += 1
            if report_count > 100:
                break
        print(report, end='')
        all_consensus_distances.extend(consensus_distances)
        for key, messages in scenarios.items():
            scenario_summaries[key] += scenarios[key]
    for key, messages in sorted(scenario_summaries.items()):
        if messages:
            sample_names = {message.split()[0] for message in messages}
            summary = [key, len(messages), 'changes']
            body = ''.join(messages).rstrip('.')
            if body:
                summary.extend(['in', len(sample_names), 'samples'])
            print(*summary, end='.\n')
            print(body, end='')

    distance_data = pd.DataFrame(all_consensus_distances)
    plot_distances(distance_data,
                   'consensus_distances.svg',
                   'Consensus Distances Between v7.7 and v7.8')
    plot_distances(distance_data[distance_data['distance'] != 0],
                   'consensus_distances_nonzero.svg',
                   'Non-zero Consensus Distances Between v7.7 and v7.8')
    print('Finished {} samples.'.format(i))


if __name__ == '__main__':
    main()
