""" Compare result files in shared folder with previous release. """
from argparse import ArgumentParser
import csv
from collections import namedtuple, defaultdict, Counter
from difflib import Differ
from enum import IntEnum
from functools import partial
from itertools import groupby, zip_longest, chain
from glob import glob
from multiprocessing.pool import Pool
from operator import itemgetter
import os

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import Levenshtein

MICALL_VERSION = '7.11'

MiseqRun = namedtuple('MiseqRun', 'source_path target_path is_done')
MiseqRun.__new__.__defaults__ = (None,) * 3
SampleFiles = namedtuple(
    'SampleFiles',
    'cascade coverage_scores g2p_summary consensus nuc_limits remap_counts')
SampleFiles.__new__.__defaults__ = (None,) * 6
Sample = namedtuple('Sample', 'run name source_files target_files')
ConsensusDistance = namedtuple('ConsensusDistance',
                               'region cutoff distance pct_diff')
ConsensusDistance.__new__.__defaults__ = (None,) * 4
SampleComparison = namedtuple('SampleComparison',
                              ['diffs',  # multi-line string
                               'scenarios',  # {Scenarios: [description]}
                               'consensus_distances'])  # [ConsensusDistance]


class Scenarios(IntEnum):
    NONE = 0
    REMAP_COUNTS_CHANGED = 1
    MAIN_CONSENSUS_CHANGED = 2
    OTHER_CONSENSUS_CHANGED = 4


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
                                   'version_' + MICALL_VERSION)
        done_path = os.path.join(target_path, 'doneprocessing')
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
            source_nuc_limits = group_nucs(os.path.join(run.source_path,
                                                        'nuc.csv'))
            target_nuc_limits = group_nucs(os.path.join(run.target_path,
                                                        'nuc.csv'))
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
                                     source_nuc_limits.get(sample_name),
                                     source_counts.get(sample_name)),
                         SampleFiles(target_cascades.get(sample_name),
                                     target_coverages.get(sample_name),
                                     target_g2ps.get(sample_name),
                                     target_consensus.get(sample_name),
                                     target_nuc_limits.get(sample_name),
                                     target_counts.get(sample_name)))
    if missing_targets:
        print('Missing targets: ', missing_targets)
    if missing_sources:
        print('Missing sources: ', missing_sources)
    print('\n'.join(missing_files))
    if bad_files:
        print(bad_files)


def group_samples(output_file_name):
    with open(output_file_name) as output_file:
        return group_samples_file(output_file)


def group_samples_file(output_file):
    reader = csv.DictReader(output_file)
    return {key: list(rows)
            for key, rows in groupby(reader, itemgetter('sample'))}


def group_nucs(output_file_name):
    with open(output_file_name) as output_file:
        return group_nucs_file(output_file)


def group_nucs_file(output_file):
    reader = csv.DictReader(output_file)
    groups = {}
    for sample, sample_rows in groupby(reader, itemgetter('sample')):
        sample_limits = {}
        for seed, seed_rows in groupby(sample_rows, itemgetter('seed')):
            seed_limits = []
            for region, region_rows in groupby(seed_rows, itemgetter('region')):
                positions = [int(row['query.nuc.pos'])
                             for row in region_rows
                             if row['query.nuc.pos'] != '']
                seed_limits.append((region, min(positions), max(positions)))
            sample_limits[seed] = seed_limits
        groups[sample] = sample_limits
    return groups


def get_run_name(sample):
    return os.path.basename(
        os.path.dirname(os.path.dirname(sample.run.target_path)))


def compare_g2p(sample, diffs):
    source_fields = sample.source_files.g2p_summary
    target_fields = sample.target_files.g2p_summary
    if source_fields == target_fields:
        return
    if MICALL_VERSION == '7.11' and source_fields is None and sample.name.startswith('90308A'):
        return
    assert source_fields is not None, (sample.run, sample.name)
    assert target_fields is not None, (sample.run, sample.name)
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
        if (source_compare == '-' and
                sample.name.startswith('90308A') and
                MICALL_VERSION == '7.11'):
            # One sample failed in 7.10.
            pass
        elif source_compare != target_compare:
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
            else:
                diffs.append(message)


def adjust_region(region):
    if region.startswith('HCV'):
        parts = region.split('-')
        return 'HCV-' + parts[-1]
    return region


def is_consensus_interesting(region, cutoff):
    if region.startswith('HLA-'):
        return cutoff == '0.250'
    return cutoff == 'MAX'


def map_consensus_sequences(files):
    consensus_sequences = {}
    if None in (files.coverage_scores, files.nuc_limits, files.consensus):
        return consensus_sequences
    covered_regions = {(row.get('seed'), row.get('region'))
                       for row in files.coverage_scores
                       if row['on.score'] == '4'}
    # noinspection PyArgumentList
    region_limits = defaultdict(list,
                                {seed: [(adjust_region(region), first, last)
                                        for (region, first, last) in regions
                                        if (seed, region) in covered_regions]
                                 for seed, regions in files.nuc_limits.items()})

    for row in files.consensus:
        seed = row['region']
        cutoff = row['consensus-percent-cutoff']
        sequence = row['sequence'].rstrip('-')
        offset = int(row['offset'])
        for region, first, last in region_limits[seed]:
            if first > offset:
                subsequence = sequence[first-offset-1:]
            else:
                subsequence = '-' * (offset-first+1) + sequence
            subsequence = subsequence[:last-first+1]
            consensus_sequences[(seed, region, cutoff)] = subsequence
    return consensus_sequences


def display_consensus(sequence):
    if sequence is None:
        return []
    return [sequence + '\n']


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


def calculate_distance(region, cutoff, sequence1, sequence2):
    if sequence1 is None or sequence2 is None:
        return
    if len(sequence1) > len(sequence2):
        sequence2 += '-' * (len(sequence1) - len(sequence2))
    elif len(sequence2) > len(sequence1):
        sequence1 += '-' * (len(sequence2) - len(sequence1))
    distance = Levenshtein.distance(sequence1, sequence2)
    return ConsensusDistance(region=region,
                             cutoff=cutoff,
                             distance=distance,
                             pct_diff=distance/len(sequence1)*100)


def compare_consensus(sample,
                      source_seqs,
                      target_seqs,
                      diffs,
                      scenarios_reported,
                      scenarios):
    consensus_distances = []
    duplicates = find_duplicates(sample)
    if duplicates:
        diffs.extend(duplicates)
        return consensus_distances
    run_name = get_run_name(sample)
    keys = sorted(set(source_seqs.keys()) | target_seqs.keys())
    for key in keys:
        seed, region, cutoff = key
        source_fields = source_seqs.get(key)
        target_fields = target_seqs.get(key)
        consensus_distance = calculate_distance(region,
                                                cutoff,
                                                source_fields,
                                                target_fields)
        if False and consensus_distance.pct_diff > 5:
            print(sample.run.source_path, sample.name, consensus_distance)
            print(source_fields)
            print(target_fields)
        if consensus_distance is not None:
            consensus_distances.append(consensus_distance)
        if source_fields == target_fields:
            continue
        if (MICALL_VERSION == '7.11' and source_fields is None and
                sample.name.startswith('90308A')):
            continue
        is_main = is_consensus_interesting(region, cutoff)
        if is_main and scenarios_reported & Scenarios.MAIN_CONSENSUS_CHANGED:
            scenarios[Scenarios.MAIN_CONSENSUS_CHANGED].append('.')
        elif (not is_main and
                scenarios_reported & Scenarios.OTHER_CONSENSUS_CHANGED):
            scenarios[Scenarios.OTHER_CONSENSUS_CHANGED].append('.')
        else:
            diffs.append('{}:{} consensus: {} {} {}'.format(run_name,
                                                            sample.name,
                                                            seed,
                                                            region,
                                                            cutoff))
            diff = list(differ.compare(display_consensus(source_fields),
                                       display_consensus(target_fields)))
            diffs.extend(line.rstrip() for line in diff)
    return consensus_distances


def trim_consensus_sequences(target_seqs):
    old_sections = {('HCV-2', 'HCV-NS2'): slice(36, None),
                    ('HCV-2', 'HCV-NS5a'): slice(None, -114),
                    ('HCV-3', 'HCV-E1'): slice(9, None),
                    ('HCV-4', 'HCV-E1'): slice(6, None),
                    ('HCV-5', 'HCV-E1'): slice(6, None),
                    ('HCV-6', 'HCV-E1'): slice(6, None),
                    ('HIV1-', 'HIV1B-vpr'): slice(None, 234)}
    keys = list(target_seqs)
    for key in keys:
        seed_name, region_name, cutoff = key
        genotype = seed_name[:5]
        old_section_slice = old_sections.get((genotype, region_name))
        if old_section_slice is not None:
            target_seqs[key] = target_seqs[key][old_section_slice]


def compare_sample(sample,
                   scenarios_reported=Scenarios.NONE):
    scenarios = defaultdict(list)
    diffs = []
    compare_g2p(sample, diffs)
    compare_coverage(sample, diffs, scenarios_reported, scenarios)
    source_seqs = map_consensus_sequences(sample.source_files)
    target_seqs = map_consensus_sequences(sample.target_files)
    if MICALL_VERSION == '7.10':
        trim_consensus_sequences(target_seqs)
    consensus_distances = compare_consensus(sample,
                                            source_seqs,
                                            target_seqs,
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


def plot_distances(distance_data, filename, title, plot_variable='distance'):
    seeds = sorted(set(distance_data['region']))
    distance_data = distance_data.sort_values(['region', 'cutoff'])
    sns.set()
    num_plots = len(seeds)
    figure, axes_sets = plt.subplots(nrows=num_plots, ncols=1, squeeze=False)
    axes_sets = list(chain(*axes_sets))  # 2-dim array -> 1-dim list
    for ax, seed in zip(axes_sets, seeds):
        seed_data = distance_data[distance_data['region'] == seed]
        seed_data = seed_data.assign(
            count=lambda df: df['cutoff'].map(
                df.groupby(by=['cutoff'])[plot_variable].count()))
        seed_data['cutoff_n'] = seed_data.apply(format_cutoff, 'columns')

        sns.violinplot(x='cutoff_n',
                       y=plot_variable,
                       data=seed_data,
                       cut=0,
                       alpha=0.7,
                       ax=ax)
        plt.setp(ax.lines, zorder=100)
        plt.setp(ax.collections, zorder=100)
        sns.swarmplot(x='cutoff_n',
                      y=plot_variable,
                      data=seed_data,
                      color='k',
                      ax=ax)
        ax.set_ylabel(seed + '\n' + plot_variable)
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
                                scenarios_reported=Scenarios.OTHER_CONSENSUS_CHANGED),
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
    non_zero_distances = distance_data[distance_data['distance'] != 0]
    region_names = sorted(non_zero_distances['region'].unique())
    names_iter = iter(region_names)
    for page_num, region_group in enumerate(zip_longest(names_iter, names_iter, names_iter), 1):
        group_distances = distance_data[distance_data['region'].isin(region_group)]
        plot_distances(group_distances,
                       'consensus_distances_{}.svg'.format(page_num),
                       'Consensus Distances Between Previous and v' + MICALL_VERSION)
        plot_distances(group_distances,
                       'consensus_diffs_{}.svg'.format(page_num),
                       'Consensus Differences Between Previous and v' + MICALL_VERSION,
                       'pct_diff')
    print('Finished {} samples.'.format(i))


if __name__ == '__main__':
    main()
