""" Compare result files in shared folder with previous release. """
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import csv
from collections import namedtuple, defaultdict
from concurrent.futures.process import ProcessPoolExecutor, BrokenProcessPool
from difflib import Differ
from enum import IntEnum
from functools import partial
from itertools import groupby, zip_longest, chain
from glob import glob
from operator import itemgetter
import os
import logging

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
# noinspection PyPackageRequirements
import Levenshtein

from micall.utils.primer_tracker import PrimerTracker
from micall.utils.report_amino import SeedNucleotide, MAX_CUTOFF
from micall.utils.translation import translate
from micall.utils.analyze import get_available_memory

MICALL_VERSION = '7.15'
#                ^^^^^^ Version of the MiCall release being tested.
#                This is the new version against which older versions are compared.
# The version for the older revision is determined dynamically in the `find_runs` function.
# The source folder is inspected to find all previous result versions for each run.
# These versions are then sorted and the latest one is selected for comparison.

MiseqRun = namedtuple('MiseqRun', 'source_path target_path is_done')
MiseqRun.__new__.__defaults__ = (None,) * 3
SampleFiles = namedtuple(
    'SampleFiles',
    'cascade coverage_scores g2p_summary region_consensus remap_counts conseq')
SampleFiles.__new__.__defaults__ = (None,) * 6
Sample = namedtuple('Sample', 'run name source_files target_files')
ConsensusDistance = namedtuple('ConsensusDistance',
                               'region cutoff distance pct_diff')
ConsensusDistance.__new__.__defaults__ = (None,) * 4
SampleComparison = namedtuple('SampleComparison',
                              ['diffs',  # multi-line string
                               'scenarios',  # {Scenarios: [description]}
                               'consensus_distances'])  # [ConsensusDistance]

logger = logging.getLogger(__name__)


class Scenarios(IntEnum):
    NONE = 0
    REMAP_COUNTS_CHANGED = 1
    MAIN_CONSENSUS_CHANGED = 2
    OTHER_CONSENSUS_CHANGED = 4
    CONSENSUS_DELETIONS_CHANGED = 8
    VPR_FRAME_SHIFT_FIXED = 16
    CONSENSUS_EXTENDED = 32
    INSERTIONS_ADDED = 64


differ = Differ()


def parse_args(default_max_active):
    # noinspection PyTypeChecker
    parser = ArgumentParser(
        description='Compare sample results for testing a new release.',
        formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('--denovo',
                        action='store_true',
                        help='Compare old remapped results to new assembled results.')
    parser.add_argument('source_folder',
                        help='Main RAWDATA folder with results from previous version.')
    parser.add_argument('target_folder',
                        help='Testing RAWDATA folder to compare with.')
    parser.add_argument('--workers',
                        default=default_max_active,
                        type=int,
                        help='Number of parallel workers to process the samples.')

    verbosity_group = parser.add_mutually_exclusive_group()
    verbosity_group.add_argument('--verbose', action='store_true', help='Increase output verbosity.')
    verbosity_group.add_argument('--no-verbose', action='store_true', help='Normal output verbosity.', default=True)
    verbosity_group.add_argument('--debug', action='store_true', help='Maximum output verbosity.')
    verbosity_group.add_argument('--quiet', action='store_true', help='Minimize output verbosity.')

    return parser.parse_args()


def find_runs(source_folder, target_folder, use_denovo):
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
        try:
            source_versions.sort(key=parse_version)
        except ValueError as ex:
            message = f'Unexpected results file name in {run_name}.'
            raise ValueError(message) from ex
        source_path = os.path.join(source_results_path, source_versions[-1])

        if use_denovo:
            target_path = os.path.join(target_path, 'denovo')
            source_path = os.path.join(source_path, 'denovo')

        logger.debug("Comparing %r with %r.", source_path, target_path)
        yield MiseqRun(source_path, target_path, is_done)


def parse_version(version_name):
    version_text = version_name.split('_')[-1]
    if version_text.endswith('.zip'):
        version_text = version_text[:-4]
    version_text, possible_dash, possible_modifiers = version_text.partition("-")
    version_numbers = tuple(map(int, version_text.split('.')))
    return (version_numbers, possible_modifiers)


def report_source_versions(runs):
    version_runs = defaultdict(list)  # {version: [source_path]}
    for run in runs:
        version = os.path.basename(run.source_path)
        version_runs[version].append(run.source_path)
        yield run
    if version_runs:
        max_count = max(len(paths) for paths in version_runs.values())
    else:
        max_count = 0
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
            source_conseq = group_samples(os.path.join(run.source_path,
                                                       'conseq.csv'))
            target_conseq = group_samples(os.path.join(run.target_path,
                                                       'conseq.csv'))
            source_region_consensus = group_nucs(os.path.join(run.source_path,
                                                              'nuc.csv'))
            target_region_consensus = group_nucs(os.path.join(run.target_path,
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
                                     source_region_consensus.get(sample_name),
                                     source_counts.get(sample_name),
                                     source_conseq.get(sample_name)),
                         SampleFiles(target_cascades.get(sample_name),
                                     target_coverages.get(sample_name),
                                     target_g2ps.get(sample_name),
                                     target_region_consensus.get(sample_name),
                                     target_counts.get(sample_name),
                                     target_conseq.get(sample_name)))
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
        sample_conseqs = {}
        for seed, seed_rows in groupby(sample_rows, itemgetter('seed')):
            for region, region_rows in groupby(seed_rows, itemgetter('region')):
                consensus = [(choose_consensus(row), row)
                             for row in region_rows]
                sample_conseqs[(seed, region)] = consensus
        groups[sample] = sample_conseqs
    return groups


def choose_consensus(nuc_row: dict) -> str:
    coverage = int(nuc_row['coverage'])
    if coverage < 100:
        return 'x'
    nuc = SeedNucleotide()
    for nuc_seq in nuc.COUNTED_NUCS:
        source_nuc = 'del' if nuc_seq == '-' else nuc_seq
        nuc.count_nucleotides(nuc_seq, int(nuc_row[source_nuc]))
    consensus = nuc.get_consensus(MAX_CUTOFF)
    if int(nuc_row['ins']) > coverage / 2:
        consensus += 'i'
    return consensus


def get_run_name(sample):
    dirname = os.path.dirname(os.path.dirname(sample.run.target_path))
    run_name = os.path.basename(dirname)
    if run_name == 'Results':
        run_name = os.path.basename(os.path.dirname(dirname))
    return run_name


def compare_g2p(sample, diffs):
    source_fields = sample.source_files.g2p_summary
    target_fields = sample.target_files.g2p_summary
    if source_fields == target_fields:
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
    filtered_scores = {(score['project'],
                        score['region']): (score['on.score'],
                                           score.get('seed'),
                                           score.get('which.key.pos'))
                       for score in coverage_scores}
    if MICALL_VERSION == '7.15':
        filtered_scores = {
            (project, region): scores
            for (project, region), scores in filtered_scores.items()
            if project != 'HIVB' and region not in ('HIV1B-tat', 'SARS-CoV-2-nsp11')}

    return filtered_scores


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
        if MICALL_VERSION == '7.15':
            if key[0] == 'HIVGHA':
                # ignore HIVGHA project code
                continue
            elif target_seed is None:
                project, region = key
                hivgha_result = target_scores.get(('HIVGHA', region))
                if hivgha_result is not None:
                    # retrieve HIVGHA project code results for comparison
                    (target_score, target_seed, target_key_pos) = hivgha_result
                    if target_seed in ('HIV1-CRF02_AG-GH-AB286855-seed', 'HIV1-CRF06_CPX-GH-AB286851-seed'):
                        # ignore samples that switched to the new HIVGHA seeds
                        continue
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


def display_consensus(sequence):
    if sequence is None:
        return []
    return [sequence + '\n']


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


def has_big_prevalence_change(nuc: str, old_counts: dict, new_counts: dict):
    """ Check if a nucleotide's prevalence has changed significantly.
    
    True if it has changed by more than 5 percentage points.
    Ignore deletions and low coverage. Ties in the MAX consensus are always
    significant.
    """
    if nuc in 'x-':
        return False
    if nuc not in 'ACTG':
        return True
    if old_counts['coverage'] == '0' or new_counts['coverage'] == '0':
        return True
    old_prevalence = int(old_counts[nuc]) / int(old_counts['coverage'])
    new_prevalence = int(new_counts[nuc]) / int(new_counts['coverage'])
    return 0.05 < abs(new_prevalence - old_prevalence)


def compare_consensus(sample: Sample,
                      diffs,
                      scenarios_reported,
                      scenarios,
                      use_denovo=False):
    consensus_distances = []
    source_seqs = filter_consensus_sequences(sample.source_files, use_denovo)
    target_seqs = filter_consensus_sequences(sample.target_files, use_denovo)
    run_name = get_run_name(sample)
    source_counts = map_remap_counts(sample.source_files.remap_counts)
    target_counts = map_remap_counts(sample.target_files.remap_counts)
    keys = sorted(set(source_seqs.keys()) | target_seqs.keys())
    primer_trackers = {}
    for row in sample.target_files.conseq or ():
        if row['consensus-percent-cutoff'] == 'MAX':
            conseq = row['sequence']
            offset = int(row['offset'])
            conseq = 'x' * offset + conseq
            seed_name = row['region']
            primer_trackers[seed_name] = PrimerTracker(conseq, seed_name)
    cutoff = MAX_CUTOFF
    for key in keys:
        seed, region = key
        if MICALL_VERSION == '7.15':
            if region in ('HIV1B-tat',
                          'SARS-CoV-2-nsp11',
                          "SARS-CoV-2-TRS-B-1",
                          "SARS-CoV-2-TRS-B-2",
                          "SARS-CoV-2-TRS-B-3",
                          "SARS-CoV-2-TRS-B-4",
                          "SARS-CoV-2-TRS-B-5",
                          "SARS-CoV-2-TRS-B-6",
                          "SARS-CoV-2-TRS-B-7",
                          "SARS-CoV-2-TRS-B-8",
                          "SARS-CoV-2-TRS-B-9"):
                continue
            if seed in ('HIV1-CRF02_AG-GH-AB286855-seed', 'HIV1-CRF06_CPX-GH-AB286851-seed'):
                # ignore new HIVGHA seeds
                continue
        primer_tracker = primer_trackers.get(seed)
        source_details = source_seqs.get(key)
        target_details = target_seqs.get(key)
        if MICALL_VERSION == '7.15' and target_details is None:
            # ignore samples that map to the new HIVGHA seeds
            target_details = target_seqs.get(('HIV1-CRF02_AG-GH-AB286855-seed', region)) or \
                target_seqs.get(('HIV1-CRF06_CPX-GH-AB286851-seed', region))
            if target_details is not None:
                continue
        source_nucs = []
        target_nucs = []

        # Note: if either source or target region is missing, it might be because its coverage score is below 4.
        if source_details is None:
            has_big_change = True
            target_nucs = [nuc for nuc, row in target_details]
        elif target_details is None:
            has_big_change = True
            source_nucs = [nuc for nuc, row in source_details]
        else:
            has_big_change = False
            dummy_row = {'refseq.nuc.pos': '-1',
                         'coverage': '0'}
            if MICALL_VERSION == '7.15':
                if region.split('-')[-1] == 'NS5b':
                    NS5b_nucs = [nuc for nuc, row in target_details]
                    new_positions = [int(row['refseq.nuc.pos']) for nuc, row in target_details]
                    old_positions = [int(row['refseq.nuc.pos']) for nuc, row in source_details]
                    for position in old_positions:
                        if position not in new_positions:
                            target_details.insert(position-1, ('x', {'refseq.nuc.pos': position, 'coverage': '0'}))
                    last_codon = ''
                    for nuc in NS5b_nucs[-3:]:
                        last_codon = last_codon + nuc
                    try:
                        last_amino = translate(last_codon)
                    except KeyError:
                        last_amino = 'X'
                    if last_amino == '*':
                        target_details = target_details[:-3]
            for source_item, target_item in zip_longest(source_details,
                                                        target_details,
                                                        fillvalue=('', dummy_row)):
                source_consensus, source_row = source_item
                target_consensus, target_row = target_item
                if source_consensus != target_consensus:
                    if has_big_prevalence_change(source_consensus,
                                                 source_row,
                                                 target_row):
                        has_big_change = True
                    else:
                        source_consensus = source_consensus.lower()
                    if has_big_prevalence_change(target_consensus,
                                                 target_row,
                                                 source_row):
                        has_big_change = True
                    else:
                        target_consensus = target_consensus.lower()
                    if not has_big_change and 'x' in (source_consensus,
                                                      target_consensus):
                        old_coverage = int(source_row['coverage'])
                        new_coverage = int(target_row['coverage'])
                        has_big_change = (old_coverage < new_coverage/2 or
                                          new_coverage < old_coverage/2)

                source_nucs.append(source_consensus)
                target_nucs.append(target_consensus)

        source_seq = ''.join(source_nucs) or None
        target_seq = ''.join(target_nucs) or None
        consensus_distance = calculate_distance(region,
                                                cutoff,
                                                source_seq,
                                                target_seq)
        if False and consensus_distance.pct_diff > 5:
            print(sample.run.source_path, sample.name, consensus_distance)
            print(source_seq)
            print(target_seq)
        if consensus_distance is not None:
            consensus_distances.append(consensus_distance)
        if not has_big_change:
            continue
        is_main = is_consensus_interesting(region, cutoff)

        trimmed_source_seq = (source_seq and
                              source_seq.replace('-', '').replace('x', ''))
        trimmed_target_seq = (target_seq and
                              target_seq.replace('-', '').replace('x', ''))
        if (trimmed_source_seq == trimmed_target_seq and
                (scenarios_reported & Scenarios.CONSENSUS_DELETIONS_CHANGED)):
            scenarios[Scenarios.CONSENSUS_DELETIONS_CHANGED].append('.')
            continue
        if source_seq and target_seq:
            stripped_source_seq = source_seq.strip('x')
            stripped_target_seq = target_seq.strip('x')
            if stripped_source_seq in stripped_target_seq:
                scenarios[Scenarios.CONSENSUS_EXTENDED].append('.')
                continue
        if (use_denovo and
                is_main and
                seed == 'HIV1-B' and
                region == 'HIV1B-vpr' and
                (scenarios_reported & Scenarios.VPR_FRAME_SHIFT_FIXED)):
            if source_seq and source_seq[212] == '-':
                source_seq = source_seq[:212] + source_seq[213:273]
                target_seq = target_seq and target_seq[:272]
                if source_seq == target_seq:
                    scenarios[Scenarios.VPR_FRAME_SHIFT_FIXED].append('.')
                    continue
        if (source_counts != target_counts and
                is_main and
                scenarios_reported & Scenarios.REMAP_COUNTS_CHANGED):
            scenario = f'  {run_name}:{sample.name} consensus {region}\n'
            scenarios[Scenarios.REMAP_COUNTS_CHANGED].append(scenario)
        elif is_main and scenarios_reported & Scenarios.MAIN_CONSENSUS_CHANGED:
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
            diff = list(differ.compare(display_consensus(source_seq),
                                       display_consensus(target_seq)))
            diffs.extend(line.rstrip() for line in diff)
    return consensus_distances


def filter_consensus_sequences(files, use_denovo):
    region_consensus = files.region_consensus
    coverage_scores = files.coverage_scores
    if not (region_consensus and coverage_scores):
        return {}

    if MICALL_VERSION == '7.15':
        coverage_scores = [row
                           for row in coverage_scores
                           if row['project'] != 'HIVB']

    covered_regions = {(row.get('seed'), row.get('region'))
                       for row in coverage_scores
                       if row['on.score'] == '4'}
    return {adjust_seed(seed, region, use_denovo): consensus
            for (seed, region), consensus in region_consensus.items()
            if (seed, region) in covered_regions}


def adjust_seed(seed: str, region: str, use_denovo: bool):
    if use_denovo:
        if region == 'V3LOOP':
            seed = 'some-HIV-seed'
        elif seed.startswith('HIV'):
            seed = '-'.join(seed.split('-')[:2]) + '-?-seed'
    return seed, region


def compare_sample(sample,
                   scenarios_reported=Scenarios.NONE,
                   use_denovo=False):
    scenarios = defaultdict(list)
    diffs = []
    compare_g2p(sample, diffs)
    compare_coverage(sample, diffs, scenarios_reported, scenarios)
    consensus_distances = compare_consensus(sample,
                                            diffs,
                                            scenarios_reported,
                                            scenarios,
                                            use_denovo)
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
    available_memory = get_available_memory()
    recommended_memory = int((1 << 30) * 1.5)  # 1.5GB
    default_max_active = max(1, available_memory // recommended_memory)
    args = parse_args(default_max_active)

    if args.quiet:
        logger.setLevel(logging.ERROR)
    elif args.verbose:
        logger.setLevel(logging.INFO)
    elif args.debug:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.WARN)

    with ProcessPoolExecutor() as pool:
        runs = find_runs(args.source_folder, args.target_folder, args.denovo)
        runs = report_source_versions(runs)
        samples = read_samples(runs)
        # noinspection PyTypeChecker
        scenarios_reported = (Scenarios.OTHER_CONSENSUS_CHANGED |
                              Scenarios.CONSENSUS_DELETIONS_CHANGED |
                              Scenarios.VPR_FRAME_SHIFT_FIXED |
                              Scenarios.CONSENSUS_EXTENDED |
                              Scenarios.REMAP_COUNTS_CHANGED)
        try:
            results = pool.map(partial(compare_sample,
                                       scenarios_reported=scenarios_reported,
                                       use_denovo=args.denovo),
                               samples,
                               chunksize=args.workers)
        except BrokenProcessPool:
            print("Broken Process Pool - probably the memory usage is too high. Try again with fewer workers!")
            print("Current number of workers: {args.workers}")
            raise
        scenario_summaries = defaultdict(list)
        i = 0
        all_consensus_distances = []
        report_limit = 1000
        report_count = 0
        for i, (report, scenarios, consensus_distances) in enumerate(results, 1):
            if report:
                report_count += 1
                if report_count > report_limit:
                    print(f'More than {report_limit} reports, skipping the rest.')
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

    report_distances(all_consensus_distances)
    print('Finished {} samples.'.format(i))


def report_distances(all_consensus_distances):
    if not all_consensus_distances:
        return
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


if __name__ == '__main__':
    main()
