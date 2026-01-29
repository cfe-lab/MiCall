import logging
import os
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from csv import DictReader, DictWriter
from itertools import groupby
from operator import itemgetter
from pathlib import Path
from random import shuffle
from statistics import median
import sys

import matplotlib.pyplot as plt

logger = logging.getLogger(__name__)


def parse_args():
    default_runs_path = Path('~/data/RAW_DATA/MiSeq/runs').expanduser()
    parser = ArgumentParser(
        description='Report median coverage for all samples.',
        formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('runs',
                        nargs='?',
                        type=Path,
                        default=default_runs_path,
                        help='Runs folder to search for samples.')
    parser.add_argument('report',
                        nargs='?',
                        default='-',
                        help='File to write coverage report or read coverage from.')
    return parser.parse_args()


def load_coverage(report_file, runs_folder):
    writer = DictWriter(report_file,
                        ['sample', 'reads', 'coverage', 'is_good'],
                        lineterminator=os.linesep)
    writer.writeheader()
    result_folders = sorted(runs_folder.glob('*/Results'), reverse=True)
    version_name = None
    for result_folder in result_folders:
        if version_name is None:
            version_names = sorted(result_folder.iterdir())
            version_name = version_names[-1].name
            logging.info('Loading version %s.', version_name)
        cascade_path = result_folder / version_name / 'cascade.csv'
        amino_path = result_folder / version_name / 'amino.csv'
        coverage_path = result_folder / version_name / 'coverage_scores.csv'
        missing_paths = [path.name
                         for path in (amino_path, cascade_path, coverage_path)
                         if not path.exists()]
        if missing_paths:
            logging.warning('Folder %s is missing files: %s.',
                            result_folder.parent.name,
                            missing_paths)
            continue
        logging.debug('Loading folder %s.', result_folder.parent.name)
        with coverage_path.open() as coverage_file:
            reader = DictReader(coverage_file)
            good_coverages = {row['sample']
                              for row in reader
                              if row['on.score'] == '4'}
        with cascade_path.open() as cascade_file:
            reader = DictReader(cascade_file)
            read_counts = {row['sample']: int(row['demultiplexed'])
                           for row in reader}
        unseen_samples = set(read_counts)
        with amino_path.open() as amino_file:
            reader = DictReader(amino_file)
            for sample, sample_rows in groupby(reader, itemgetter('sample')):
                is_good = sample in good_coverages
                read_count = read_counts[sample]
                unseen_samples.discard(sample)
                coverage = median(int(row['coverage'])
                                  for row in sample_rows
                                  if row['coverage'] != '0')
                writer.writerow(dict(sample=sample,
                                     reads=read_count,
                                     coverage=coverage,
                                     is_good=is_good))
        for sample in sorted(unseen_samples):
            writer.writerow(dict(sample=sample,
                                 reads=read_counts[sample],
                                 coverage=0,
                                 is_good=False))


def display_coverage(report_file):
    # TODO: Try KDE plots.
    expected_widths = {'HCV': 9600,
                       'MidHCV': 1200,
                       'NS5a': 1200,
                       'V3LOOP': 312}
    group_size = 1000
    alpha = 100/group_size
    reader = DictReader(report_file)
    logger.info('Displaying.')
    rows = list(reader)
    for row in rows:
        row['project'] = row['sample'].split('_')[0].split('-')[-1]
        row['reads'] = int(row['reads'])
        row['coverage'] = float(row['coverage'])
        expected_width = expected_widths.get(row['project'], 1000000)
        row['expected_depth'] = row['reads'] / expected_width
    rows = [row
            for row in rows
            if row['is_good'] == 'True' and
            row['project'] in ('HCV', 'MidHCV', 'V3LOOP')]
    shuffle(rows)
    rows.sort(key=itemgetter('project'), reverse=True)

    plt.subplot(211)
    for project, group_rows in groupby(rows, itemgetter('project')):
        group_rows = list(group_rows)[:group_size]

        read_counts = [row['reads'] for row in group_rows]
        coverage = [row['coverage'] for row in group_rows]
        plt.plot(read_counts, coverage, 'o', label=project, alpha=alpha)
    plt.xlim([-50_000, 1_500_000])
    plt.ylim([-10_000, 100_000])
    plt.title('Coverage vs. Read Count')
    plt.xlabel('read count')
    plt.ylabel('median coverage')
    add_legend()
    plt.subplot(212)
    for project, group_rows in groupby(rows, itemgetter('project')):
        group_rows = list(group_rows)[:group_size]
        expected_depth = [row['expected_depth'] for row in group_rows]
        coverage = [row['coverage'] for row in group_rows]
        plt.plot(expected_depth, coverage, 'o', label=project, alpha=alpha)
    plt.xlim([-10, 300])
    plt.ylim([-10_000, 100_000])
    plt.title('Coverage vs. Expected Depth')
    plt.xlabel('read count / expected width')
    plt.ylabel('median coverage')
    add_legend()
    plt.tight_layout()
    plt.show()


def add_legend():
    legend = plt.legend(ncol=2)
    for handle in legend.legendHandles:
        # noinspection PyProtectedMember
        handle._legmarker.set_alpha(1.0)


def main():
    logging.basicConfig(level=logging.DEBUG)
    args = parse_args()
    runs_folder: Path = args.runs
    report_path = args.report
    if not os.path.exists(report_path):
        if report_path == '-':
            report_file = sys.stdout
        else:
            report_file = open(report_path, 'w')
        with report_file:
            load_coverage(report_file, runs_folder)
    if os.path.exists(report_path):
        with open(report_path, 'r') as report_file:
            display_coverage(report_file)


if __name__ == '__main__':
    main()
