""" Copy result files to shared folder so they can be compared.

They should also be processed by the report scripts.
"""
from __future__ import print_function
from argparse import ArgumentParser
import csv
from difflib import Differ
from glob import glob
from operator import itemgetter
import os


from micall.settings import rawdata_mount, pipeline_version, DONE_PROCESSING


def parse_args():
    parser = ArgumentParser(description='Compare sample results for testing a new release.')
    parser.add_argument('source_folder',
                        help='Main RAWDATA folder with results from previous version.')
    parser.add_argument('target_folder',
                        help='Testing RAWDATA folder to compare with.')
    return parser.parse_args()


def select_columns(csv_path,
                   column_names,
                   all_sample_names=None,
                   filter_sample_names=None):
    """ Return lines from a CSV file with only some of the columns.

    :param csv_path: path to the CSV file
    :param column_names: columns to include
    :param all_sample_names: a set to receive all values for the 'sample'
        column, or None
    :param filter_sample_names: a set of sample names to include in the
        output, or None to not filter
    """
    lines = []
    with open(csv_path, 'rU') as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row.get('on.score', None) == '1':
                continue
            if all_sample_names is not None:
                all_sample_names.add(row['sample'])
            if filter_sample_names is None or row['sample'] in filter_sample_names:
                lines.append(',   '.join(itemgetter(*column_names)(row)) + '\n')
    return lines


def compare_coverage_scores(source_path, target_path):
    columns = ['sample',
               'project',
               'region',
               'seed',
               'on.score']
    filename = 'coverage_scores.csv'
    sample_names = set()
    target_columns = select_columns(os.path.join(target_path, filename),
                                    columns,
                                    all_sample_names=sample_names)
    source_columns = select_columns(os.path.join(source_path, filename),
                                    columns,
                                    filter_sample_names=sample_names)
    differ = Differ()
    for i, diff in enumerate(differ.compare(source_columns, target_columns)):
        if i == 0:
            print('Coverage score changes in {}:'.format(target_path))
        print(diff, end='')


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
        compare_coverage_scores(source_path, target_path)
        print('Done: ' + run_name)
    print('Done.')


if __name__ == '__main__':
    main()
