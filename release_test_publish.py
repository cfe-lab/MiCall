""" Copy result files to shared folder so they can be compared.

They should also be processed by the report scripts.
"""
import os
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from glob import glob

from os import makedirs
from pathlib import Path
from subprocess import run


def parse_args():
    # noinspection PyTypeChecker
    parser = ArgumentParser(
        description='Publish sample results for testing a new release.',
        formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        '--test_folder',
        type=Path,
        default=os.environ.get('MICALL_RAW_DATA',
                               Path.home() / "data/RAW_DATA"),
        help='Local RAWDATA folder that was used to run tests.')
    parser.add_argument('publish_folder',
                        help='Testing RAWDATA folder to receive copies from local RAWDATA.')
    parser.add_argument(
        '--pipeline_version',
        default='0-dev',
        help='version suffix for folder names')
    parser.add_argument(
        '--delete',
        action='store_true',
        help='delete old files from publish_folder')
    return parser.parse_args()


def main():
    args = parse_args()
    run_paths = glob(os.path.join(args.test_folder, 'MiSeq', 'runs', '*'))
    run_paths.sort()
    for run_path in run_paths:
        run_name = os.path.basename(run_path)
        if run_name == 'suspended':
            continue
        results_path = os.path.join(run_path,
                                    'Results',
                                    'version_' + args.pipeline_version)
        done_path = os.path.join(results_path, 'done_all_processing')
        target_path = os.path.join(args.publish_folder,
                                   'MiSeq',
                                   'runs',
                                   run_name,
                                   'Results')
        if not os.path.exists(done_path):
            print('Not done: ' + run_name)
            continue
        makedirs(target_path, exist_ok=True)
        rsync_args = ['rsync', '-a', results_path, target_path]
        if args.delete:
            rsync_args.append('--delete')
        run(rsync_args)
        print('Done: ' + run_name)
    print('Done.')


if __name__ == '__main__':
    main()
