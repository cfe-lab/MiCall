""" Copy result files to shared folder so they can be compared.

They should also be processed by the report scripts.
"""
import os
import shutil
from argparse import ArgumentParser
from glob import glob

import errno

from micall.settings import rawdata_mount, pipeline_version, DONE_PROCESSING


def parse_args():
    parser = ArgumentParser(description='Publish sample results for testing a new release.')
    parser.add_argument('target_folder',
                        help='Testing RAWDATA folder to receive copies from local RAWDATA.')
    return parser.parse_args()


def main():
    args = parse_args()
    run_paths = glob(os.path.join(rawdata_mount, 'MiSeq', 'runs', '*'))
    run_paths.sort()
    for run_path in run_paths:
        run_name = os.path.basename(run_path)
        if run_name == 'suspended':
            continue
        results_path = os.path.join(run_path,
                                    'Results',
                                    'version_' + pipeline_version)
        done_path = os.path.join(results_path, DONE_PROCESSING)
        target_path = os.path.join(args.target_folder,
                                   run_name,
                                   'Results',
                                   'version_' + pipeline_version)
        try:
            shutil.rmtree(target_path)
        except OSError as ex:
            if ex.errno != errno.ENOENT:
                raise
        if not os.path.exists(done_path):
            print('Not done: ' + run_name)
            continue
        shutil.copytree(results_path, target_path)
        print('Done: ' + run_name)
    print('Done.')


if __name__ == '__main__':
    main()
