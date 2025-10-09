""" Set up raw data for the release test.

Copies a set of run folders and selected FASTQ files from the main raw data
folder to your developer workstation. Configure the set of samples to copy
in test_samples.csv.
You'll need to upload the QC data from the InterOp folders before you can
start MISEQ_MONITOR.py.
"""

from argparse import ArgumentParser
from csv import DictReader
from datetime import datetime
import errno
from glob import glob
from itertools import groupby
from operator import attrgetter
import os
from pathlib import Path
from shutil import copy, rmtree, copytree

from micall.utils.sample_sheet_parser import read_sample_sheet_and_overrides
from micall.utils.list_fastq_files import list_fastq_files, find_fastq_source_folder

NEEDS_PROCESSING = 'needsprocessing'
ERROR_PROCESSING = 'errorprocessing'


def parse_args():
    parser = ArgumentParser(description='Prepare sample runs for testing a new release.')
    parser.add_argument(
        'old_folder',
        help='Main RAWDATA folder that will be copied to local RAWDATA.')
    parser.add_argument(
        'test_folder',
        type=Path,
        default=os.environ.get('MICALL_RAW_DATA',
                               Path.home() / "data/RAW_DATA"),
        help='Local RAWDATA folder that will be used to run tests.')
    parser.add_argument(
        '-s',
        '--samples_csv',
        default='test_samples.csv',
        help='CSV file with sample and run names.')
    parser.add_argument(
        '-m',
        '--min_run_name',
        help='Select all later runs, (e.g., "161201"). Overrides samples_csv.')
    parser.add_argument(
        '--pipeline_version',
        default='0-dev',
        help='version suffix for folder names')
    parser.add_argument(
        '-n',
        '--no_links',
        action='store_true',
        help='Copy data files instead of using symlinks.')
    return parser.parse_args()


class Sample(object):
    def __init__(self, run_name, extract_num, config):
        self.run_name = run_name
        self.extract_num = extract_num
        self.config = config
        self.fastq_paths = None
        self.source_fastq_folder = None  # Track where FASTQ files were found

    def find(self, source_folder, qai_run_names=None):
        """ Find matching samples in the source folder.

        Puts all the sample names in self.fastq_paths.
        :param str source_folder: the folder to search for samples that match
            self.extract_num
        :param set qai_run_names: a set to add the new run name to as it would
            be formatted on QAI.
        """
        run_path = os.path.join(source_folder, 'MiSeq', 'runs', self.run_name)
        if os.path.exists(run_path):
            self.run_name = run_path
        else:
            run_parts = str(self.run_name).split('.')
            if len(run_parts) < 2:
                run_parts.append('M01841')
            format_string = '%d-%b-%y' if len(run_parts[0]) == 9 else '%d-%b-%Y'
            run_date = datetime.strptime(run_parts[0], format_string)
            base_run_name = run_date.strftime('%y%m%d') + '_' + run_parts[1]
            pattern = os.path.join(source_folder,
                                   'MiSeq',
                                   'runs',
                                   base_run_name+'*',
                                   NEEDS_PROCESSING)
            matches = glob(pattern)
            if len(matches) != 1:
                raise RuntimeError('Expected one match for {}, but found: {}'.format(
                    pattern,
                    matches))
            self.run_name = os.path.dirname(matches[0])

        # Use list_fastq_files to search in BaseCalls, Alignment_*/*/Fastq, or run path
        pattern = self.extract_num + '*_R1_*'
        fastq_files = list_fastq_files(self.run_name, pattern, fallback_to_run_path=False)
        if not fastq_files:
            raise RuntimeError('No matches found for run {}'.format(self.run_name))
        matches = sorted([str(f) for f in fastq_files])
        self.fastq_paths = matches
        
        # Store the source folder where FASTQ files were found
        self.source_fastq_folder = find_fastq_source_folder(self.run_name, pattern)
        
        if qai_run_names is not None:
            sample_sheet_path = Path(self.run_name) / 'SampleSheet.csv'
            try:
                sample_sheet = read_sample_sheet_and_overrides(sample_sheet_path)
            except ValueError:
                print(f'Bad sample sheet for {self.run_name}.')
            else:
                qai_run_name = sample_sheet['Project Name']
                qai_run_names.add(qai_run_name)

    def setup_samples(self):
        base_run_name = os.path.basename(self.run_name)
        
        # Compute the relative path from run folder to source FASTQ folder
        # This could be 'Data/Intensities/BaseCalls' or 'Alignment_20/L001/Fastq' etc.
        assert self.source_fastq_folder is not None, "source_fastq_folder is not set"
        source_relative_path = self.source_fastq_folder.relative_to(self.run_name)

        # Replicate the same directory structure in the target
        target_data_path = os.path.join(self.config.test_folder,
                                        'MiSeq',
                                        'runs',
                                        base_run_name,
                                        str(source_relative_path))
        suspended_path = os.path.join(target_data_path, 'suspended')
        sample_paths = []
        assert self.fastq_paths is not None, "fastq_paths is not set"
        for fastq1_path in self.fastq_paths:
            fastq2_path = fastq1_path.replace('_R1_', '_R2_')
            for fastq_path in (fastq1_path, fastq2_path):
                base_name = os.path.basename(fastq_path)
                target_path = os.path.join(target_data_path, base_name)
                sample_paths.append(target_path)
                if not os.path.exists(target_path):
                    suspended_fastq_path = os.path.join(suspended_path,
                                                        base_name)
                    if os.path.exists(suspended_fastq_path):
                        os.rename(suspended_fastq_path,
                                  target_path)
                    elif self.config.no_links:
                        copy(fastq_path, target_path)
                    else:
                        os.symlink(fastq_path, target_path)
        return sample_paths

    def setup_run(self):
        base_run_name = os.path.basename(self.run_name)
        target_run_path = os.path.join(self.config.test_folder,
                                       'MiSeq',
                                       'runs',
                                       base_run_name)
        suspended_run_path = os.path.join(self.config.test_folder,
                                          'MiSeq',
                                          'runs',
                                          'suspended',
                                          base_run_name)
        if os.path.exists(target_run_path):
            pass
        elif os.path.exists(suspended_run_path):
            os.rename(suspended_run_path, target_run_path)
        else:
            # Create directory structure matching the source folder
            assert self.source_fastq_folder is not None, "source_fastq_folder is not set"
            source_relative_path = self.source_fastq_folder.relative_to(self.run_name)
            fastq_folder_path = os.path.join(target_run_path, str(source_relative_path))
            os.makedirs(fastq_folder_path)

            interop_source = os.path.join(self.run_name, 'InterOp')
            interop_target = os.path.join(target_run_path, 'InterOp')
            if self.config.no_links:
                copytree(interop_source, interop_target)
            else:
                os.symlink(interop_source, interop_target)
            for filename in ('RunInfo.xml',
                             'SampleSheet.csv',
                             'SampleSheetOverrides.csv',
                             'needsprocessing'):

                source = os.path.join(self.run_name, filename)
                target = os.path.join(target_run_path, filename)
                if os.path.exists(source):
                    copy(source, target)

        results_path = os.path.join(target_run_path,
                                    'Results',
                                    'version_' + self.config.pipeline_version)
        try:
            rmtree(results_path)
        except OSError as ex:
            if ex.errno != errno.ENOENT:
                raise
            os.makedirs(os.path.dirname(results_path), exist_ok=True)
        try:
            os.remove(os.path.join(target_run_path,
                                   ERROR_PROCESSING))
        except OSError as ex:
            if ex.errno != errno.ENOENT:
                raise
        return target_run_path


def suspend_inactive_runs(active_runs, rawdata_mount):
    active_run_names = set(map(os.path.basename, active_runs))
    local_runs = glob(os.path.join(rawdata_mount, 'MiSeq', 'runs', '*'))
    for run_path in local_runs:
        base_run_name = os.path.basename(run_path)
        if base_run_name not in active_run_names and base_run_name != 'suspended':
            suspended_folder_path = os.path.join(rawdata_mount,
                                                 'MiSeq',
                                                 'runs',
                                                 'suspended')
            try:
                os.mkdir(suspended_folder_path)
            except FileExistsError:
                pass
            suspended_path = os.path.join(suspended_folder_path,
                                          base_run_name)
            os.rename(run_path, suspended_path)


def suspend_inactive_samples(run_path, active_sample_paths):
    # Find the actual FASTQ folder (could be BaseCalls or Alignment_*/*/Fastq)
    # Use a generic pattern to find any FASTQ files
    fastq_folder = find_fastq_source_folder(run_path, '*_R1_*')
    if fastq_folder is None:
        # No FASTQ files found, nothing to suspend
        return
        
    local_samples = glob(os.path.join(str(fastq_folder), '*'))
    excluded = list(active_sample_paths)
    suspended_path = os.path.join(str(fastq_folder), 'suspended')
    excluded.append(suspended_path)
    for sample_path in local_samples:
        if sample_path not in excluded:
            if not os.path.exists(suspended_path):
                os.makedirs(suspended_path)
            os.rename(sample_path,
                      os.path.join(suspended_path, os.path.basename(sample_path)))
    try:
        os.rmdir(suspended_path)
    except OSError as ex:
        if ex.errno not in (errno.ENOTEMPTY, errno.ENOENT):
            raise


def find_run_folders(source_folder, min_run_name):
    pattern = os.path.join(source_folder, 'MiSeq', 'runs', '*')
    matches = glob(pattern)
    folder_names = [os.path.basename(match) for match in matches]
    recent_folders = [folder
                      for folder in folder_names
                      if folder >= min_run_name and folder[0].isnumeric()]
    folders_with_data = [folder
                         for folder in recent_folders
                         if os.path.exists(os.path.join(source_folder,
                                                        'MiSeq',
                                                        'runs',
                                                        folder,
                                                        NEEDS_PROCESSING))]
    return folders_with_data


def main():
    args = parse_args()
    if args.min_run_name is not None:
        print('Setting up',
              args.pipeline_version,
              'for runs after',
              args.min_run_name)
        run_folders = find_run_folders(args.old_folder, args.min_run_name)
        samples = [Sample(folder, '*', args) for folder in run_folders]
    else:
        print('Setting up',
              args.pipeline_version,
              'for runs in',
              args.samples_csv)
        with open(args.samples_csv, 'r') as samples_csv:
            samples = [Sample(row['run'], row['enum'], args)
                       for row in DictReader(samples_csv)]
    qai_run_names = set()
    for sample in samples:
        sample.find(args.old_folder, qai_run_names)
    for qai_run_name in sorted(qai_run_names):
        print('LabMiseqRun.import("{}")'.format(qai_run_name))
    samples.sort(key=lambda s: (s.run_name, s.extract_num))
    runs = []
    for run, run_samples in groupby(samples, key=attrgetter('run_name')):
        runs.append(run)
        run_name = os.path.basename(run)
        print('--' + run_name)
        sample_paths = []
        run_path = None
        for sample in run_samples:
            run_path = sample.setup_run()
            sample_paths.extend(sample.setup_samples())
            print('  ' + ', '.join(map(os.path.basename, sample.fastq_paths)))
        suspend_inactive_samples(run_path, active_sample_paths=sample_paths)
    suspend_inactive_runs(active_runs=runs, rawdata_mount=args.test_folder)


if __name__ == '__main__':
    main()
