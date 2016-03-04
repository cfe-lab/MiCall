from collections import namedtuple, defaultdict, Counter
import csv
from datetime import datetime, timedelta
from glob import glob
import hashlib
import itertools
import logging
from operator import itemgetter
import os
import shutil
import subprocess
import sys
from xml.etree import ElementTree

from kiveapi.errors import KiveRunFailedException
from micall import settings
from micall.monitor import qai_helper, update_qai
from micall.monitor.kive_download import kive_login, download_results

MAX_RUN_NAME_LENGTH = 60
logger = logging.getLogger("kive_loader")


class KiveLoader(object):
    RunInfo = namedtuple('RunInfo', 'miseq_run_id read_lengths indexes_length')

    def __init__(self,
                 launch_limit=sys.maxint,
                 status_delay=30,
                 folder_delay=3600,
                 retry_delay=3600):
        """ Prepare an instance.

        @param launch_limit: the maximum number active runs
        @param status_delay: seconds delay between checking status of all runs
        @param folder_delay: seconds delay between scanning for new folders
        @param retry_delay: seconds to continue retrying after an error
        """
        self.launch_limit = launch_limit
        self.status_delay = status_delay
        self.folder_delay = folder_delay
        self.retry_delay = retry_delay
        self.folders = None
        self.latest_folder = None
        self.kive = None
        self.preexisting_runs = None
        self.active_runs = []
        self.batches = defaultdict(list)  # {folder: [run]}
        self.retry_counts = Counter()  # {folder: count}
        self.is_status_available = True

    def poll(self):
        try:
            self.downloading_folder = None
            self.check_folders()
            if self.file_count >= len(self.files):
                if self.folder_count >= len(self.folders):
                    self.folder = None
                else:
                    self.folder = self.folders[self.folder_count]
                    self.folder_count += 1
                    self.files = self.find_files(self.folder)
                    self.file_count = 0
                    self.trimmed_folder = self.trim_folder(self.folder)
                    quality_file = self.download_quality(self.folder)
                    self.check_kive_connection()
                    self.quality_dataset = self.prepare_kive_dataset(
                        quality_file,
                        'phiX174 quality from MiSeq run ' + self.trimmed_folder,
                        self.quality_cdt)
            if not self.can_launch():
                self.check_run_status()
            if not self.can_launch():
                return self.status_delay
            file1 = self.files[self.file_count]
            file2 = file1.replace('_R1_', '_R2_')
            self.file_count += 1
            description = 'read from MiSeq run ' + self.trimmed_folder
            dataset1 = self.prepare_kive_dataset(file1,
                                                 'forward ' + description,
                                                 None)
            dataset2 = self.prepare_kive_dataset(file2,
                                                 'reverse ' + description,
                                                 None)
            if self.preexisting_runs is None:
                self.preexisting_runs = self.find_preexisting_runs()
            run_key = self.get_run_key(self.quality_dataset, dataset1, dataset2)
            run = self.preexisting_runs.pop(run_key, None)
            if run is None:
                run = self.launch_run(self.quality_dataset, dataset1, dataset2)
            self.active_runs.append(run)
            self.batches[self.folder].append(run)
            return 0
        except StandardError as ex:
            failed_folder = self.downloading_folder or self.folder
            is_reset_needed = True
            delay_fractions = [1.0/60, 5.0/60]
            retry_count = self.retry_counts[failed_folder]
            self.retry_counts[failed_folder] += 1
            if retry_count > len(delay_fractions):
                del self.retry_counts[failed_folder]
                message = str(ex)
                self.mark_folder_disabled(failed_folder,
                                          message,
                                          exc_info=True)
                delay = 0
            else:
                self.log_retry(failed_folder)
                if retry_count < len(delay_fractions):
                    delay = self.retry_delay * delay_fractions[retry_count]
                else:
                    delay = self.retry_delay * (1-sum(delay_fractions))
                is_reset_needed = self.downloading_folder is None
            if is_reset_needed:
                self.folders = None
                self.reset_folders()
            return delay

    def check_folders(self):
        now = self.get_time()
        if self.folders is None or now >= self.folder_scan_time:
            self.folder_scan_time = now + timedelta(seconds=self.folder_delay)
            new_folders = self.find_folders()
            if self.folders is None or set(new_folders).difference(self.folders):
                # First time or we found a new folder
                self.folders = new_folders
                self.reset_folders()
            elif not self.active_runs:
                logger.info('No folders need processing.')

    def reset_folders(self):
        self.folder = None
        self.folder_count = 0
        self.files = []
        self.file_count = 0
        self.active_runs = []
        self.preexisting_runs = None
        self.batches.clear()

    def can_launch(self):
        if self.folder is None:
            return False
        return (self.folder == self.latest_folder or
                len(self.active_runs) < self.launch_limit)

    def check_run_status(self):
        for i in reversed(range(len(self.active_runs))):
            run = self.active_runs[i]
            try:
                is_complete = self.is_run_complete(run)
                self.is_status_available = True
            except KiveRunFailedException:
                is_complete = True
            except:
                is_complete = False
                if self.is_status_available:
                    # First failure, so log it.
                    logger.warn('Unable to check run status.', exc_info=True)
                    self.is_status_available = False
            if is_complete:
                del self.active_runs[i]
        for folder, runs in self.batches.items():
            if folder != self.folder and all(run not in self.active_runs
                                             for run in runs):
                self.downloading_folder = folder
                self.download_results(folder, runs)
                del self.batches[folder]
                self.downloading_folder = None

    def find_folders(self):
        """ Find run folders ready to be processed.

        Also set self.latest_folder to be the latest folder that is marked
        with a needsprocessing file, even if it's done processing or failed.

        @return: a list of paths to the folders in the order they should be
            processed
        """
        # flag indicates that Illumina MiseqReporter has completed pre-processing, files available on NAS
        flag_files = glob(settings.rawdata_mount +
                          'MiSeq/runs/*/{}'.format(settings.NEEDS_PROCESSING))
        flag_files.sort(reverse=True)
        folders = []
        for i, flag_file in enumerate(flag_files):
            folder = os.path.dirname(flag_file)
            result_path = os.path.join(folder,
                                       'Results',
                                       'version_{}'.format(settings.pipeline_version))
            done_path = os.path.join(result_path, settings.DONE_PROCESSING)
            if i == 0:
                self.latest_folder = folder
            if self.is_marked_as_disabled(folder):
                continue
            if not self.is_quality_control_uploaded(folder):
                continue
            # if doneprocessing file already exists, then do not re-process
            if settings.production:
                if os.path.exists(done_path):
                    continue
                if os.path.exists(result_path):
                    # Not done, but results folder exists. Assume it's incomplete,
                    # delete it, and rerun.
                    shutil.rmtree(result_path)
                # running in development mode - do all runs even if already processed
                # note that results will not be uploaded
            folders.append(folder)

        return folders

    def is_marked_as_disabled(self, folder):
        return os.path.exists(os.path.join(folder, settings.ERROR_PROCESSING))

    def is_quality_control_uploaded(self, folder):
        return os.path.exists(os.path.join(folder, settings.QC_UPLOADED))

    def find_files(self, folder):
        """ Find FASTQ files within a folder.

        @return: a list of paths to the files within the folder.
        """
        gz_files = []
        for filepath in glob(os.path.join(
                folder,
                'Data/Intensities/BaseCalls/*_R1_001.fastq.gz')):
            # Report number of reads failing to demultiplex to the log
            if not os.path.basename(filepath).startswith('Undetermined'):
                gz_files.append(filepath)
            else:
                # Do word count directly on stream redirected from gunzip
                p1 = subprocess.Popen(['gunzip', '-c', filepath], stdout=subprocess.PIPE)
                p2 = subprocess.Popen(['wc', '-l'], stdin=p1.stdout, stdout=subprocess.PIPE)
                output = p2.communicate()[0]
                failed_demultiplexing = output.strip(' \n')
                logger.info("%s reads failed to demultiplex in %s (skipping file)",
                            failed_demultiplexing,
                            filepath)

        return sorted(gz_files)

    def prepare_kive_dataset(self, filename, description, cdt):
        """ Upload a dataset to Kive, if it's not already.

        @return: the dataset object from the Kive API wrapper
        """
        dataset = self.find_kive_dataset(filename, cdt)
        return dataset or self.upload_kive_dataset(filename, description, cdt)

    def find_kive_dataset(self, filename, cdt):
        """ Search for a dataset in Kive by name and checksum.

        @return: the dataset object from the Kive API wrapper, or None
        """
        self.check_kive_connection()
        dataset_name = os.path.basename(filename)
        CHUNK_SIZE = 4096
        hash = hashlib.md5()
        with open(filename, 'rb') as f:
            for chunk in iter(lambda: f.read(CHUNK_SIZE), b""):
                hash.update(chunk)
            checksum = hash.hexdigest()
            datasets = self.kive.find_datasets(name=dataset_name,
                                               md5=checksum,
                                               cdt=cdt)
            needed_groups = set(settings.kive_groups_allowed)
            for dataset in datasets:
                missing_groups = needed_groups - set(dataset.groups_allowed)
                if not missing_groups:
                    logger.info('dataset already in Kive: %r', dataset_name)
                    return dataset
        return None

    def check_kive_connection(self):
        if self.kive is None:
            self.kive = kive_login(settings.kive_server_url,
                                   settings.kive_user,
                                   settings.kive_password)
            # retrieve Pipeline object based on version
            self.pipeline = self.kive.get_pipeline(
                settings.pipeline_version_kive_id)

            # retrieve quality.csv compound data type
            self.quality_cdt = self.kive.get_cdt(settings.quality_cdt_kive_id)

    def upload_kive_dataset(self, filename, description, cdt):
        """ Upload a dataset to Kive.

        @return: the dataset object from the Kive API wrapper
        """
        self.check_kive_connection()
        dataset_name = os.path.basename(filename)
        logger.info('uploading dataset %r', dataset_name)
        with open(filename, 'rb') as f:
            dataset = self.kive.add_dataset(name=dataset_name,
                                            description=description,
                                            handle=f,
                                            cdt=cdt,
                                            groups=settings.kive_groups_allowed)
        return dataset

    def trim_folder(self, folder):
        """ Trim folder path to minimum unique name.

        Folder name is /some/path/YYMMDD_MACHINEID_SEQ_RANDOMNUMBER.
        Date and machine id should be unique together.
        """
        run_name = os.path.basename(folder)
        return '_'.join(run_name.split('_')[:2])

    def parse_run_info(self, run_info_path):
        runInfoTree = ElementTree.parse(run_info_path)
        runInfoRoot = runInfoTree.getroot()
        run = runInfoRoot[0]
        miseq_run_id = run.attrib['Id']
        indexes_length = 0
        read_lengths = []
        for read in run.iter('Read'):
            num_cycles = int(read.attrib['NumCycles'])
            if read.attrib['IsIndexedRead'] == 'Y':
                indexes_length += num_cycles
            else:
                read_lengths.append(num_cycles)
        return KiveLoader.RunInfo(miseq_run_id=miseq_run_id,
                                  read_lengths=read_lengths,
                                  indexes_length=indexes_length)

    def download_quality(self, folder):
        """ Download quality control data for the run.

        @return path for the quality CSV file
        """
        trimmed_folder = self.trim_folder(folder)
        destination_folder = os.path.join(settings.home,
                                          os.path.basename(folder))
        destination = os.path.join(destination_folder,
                                   '{}_quality.csv'.format(trimmed_folder))
        if not os.path.exists(destination_folder):
            os.makedirs(destination_folder)

        run_info_path = os.path.join(folder, 'RunInfo.xml')
        run_info = self.parse_run_info(run_info_path)
        direction_params = [(run_info.read_lengths[0], 1),
                            (run_info.read_lengths[1], -1)]
        with qai_helper.Session() as session:
            session.login(settings.qai_path, settings.qai_user, settings.qai_password)
            metrics = session.get_json(
                '/miseqqc_errormetrics?runid=' + run_info.miseq_run_id)
            if not metrics:
                raise StandardError(
                    'No quality control metrics found for run ' + run_info.miseq_run_id)

        with open(destination, 'w') as f:
            writer = csv.DictWriter(f,
                                    ['tile', 'cycle', 'errorrate'],
                                    lineterminator=os.linesep)
            writer.writeheader()
            for tile, tile_metrics in itertools.groupby(metrics,
                                                        itemgetter('tile')):
                expected_cycle = 0
                metric = next(tile_metrics)
                for cycles, sign in direction_params:
                    report_cycle = sign
                    for _ in range(1, cycles+1):
                        expected_cycle += 1
                        cycle = int(metric['cycle'])
                        if cycle == expected_cycle:
                            errorrate = metric['errorrate']
                            try:
                                metric = next(tile_metrics)
                            except StopIteration:
                                metric = dict(cycle=-1)
                        else:
                            errorrate = ''
                        writer.writerow(dict(tile=tile,
                                             cycle=report_cycle,
                                             errorrate=errorrate))
                        report_cycle += sign
                    expected_cycle += run_info.indexes_length
        return destination

    def get_sample_name(self, fastq1):
        """ Format sample name from a fastq1 Dataset object. """
        short_name, sample_num = fastq1.name.split('_')[:2]
        sample_name = short_name + '_' + sample_num
        return sample_name

    def launch_run(self, quality, fastq1, fastq2):
        """ Launch a run on Kive.

        @param folder: the path to the run folder
        @param quality: a Dataset object for the quality file
        @param fastq1: a Dataset object for the forward reads
        @param fastq2: a Dataset object for the reverse reads
        @return: (sample_name, run_status)
        """
        self.check_kive_connection()
        sample_name = self.get_sample_name(fastq1)
        name = '{} - {} ({})'.format(self.pipeline.family,
                                     sample_name,
                                     self.trimmed_folder)
        name = name[:MAX_RUN_NAME_LENGTH]

        logger.info('launching %s', name)
        # Note: order of inputs is critical
        status = self.kive.run_pipeline(pipeline=self.pipeline,
                                        inputs=[quality, fastq1, fastq2],
                                        name=name,
                                        groups=settings.kive_groups_allowed)
        return sample_name, status

    def find_preexisting_runs(self):
        """ Query Kive for all active runs, filter by pipeline.

        @return: {(quality_id, fastq1_id, fastq2_id): run}
        """
        runs = self.kive.find_runs(active=True)
        map = {}
        for run in runs:
            if run.pipeline_id == self.pipeline.pipeline_id:
                inputs = sorted(run.raw['inputs'], key=itemgetter('index'))
                fastq1 = self.kive.get_dataset(inputs[1]['dataset'])
                sample_name = self.get_sample_name(fastq1)
                input_ids = tuple(input['dataset'] for input in inputs)
                map[input_ids] = (sample_name, run)
        return map

    def get_run_key(self, quality, fastq1, fastq2):
        """ Calculate the key to look up preexisting runs for the given inputs.
        """
        return (quality.dataset_id, fastq1.dataset_id, fastq2.dataset_id)

    def is_run_complete(self, run):
        """ Check if a Kive run is complete.

        @param run: a pair of (sample_name, run_status)
        """
        _sample_name, run_status = run
        return run_status.is_complete()

    def mark_folder_disabled(self, folder, message, exc_info=None):
        """ Mark a run that failed, so it won't be processed again.

        @param folder: path to the run folder that had an error
        @param message: a description of the error
        @param exc_info: details about the error's exception in the standard tuple,
            True to look up the current exception, or None if there is no exception
            to report
        """
        failure_message = message + " - skipping run " + folder
        logger.error(failure_message, exc_info=exc_info)
        if settings.production:
            with open(os.path.join(folder, settings.ERROR_PROCESSING), 'w') as f:
                f.write(message)
        else:
            # in development mode - exit the monitor if a run fails
            sys.exit()
        return failure_message

    def log_retry(self, folder):
        logger.warn('Retrying folder %r.', folder, exc_info=True)

    def download_results(self, folder, runs):
        """ Download results from Kive.

        @param folder: the run folder
        @param runs: [(sample_name, run_status)] a sequence of pairs
        """
        # First, check that all runs in the batch were successful.
        for sample_name, run_status in runs:
            try:
                run_status.is_successful()
            except KiveRunFailedException:
                message = 'Sample {} failed in Kive.'.format(sample_name)
                self.mark_folder_disabled(folder, message, exc_info=True)
                return

        results_parent = os.path.join(folder, 'Results')
        if not os.path.exists(results_parent):
            os.mkdir(results_parent)
        results_folder = os.path.join(results_parent,
                                      'version_' + settings.pipeline_version)
        if not os.path.exists(results_folder):
            os.mkdir(results_folder)
        run_folder = os.path.join(settings.home, os.path.basename(folder))
        logger.info('downloading results for %r', folder)
        download_results(runs, results_folder, run_folder)
        update_qai.process_folder(results_folder)
        with open(os.path.join(results_folder, settings.DONE_PROCESSING), 'w'):
            pass  # Leave the file empty
        logger.info('completed folder %r', folder)

    def get_time(self):
        """ Get the current system time.

        Wrapped in a method to make it easier to mock when testing."""
        return datetime.now()

if __name__ == '__live_coding__':
    import unittest

    def setUp(self):
        logging.disable(logging.CRITICAL)  # avoid polluting test output
        self.loader = KiveLoader()
        self.existing_datasets = []
        self.existing_runs = {}
        self.uploaded = []
        self.launched = []
        self.completed = []
        self.downloaded = []
        self.now = datetime(2000, 1, 1)
        self.quality_cdt = 'quality CDT'
        self.disabled_folders = []
        self.folder_retries = []

        def check_kive_connection():
            self.loader.quality_cdt = self.quality_cdt
        self.loader.check_kive_connection = check_kive_connection
        self.loader.find_folders = lambda: []
        self.loader.find_files = lambda folder: []
        self.loader.find_preexisting_runs = lambda: self.existing_runs
        self.loader.get_run_key = lambda quality, fastq1, fastq2: (quality,
                                                                   fastq1,
                                                                   fastq2)
        self.loader.upload_kive_dataset = lambda filename, description, cdt: (
            self.uploaded.append((filename, description, cdt)) or filename)
        self.loader.download_quality = lambda folder: folder + '/quality.csv'
        self.loader.find_kive_dataset = lambda filename, cdt: (
            filename if filename in self.existing_datasets else None)
        self.loader.launch_run = lambda quality, fastq1, fastq2: (
            self.launched.append((quality, fastq1, fastq2)) or
            (quality, fastq1, fastq2))
        self.loader.is_run_complete = lambda run: run in self.completed
        self.loader.download_results = lambda folder, runs: (
            self.downloaded.append((folder, runs)))
        self.loader.get_time = lambda: self.now
        self.loader.mark_folder_disabled = lambda folder, message, exc_info: (
            self.disabled_folders.append(folder))
        self.loader.log_retry = lambda folder: self.folder_retries.append(folder)

    def test_something(self):
        # def test_failed_results_download(self):
        def download_results(folder):
            raise StandardError('Mock Kive failure.')

        self.loader.download_results = download_results

        self.loader.find_folders = lambda: ['run1']
        self.loader.find_files = lambda folder: [folder + '/sample1_R1_x.fastq']
        retry_delays = []

        first_delay = self.loader.poll()  # launch 1
        self.completed = self.launched[:]
        retry_delays.append(self.loader.poll())
        retry_delays.append(self.loader.poll())
        retry_delays.append(self.loader.poll())
        final_delay = self.loader.poll()

        self.assertEqual(0, first_delay)
        self.assertEqual(self.loader.retry_delay, sum(retry_delays))
        self.assertEqual(0, final_delay)

    class DummyTest(unittest.TestCase):
        def test_delegation(self):
            setUp(self)
            test_something(self)

    suite = unittest.TestSuite()
    suite.addTest(DummyTest("test_delegation"))
    test_results = unittest.TextTestRunner().run(suite)

    print(test_results.errors)
    print(test_results.failures)
