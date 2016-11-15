from collections import namedtuple, defaultdict, Counter
import csv
from datetime import datetime, timedelta
from glob import glob
import hashlib
import itertools
import logging
from operator import itemgetter
import os
import re
import shutil
import subprocess
import sys
from xml.etree import ElementTree

from kiveapi.errors import KiveRunFailedException, KiveClientException
from micall import settings
from micall.monitor import qai_helper, update_qai
from micall.monitor.kive_download import kive_login, download_results

MAX_RUN_NAME_LENGTH = 60
RUN_ACTIVE = 'active'
RUN_CANCELLED = 'cancelled'
RUN_COMPLETED = 'completed'
RUN_FAILED = 'failed'
RUN_PURGED = 'purged'
logger = logging.getLogger("kive_loader")


def kive_retries(target):
    def wrapper(self, *args, **kwargs):
        try:
            target(self, *args, **kwargs)
        except KiveClientException:
            logger.warn('Retrying with a fresh Kive login.', exc_info=True)
            self.refresh_login()
            target(self, *args, **kwargs)

    return wrapper


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
        self.latest_folder = self.folder = None
        self.kive = None
        self.preexisting_runs = None
        self.active_runs = []  # [(sample_name, run_status)]
        self.active_inputs = {}  # {run_id: (pipeline_id, dataset1, dataset2, quality, batch_id)}
        self.batches = defaultdict(list)  # {folder: [run]}
        self.retry_counts = Counter()  # {folder: count}
        self.is_status_available = True
        self.pipelines = {}
        self.miseq_runs_path = os.path.join(settings.rawdata_mount,
                                            'MiSeq/runs')

    def add_pipeline(self, id, inputs, format='', pattern=''):
        """ Add another pipeline definition for launching.

        :param int id: the pipeline id in Kive
        :param list(str) inputs: input names, in order
        :param str format: format string for run names, for example:
            'My Pipeline - {sample} ({folder})'
        :param str pattern: regular expression to match sample names that this
            pipeline should be launched for. Can match any part of the sample
            name.
        """
        self.pipelines[id] = dict(inputs=inputs, format=format, pattern=pattern)

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
                    self.batch_id = self.create_batch(self.trimmed_folder)
            if not self.can_launch():
                self.check_run_status()
            if not self.can_launch():
                return (self.status_delay
                        if self.active_runs
                        else self.folder_delay)
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
            for pipeline_id in self.pipelines:
                if self.preexisting_runs is None:
                    self.preexisting_runs = self.find_preexisting_runs()
                run_key = self.get_run_key(pipeline_id, self.quality_dataset, dataset1, dataset2)
                run = self.preexisting_runs.pop(run_key, None)
                if run is None:
                    run = self.launch(pipeline_id,
                                      dataset1,
                                      dataset2,
                                      self.quality_dataset,
                                      self.batch_id)
                if run is not None:
                    self.active_runs.append(run)
                    self.batches[self.folder].append(run)
                    run_id = self.get_run_id(run)
                    self.active_inputs[run_id] = (pipeline_id,
                                                  dataset1,
                                                  dataset2,
                                                  self.quality_dataset,
                                                  self.batch_id)
            return 0
        except Exception as ex:
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

    def create_batch(self, trimmed_folder):
        name = trimmed_folder + ' ' + settings.pipeline_version
        description = 'MiCall batch for folder {}, pipeline version {}.'.format(
            trimmed_folder,
            settings.pipeline_version)
        batch = self.kive.create_run_batch(name,
                                           description=description,
                                           users=[],
                                           groups=settings.kive_groups_allowed)
        return batch.id

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
        self.active_inputs = {}
        self.preexisting_runs = None
        self.batches.clear()

    def can_launch(self):
        if self.folder is None:
            return False
        return (self.folder == self.latest_folder or
                len(self.active_runs) < self.launch_limit)

    def get_run_id(self, run):
        _sample_name, run_status = run
        return run_status.run_id

    def check_run_status(self):
        for folder, runs in self.batches.items():
            is_batch_active = folder == self.folder
            for run in runs:
                try:
                    active_index = self.active_runs.index(run)
                except ValueError:
                    active_index = None
                if active_index is None:
                    try:
                        output_status = self.fetch_output_status(run)
                        self.is_status_available = True
                        if output_status == RUN_PURGED:
                            self.relaunch(run, folder)
                            is_batch_active = True
                    except:
                        if self.is_status_available:
                            logger.warn('Unable to check output status.',
                                        exc_info=True)
                            self.is_status_available = False
                else:
                    try:
                        run_status = self.fetch_run_status(run)
                        self.is_status_available = True
                        is_complete = run_status == RUN_COMPLETED
                        if is_complete:
                            del self.active_runs[active_index]
                        else:
                            is_batch_active = True
                            if run_status == RUN_CANCELLED:
                                # Rerun the cancelled run.
                                self.relaunch(run, folder)
                    except:
                        is_batch_active = True
                        if self.is_status_available:
                            # First failure, so log it.
                            logger.warn('Unable to check run status.',
                                        exc_info=True)
                            self.is_status_available = False
            if not is_batch_active:
                self.downloading_folder = folder
                self.download_results(folder, runs)
                for run in runs:
                    del self.active_inputs[self.get_run_id(run)]
                del self.batches[folder]
                self.downloading_folder = None

    def relaunch(self, run, batch_folder):
        run_id = self.get_run_id(run)
        inputs = self.active_inputs[run_id]
        new_run = self.launch(*inputs)
        new_run_id = self.get_run_id(new_run)
        del self.active_inputs[run_id]
        self.active_inputs[new_run_id] = inputs
        batch_runs = self.batches[batch_folder]
        batch_index = batch_runs.index(run)
        batch_runs[batch_index] = new_run
        try:
            active_index = self.active_runs.index(run)
            self.active_runs[active_index] = new_run
        except ValueError:
            self.active_runs.append(new_run)

    def find_folders(self):
        """ Find run folders ready to be processed.

        Also set self.latest_folder to be the latest folder that is marked
        with a needsprocessing file, even if it's done processing or failed.

        @return: a list of paths to the folders in the order they should be
            processed
        """
        # flag indicates that Illumina MiseqReporter has completed pre-processing, files available on NAS
        flag_files = glob(os.path.join(self.miseq_runs_path,
                                       '*',
                                       settings.NEEDS_PROCESSING))
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
            self.refresh_login()
            # retrieve Pipeline object based on version
            for pipeline_id, pipeline in self.pipelines.iteritems():
                pipeline['kive'] = self.kive.get_pipeline(pipeline_id)

            # retrieve quality.csv compound data type
            self.quality_cdt = self.kive.get_cdt(settings.quality_cdt_kive_id)

            # retrieve external file directory
            directories = self.kive.get('/api/externalfiledirectories',
                                        is_json=True).json()
            self.external_directory_name = self.external_directory_path = None
            for directory in directories:
                if self.miseq_runs_path.startswith(directory['path']):
                    self.external_directory_name = directory['name']
                    self.external_directory_path = directory['path']
                    break

    def refresh_login(self):
        self.kive = kive_login(settings.kive_server_url,
                               settings.kive_user,
                               settings.kive_password)

    def upload_kive_dataset(self, filename, description, cdt):
        """ Upload a dataset to Kive.

        @return: the dataset object from the Kive API wrapper
        """
        self.check_kive_connection()
        dataset_name = os.path.basename(filename)
        logger.info('uploading dataset %r', dataset_name)
        if (self.external_directory_name is None or
                not filename.startswith(self.external_directory_name)):
            with open(filename, 'rb') as f:
                dataset = self.kive.add_dataset(
                    name=dataset_name,
                    description=description,
                    handle=f,
                    cdt=cdt,
                    groups=settings.kive_groups_allowed)
        else:
            external_path = os.path.relpath(filename,
                                            self.external_directory_path)
            dataset = self.kive.add_dataset(
                name=dataset_name,
                description=description,
                handle=None,
                externalfiledirectory=self.external_directory_name,
                external_path=external_path,
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

    def get_run_name(self, pipeline_id, sample_name):
        pipeline = self.pipelines[pipeline_id]
        name = pipeline['format'].format(sample=sample_name,
                                         folder=self.trimmed_folder)
        name = name[:MAX_RUN_NAME_LENGTH]
        return name

    def launch(self, pipeline_id, fastq1, fastq2, quality, batch_id):
        """ Prepare and launch a run on Kive.

        :param int pipeline_id: the pipeline to launch with
        :param fastq1: a Dataset object for the forward reads
        :param fastq2: a Dataset object for the reverse reads
        :param quality: a Dataset object for the quality file
        :return: (sample_name, run_status)
        """
        pipeline = self.pipelines[pipeline_id]
        sample_name = self.get_sample_name(fastq1)
        is_match = re.search(pipeline['pattern'], sample_name) is not None
        if not is_match:
            return
        name = self.get_run_name(pipeline_id, sample_name)
        inputs = self.build_inputs(pipeline_id, fastq1, fastq2, quality)
        return sample_name, self.launch_run(pipeline_id, name, inputs, batch_id)

    def build_inputs(self, pipeline_id, fastq1, fastq2, quality):
        # Note: order of inputs is critical
        input_dict = dict(fastq1=fastq1, fastq2=fastq2, quality=quality)
        pipeline = self.pipelines[pipeline_id]
        return [input_dict[s] for s in pipeline['inputs']]

    def launch_run(self, pipeline_id, name, inputs, batch_id):
        self.check_kive_connection()

        logger.info('launching %s', name)
        kive_pipeline = self.pipelines[pipeline_id]['kive']
        status = self.kive.run_pipeline(pipeline=kive_pipeline,
                                        inputs=inputs,
                                        name=name,
                                        runbatch=batch_id,
                                        groups=settings.kive_groups_allowed)
        return status

    def find_preexisting_runs(self):
        """ Query Kive for all active runs, filter by pipeline.

        @return: {(quality_id, fastq1_id, fastq2_id): run}
        """
        runs = self.kive.find_runs(active=True)
        map = {}
        for run in runs:
            pipeline_config = self.pipelines.get(run.pipeline_id)
            if pipeline_config is not None:
                try:
                    run.is_complete()
                    fastq1_index = pipeline_config['inputs'].index('fastq1')
                    inputs = sorted(run.raw['inputs'], key=itemgetter('index'))
                    fastq1 = self.kive.get_dataset(inputs[fastq1_index]['dataset'])
                    sample_name = self.get_sample_name(fastq1)
                    input_ids = tuple(input['dataset'] for input in inputs)
                    key = (run.pipeline_id, ) + input_ids
                    map[key] = (sample_name, run)
                except KiveRunFailedException:
                    # Failed or cancelled, rerun.
                    pass
        return map

    def get_run_key(self, pipeline_id, quality, fastq1, fastq2):
        """ Calculate the key to look up preexisting runs for the given inputs.
        """
        inputs = self.build_inputs(pipeline_id, fastq1, fastq2, quality)
        input_ids = tuple(inp.dataset_id for inp in inputs)
        return (pipeline_id, ) + input_ids

    @kive_retries
    def fetch_run_status(self, run):
        """ Fetch the run status from Kive.

        :return: RUN_ACTIVE, RUN_COMPLETED, or RUN_CANCELLED
        """
        _sample_name, run_status = run
        details = self.kive.get_run(run_status.run_id).raw
        if details['stopped_by'] is not None:
            return RUN_CANCELLED
        if details['end_time'] is not None:
            return RUN_COMPLETED
        return RUN_ACTIVE

    @kive_retries
    def fetch_output_status(self, run):
        """ Fetch the output status for a completed run from Kive.

        :return: RUN_COMPLETED, RUN_FAILED, or RUN_PURGED
        """
        _sample_name, run_status = run
        outputs = run_status.get_results()
        for output in outputs.itervalues():
            if not output.raw['is_ok']:
                return RUN_FAILED
            if output.raw['id'] is None:
                return RUN_PURGED
        return RUN_COMPLETED

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
        try:
            try:
                # First, check that all runs in the batch were successful.
                for sample_name, run_status in runs:
                    run_status.is_successful()
            except KiveRunFailedException:
                _, _, traceback = sys.exc_info()
                message = 'Sample {} failed in Kive.'.format(sample_name)
                raise KiveRunFailedException, message, traceback

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
        except Exception:
            message = 'Downloading results failed.'
            self.mark_folder_disabled(folder, message, exc_info=True)

    def get_time(self):
        """ Get the current system time.

        Wrapped in a method to make it easier to mock when testing."""
        return datetime.now()

if __name__ == '__live_coding__':

    import unittest
    from micall.tests.kive_loader_test import KiveLoaderTest

    suite = unittest.TestSuite()
    suite.addTest(KiveLoaderTest("test_download_one_sample"))
    test_results = unittest.TextTestRunner().run(suite)

    print(test_results.errors)
    print(test_results.failures)
