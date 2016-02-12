from collections import namedtuple, defaultdict
import csv
from glob import glob
import hashlib
import itertools
import logging
import operator
import os
import shutil
from xml.etree import ElementTree

from micall import settings
from micall.monitor import qai_helper, update_qai
from micall.monitor.kive_download import kive_login, download_results

MAX_RUN_NAME_LENGTH = 60
logger = logging.getLogger("kive_loader")


class KiveLoader(object):
    RunInfo = namedtuple('RunInfo', 'miseq_run_id read_lengths indexes_length')

    def __init__(self):
        self.folders = None
        self.kive = None
        self.active_runs = []
        self.batches = defaultdict(list)  # {folder: [run]}
        self.status_check_index = 0

    def poll(self):
        if self.folders is None:
            self.folders = self.find_folders()
            self.folder_count = 0
            self.files = []
            self.file_count = 0
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
        self.check_run_status()
        if self.folder is None:
            return
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
        run = self.launch_run(self.quality_dataset, dataset1, dataset2)
        self.active_runs.append(run)
        self.batches[self.folder].append(run)

    def check_run_status(self):
        if not self.active_runs:
            return
        self.status_check_index = (self.status_check_index + 1) % len(self.active_runs)
        run = self.active_runs[self.status_check_index]
        if self.is_run_complete(run):
            self.active_runs.remove(run)
        for folder, runs in self.batches.items():
            if folder != self.folder and all(run not in self.active_runs
                                             for run in runs):
                self.download_results(folder, runs)
                del self.batches[folder]

    def find_folders(self):
        """ Find run folders ready to be processed.

        @return: a list of paths to the folders in the order they should be
            processed
        """
        # flag indicates that Illumina MiseqReporter has completed pre-processing, files available on NAS
        flag_files = glob(settings.rawdata_mount +
                          'MiSeq/runs/*/{}'.format(settings.NEEDS_PROCESSING))
        folders = []
        for flag_file in flag_files:
            folder = os.path.dirname(flag_file)
            result_path = os.path.join(folder,
                                       'Results',
                                       'version_{}'.format(settings.pipeline_version))
            done_path = os.path.join(result_path, settings.DONE_PROCESSING)
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

        folders.sort(reverse=True)
        return folders

    def is_marked_as_disabled(self, folder):
        return os.path.exists(os.path.join(folder, settings.ERROR_PROCESSING))

    def is_quality_control_uploaded(self, folder):
        return os.path.exists(os.path.join(folder, settings.QC_UPLOADED))

    def find_files(self, folder):
        """ Find FASTQ files within a folder.

        @return: a list of paths to the files within the folder.
        """
        gz_files = glob(os.path.join(folder,
                                     'Data/Intensities/BaseCalls/*_R1_001.fastq.gz'))

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
            datasets = self.kive.find_datasets(dataset_name=dataset_name,
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
        destination = os.path.join(settings.home,
                                   os.path.basename(folder),
                                   '{}_quality.csv'.format(trimmed_folder))
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
                                                        operator.itemgetter('tile')):
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

    def launch_run(self, quality, fastq1, fastq2):
        """ Launch a run on Kive.

        @param folder: the path to the run folder
        @param quality: a Dataset object for the quality file
        @param fastq1: a Dataset object for the forward reads
        @param fastq2: a Dataset object for the reverse reads
        @return: (sample_name, run_status)
        """
        self.check_kive_connection()
        short_name, sample_num = fastq1.name.split('_')[:2]
        sample_name = short_name + '_' + sample_num
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

    def is_run_complete(self, run):
        """ Check if a Kive run is complete.

        @param run: a pair of (sample_name, run_status)
        """
        _sample_name, run_status = run
        return run_status.is_complete()

    def download_results(self, folder, runs):
        """ Download results from Kive.

        @param folder: the run folder
        @param runs: [(sample_name, run_status)] a sequence of pairs
        """
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

if __name__ == '__live_coding__':
    import unittest

    def setUp(self):
        self.loader = KiveLoader()
        self.existing_datasets = []
        self.uploaded = []
        self.launched = []
        self.completed = []
        self.downloaded = []
        self.quality_cdt = 'quality CDT'

        def check_kive_connection():
            self.loader.quality_cdt = self.quality_cdt
        self.loader.check_kive_connection = check_kive_connection
        self.loader.find_folders = lambda: []
        self.loader.find_files = lambda folder: []
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

    def test_something(self):
        # def test_download_two_folders(self):
        self.loader.find_folders = lambda: ['run1', 'run2']
        self.loader.find_files = lambda folder: [folder + '/sample1_R1_x.fastq']
        expected_downloaded1 = [('run1', [('run1/quality.csv',
                                           'run1/sample1_R1_x.fastq',
                                           'run1/sample1_R2_x.fastq')])]
        expected_downloaded2 = [('run1', [('run1/quality.csv',
                                           'run1/sample1_R1_x.fastq',
                                           'run1/sample1_R2_x.fastq')]),
                                ('run2', [('run2/quality.csv',
                                           'run2/sample1_R1_x.fastq',
                                           'run2/sample1_R2_x.fastq')])]

        self.loader.poll()  # launch 1
        self.loader.poll()  # check status 1, launch 2
        self.completed = self.launched[:1]  # finish sample 1
        self.loader.poll()  # check status 2
        self.loader.poll()  # check status 1, sample 2 irrelevant
        downloaded1 = self.downloaded[:]

        self.completed = self.launched[:]  # finish sample 2
        self.loader.poll()  # finished, so now download
        downloaded2 = self.downloaded[:]

        self.assertEqual(expected_downloaded1, downloaded1)
        self.assertEqual(expected_downloaded2, downloaded2)

    class DummyTest(unittest.TestCase):
        def test_delegation(self):
            setUp(self)
            test_something(self)

    suite = unittest.TestSuite()
    suite.addTest(DummyTest("test_delegation"))
    test_results = unittest.TextTestRunner().run(suite)

    print(test_results.errors)
    print(test_results.failures)
