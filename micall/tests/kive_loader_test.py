from datetime import timedelta, datetime
import logging
import unittest

from micall.monitor.kive_loader import KiveLoader
from kiveapi.errors import KiveRunFailedException


class KiveLoaderTest(unittest.TestCase):
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

    def test_idle(self):
        delay = self.loader.poll()

        self.assertEqual(self.loader.status_delay, delay)

    def test_upload_one_sample(self):
        self.loader.find_folders = lambda: ['run2', 'run1']
        self.loader.find_files = lambda folder: [folder + '/sample1_R1_x.fastq',
                                                 folder + '/sample2_R1_x.fastq']
        expected_uploaded = [('run2/quality.csv',
                              'phiX174 quality from MiSeq run run2',
                              self.quality_cdt),
                             ('run2/sample1_R1_x.fastq',
                              'forward read from MiSeq run run2',
                              None),
                             ('run2/sample1_R2_x.fastq',
                              'reverse read from MiSeq run run2',
                              None)]

        self.loader.poll()

        self.assertEqual(expected_uploaded, self.uploaded)

    def test_launch_one_sample(self):
        self.loader.find_folders = lambda: ['run2', 'run1']
        self.loader.find_files = lambda folder: [folder + '/sample1_R1_x.fastq',
                                                 folder + '/sample2_R1_x.fastq']
        expected_launched = [('run2/quality.csv',
                              'run2/sample1_R1_x.fastq',
                              'run2/sample1_R2_x.fastq')]

        self.loader.poll()

        self.assertEqual(expected_launched, self.launched)

    def test_launch_two_samples(self):
        self.loader.find_folders = lambda: ['run1']
        self.loader.find_files = lambda folder: [folder + '/sample1_R1_x.fastq',
                                                 folder + '/sample2_R1_x.fastq']
        expected_launched = [('run1/quality.csv',
                              'run1/sample1_R1_x.fastq',
                              'run1/sample1_R2_x.fastq'),
                             ('run1/quality.csv',
                              'run1/sample2_R1_x.fastq',
                              'run1/sample2_R2_x.fastq')]

        self.loader.poll()
        self.loader.poll()

        self.assertEqual(expected_launched, self.launched)

    def test_launch_two_folders(self):
        self.loader.find_folders = lambda: ['run2', 'run1']
        self.loader.find_files = lambda folder: [folder + '/sample1_R1_x.fastq']
        expected_launched = [('run2/quality.csv',
                              'run2/sample1_R1_x.fastq',
                              'run2/sample1_R2_x.fastq'),
                             ('run1/quality.csv',
                              'run1/sample1_R1_x.fastq',
                              'run1/sample1_R2_x.fastq')]

        self.loader.poll()
        self.loader.poll()

        self.assertEqual(expected_launched, self.launched)

    def test_launch_finished(self):
        self.loader.find_folders = lambda: ['run1']
        self.loader.find_files = lambda folder: [folder + '/sample1_R1_x.fastq']
        expected_launched = [('run1/quality.csv',
                              'run1/sample1_R1_x.fastq',
                              'run1/sample1_R2_x.fastq')]

        self.loader.poll()
        self.loader.poll()

        self.assertEqual(expected_launched, self.launched)

    def test_trimmed_folder(self):
        self.loader.find_folders = lambda: ['path/160214_M01234_AX23']
        self.loader.find_files = lambda folder: [folder + '/sample1_R1_x.fastq']
        expected_uploaded = [('path/160214_M01234_AX23/quality.csv',
                              'phiX174 quality from MiSeq run 160214_M01234',
                              'quality CDT'),
                             ('path/160214_M01234_AX23/sample1_R1_x.fastq',
                              'forward read from MiSeq run 160214_M01234',
                              None),
                             ('path/160214_M01234_AX23/sample1_R2_x.fastq',
                              'reverse read from MiSeq run 160214_M01234',
                              None)]

        self.loader.poll()

        self.assertEqual(expected_uploaded, self.uploaded)

    def test_datasets_exist(self):
        self.existing_datasets = ['run1/quality.csv',
                                  'run1/sample1_R1_x.fastq',
                                  'run1/sample1_R2_x.fastq']
        self.loader.find_folders = lambda: ['run1']
        self.loader.find_files = lambda folder: [folder + '/sample1_R1_x.fastq']
        expected_uploaded = []

        self.loader.poll()

        self.assertEqual(expected_uploaded, self.uploaded)

    def test_some_datasets_exist(self):
        self.existing_datasets = ['run1/quality.csv',
                                  'run1/sample1_R2_x.fastq']
        self.loader.find_folders = lambda: ['run1']
        self.loader.find_files = lambda folder: [folder + '/sample1_R1_x.fastq']
        expected_uploaded = [('run1/sample1_R1_x.fastq',
                              'forward read from MiSeq run run1',
                              None)]

        self.loader.poll()

        self.assertEqual(expected_uploaded, self.uploaded)

    def test_download_one_sample(self):
        self.loader.find_folders = lambda: ['run1']
        self.loader.find_files = lambda folder: [folder + '/sample1_R1_x.fastq']
        expected_downloaded = [('run1', [('run1/quality.csv',
                                          'run1/sample1_R1_x.fastq',
                                          'run1/sample1_R2_x.fastq')])]

        self.loader.poll()  # launches
        self.loader.poll()  # checks status, not finished
        downloaded1 = self.downloaded[:]

        self.completed = self.launched[:]  # finish run
        self.loader.poll()  # finished, so now download
        downloaded2 = self.downloaded[:]

        self.assertEqual([], downloaded1)
        self.assertEqual(expected_downloaded, downloaded2)

    def test_download_two_samples(self):
        self.loader.find_folders = lambda: ['run1']
        self.loader.find_files = lambda folder: [folder + '/sample1_R1_x.fastq',
                                                 folder + '/sample2_R1_x.fastq']
        expected_downloaded = [('run1', [('run1/quality.csv',
                                          'run1/sample1_R1_x.fastq',
                                          'run1/sample1_R2_x.fastq'),
                                         ('run1/quality.csv',
                                          'run1/sample2_R1_x.fastq',
                                          'run1/sample2_R2_x.fastq')])]

        self.loader.poll()  # launch 1
        self.loader.poll()  # launch 2
        self.completed = self.launched[:1]  # finish sample 1
        self.loader.poll()  # check status, sample 2 not finished
        downloaded1 = self.downloaded[:]

        self.completed = self.launched[:]  # finish sample 2
        self.loader.poll()  # finished, so now download
        downloaded2 = self.downloaded[:]

        self.assertEqual([], downloaded1)
        self.assertEqual(expected_downloaded, downloaded2)

    def test_download_two_folders(self):
        self.loader.find_folders = lambda: ['run2', 'run1']
        self.loader.find_files = lambda folder: [folder + '/sample1_R1_x.fastq']
        expected_downloaded1 = [('run2', [('run2/quality.csv',
                                           'run2/sample1_R1_x.fastq',
                                           'run2/sample1_R2_x.fastq')])]
        expected_downloaded2 = [('run2', [('run2/quality.csv',
                                           'run2/sample1_R1_x.fastq',
                                           'run2/sample1_R2_x.fastq')]),
                                ('run1', [('run1/quality.csv',
                                           'run1/sample1_R1_x.fastq',
                                           'run1/sample1_R2_x.fastq')])]

        self.loader.poll()  # launch 1
        self.loader.poll()  # launch 2
        self.completed = self.launched[:1]  # finish sample 1
        self.loader.poll()  # check status, download 1
        downloaded1 = self.downloaded[:]

        self.completed = self.launched[:]  # finish sample 2
        self.loader.poll()  # check status, download 2
        downloaded2 = self.downloaded[:]

        self.assertEqual(expected_downloaded1, downloaded1)
        self.assertEqual(expected_downloaded2, downloaded2)

    def test_limit_launches(self):
        self.loader.launch_limit = 3
        self.loader.find_folders = lambda: ['run2', 'run1']
        self.loader.find_files = lambda folder: [folder + '/sample1_R1_x.fastq',
                                                 folder + '/sample2_R1_x.fastq']
        expected_launched1 = [('run2/quality.csv',
                               'run2/sample1_R1_x.fastq',
                               'run2/sample1_R2_x.fastq'),
                              ('run2/quality.csv',
                               'run2/sample2_R1_x.fastq',
                               'run2/sample2_R2_x.fastq'),
                              ('run1/quality.csv',
                               'run1/sample1_R1_x.fastq',
                               'run1/sample1_R2_x.fastq')]
        expected_launched2 = [('run2/quality.csv',
                               'run2/sample1_R1_x.fastq',
                               'run2/sample1_R2_x.fastq'),
                              ('run2/quality.csv',
                               'run2/sample2_R1_x.fastq',
                               'run2/sample2_R2_x.fastq'),
                              ('run1/quality.csv',
                               'run1/sample1_R1_x.fastq',
                               'run1/sample1_R2_x.fastq'),
                              ('run1/quality.csv',
                               'run1/sample2_R1_x.fastq',
                               'run1/sample2_R2_x.fastq')]

        self.loader.poll()  # launch run2, sample1
        self.loader.poll()  # launch run2, sample2
        self.loader.poll()  # launch run1, sample1
        self.loader.poll()  # don't launch anything, limit reached
        launched1 = self.launched[:]

        self.completed = self.launched[:1]
        self.loader.poll()  # launch run1, sample2
        launched2 = self.launched[:]

        self.assertEqual(expected_launched1, launched1)
        self.assertEqual(expected_launched2, launched2)

    def test_no_limit_on_latest_run(self):
        self.loader.launch_limit = 1
        self.loader.latest_folder = 'run2'
        self.loader.find_folders = lambda: ['run2', 'run1']
        self.loader.find_files = lambda folder: [folder + '/sample1_R1_x.fastq',
                                                 folder + '/sample2_R1_x.fastq']
        expected_launched = [('run2/quality.csv',
                              'run2/sample1_R1_x.fastq',
                              'run2/sample1_R2_x.fastq'),
                             ('run2/quality.csv',
                              'run2/sample2_R1_x.fastq',
                              'run2/sample2_R2_x.fastq')]

        self.loader.poll()  # launch run1, sample1
        self.loader.poll()  # launch run1, sample2
        self.loader.poll()  # don't launch anything, limit reached
        launched = self.launched[:]

        self.assertEqual(expected_launched, launched)

    def test_timing_wait(self):
        self.loader.find_folders = lambda: ['run1']
        self.loader.find_files = lambda folder: [folder + '/sample1_R1_x.fastq']

        sleep1 = self.loader.poll()  # upload, no wait before next upload
        sleep2 = self.loader.poll()  # check status, wait before next check

        self.assertEqual(0, sleep1)
        self.assertEqual(self.loader.status_delay, sleep2)

    def test_new_folder(self):
        self.loader.find_folders = lambda: ['run1']
        self.loader.find_files = lambda folder: [folder + '/sample1_R1_x.fastq']
        expected_launched1 = [('run1/quality.csv',
                               'run1/sample1_R1_x.fastq',
                               'run1/sample1_R2_x.fastq')]
        expected_launched2 = [('run1/quality.csv',
                               'run1/sample1_R1_x.fastq',
                               'run1/sample1_R2_x.fastq'),
                              ('run2/quality.csv',
                               'run2/sample1_R1_x.fastq',
                               'run2/sample1_R2_x.fastq')]

        self.loader.poll()  # launch run1, sample1
        self.loader.find_folders = lambda: ['run2', 'run1']
        self.now += timedelta(seconds=self.loader.folder_delay-1)

        self.loader.poll()  # not ready to scan for new folder

        launched1 = self.launched[:]
        self.now += timedelta(seconds=1)

        self.loader.poll()  # launch run2, sample1

        launched2 = self.launched[:]

        self.assertEqual(expected_launched1, launched1)
        self.assertEqual(expected_launched2, launched2)

    def test_completed_folder_no_reset(self):
        self.loader.find_folders = lambda: ['run2', 'run1']
        self.loader.find_files = lambda folder: [folder + '/sample1_R1_x.fastq']
        expected_launched = [('run2/quality.csv',
                              'run2/sample1_R1_x.fastq',
                              'run2/sample1_R2_x.fastq'),
                             ('run1/quality.csv',
                              'run1/sample1_R1_x.fastq',
                              'run1/sample1_R2_x.fastq')]

        self.loader.poll()  # launch 2
        self.loader.poll()  # launch 1
        self.completed = self.launched[:1]  # finish run 2
        self.loader.poll()  # check status, download 2

        self.loader.find_folders = lambda: ['run1']
        self.now += timedelta(seconds=self.loader.folder_delay)

        self.loader.poll()  # check status
        launched = self.launched[:]

        self.assertEqual(expected_launched, launched)

    def test_preexisting_runs(self):
        self.loader.find_folders = lambda: ['run1']
        self.loader.find_files = lambda folder: [folder + '/sample1_R1_x.fastq',
                                                 folder + '/sample2_R1_x.fastq']
        self.existing_runs = {('run1/quality.csv',
                               'run1/sample1_R1_x.fastq',
                               'run1/sample1_R2_x.fastq'): 'dummy run status'}
        expected_launched = [('run1/quality.csv',
                              'run1/sample2_R1_x.fastq',
                              'run1/sample2_R2_x.fastq')]

        self.loader.poll()
        self.loader.poll()

        self.assertEqual(expected_launched, self.launched)

    def test_failed_run(self):
        def is_run_complete(run):
            if run == ('run1/quality.csv',
                       'run1/sample2_R1_x.fastq',
                       'run1/sample2_R2_x.fastq'):
                raise KiveRunFailedException('Mock failure.')
            return True

        self.loader.is_run_complete = is_run_complete

        self.loader.find_folders = lambda: ['run1']
        self.loader.find_files = lambda folder: [folder + '/sample1_R1_x.fastq',
                                                 folder + '/sample2_R1_x.fastq']
        expected_downloaded = [('run1', [('run1/quality.csv',
                                          'run1/sample1_R1_x.fastq',
                                          'run1/sample1_R2_x.fastq'),
                                         ('run1/quality.csv',
                                          'run1/sample2_R1_x.fastq',
                                          'run1/sample2_R2_x.fastq')])]

        self.loader.poll()  # launch 1
        self.loader.poll()  # launch 2
        self.loader.poll()  # check status, catch exception
        downloaded = self.downloaded[:]

        self.assertEqual(expected_downloaded, downloaded)

    def test_unable_to_check_status(self):
        is_kive_running = False

        def is_run_complete(run):
            if not is_kive_running:
                raise StandardError('Kive connection failed.')
            return True

        self.loader.is_run_complete = is_run_complete

        self.loader.find_folders = lambda: ['run1']
        self.loader.find_files = lambda folder: [folder + '/sample1_R1_x.fastq']
        expected_downloaded = [('run1', [('run1/quality.csv',
                                          'run1/sample1_R1_x.fastq',
                                          'run1/sample1_R2_x.fastq')])]

        self.loader.poll()  # launch 1
        self.loader.poll()  # check status, catch exception
        downloaded1 = self.downloaded[:]

        is_kive_running = True
        self.loader.poll()  # check status successfully, and download
        downloaded2 = self.downloaded[:]

        self.assertEqual([], downloaded1)
        self.assertEqual(expected_downloaded, downloaded2)

    def test_failed_quality_download(self):
        def download_quality(folder):
            raise StandardError('Mock quality failure.')

        self.loader.download_quality = download_quality

        self.loader.find_folders = lambda: ['run1']
        self.loader.find_files = lambda folder: [folder + '/sample1_R1_x.fastq',
                                                 folder + '/sample2_R1_x.fastq']
        retry_delays = []

        retry_delays.append(self.loader.poll())
        retry_delays.append(self.loader.poll())
        retry_delays.append(self.loader.poll())
        final_delay = self.loader.poll()

        self.assertEqual(self.loader.retry_delay, sum(retry_delays))
        self.assertEqual(0, final_delay)

    def test_failed_results_download(self):
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
