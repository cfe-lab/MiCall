from datetime import timedelta, datetime
import logging
import os
import unittest

from micall.monitor.kive_loader import KiveLoader, RUN_COMPLETED, RUN_ACTIVE,\
    RUN_CANCELLED, RUN_PURGED, RUN_FAILED


DEFAULT_PIPELINE_ID = 1234
EXTRA_PIPELINE_ID = 9999


class KiveLoaderTest(unittest.TestCase):
    def setUp(self):
        logging.disable(logging.CRITICAL)  # avoid polluting test output
        self.loader = KiveLoader()
        self.loader.add_pipeline(pipeline_id=DEFAULT_PIPELINE_ID,
                                 inputs=['quality', 'fastq1', 'fastq2'],
                                 name_format='Default - {sample} ({folder})')
        self.existing_datasets = []
        self.existing_runs = {}
        self.uploaded = []
        self.launched = []
        self.batches = []
        self.cancelled = []
        self.completed = []
        self.failed = []
        self.purged = []
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
        self.loader.upload_kive_dataset = lambda filename, description, cdt: (
            self.uploaded.append((filename, description, cdt)) or filename)
        self.loader.download_quality = lambda folder: folder + '/quality.csv'
        self.loader.find_kive_dataset = lambda filename, cdt: (
            filename if filename in self.existing_datasets else None)
        self.loader.get_sample_name = lambda fastq1: os.path.basename(fastq1).split('_')[0]
        self.loader.fetch_run_status = lambda run: (RUN_COMPLETED
                                                    if run[1] in self.completed
                                                    else RUN_CANCELLED
                                                    if run[1] in self.cancelled
                                                    else RUN_ACTIVE)
        self.loader.fetch_output_status = lambda run: (RUN_PURGED
                                                       if run[1] in self.purged
                                                       else RUN_FAILED
                                                       if run[1] in self.failed
                                                       else RUN_COMPLETED)
        self.loader.get_run_id = lambda run: run
        self.loader.download_results = lambda folder, runs: (
            self.downloaded.append((folder, runs)))
        self.loader.get_time = lambda: self.now
        self.loader.mark_folder_disabled = lambda folder, message, exc_info: (
            self.disabled_folders.append(folder))
        self.loader.log_retry = lambda folder: self.folder_retries.append(folder)

        # Note: get_run_key must match return value of launch_run
        self.loader.get_run_key = lambda pipeline_id, quality, fastq1, fastq2: (
            self.loader.get_run_name(pipeline_id, self.loader.get_sample_name(fastq1)))
        self.loader.launch_run = lambda pipeline_id, run_name, inputs, batch_id: (
            self.launched.append(run_name) or run_name)
        self.loader.create_batch = lambda batch_name: (
            self.batches.append(batch_name) or len(self.batches))

    def test_idle(self):
        delay = self.loader.poll()

        self.assertEqual(self.loader.folder_delay, delay)

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
        expected_launched = ['Default - sample1 (run2)']

        self.loader.poll()

        self.assertEqual(expected_launched, self.launched)

    def test_launch_two_samples(self):
        self.loader.find_folders = lambda: ['run1']
        self.loader.find_files = lambda folder: [folder + '/sample1_R1_x.fastq',
                                                 folder + '/sample2_R1_x.fastq']
        expected_launched = ['Default - sample1 (run1)',
                             'Default - sample2 (run1)']
        expected_batches = ['run1']

        self.loader.poll()
        self.loader.poll()

        self.assertEqual(expected_launched, self.launched)
        self.assertEqual(expected_batches, self.batches)

    def test_launch_two_folders(self):
        self.loader.find_folders = lambda: ['run2', 'run1']
        self.loader.find_files = lambda folder: [folder + '/sample1_R1_x.fastq']
        expected_launched = ['Default - sample1 (run2)',
                             'Default - sample1 (run1)']
        expected_batches = ['run2', 'run1']

        self.loader.poll()
        self.loader.poll()

        self.assertEqual(expected_launched, self.launched)
        self.assertEqual(expected_batches, self.batches)

    def test_launch_two_pipelines(self):
        self.loader.add_pipeline(pipeline_id=EXTRA_PIPELINE_ID,
                                 inputs=['quality', 'fastq1', 'fastq2'],
                                 name_format='Extra - {sample} ({folder})')
        self.loader.find_folders = lambda: ['run2']
        self.loader.find_files = lambda folder: [folder + '/sample1_R1_x.fastq',
                                                 folder + '/sample2_R1_x.fastq']
        expected_launched = ['Default - sample1 (run2)',
                             'Extra - sample1 (run2)',
                             'Default - sample2 (run2)',
                             'Extra - sample2 (run2)']
        expected_downloaded = [('run2', [('sample1', 'Default - sample1 (run2)'),
                                         ('sample1', 'Extra - sample1 (run2)'),
                                         ('sample2', 'Default - sample2 (run2)'),
                                         ('sample2', 'Extra - sample2 (run2)')])]

        self.loader.poll()
        self.loader.poll()

        self.completed = expected_launched[:]
        self.loader.poll()

        self.assertEqual(expected_launched, self.launched)
        self.assertEqual(expected_downloaded, self.downloaded)

    def test_pipelines_diff_inputs(self):
        self.loader.launch_run = lambda pipeline_id, run_name, inputs, batch_id: (
            self.launched.append((pipeline_id, inputs)))
        self.loader.add_pipeline(pipeline_id=EXTRA_PIPELINE_ID,
                                 inputs=['fastq1'],
                                 name_format='Extra - {sample} ({folder})')
        self.loader.find_folders = lambda: ['run2', 'run1']
        self.loader.find_files = lambda folder: [folder + '/sample1_R1_x.fastq',
                                                 folder + '/sample2_R1_x.fastq']
        expected_launched = [(DEFAULT_PIPELINE_ID, ['run2/quality.csv',
                                                    'run2/sample1_R1_x.fastq',
                                                    'run2/sample1_R2_x.fastq']),
                             (EXTRA_PIPELINE_ID, ['run2/sample1_R1_x.fastq'])]

        self.loader.poll()

        self.assertEqual(expected_launched, self.launched)

    def test_pipelines_diff_patterns(self):
        self.loader.add_pipeline(pipeline_id=EXTRA_PIPELINE_ID,
                                 inputs=['quality', 'fastq1', 'fastq2'],
                                 pattern='HCV',
                                 name_format='Extra - {sample} ({folder})')
        self.loader.find_folders = lambda: ['run2']
        self.loader.find_files = lambda folder: [folder + '/sample1-HIV_R1_x.fastq',
                                                 folder + '/sample2-HCV_R1_x.fastq']
        expected_launched = ['Default - sample1-HIV (run2)',
                             'Default - sample2-HCV (run2)',
                             'Extra - sample2-HCV (run2)']
        expected_downloaded = [('run2', [('sample1-HIV', 'Default - sample1-HIV (run2)'),
                                         ('sample2-HCV', 'Default - sample2-HCV (run2)'),
                                         ('sample2-HCV', 'Extra - sample2-HCV (run2)')])]

        self.loader.poll()
        self.loader.poll()

        self.completed = expected_launched[:]
        self.loader.poll()

        self.assertEqual(expected_launched, self.launched)
        self.assertEqual(expected_downloaded, self.downloaded)

    def test_launch_finished(self):
        self.loader.find_folders = lambda: ['run1']
        self.loader.find_files = lambda folder: [folder + '/sample1_R1_x.fastq']
        expected_launched = ['Default - sample1 (run1)']

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
        expected_downloaded = [('run1', [('sample1', 'Default - sample1 (run1)')])]

        self.loader.poll()  # launches
        self.loader.poll()  # checks status, not finished
        downloaded1 = self.downloaded[:]

        self.completed = self.launched[:]  # finish run
        self.loader.poll()  # finished, so now download
        downloaded2 = self.downloaded[:]

        self.assertEqual([], downloaded1)
        self.assertEqual(expected_downloaded, downloaded2)
        self.assertEqual({}, self.loader.batches)

    def test_download_two_samples(self):
        self.loader.find_folders = lambda: ['run1']
        self.loader.find_files = lambda folder: [folder + '/sample1_R1_x.fastq',
                                                 folder + '/sample2_R1_x.fastq']
        expected_downloaded = [('run1', [('sample1', 'Default - sample1 (run1)'),
                                         ('sample2', 'Default - sample2 (run1)')])]

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
        expected_downloaded1 = [('run2', [('sample1', 'Default - sample1 (run2)')])]
        expected_downloaded2 = [('run2', [('sample1', 'Default - sample1 (run2)')]),
                                ('run1', [('sample1', 'Default - sample1 (run1)')])]

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
        expected_launched1 = ['Default - sample1 (run2)',
                              'Default - sample2 (run2)',
                              'Default - sample1 (run1)']
        expected_launched2 = ['Default - sample1 (run2)',
                              'Default - sample2 (run2)',
                              'Default - sample1 (run1)',
                              'Default - sample2 (run1)']

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
        expected_launched = ['Default - sample1 (run2)',
                             'Default - sample2 (run2)']

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
        expected_launched1 = ['Default - sample1 (run1)']
        expected_launched2 = ['Default - sample1 (run1)',
                              'Default - sample1 (run2)']

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
        expected_launched = ['Default - sample1 (run2)',
                             'Default - sample1 (run1)']

        self.loader.poll()  # launch 2
        self.loader.poll()  # launch 1
        self.completed = self.launched[:1]  # finish run 2
        self.loader.poll()  # check status, download 2

        self.loader.find_folders = lambda: ['run1']
        self.now += timedelta(seconds=self.loader.folder_delay)

        self.loader.poll()  # check status
        launched = self.launched[:]

        self.assertEqual(expected_launched, launched)

    def test_cancelled_run(self):
        self.loader.find_folders = lambda: ['run1']
        self.loader.find_files = lambda folder: [folder + '/sample1_R1_x.fastq',
                                                 folder + '/sample2_R1_x.fastq']
        expected_launched = ['Default - sample1 (run1)',
                             'Default - sample2 (run1)',
                             'Default - sample1 (run1)']
        expected_downloaded = []

        self.loader.poll()  # launch 1
        self.loader.poll()  # launch 2
        self.cancelled = self.launched[:1]  # cancel run 1
        self.loader.poll()  # check status, rerun 1

        launched = self.launched[:]
        downloaded = self.downloaded[:]

        self.assertEqual(expected_launched, launched)
        self.assertEqual(expected_downloaded, downloaded)

    def test_purged_run(self):
        self.loader.find_folders = lambda: ['run1']
        self.loader.find_files = lambda folder: [folder + '/sample1_R1_x.fastq',
                                                 folder + '/sample2_R1_x.fastq']
        expected_launched = ['Default - sample1 (run1)',
                             'Default - sample2 (run1)',
                             'Default - sample1 (run1)']
        expected_downloaded = []

        self.loader.poll()  # launch 1
        self.loader.poll()  # launch 2
        self.completed = self.launched[:1]  # complete run 1
        self.loader.poll()
        self.purged = self.launched[:1]  # purge some results from run 1
        self.loader.poll()  # check status, rerun 1

        launched = self.launched[:]
        downloaded = self.downloaded[:]

        self.assertEqual(expected_launched, launched)
        self.assertEqual(expected_downloaded, downloaded)

    def test_preexisting_runs(self):
        self.loader.find_folders = lambda: ['run1']
        self.loader.find_files = lambda folder: [folder + '/sample1_R1_x.fastq',
                                                 folder + '/sample2_R1_x.fastq']
        self.existing_runs = {'Default - sample1 (run1)': 'dummy run status'}
        expected_launched = ['Default - sample2 (run1)']

        self.loader.poll()
        self.loader.poll()

        self.assertEqual(expected_launched, self.launched)

    def test_unable_to_check_status(self):
        is_kive_running = False

        # noinspection PyUnusedLocal
        def fetch_run_status(run):
            if not is_kive_running:
                raise StandardError('Kive connection failed.')
            return RUN_COMPLETED

        self.loader.fetch_run_status = fetch_run_status

        self.loader.find_folders = lambda: ['run1']
        self.loader.find_files = lambda folder: [folder + '/sample1_R1_x.fastq']
        expected_downloaded = [('run1', [('sample1', 'Default - sample1 (run1)')])]

        self.loader.poll()  # launch 1
        self.loader.poll()  # check status, catch exception
        downloaded1 = self.downloaded[:]

        is_kive_running = True
        self.loader.poll()  # check status successfully, and download
        downloaded2 = self.downloaded[:]

        self.assertEqual([], downloaded1)
        self.assertEqual(expected_downloaded, downloaded2)

    def test_failed_quality_download(self):
        # noinspection PyUnusedLocal
        def download_quality(folder):
            raise StandardError('Mock quality failure.')

        self.loader.download_quality = download_quality

        self.loader.find_folders = lambda: ['run1']
        self.loader.find_files = lambda folder: [folder + '/sample1_R1_x.fastq',
                                                 folder + '/sample2_R1_x.fastq']

        retry_delays = [self.loader.poll(), self.loader.poll(), self.loader.poll()]
        final_delay = self.loader.poll()

        self.assertEqual(self.loader.retry_delay, sum(retry_delays))
        self.assertEqual(0, final_delay)

    def test_failed_results_download(self):
        self.loader.download_results = lambda folder, runs: 1/0  # fails
        expected_disabled = ['run1']

        self.loader.find_folders = lambda: ['run1']
        self.loader.find_files = lambda folder: [folder + '/sample1_R1_x.fastq']

        self.loader.poll()  # launch 1
        self.completed = self.launched[:]
        self.loader.poll()

        self.assertEqual(expected_disabled, self.disabled_folders)

    def test_removing_failure_will_retry_folder(self):
        self.loader.find_folders = lambda: ['run2', 'run1']
        self.loader.find_files = lambda folder: [folder + '/sample1_R1_x.fastq']
        expected_launched1 = ['Default - sample1 (run2)',
                              'Default - sample1 (run1)']
        expected_launched2 = ['Default - sample1 (run2)',
                              'Default - sample1 (run1)',
                              'Default - sample1 (run2)']

        self.loader.poll()  # launch 2
        self.loader.poll()  # launch 1
        self.completed = self.launched[:1]  # fail run 2
        self.loader.download_results = lambda folder, runs: 1/0  # fails
        self.loader.poll()  # check status, mark 2 as failed

        self.loader.find_folders = lambda: ['run1']  # run2 marked as failed
        self.now += timedelta(seconds=self.loader.folder_delay)

        self.loader.poll()  # Nothing to launch
        launched1 = self.launched[:]

        self.loader.find_folders = lambda: ['run2', 'run1']  # run2 unmarked
        self.now += timedelta(seconds=self.loader.folder_delay)

        self.loader.poll()  # Should relaunch run2
        launched2 = self.launched[:]

        self.assertEqual(expected_launched1, launched1)
        self.assertEqual(expected_launched2, launched2)
