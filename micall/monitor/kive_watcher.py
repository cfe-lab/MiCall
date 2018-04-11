import hashlib
import logging
import os
import shutil
import tarfile
from collections import namedtuple
from datetime import datetime, timedelta
from enum import Enum
from pathlib import Path
from queue import Full

from io import StringIO, BytesIO
from time import sleep

from micall.drivers.run_info import parse_read_sizes
from micall.monitor import error_metrics_parser
from micall.monitor.sample_watcher import FolderWatcher, ALLOWED_GROUPS, SampleWatcher, PipelineType
from micall.resistance.resistance import find_groups

try:
    from kiveapi import KiveAPI
    from kiveapi.runstatus import RunStatus
except ImportError:
    # Ignore import errors during testing.
    KiveAPI = RunStatus = None

try:
    from requests.adapters import HTTPAdapter
except ImportError:
    # Ignore import errors during testing.
    HTTPAdapter = None

logger = logging.getLogger(__name__)
FOLDER_SCAN_INTERVAL = timedelta(minutes=2)
SLEEP_SECONDS = 60
DOWNLOADED_RESULTS = ['remap_counts_csv',
                      'conseq_csv',
                      'conseq_ins_csv',
                      'failed_csv',
                      'nuc_csv',
                      'amino_csv',
                      'coord_ins_csv',
                      'failed_align_csv',
                      'g2p_csv',
                      'g2p_summary_csv',
                      'coverage_scores_csv',
                      'coverage_maps_tar',
                      'cascade_csv',
                      'mixed_counts_csv',
                      'mixed_amino_csv',
                      'mixed_amino_merged_csv',
                      'resistance_csv',
                      'mutations_csv']

# noinspection PyArgumentList
FolderEventType = Enum('FolderEventType', 'ADD_SAMPLE FINISH_FOLDER')
FolderEvent = namedtuple('FolderEvent', 'base_calls type sample_group')


def open_kive(server_url):
    if KiveAPI is None:
        raise ImportError('Kive API failed to import. Is it installed?')
    if HTTPAdapter is None:
        raise ImportError('requests module failed to import. Is it installed?')
    session = KiveAPI(server_url)
    session.mount('https://', HTTPAdapter(max_retries=20))
    return session


def now():
    """ Get the current date/time.

    Wrapped in a local function so it can be mocked during testing.
    """
    return datetime.now()


def find_samples(raw_data_folder, pipeline_version, sample_queue, wait=True):
    while True:
        is_complete = scan_samples(raw_data_folder,
                                   pipeline_version,
                                   sample_queue,
                                   wait)
        if is_complete and not wait:
            break


def scan_samples(raw_data_folder, pipeline_version, sample_queue, wait):
    next_scan = now() + FOLDER_SCAN_INTERVAL
    flag_paths = sorted(
        raw_data_folder.glob("MiSeq/runs/*/needsprocessing"),
        reverse=True)
    is_found = False
    for flag_path in flag_paths:
        run_path = flag_path.parent
        done_path = (run_path /
                     f"Results/version_{pipeline_version}/doneprocessing")
        if done_path.exists():
            continue
        base_calls_path = run_path / "Data/Intensities/BaseCalls"
        fastq_files = base_calls_path.glob("*_R1_*.fastq.gz")
        sample_sheet_path = run_path / "SampleSheet.csv"
        file_names = [f.name for f in fastq_files]
        for sample_group in find_groups(file_names, sample_sheet_path):
            is_found = True
            is_sent = send_event(sample_queue,
                                 FolderEvent(base_calls_path,
                                             FolderEventType.ADD_SAMPLE,
                                             sample_group),
                                 next_scan)
            if not is_sent:
                return False
        is_sent = send_event(sample_queue,
                             FolderEvent(base_calls_path,
                                         FolderEventType.FINISH_FOLDER,
                                         None),
                             next_scan)
        if not is_sent:
            return False
    if not is_found:
        logger.info('No folders need processing.')
    while wait and now() < next_scan:
        sleep(SLEEP_SECONDS)
    return True


def send_event(sample_queue, folder_event, next_scan):
    is_sent = False
    while not is_sent and now() < next_scan:
        try:
            sample_queue.put(folder_event,
                             timeout=SLEEP_SECONDS)
            logger.debug('Put %s in queue.', folder_event)
            is_sent = True
        except Full:
            pass
    return is_sent


def trim_name(sample_name):
    return '_'.join(sample_name.split('_')[:2])


def get_output_filename(output_name):
    return '.'.join(output_name.rsplit('_', 1))


class KiveWatcher:
    def __init__(self,
                 config=None,
                 result_handler=lambda result_folder: None):
        """ Initialize.

        :param config: command line arguments
        :param result_handler: called when a run folder has collated all the
            results into a result folder"""
        self.config = config
        self.result_handler = result_handler
        self.session = None
        self.loaded_folders = set()  # base_calls folders with all samples loaded
        self.folder_watchers = {}  # {base_calls_folder: FolderWatcher}
        self.pipelines = {}  # {pipeline_id: pipeline}
        self.external_directory_path = self.external_directory_name = None
        self.input_pipeline_ids = config and dict(
            quality_csv=config.micall_filter_quality_pipeline_id,
            bad_cycles_csv=config.micall_main_pipeline_id,
            main_amino_csv=config.micall_resistance_pipeline_id,
            midi_amino_csv=config.micall_resistance_pipeline_id) or {}

    def is_full(self):
        active_count = sum(len(folder_watcher.active_samples)
                           for folder_watcher in self.folder_watchers.values())
        return active_count >= self.config.max_active

    def get_kive_pipeline(self, pipeline_id):
        self.check_session()
        kive_pipeline = self.pipelines.get(pipeline_id)
        if kive_pipeline is None:
            kive_pipeline = self.session.get_pipeline(pipeline_id)
            self.pipelines[pipeline_id] = kive_pipeline
        return kive_pipeline

    def get_kive_input(self, input_name):
        pipeline_id = self.input_pipeline_ids[input_name]
        kive_pipeline = self.get_kive_pipeline(pipeline_id)
        for kive_input in kive_pipeline.inputs:
            if kive_input.dataset_name == input_name:
                return kive_input
        raise ValueError('Input {} not found on pipeline id {}.'.format(
            input_name,
            pipeline_id))

    def create_batch(self, folder_watcher):
        batch_name = folder_watcher.run_name + ' v' + self.config.pipeline_version
        description = 'MiCall batch for folder {}, pipeline version {}.'.format(
            folder_watcher.run_name,
            self.config.pipeline_version)
        batch = self.session.create_run_batch(batch_name,
                                              description=description,
                                              users=[],
                                              groups=ALLOWED_GROUPS)
        folder_watcher.batch = batch

    def find_kive_dataset(self, source_file, dataset_name, cdt):
        """ Search for a dataset in Kive by name and checksum.

        :param source_file: open file object to read from
        :param str dataset_name: dataset name to search for
        :param cdt: CompoundDatatype object returned by Kive API
        :return: the dataset object from the Kive API wrapper, or None
        """
        chunk_size = 4096
        digest = hashlib.md5()
        for chunk in iter(lambda: source_file.read(chunk_size), b""):
            digest.update(chunk)
        checksum = digest.hexdigest()
        datasets = self.session.find_datasets(name=dataset_name,
                                              md5=checksum,
                                              cdt=cdt,
                                              uploaded=True)
        needed_groups = set(ALLOWED_GROUPS)
        for dataset in datasets:
            missing_groups = needed_groups - set(dataset.groups_allowed)
            if dataset.name == dataset_name and not missing_groups:
                logger.info('dataset already in Kive: %r', dataset_name)
                return dataset
        return None

    def upload_kive_dataset(self, source_file, dataset_name, cdt, description):
        """ Upload a dataset to Kive.

        :param source_file: open file object to read from
        :param str dataset_name:
        :param cdt: CompoundDatatype object returned by Kive API
        :param str description:
        :return: the dataset object from the Kive API wrapper, or None
        """
        logger.info('uploading dataset %r', dataset_name)
        filepath = Path(getattr(source_file, 'name', ''))
        if (self.external_directory_name is None or
                self.external_directory_path not in filepath.parents):
            dataset = self.session.add_dataset(
                name=dataset_name,
                description=description,
                handle=source_file,
                cdt=cdt,
                groups=ALLOWED_GROUPS)
        else:
            external_path = os.path.relpath(filepath,
                                            self.external_directory_path)
            dataset = self.session.add_dataset(
                name=dataset_name,
                description=description,
                handle=None,
                externalfiledirectory=self.external_directory_name,
                external_path=external_path,
                cdt=cdt,
                groups=ALLOWED_GROUPS)
        return dataset

    def add_sample_group(self, base_calls, sample_group):
        self.check_session()
        folder_watcher = self.folder_watchers.get(base_calls)
        if folder_watcher is None:
            folder_watcher = self.add_folder(base_calls)
            self.create_batch(folder_watcher)
            self.upload_filter_quality(folder_watcher)
            shutil.rmtree(self.get_results_path(folder_watcher),
                          ignore_errors=True)

        for sample_watcher in folder_watcher.sample_watchers:
            if sample_watcher.sample_group == sample_group:
                return sample_watcher

        sample_watcher = SampleWatcher(sample_group)
        for fastq1 in filter(None, sample_group.names):
            fastq2 = fastq1.replace('_R1_', '_R2_')
            for fastq_name, direction in ((fastq1, 'forward'), (fastq2, 'reverse')):
                with (base_calls / fastq_name).open('rb') as fastq_file:
                    fastq_dataset = self.find_or_upload_dataset(
                        fastq_file,
                        fastq_name,
                        direction + ' read from MiSeq run ' + folder_watcher.run_name,
                        compounddatatype=None)
                    sample_watcher.fastq_datasets.append(fastq_dataset)
        folder_watcher.sample_watchers.append(sample_watcher)
        return sample_watcher

    def add_folder(self, base_calls):
        folder_watcher = FolderWatcher(base_calls, self)
        self.folder_watchers[base_calls] = folder_watcher
        return folder_watcher

    def finish_folder(self, base_calls):
        """ Record that all samples have been loaded for a folder.

        Processing isn't finished yet, but all samples have been loaded.
        """
        self.loaded_folders.add(base_calls)

    def poll_runs(self):
        completed_folders = []
        for folder, folder_watcher in self.folder_watchers.items():
            folder_watcher.poll_runs()
            if folder in self.loaded_folders and folder_watcher.is_complete:
                completed_folders.append(folder)
        for folder in completed_folders:
            folder_watcher = self.folder_watchers.pop(folder)
            results_path = self.collate_folder(folder_watcher)
            if (results_path / "coverage_scores.csv").exists():
                self.result_handler(results_path)
            (results_path / "doneprocessing").touch()
            if not self.folder_watchers:
                logger.info('No more folders to process.')

    def collate_folder(self, folder_watcher):
        results_path = self.get_results_path(folder_watcher)
        logger.info('Collating results in %s', results_path)
        scratch_path = results_path / "scratch"
        for output_name in DOWNLOADED_RESULTS:
            if output_name == 'coverage_maps_tar':
                self.extract_coverage_maps(folder_watcher)
                continue
            source_count = 0
            filename = get_output_filename(output_name)
            target_path = results_path / filename
            with target_path.open('w') as target:
                for sample_name in folder_watcher.all_samples:
                    sample_name = trim_name(sample_name)
                    source_path = scratch_path / sample_name / filename
                    try:
                        with source_path.open() as source:
                            for i, line in enumerate(source):
                                if i != 0:
                                    prefix = sample_name
                                elif source_count == 0:
                                    prefix = 'sample'
                                else:
                                    continue
                                target.write(prefix + ',' + line)
                            source_count += 1
                    except FileNotFoundError:
                        # Skip the file.
                        pass
            if not source_count:
                target_path.unlink()
        shutil.rmtree(scratch_path)
        return results_path

    def extract_coverage_maps(self, folder_watcher):
        results_path = self.get_results_path(folder_watcher)
        coverage_path = results_path / "coverage_maps"
        coverage_path.mkdir()
        scratch_path = results_path / "scratch"
        for sample_name in folder_watcher.all_samples:
            sample_name = trim_name(sample_name)
            source_path = scratch_path / sample_name / 'coverage_maps.tar'
            try:
                with tarfile.open(source_path) as f:
                    for source_info in f:
                        filename = os.path.basename(source_info.name)
                        target_path = coverage_path / (sample_name + '.' + filename)
                        with f.extractfile(source_info) as source, \
                                open(target_path, 'wb') as target:
                            shutil.copyfileobj(source, target)
            except FileNotFoundError:
                pass

    def run_pipeline(self, folder_watcher, pipeline_type, sample_watcher):
        if pipeline_type == PipelineType.FILTER_QUALITY:
            quality_pipeline = self.get_kive_pipeline(
                self.config.micall_filter_quality_pipeline_id)
            return self.session.run_pipeline(
                quality_pipeline,
                [folder_watcher.quality_dataset],
                'MiCall filter quality on ' + folder_watcher.run_name,
                runbatch=folder_watcher.batch,
                groups=ALLOWED_GROUPS)
        if pipeline_type == PipelineType.RESISTANCE:
            resistance_pipeline = self.get_kive_pipeline(
                self.config.micall_resistance_pipeline_id)
            main_runs = filter(None,
                               (sample_watcher.runs.get(pipeline_type)
                                for pipeline_type in (PipelineType.MAIN,
                                                      PipelineType.MIDI)))
            input_datasets = [run.get_results()['amino_csv']
                              for run in main_runs]
            main_aminos_input = self.get_kive_input('main_amino_csv')
            for input_dataset in input_datasets:
                # TODO: remove this when Kive API sets CDT properly (issue #729)
                input_dataset.cdt = main_aminos_input.compounddatatype
            if len(input_datasets) == 1:
                input_datasets *= 2
            return self.session.run_pipeline(
                resistance_pipeline,
                input_datasets,
                'MiCall resistance on ' + sample_watcher.sample_group.enum,
                runbatch=folder_watcher.batch,
                groups=ALLOWED_GROUPS)
        if pipeline_type in (PipelineType.MIXED_HCV_MAIN,
                             PipelineType.MIXED_HCV_MIDI):
            if self.config.mixed_hcv_pipeline_id is None:
                return None
            mixed_hcv_pipeline = self.get_kive_pipeline(
                self.config.mixed_hcv_pipeline_id)
            if pipeline_type == PipelineType.MIXED_HCV_MAIN:
                input_datasets = sample_watcher.fastq_datasets[:2]
                sample_name = sample_watcher.sample_group.names[0]
            else:
                input_datasets = sample_watcher.fastq_datasets[2:]
                sample_name = sample_watcher.sample_group.names[1]
            sample_name = trim_name(sample_name)
            return self.session.run_pipeline(
                mixed_hcv_pipeline,
                input_datasets,
                'Mixed HCV on ' + trim_name(sample_name),
                runbatch=folder_watcher.batch,
                groups=ALLOWED_GROUPS)
        if pipeline_type == PipelineType.MAIN:
            fastq1, fastq2 = sample_watcher.fastq_datasets[:2]
            sample_name = sample_watcher.sample_group.names[0]
        else:
            assert pipeline_type == PipelineType.MIDI
            fastq1, fastq2 = sample_watcher.fastq_datasets[2:]
            sample_name = sample_watcher.sample_group.names[1]
        main_pipeline = self.get_kive_pipeline(
            self.config.micall_main_pipeline_id)
        if folder_watcher.bad_cycles_dataset is None:
            results = folder_watcher.filter_quality_run.get_results()
            folder_watcher.bad_cycles_dataset = results['bad_cycles_csv']
            bad_cycles_input = self.get_kive_input('bad_cycles_csv')
            # TODO: remove this when Kive API sets CDT properly (issue #729)
            folder_watcher.bad_cycles_dataset.cdt = bad_cycles_input.compounddatatype

        run_name = 'MiCall main on ' + trim_name(sample_name)
        return self.session.run_pipeline(
            main_pipeline,
            [fastq1, fastq2, folder_watcher.bad_cycles_dataset],
            run_name,
            runbatch=folder_watcher.batch,
            groups=ALLOWED_GROUPS)

    def fetch_run_status(self, run, folder_watcher, pipeline_type, sample_watcher):
        is_complete = run.is_complete()
        if is_complete and pipeline_type != PipelineType.FILTER_QUALITY:
            sample_name = (sample_watcher.sample_group.names[1]
                           if pipeline_type in (PipelineType.MIDI,
                                                PipelineType.MIXED_HCV_MIDI)
                           else sample_watcher.sample_group.names[0])
            results_path = self.get_results_path(folder_watcher)
            scratch_path = results_path / "scratch" / trim_name(sample_name)
            scratch_path.mkdir(parents=True, exist_ok=True)
            results = run.get_results()
            for output_name in DOWNLOADED_RESULTS:
                dataset = results.get(output_name)
                if dataset is None:
                    continue
                filename = get_output_filename(output_name)
                with (scratch_path / filename).open('wb') as f:
                    dataset.download(f)

        return is_complete

    def get_results_path(self, folder_watcher):
        version_name = f'version_{self.config.pipeline_version}'
        results_path = folder_watcher.run_folder / "Results" / version_name
        return results_path

    def upload_filter_quality(self, folder_watcher):
        read_sizes = parse_read_sizes(folder_watcher.run_folder / "RunInfo.xml")
        read_lengths = [read_sizes.read1,
                        read_sizes.index1,
                        read_sizes.index2,
                        read_sizes.read2]
        quality_input = self.get_kive_input('quality_csv')
        error_path = folder_watcher.run_folder / "InterOp/ErrorMetricsOut.bin"
        quality_csv = StringIO()
        with error_path.open('rb') as error_file:
            records = error_metrics_parser.read_errors(error_file)
            error_metrics_parser.write_phix_csv(quality_csv,
                                                records,
                                                read_lengths)
        quality_csv_bytes = BytesIO()
        quality_csv_bytes.write(quality_csv.getvalue().encode('utf8'))
        quality_csv_bytes.seek(0)
        folder_watcher.quality_dataset = self.find_or_upload_dataset(
            quality_csv_bytes,
            folder_watcher.run_name + '_quality.csv',
            'Error rates for {} run.'.format(folder_watcher.run_name),
            quality_input.compounddatatype)

    def check_session(self):
        if self.session is None:
            self.session = open_kive(self.config.kive_server)
            self.session.login(self.config.kive_user, self.config.kive_password)

            # retrieve external file directory
            directories = self.session.get('/api/externalfiledirectories',
                                           is_json=True).json()
            self.external_directory_name = self.external_directory_path = None
            for directory in directories:
                directory_path = Path(directory['path'])
                if (directory_path == self.config.raw_data or
                        directory_path in self.config.raw_data.parents):
                    self.external_directory_name = directory['name']
                    self.external_directory_path = directory_path
                    break

    def find_or_upload_dataset(self,
                               dataset_file,
                               dataset_name,
                               description,
                               compounddatatype):
        dataset = self.find_kive_dataset(dataset_file,
                                         dataset_name,
                                         compounddatatype)
        if dataset is None:
            dataset_file.seek(0)
            dataset = self.upload_kive_dataset(
                dataset_file,
                dataset_name,
                compounddatatype,
                description)
        return dataset
