import hashlib
import logging
import os
import re
import shutil
import tarfile
from collections import namedtuple
from datetime import datetime, timedelta
from enum import Enum
from itertools import count
from operator import itemgetter
from pathlib import Path
from queue import Full

from io import StringIO, BytesIO
from time import sleep

from requests.adapters import HTTPAdapter
from kiveapi import KiveAPI, KiveClientException, KiveRunFailedException

from micall.drivers.run_info import parse_read_sizes
from micall.monitor import error_metrics_parser
from micall.monitor.sample_watcher import FolderWatcher, ALLOWED_GROUPS, SampleWatcher, PipelineType
from micall.monitor.find_groups import find_groups

logger = logging.getLogger(__name__)
FOLDER_SCAN_INTERVAL = timedelta(hours=1)
SLEEP_SECONDS = 60
MINIMUM_RETRY_WAIT = timedelta(seconds=5)
MAXIMUM_RETRY_WAIT = timedelta(days=1)
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
                      'mutations_csv',
                      'resistance_fail_csv',
                      'resistance_consensus_csv']

# noinspection PyArgumentList
FolderEventType = Enum('FolderEventType', 'ADD_SAMPLE FINISH_FOLDER')
FolderEvent = namedtuple('FolderEvent', 'base_calls type sample_group')


def open_kive(server_url):
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
        try:
            is_complete = scan_samples(raw_data_folder,
                                       pipeline_version,
                                       sample_queue,
                                       wait)
            if is_complete and not wait:
                break
        except Exception:
            logger.error("Failed while finding samples.", exc_info=True)
            raise


def scan_samples(raw_data_folder, pipeline_version, sample_queue, wait):
    next_scan = now() + FOLDER_SCAN_INTERVAL
    flag_paths = sorted(
        raw_data_folder.glob("MiSeq/runs/*/needsprocessing"),
        reverse=True)
    is_found = False
    for flag_path in flag_paths:
        run_path = flag_path.parent
        error_path = run_path / "errorprocessing"
        if error_path.exists():
            continue
        done_path = (run_path /
                     f"Results/version_{pipeline_version}/doneprocessing")
        if done_path.exists():
            continue
        base_calls_path = run_path / "Data/Intensities/BaseCalls"
        sample_groups = find_sample_groups(run_path, base_calls_path)
        for sample_group in sample_groups:
            is_found = True
            is_sent = send_event(sample_queue,
                                 FolderEvent(base_calls_path,
                                             FolderEventType.ADD_SAMPLE,
                                             sample_group),
                                 next_scan)
            if not is_sent:
                return False
        if sample_groups:
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


def find_sample_groups(run_path, base_calls_path):
    # noinspection PyBroadException
    try:
        fastq_files = base_calls_path.glob("*_R1_*.fastq.gz")
        sample_sheet_path = run_path / "SampleSheet.csv"
        file_names = [f.name for f in fastq_files]
        sample_groups = list(find_groups(file_names, sample_sheet_path))
        sample_groups.sort(key=lambda group: get_sample_number(group.names[0]),
                           reverse=True)
    except Exception:
        logger.error("Finding sample groups in %s", run_path, exc_info=True)
        (run_path / "errorprocessing").write_text(
            "Finding sample groups failed.\n")
        sample_groups = []
    return sample_groups


def get_sample_number(fastq_name):
    match = re.match(r'.*_S(\d+)_', fastq_name)
    return int(match.group(1))


def send_event(sample_queue, folder_event, next_scan):
    is_sent = False
    while not is_sent and now() < next_scan:
        try:
            sample_queue.put(folder_event,
                             timeout=SLEEP_SECONDS)
            is_sent = True
        except Full:
            pass
    return is_sent


def trim_name(sample_name):
    return '_'.join(sample_name.split('_')[:2])


def get_output_filename(output_name):
    return '.'.join(output_name.rsplit('_', 1))


def wait_for_retry(attempt_count):
    delay = calculate_retry_wait(MINIMUM_RETRY_WAIT,
                                 MAXIMUM_RETRY_WAIT,
                                 attempt_count)
    logger.error('Waiting %s before retrying.', delay, exc_info=True)
    sleep(delay.total_seconds())


def calculate_retry_wait(min_wait, max_wait, attempt_count):
    min_seconds = int(min_wait.total_seconds())
    seconds = min_seconds * (2 ** (attempt_count - 1))
    seconds = min(seconds, max_wait.total_seconds())
    return timedelta(seconds=seconds)


class KiveWatcher:
    def __init__(self,
                 config=None,
                 result_handler=lambda result_folder: None,
                 retry=False):
        """ Initialize.

        :param config: command line arguments
        :param result_handler: called when a run folder has collated all the
            results into a result folder
        :param bool retry: should the main methods retry forever?
        """
        self.config = config
        self.result_handler = result_handler
        self.retry = retry
        self.session = None
        self.loaded_folders = set()  # base_calls folders with all samples loaded
        self.folder_watchers = {}  # {base_calls_folder: FolderWatcher}
        self.pipelines = {}  # {pipeline_id: pipeline}
        self.external_directory_path = self.external_directory_name = None
        self.find_default_pipelines()
        self.input_pipeline_ids = config and dict(
            quality_csv=config.micall_filter_quality_pipeline_id,
            bad_cycles_csv=config.micall_main_pipeline_id,
            main_amino_csv=config.micall_resistance_pipeline_id,
            midi_amino_csv=config.micall_resistance_pipeline_id) or {}

        # Active runs started by other users.
        self.other_runs = None  # {(pipeline_id, dataset_id, dataset_id, ...): run}

    def find_default_pipelines(self):
        default_family_names = ['filter quality',
                                'main',
                                'resistance']
        family_data = None
        for family_name in default_family_names:
            attribute_name = (
                    'micall_' + family_name.replace(' ', '_') + '_pipeline_id')
            if getattr(self.config, attribute_name) is not None:
                continue
            if family_data is None:
                self.check_session()
                filter_text = 'filters[0][key]=smart&filters[0][val]=micall'
                family_data = self.session.get(
                    '@api_pipeline_families',
                    context={'filters': filter_text}).json()

            search_text = 'micall ' + family_name
            for family in family_data:
                if search_text in family['name'].lower():
                    for member in family['members']:
                        setattr(self.config, attribute_name, member['id'])
                        break
                    break
            else:
                raise RuntimeError(f'Argument {attribute_name} not set, and '
                                   f'no pipeline found named {search_text!r}.')

    def is_full(self):
        active_count = sum(len(folder_watcher.active_samples)
                           for folder_watcher in self.folder_watchers.values())
        return active_count >= self.config.max_active

    def is_idle(self):
        return not self.folder_watchers

    def get_kive_pipeline(self, pipeline_id):
        self.check_session()
        kive_pipeline = self.pipelines.get(pipeline_id)
        if kive_pipeline is None:
            kive_pipeline = self.kive_retry(
                lambda: self.session.get_pipeline(pipeline_id))
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
        batch = self.kive_retry(
            lambda: self.session.create_run_batch(batch_name,
                                                  description=description,
                                                  users=[],
                                                  groups=ALLOWED_GROUPS))
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
        datasets = self.kive_retry(
            lambda: self.session.find_datasets(name=dataset_name,
                                               md5=checksum,
                                               cdt=cdt,
                                               uploaded=True))
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
        """ Add a sample group (main and optional midi sample) to process.

        Also checks to see whether the folder finished processing since the
        last folder scan.
        :param base_calls: path to the BaseCalls folder with FASTQ files in it
        :param SampleGroup sample_group: the sample(s) to add
        :return: SampleWatcher for the sample group, or None if that folder has
            already finished processing
        """
        for attempt_count in count(1):
            # noinspection PyBroadException
            try:
                self.check_session()
                folder_watcher = self.folder_watchers.get(base_calls)
                if folder_watcher is None:
                    folder_watcher = FolderWatcher(base_calls, self)

                    # Check if folder has finished since it was scanned.
                    results_path = self.get_results_path(folder_watcher)
                    done_path = results_path / "doneprocessing"
                    if done_path.exists():
                        return None
                    error_path = folder_watcher.run_folder / "errorprocessing"
                    if error_path.exists():
                        return None

                    self.create_batch(folder_watcher)
                    self.upload_filter_quality(folder_watcher)
                    shutil.rmtree(results_path, ignore_errors=True)
                    self.folder_watchers[base_calls] = folder_watcher

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
            except Exception:
                if not self.retry:
                    raise
                wait_for_retry(attempt_count)

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
        for attempt_count in count(1):
            # noinspection PyBroadException
            try:
                self.check_session()
                self.find_other_runs()
                completed_folders = []
                for folder, folder_watcher in self.folder_watchers.items():
                    folder_watcher.poll_runs()
                    if folder in self.loaded_folders and folder_watcher.is_complete:
                        completed_folders.append(folder)
                for folder in completed_folders:
                    folder_watcher = self.folder_watchers.pop(folder)
                    results_path = self.collate_folder(folder_watcher)
                    if results_path is None:
                        continue
                    if (results_path / "coverage_scores.csv").exists():
                        self.result_handler(results_path)
                    (results_path / "doneprocessing").touch()
                    if not self.folder_watchers:
                        logger.info('No more folders to process.')
                return
            except Exception:
                if not self.retry:
                    raise
                self.other_runs = None  # Scan again.
                wait_for_retry(attempt_count)

    def find_other_runs(self):
        """ Check for runs started by other users. """
        if self.other_runs is not None:
            # Already checked.
            return
        runs = self.kive_retry(lambda: self.session.find_runs(active=True))
        run_map = {}
        pipeline_ids = set(self.input_pipeline_ids.values())
        for run in runs:
            if run.pipeline_id not in pipeline_ids:
                continue
            try:
                if self.kive_retry(run.is_complete):
                    outputs = self.kive_retry(run.get_results)
                    if any(output.dataset_id is None
                           for output in outputs.values()):
                        # Output has been purged, can't reuse the run.
                        continue
                input_list = run.raw['inputs']
                inputs = sorted(input_list, key=itemgetter('index'))
                input_ids = tuple(inp['dataset'] for inp in inputs)
                key = (run.pipeline_id, ) + input_ids
                run_map[key] = run
            except KiveRunFailedException:
                # Failed or cancelled, rerun.
                pass
        self.other_runs = run_map

    def collate_folder(self, folder_watcher):
        """ Collate scratch files for a run folder.

        :param FolderWatcher folder_watcher: holds details about the run folder
        """
        results_path = self.get_results_path(folder_watcher)
        failed_sample_names = [
            sample_watcher.sample_group.enum
            for sample_watcher in folder_watcher.sample_watchers
            if sample_watcher.is_failed]
        if failed_sample_names:
            run_path = (results_path / "../..").resolve()
            error_message = 'Samples failed in Kive: {}.'.format(
                ', '.join(failed_sample_names))
            (run_path / 'errorprocessing').write_text(error_message + '\n')
            logger.error('Error in folder %s: %s', run_path, error_message)
            return
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
            return self.find_or_launch_run(
                self.config.micall_filter_quality_pipeline_id,
                [folder_watcher.quality_dataset],
                'MiCall filter quality on ' + folder_watcher.run_name,
                folder_watcher.batch)
        if pipeline_type == PipelineType.RESISTANCE:
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
            return self.find_or_launch_run(
                self.config.micall_resistance_pipeline_id,
                input_datasets,
                'MiCall resistance on ' + sample_watcher.sample_group.enum,
                folder_watcher.batch)
        if pipeline_type in (PipelineType.MIXED_HCV_MAIN,
                             PipelineType.MIXED_HCV_MIDI):
            if self.config.mixed_hcv_pipeline_id is None:
                return None
            if pipeline_type == PipelineType.MIXED_HCV_MAIN:
                input_datasets = sample_watcher.fastq_datasets[:2]
                sample_name = sample_watcher.sample_group.names[0]
            else:
                input_datasets = sample_watcher.fastq_datasets[2:]
                sample_name = sample_watcher.sample_group.names[1]
            sample_name = trim_name(sample_name)
            return self.find_or_launch_run(
                self.config.mixed_hcv_pipeline_id,
                input_datasets,
                'Mixed HCV on ' + trim_name(sample_name),
                folder_watcher.batch)
        if pipeline_type == PipelineType.MAIN:
            fastq1, fastq2 = sample_watcher.fastq_datasets[:2]
            sample_name = sample_watcher.sample_group.names[0]
        else:
            assert pipeline_type == PipelineType.MIDI
            fastq1, fastq2 = sample_watcher.fastq_datasets[2:]
            sample_name = sample_watcher.sample_group.names[1]
        if folder_watcher.bad_cycles_dataset is None:
            results = folder_watcher.filter_quality_run.get_results()
            folder_watcher.bad_cycles_dataset = results['bad_cycles_csv']
            bad_cycles_input = self.get_kive_input('bad_cycles_csv')
            # TODO: remove this when Kive API sets CDT properly (issue #729)
            folder_watcher.bad_cycles_dataset.cdt = bad_cycles_input.compounddatatype

        run_name = 'MiCall main on ' + trim_name(sample_name)
        return self.find_or_launch_run(
            self.config.micall_main_pipeline_id,
            [fastq1, fastq2, folder_watcher.bad_cycles_dataset],
            run_name,
            folder_watcher.batch)

    def find_or_launch_run(self,
                           pipeline_id,
                           inputs,
                           run_name,
                           run_batch):
        run_key = (pipeline_id,) + tuple(inp.dataset_id for inp in inputs)
        if self.other_runs:
            run = self.other_runs.pop(run_key, None)
            if run is not None:
                return run
        pipeline = self.get_kive_pipeline(pipeline_id)
        return self.session.run_pipeline(pipeline,
                                         inputs,
                                         run_name,
                                         runbatch=run_batch,
                                         groups=ALLOWED_GROUPS)

    def kive_retry(self, target):
        """ Add a single retry to a Kive API call.

        Tries to call the target function, then refreshes the session login if
        the call fails and tries a second time.
        """
        try:
            return target()
        except KiveClientException:
            self.session.login(self.config.kive_user, self.config.kive_password)
            return target()

    def fetch_run_status(self, run, folder_watcher, pipeline_type, sample_watcher):
        self.check_session()
        try:
            is_complete = self.kive_retry(lambda: run.is_complete())
        except KiveRunFailedException as ex:
            if 'fail' in ex.args[0]:
                raise
            return self.run_pipeline(folder_watcher, pipeline_type, sample_watcher)
        if is_complete and pipeline_type != PipelineType.FILTER_QUALITY:
            sample_name = (sample_watcher.sample_group.names[1]
                           if pipeline_type in (PipelineType.MIDI,
                                                PipelineType.MIXED_HCV_MIDI)
                           else sample_watcher.sample_group.names[0])
            results_path = self.get_results_path(folder_watcher)
            scratch_path = results_path / "scratch" / trim_name(sample_name)
            scratch_path.mkdir(parents=True, exist_ok=True)
            results = self.kive_retry(run.get_results)
            for output_name in DOWNLOADED_RESULTS:
                dataset = results.get(output_name)
                if dataset is None:
                    continue
                filename = get_output_filename(output_name)
                with (scratch_path / filename).open('wb') as f:
                    self.kive_retry(lambda: dataset.download(f))

        if is_complete:
            return None
        return run

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
