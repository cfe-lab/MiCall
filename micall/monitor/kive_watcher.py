import hashlib
import logging
import os
import re
import shutil
import tarfile
from collections import namedtuple
from csv import DictWriter, DictReader
from datetime import datetime, timedelta
from enum import Enum
from itertools import count
from pathlib import Path
from queue import Full, Queue
from typing import IO, Callable, Mapping, Optional, Sequence, Iterable, TextIO, Tuple, TypeVar

from io import StringIO, BytesIO
from time import sleep
from zipfile import ZipFile, ZIP_DEFLATED

# noinspection PyPackageRequirements
import kiveapi
from requests.adapters import HTTPAdapter
from kiveapi import KiveAPI, KiveClientException, KiveRunFailedException
import urllib3

from micall.drivers.run_info import parse_read_sizes
from micall.monitor.sample_watcher import Batch, FolderWatcher, ALLOWED_GROUPS, Item, SampleWatcher, PipelineType, PIPELINE_GROUPS, Run, RunDataset, RunCreationDataset, ConfigInterface
from micall.monitor.find_groups import SampleGroup, find_groups
from micall.monitor import disk_operations
from micall.utils.check_sample_sheet import check_sample_name_consistency
from micall.utils.list_fastq_files import find_fastq_source_folder, list_fastq_file_names
from miseqinteropreader.error_metrics_parser import write_phix_csv
from miseqinteropreader.interop_reader import InterOpReader
from miseqinteropreader import ReadLengths4

logger = logging.getLogger(__name__)
FOLDER_SCAN_INTERVAL = timedelta(hours=1)
SLEEP_SECONDS = 60
MINIMUM_RETRY_WAIT = timedelta(seconds=5)
MAXIMUM_RETRY_WAIT = timedelta(days=1)
MAX_RUN_NAME_LENGTH = 60
DOWNLOADED_RESULTS = ['remap_counts_csv',
                      'conseq_csv',
                      'conseq_all_csv',
                      'conseq_stitched_csv',
                      'conseq_region_csv',
                      'concordance_csv',
                      'concordance_seed_csv',
                      'insertions_csv',
                      'failed_csv',
                      'nuc_csv',
                      'amino_csv',
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
                      'nuc_mutations_csv',
                      'resistance_fail_csv',
                      'resistance_consensus_csv',
                      'wg_fasta',
                      'mid_fasta',
                      'unstitched_cascade_csv',
                      'unstitched_conseq_csv',
                      'unstitched_contigs_csv',
                      'contigs_csv',
                      'alignment_svg',
                      'alignment_png',
                      'assembly_fasta',
                      'genome_coverage_csv',
                      'genome_coverage_svg',
                      'genome_concordance_svg',
                      'read_entropy_csv',
                      'outcome_summary_csv',  # proviral outputs
                      'conseqs_primers_csv',
                      'contigs_primers_csv',
                      'table_precursor_csv',
                      'proviral_landscape_csv',
                      'hivseqinr_results_tar',  # TODO: remove when proviral is finally updated.
                      'detailed_results_tar']

# noinspection PyArgumentList
FolderEventType = Enum('FolderEventType', 'ADD_SAMPLE FINISH_FOLDER')
FolderEvent = namedtuple('FolderEvent', 'base_calls type sample_group')

kiveapi.kiveapi.logger.setLevel(logging.ERROR)  # Suppress routine credential refresh noise

# Suppress urllib3 connection retry warnings
urllib3.connectionpool.log.setLevel(logging.ERROR)

T = TypeVar('T')

def open_kive(server_url: str) -> KiveAPI:
    session = KiveAPI(server_url)
    session.mount('https://', HTTPAdapter(max_retries=20))
    return session


def now() -> datetime:
    """ Get the current date/time.

    Wrapped in a local function so it can be mocked during testing.
    """
    return datetime.now()


def find_samples(raw_data_folder: Path,
                 pipeline_version: str,
                 sample_queue: Queue[FolderEvent],
                 qai_upload_queue: Queue[Tuple[Path, PipelineType]],
                 wait: bool = True,
                 retry: bool = True) -> None:
    attempt_count = 0
    start_time = None
    while True:
        # noinspection PyBroadException
        try:
            is_complete = scan_samples(raw_data_folder,
                                       pipeline_version,
                                       sample_queue,
                                       wait)

            # After scanning for new samples, enqueue QAI upload tasks for finished runs.
            qai_next_scan = now() + FOLDER_SCAN_INTERVAL
            for version_folder in scan_qai_done_flags(raw_data_folder, pipeline_version):
                is_sent = send_qai_event(qai_upload_queue,
                                         (version_folder, PipelineType.MAIN),
                                         qai_next_scan)
                if not is_sent:
                    # Let the caller loop and try again later
                    is_complete = False
                    break

            attempt_count = 0  # Reset after success
            start_time = None  # Reset start time after success
            if is_complete and not wait:
                break
        except Exception:
            if not retry:
                raise
            attempt_count += 1

            # Record start time on first failure
            if start_time is None:
                start_time = datetime.now()

            wait_for_retry(attempt_count, start_time)


def get_version_key(version_path: Path) -> tuple[int, ...]:
    version_name = version_path.name
    version_parts = tuple(int(part) for part in re.findall(r'\d+', version_name))
    return version_parts


def compress_old_versions(version_path: Path) -> None:
    results_path: Path = version_path.parent
    if not results_path.exists():
        return

    version_key = get_version_key(version_path)
    version_folders = {}
    for other_version_path in results_path.iterdir():
        if not other_version_path.is_dir():
            continue
        if not other_version_path.name.startswith('version_'):
            continue
        other_version_key = get_version_key(other_version_path)
        if other_version_key < version_key:
            version_folders[other_version_key] = other_version_path
    if not version_folders:
        return
    version_items = sorted(version_folders.items())
    version_items.pop()
    for version_key, version_folder in version_items:
        # noinspection PyUnresolvedReferences
        zip_name = version_folder.name + '.zip'
        # noinspection PyTypeChecker
        with open(results_path / zip_name, 'wb') as f:
            with ZipFile(f, 'w', ZIP_DEFLATED) as z:
                # noinspection PyTypeChecker
                for folder_path_str, folders, files in os.walk(version_folder):
                    # noinspection PyTypeChecker
                    folder_path = Path(folder_path_str)
                    folder_rel = folder_path.relative_to(results_path)
                    for file_name in files:
                        z.write(folder_path/file_name, folder_rel/file_name)
        shutil.rmtree(version_folder)


def scan_samples(raw_data_folder: Path, pipeline_version: str, sample_queue: Queue[FolderEvent], wait: bool):
    next_scan = now() + FOLDER_SCAN_INTERVAL
    flag_paths = sorted(scan_flag_paths(raw_data_folder), reverse=True)
    is_found = False
    for flag_path in flag_paths:
        run_path = flag_path.parent
        error_path = run_path / "errorprocessing"
        if error_path.exists():
            continue
        done_path = (run_path /
                     f"Results/version_{pipeline_version}/done_all_processing")
        if done_path.exists():
            continue
        fastq_folder = find_fastq_source_folder(run_path)
        if fastq_folder is None:
            logger.warning("No FASTQ files found in %s", run_path)
            disk_operations.write_text(error_path,
                                       "No FASTQ files found.\n")
            continue
        sample_groups = find_sample_groups(run_path, fastq_folder)
        for sample_group in sample_groups:
            is_found = True
            is_sent = send_event(sample_queue,
                                 FolderEvent(fastq_folder,
                                             FolderEventType.ADD_SAMPLE,
                                             sample_group),
                                 next_scan)
            if not is_sent:
                return False
        if sample_groups:
            is_sent = send_event(sample_queue,
                                 FolderEvent(fastq_folder,
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

def scan_flag_paths(raw_data_folder: Path) -> Iterable[Path]:
    return raw_data_folder.glob("MiSeq/runs/*/needsprocessing")


def scan_qai_done_flags(raw_data_folder: Path, pipeline_version: str) -> Iterable[Path]:
    """Find version folders with done_all_processing set and not yet uploaded to QAI.

    We consider a version folder eligible if:
    - Results/version_{pipeline_version}/done_all_processing exists, and
    - Results/version_{pipeline_version}/done_qai_upload does NOT exist.
    """
    pattern = f"MiSeq/runs/*/Results/version_{pipeline_version}/done_all_processing"
    for done_flag in raw_data_folder.glob(pattern):
        version_folder = done_flag.parent
        if not (version_folder / 'done_qai_upload').exists():
            yield version_folder


def find_sample_groups(run_path: Path, base_calls_path: Path) -> Sequence[SampleGroup]:
    # noinspection PyBroadException
    try:
        file_names = list_fastq_file_names(run_path, "*_R1_*.fastq.gz", fallback_to_run_path=False)
        sample_sheet_path = run_path / "SampleSheet.csv"

        # Check sample name consistency between sample sheet and FASTQ files
        check_sample_name_consistency(sample_sheet_path, file_names, run_path)

        sample_groups = list(find_groups(file_names, str(sample_sheet_path)))
        sample_groups.sort(key=lambda group: get_sample_number(group.names[0] if group.names[0] is not None else ''),
                           reverse=True)
    except Exception:
        logger.error("Finding sample groups in %s", run_path, exc_info=True)
        disk_operations.write_text(run_path / "errorprocessing",
                                 "Finding sample groups failed.\n")
        sample_groups = []
    return sample_groups


def get_sample_number(fastq_name: str) -> int:
    match = re.match(r'.*_S(\d+)_', fastq_name)
    if match is None:
        raise ValueError(f"FASTQ filename does not match expected pattern: {fastq_name}")
    return int(match.group(1))


def send_event(sample_queue: Queue[FolderEvent], folder_event: FolderEvent, next_scan: datetime) -> bool:
    is_sent = False
    while not is_sent and now() < next_scan:
        try:
            sample_queue.put(folder_event,
                             timeout=SLEEP_SECONDS)
            is_sent = True
        except Full:
            pass
    return is_sent


def send_qai_event(qai_queue: Queue[Tuple[Path, PipelineType]], qai_event: Tuple[Path, PipelineType], next_scan: datetime) -> bool:
    """Send an event to the QAI upload queue with timeout and retry logic."""
    is_sent = False
    while not is_sent and now() < next_scan:
        try:
            qai_queue.put(qai_event, timeout=SLEEP_SECONDS)
            is_sent = True
        except Full:
            pass
    return is_sent


def trim_name(sample_name: str) -> str:
    return '_'.join(sample_name.split('_')[:2])


def trim_run_name(run_name: str) -> str:
    if len(run_name) > MAX_RUN_NAME_LENGTH:
        split_index = run_name.rfind('_')
        suffix_length = len(run_name) - split_index
        if split_index == -1 or suffix_length > 10:
            suffix_length = 0
            suffix = ''
        else:
            suffix = run_name[split_index:]

        run_name = (run_name[:MAX_RUN_NAME_LENGTH - suffix_length - 3] +
                    '...' +
                    suffix)

    return run_name


def get_output_filename(output_name: str) -> str:
    return '.'.join(output_name.rsplit('_', 1))


def wait_for_retry(attempt_count: int, start_time: datetime) -> None:
    """Wait with exponential backoff, logging warnings after 1 hour, info messages before."""
    delay = calculate_retry_wait(MINIMUM_RETRY_WAIT,
                                 MAXIMUM_RETRY_WAIT,
                                 attempt_count)

    # Determine if we should log based on elapsed time
    elapsed = datetime.now() - start_time
    should_log_warning = elapsed >= timedelta(hours=1)

    if should_log_warning:
        logger.warning('Waiting %s before retrying.', delay, exc_info=True)
    else:
        logger.info('Waiting %s before retrying.', delay, exc_info=True)
    sleep(delay.total_seconds())


def calculate_retry_wait(min_wait: timedelta, max_wait: timedelta, attempt_count: int) -> timedelta:
    min_seconds = int(min_wait.total_seconds())
    seconds = min_seconds * (2 ** (attempt_count - 1))
    seconds = min(seconds, max_wait.total_seconds())
    return timedelta(seconds=seconds)


def get_scratch_path(results_path: Path, pipeline_group: PipelineType) -> Path:
    if pipeline_group == PipelineType.MAIN:
        scratch_name = "scratch"
    elif pipeline_group == PipelineType.DENOVO_MAIN:
        scratch_name = "scratch_denovo"
    elif pipeline_group == PipelineType.PROVIRAL:
        scratch_name = "scratch_proviral"
    else:
        assert pipeline_group == PipelineType.MIXED_HCV_MAIN
        scratch_name = "scratch_mixed_hcv"
    scratch_path = results_path / scratch_name
    return scratch_path


def get_collated_path(results_path: Path, pipeline_group: PipelineType) -> Path:
    if pipeline_group == PipelineType.MAIN:
        target_path = results_path
    elif pipeline_group == PipelineType.DENOVO_MAIN:
        target_path = results_path / "denovo"
    elif pipeline_group == PipelineType.PROVIRAL:
        target_path = results_path / "proviral"
    else:
        assert pipeline_group == PipelineType.MIXED_HCV_MAIN
        target_path = results_path / "mixed_hcv"
    return target_path


class KiveWatcher:
    def __init__(self,
                 config: ConfigInterface,
                 qai_upload_queue: Queue[Tuple[Path, PipelineType]],
                 retry=False):
        """ Initialize.

        :param config: command line arguments
        :param qai_upload_queue: notified when a run folder has collated all the
            results into a result folder
        :param bool retry: should the main methods retry forever?
        """
        self.config = config
        self.retry = retry
        self._session = None
        self.loaded_folders: set[str] = set()  # base_calls folders with all samples loaded
        self.folder_watchers: dict[str, FolderWatcher] = {}  # {base_calls_folder: FolderWatcher}
        self.app_urls: dict[str | int, str] = {}  # {app_id: app_url}
        self.app_args: dict[str | int, dict[str, str]] = {}  # {app_id: {arg_name: arg_url}}
        self.external_directory_path: Optional[Path] = None
        self.external_directory_name: Optional[str] = None
        self.qai_upload_queue = qai_upload_queue

    def is_full(self) -> bool:
        if self.config is None:
            raise RuntimeError("KiveWatcher config is not set.")
        active_count = sum(folder_watcher.active_run_count
                           for folder_watcher in self.folder_watchers.values())
        return active_count >= self.config.max_active

    def is_idle(self) -> bool:
        return not self.folder_watchers

    def get_kive_app(self, app_id: str | int) -> str:
        self.get_kive_arguments(app_id)
        return self.app_urls[app_id]

    def get_kive_arguments(self, app_id: str | int) -> dict[str, str]:
        """ Get a dictionary of argument URL's for a container app. """
        self.check_session()
        kive_app = self.app_args.get(app_id)
        if kive_app is None:
            arguments = self.kive_retry(
                lambda: self.session.endpoints.containerapps.get(
                    f'{app_id}/argument_list/'))
            kive_app = {argument['name']: argument['url']
                        for argument in arguments
                        if argument['type'] == 'I'}
            self.app_args[app_id] = kive_app
            self.app_urls[app_id] = arguments[0]['app']
        return kive_app

    def get_kive_container_name(self, app_id: str | int) -> str:
        """ Get the container name for a container app. """
        self.check_session()
        kive_app = self.kive_retry(
            lambda: self.session.endpoints.containerapps.get(app_id)
        )
        return kive_app['container_name']

    def get_container_version(self, app_id: str | int) -> str:
        """ Get the container version for a container app. """
        full_name = self.get_kive_container_name(app_id)
        container_name, separator, version = full_name.rpartition(':')
        if not separator:
            # Weird case: no version found. Let's just return the full name.
            logger.warning('No version found in container name %r (id %r)', full_name, app_id)
            return full_name
        return version

    def create_batch(self, folder_watcher: FolderWatcher) -> None:
        if self.config is None:
            raise RuntimeError("KiveWatcher config is not set.")

        batch_name = folder_watcher.run_name + ' v' + self.config.pipeline_version
        description = 'MiCall batch for folder {}, pipeline version {}.'.format(
            folder_watcher.run_name,
            self.config.pipeline_version)
        old_batches = self.kive_retry(
                lambda: self.session.endpoints.batches.filter(
                    'name', batch_name))
        batch = self.find_name_and_permissions_match(old_batches,
                                                     batch_name,
                                                     'batch')
        if batch is None:
            batch = self.kive_retry(
                lambda: self.session.endpoints.batches.post(json=dict(
                    name=batch_name,
                    description=description,
                    groups_allowed=ALLOWED_GROUPS)))

        ret: Batch = batch  # type: ignore
        folder_watcher.batch = ret

    def find_kive_dataset(self, source_file: IO[bytes], dataset_name: str) -> Optional[RunDataset]:
        """ Search for a dataset in Kive by name and checksum.

        :param source_file: open file object to read from
        :param str dataset_name: dataset name to search for
        :return: the dataset object from the Kive API wrapper, or None
        """
        chunk_size = 4096
        digest = hashlib.md5()
        for chunk in iter(lambda: source_file.read(chunk_size), b""):
            digest.update(chunk)
        checksum = digest.hexdigest()
        datasets = self.kive_retry(
            lambda: self.session.endpoints.datasets.filter(
                'name', dataset_name,
                'md5', checksum,
                'uploaded', True))
        found = self.find_name_and_permissions_match(datasets,
                                                     dataset_name,
                                                    'dataset')
        if found is None:
            return None

        ret: RunDataset = found  # type: ignore
        return ret

    @staticmethod
    def find_name_and_permissions_match(items: Iterable[Item], name: str, type_name: str) -> Optional[Item]:
        needed_groups = set(ALLOWED_GROUPS)
        for item in items:
            missing_groups = needed_groups - set(item['groups_allowed'])
            if item['name'] == name and not missing_groups:
                logger.info('%s already in Kive: %r', type_name, name)
                return item
        return None

    def upload_kive_dataset(self, source_file, dataset_name, description):
        """ Upload a dataset to Kive.

        :param source_file: open file object to read from
        :param str dataset_name:
        :param str description:
        :return: the dataset object from the Kive API wrapper, or None
        """
        logger.info('uploading dataset %r', dataset_name)
        filepath = Path(getattr(source_file, 'name', ''))
        if (self.external_directory_name is None or
                self.external_directory_path not in filepath.parents):
            dataset = self.session.endpoints.datasets.post(
                files=dict(dataset_file=source_file),
                data=dict(name=dataset_name,
                          description=description,
                          users_allowed=[],
                          groups_allowed=ALLOWED_GROUPS))
        else:
            external_path = os.path.relpath(filepath,
                                            self.external_directory_path)
            dataset = self.session.endpoints.datasets.post(
                json=dict(name=dataset_name,
                          description=description,
                          externalfiledirectory=self.external_directory_name,
                          external_path=external_path,
                          users_allowed=[],
                          groups_allowed=ALLOWED_GROUPS))
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
        start_time = None
        for attempt_count in count(1):
            # noinspection PyBroadException
            try:
                self.check_session()
                folder_watcher = self.folder_watchers.get(base_calls)
                if folder_watcher is None:
                    folder_watcher = FolderWatcher(base_calls, self)

                    # Check if folder has finished since it was scanned.
                    results_path = self.get_results_path(folder_watcher)
                    zip_name = results_path.name + '.zip'
                    results_zip: Path = results_path.parent / zip_name
                    done_path = results_path / "doneprocessing"
                    if done_path.exists():
                        return None
                    error_path = folder_watcher.run_folder / "errorprocessing"
                    if error_path.exists():
                        return None

                    compress_old_versions(results_path)
                    self.create_batch(folder_watcher)
                    self.upload_filter_quality(folder_watcher)
                    if folder_watcher.quality_dataset is None:
                        return None
                    disk_operations.rmtree(results_path, ignore_errors=True)
                    disk_operations.unlink(results_zip, missing_ok=True)
                    self.folder_watchers[base_calls] = folder_watcher

                for sample_watcher in folder_watcher.sample_watchers:
                    if sample_watcher.sample_group == sample_group:
                        return sample_watcher

                sample_watcher = SampleWatcher(sample_group)
                names: Iterable[Optional[str]] = sample_group.names
                for fastq1 in filter(None, names):
                    fastq2 = fastq1.replace('_R1_', '_R2_')
                    for fastq_name, direction in ((fastq1, 'forward'), (fastq2, 'reverse')):
                        with (base_calls / fastq_name).open('rb') as fastq_file:
                            fastq_dataset = self.find_or_upload_dataset(
                                fastq_file,
                                fastq_name,
                                direction + ' read from MiSeq run ' +
                                folder_watcher.run_name)
                            sample_watcher.fastq_datasets.append(fastq_dataset)

                folder_watcher.sample_watchers.append(sample_watcher)
                return sample_watcher
            except Exception:
                if not self.retry:
                    raise

                # Record start time on first failure
                if start_time is None:
                    start_time = datetime.now()

                wait_for_retry(attempt_count, start_time)

    def add_folder(self, base_calls: str) -> FolderWatcher:
        folder_watcher = FolderWatcher(base_calls, self)
        self.folder_watchers[base_calls] = folder_watcher
        return folder_watcher

    def finish_folder(self, base_calls: str) -> None:
        """ Record that all samples have been loaded for a folder.

        Processing isn't finished yet, but all samples have been loaded.
        """
        self.loaded_folders.add(base_calls)

    def poll_runs(self) -> None:
        start_time = None
        for attempt_count in count(1):
            # noinspection PyBroadException
            try:
                self.check_session()
                poll_only_new_runs = any(
                    folder_watcher.has_new_runs
                    for folder_watcher in self.folder_watchers.values())
                for folder, folder_watcher in self.folder_watchers.items():
                    folder_watcher.poll_only_new_runs = poll_only_new_runs
                    folder_watcher.poll_runs()

                self.check_completed_folders()
                return
            except Exception:
                if not self.retry:
                    raise

                # Record start time on first failure
                if start_time is None:
                    start_time = datetime.now()

                wait_for_retry(attempt_count, start_time)

    def check_completed_folders(self) -> None:
        for folder, folder_watcher in list(self.folder_watchers.items()):
            if folder not in self.loaded_folders:
                # Still loading samples, can't be completed.
                continue
            for pipeline_group in list(folder_watcher.active_pipeline_groups):
                if not folder_watcher.is_pipeline_group_finished(pipeline_group):
                    continue
                results_path = self.collate_folder(folder_watcher,
                                                   pipeline_group)
                folder_watcher.active_pipeline_groups.remove(pipeline_group)
                if results_path is not None:
                    if (results_path / "coverage_scores.csv").exists():
                        self.qai_upload_queue.put(
                            (results_path, pipeline_group))
                    if not folder_watcher.active_pipeline_groups:
                        disk_operations.touch(results_path / "done_all_processing")
                        self.folder_watchers.pop(folder)
                if not self.folder_watchers:
                    logger.info('No more folders to process.')

    def collate_folder(self, folder_watcher: FolderWatcher, pipeline_group: PipelineType) -> Optional[Path]:
        """ Collate scratch files for a run folder.

        :param FolderWatcher folder_watcher: holds details about the run folder
        :param PipelineType pipeline_group: the group of runs to collate
        """
        results_path = self.get_results_path(folder_watcher)

        error_message = None
        if folder_watcher.is_folder_failed:
            error_message = 'Filter quality failed in Kive.'
        else:
            failed_sample_names = [
                sample_watcher.sample_group.enum
                for sample_watcher in folder_watcher.sample_watchers
                if sample_watcher.is_failed]
            if failed_sample_names:
                error_message = 'Samples failed in Kive: {}.'.format(
                    ', '.join(failed_sample_names))
        if error_message is not None:
            run_path = (results_path / "../..").resolve()
            disk_operations.write_text(run_path / 'errorprocessing', error_message + '\n')
            logger.error('Error in folder %s: %s', run_path, error_message)
            return None
        if pipeline_group == PipelineType.FILTER_QUALITY:
            return results_path
        scratch_path = get_scratch_path(results_path, pipeline_group)
        target_path = get_collated_path(results_path, pipeline_group)
        logger.info('Collating results in %s', target_path)
        self.copy_outputs(folder_watcher, scratch_path, target_path)
        disk_operations.rmtree(scratch_path)
        disk_operations.touch(target_path / 'doneprocessing')
        return results_path

    def copy_outputs(self,
                     folder_watcher: FolderWatcher,
                     scratch_path: Path,
                     results_path: Path) -> None:
        disk_operations.mkdir_p(results_path, exist_ok=True)
        for output_name in DOWNLOADED_RESULTS:
            if output_name == 'coverage_maps_tar':
                self.extract_coverage_maps(folder_watcher,
                                           scratch_path,
                                           results_path)
                continue
            if output_name.endswith('_tar'):
                self.extract_archive(folder_watcher,
                                     scratch_path,
                                     results_path,
                                     output_name)
                continue
            if output_name == 'alignment_svg':
                self.move_alignment_plot(folder_watcher,
                                         '.svg',
                                         scratch_path,
                                         results_path)
                continue
            if output_name == 'alignment_png':
                self.move_alignment_plot(folder_watcher,
                                         '.png',
                                         scratch_path,
                                         results_path)
                continue
            if output_name == 'genome_coverage_svg':
                self.move_genome_coverage(folder_watcher,
                                          scratch_path,
                                          results_path)
            source_count = 0
            filename = get_output_filename(output_name)
            target_path = results_path / filename
            with disk_operations.disk_file_operation(target_path, 'w') as target:
                for sample_name in folder_watcher.all_samples:
                    sample_name = trim_name(sample_name)
                    source_path = scratch_path / sample_name / filename
                    if not source_path.exists():
                        logger.debug('Source file %s does not exist, skipping.',
                                       source_path)
                        continue

                    with disk_operations.disk_file_operation(source_path, 'r') as source:
                        if output_name.endswith('_fasta'):
                            source_count += self.extract_fasta(source, target, sample_name)
                        else:
                            source_count += self.extract_csv(source,
                                                             target,
                                                             sample_name,
                                                             source_count)

            if not source_count:
                disk_operations.unlink(target_path)

    @staticmethod
    def extract_csv(source: TextIO, target: TextIO, sample_name: str, source_count: int) -> int:
        reader = DictReader(source)
        fieldnames = reader.fieldnames
        if fieldnames is None:
            # Empty file, nothing to copy.
            return 0
        fieldnames = list(fieldnames)
        has_sample = 'sample' in fieldnames
        if not has_sample:
            fieldnames.insert(0, 'sample')
        writer = DictWriter(target, fieldnames, lineterminator=os.linesep)
        if source_count == 0:
            # First source file, copy header.
            writer.writeheader()
        for row in reader:
            if not has_sample:
                row['sample'] = sample_name
            writer.writerow(row)
        return 1

    @staticmethod
    def extract_fasta(source: TextIO, target: TextIO, sample_name: str) -> int:
        for line in source:
            if line.startswith('>'):
                target.write(f'>{sample_name},{line[1:]}')
            else:
                target.write(line)
        return 1

    @staticmethod
    def extract_coverage_maps(folder_watcher: FolderWatcher, scratch_path: Path, results_path: Path) -> None:
        coverage_path: Path = results_path / "coverage_maps"
        disk_operations.mkdir_p(coverage_path, exist_ok=True)
        for sample_name in folder_watcher.all_samples:
            sample_name = trim_name(sample_name)
            source_path = scratch_path / sample_name / 'coverage_maps.tar'
            try:
                with tarfile.open(source_path) as f:
                    for source_info in f:
                        filename = os.path.basename(source_info.name)
                        target_path = coverage_path / (sample_name + '.' + filename)
                        source = f.extractfile(source_info)
                        if source is None:
                            raise RuntimeError(f"Failed to extract {source_info.name} from {source_path}")
                        with source, open(target_path, 'wb') as target:
                            shutil.copyfileobj(source, target)
            except FileNotFoundError:
                pass
        disk_operations.remove_empty_directory(coverage_path)

    @staticmethod
    def extract_archive(folder_watcher: FolderWatcher,
                        scratch_path: Path,
                        results_path: Path,
                        output_name: str) -> None:
        """ Extract contents of tar files.

        There will be a folder named after the output name, with a subfolder
        for each sample.
        :param folder_watcher: holds a list of all samples to extract from
        :param scratch_path: parent folder of all sample working files
        :param results_path: parent folder to extract into
        :param output_name: the name of tar files to look for with "_tar"
            instead of ".tar".
        """
        assert output_name.endswith('_tar'), output_name
        archive_name = output_name[:-4]
        output_path: Path = results_path / archive_name
        disk_operations.mkdir_p(output_path, exist_ok=True)
        for sample_name in folder_watcher.all_samples:
            sample_name = trim_name(sample_name)
            source_path = scratch_path / sample_name / (archive_name + '.tar')
            try:
                with tarfile.open(source_path) as f:
                    sample_target_path = output_path / sample_name
                    disk_operations.mkdir_p(sample_target_path, exist_ok=True)
                    for source_info in f:
                        filename = os.path.basename(source_info.name)
                        target_path = sample_target_path / filename
                        assert not target_path.exists(), target_path
                        source = f.extractfile(source_info)
                        if source is None:
                            raise RuntimeError(f"Failed to extract {source_info.name} from {source_path}")
                        with source, open(target_path, 'wb') as target:
                            shutil.copyfileobj(source, target)
                disk_operations.remove_empty_directory(sample_target_path)
            except FileNotFoundError:
                pass
        disk_operations.remove_empty_directory(output_path)

    @staticmethod
    def move_alignment_plot(folder_watcher: FolderWatcher,
                            extension: str,
                            scratch_path: Path,
                            results_path: Path) -> None:
        alignment_path: Path = results_path / "alignment"
        disk_operations.mkdir_p(alignment_path, exist_ok=True)
        for sample_name in folder_watcher.all_samples:
            sample_name = trim_name(sample_name)
            source_path = scratch_path / sample_name / f'alignment{extension}'
            target_path = alignment_path / f"{sample_name}_alignment{extension}"
            if source_path.exists():
                disk_operations.rename(source_path, target_path)
        disk_operations.remove_empty_directory(alignment_path)

    @staticmethod
    def move_genome_coverage(folder_watcher: FolderWatcher, scratch_path: Path, results_path: Path) -> None:
        plots_path = results_path / "genome_coverage"
        disk_operations.mkdir_p(plots_path, exist_ok=True)
        for sample_name in folder_watcher.all_samples:
            sample_name = trim_name(sample_name)
            source_path = scratch_path / sample_name / 'genome_coverage.svg'
            target_path = plots_path / f"{sample_name}_genome_coverage.svg"
            if source_path.exists():
                disk_operations.rename(source_path, target_path)
            concordance_path = scratch_path / sample_name / 'genome_concordance.svg'
            target_concordance_path = plots_path / f"{sample_name}_genome_concordance.svg"
            if concordance_path.exists():
                disk_operations.rename(concordance_path, target_concordance_path)
        disk_operations.remove_empty_directory(plots_path)

    def run_filter_quality_pipeline(self, folder_watcher: FolderWatcher) -> Optional[Run]:
        if folder_watcher.quality_dataset is None:
            raise RuntimeError('Quality dataset not available')
        if self.config.micall_filter_quality_pipeline_id is None:
            raise RuntimeError('Filter quality pipeline ID not configured')

        return self.find_or_launch_run(
            self.config.micall_filter_quality_pipeline_id,
            dict(quality_csv=folder_watcher.quality_dataset),
            'MiCall filter quality on ' + folder_watcher.run_name,
            folder_watcher.batch)

    def run_pipeline(self,
                     folder_watcher: FolderWatcher,
                     pipeline_type: PipelineType,
                     sample_watcher: SampleWatcher) -> Optional[Run]:
        if pipeline_type == PipelineType.FILTER_QUALITY:
            run = self.run_filter_quality_pipeline(folder_watcher)
            return run
        if pipeline_type == PipelineType.PROVIRAL:
            run = self.run_proviral_pipeline(
                sample_watcher,
                folder_watcher,
                'Proviral')
            return run
        if pipeline_type == PipelineType.RESISTANCE:
            run = self.run_resistance_pipeline(
                sample_watcher,
                folder_watcher,
                (PipelineType.MAIN, PipelineType.MIDI),
                'MiCall resistance')
            return run
        if pipeline_type == PipelineType.DENOVO_RESISTANCE:
            run = self.run_resistance_pipeline(
                sample_watcher,
                folder_watcher,
                (PipelineType.DENOVO_MAIN, PipelineType.DENOVO_MIDI),
                'MiCall denovo resistance')
            return run
        if pipeline_type in (PipelineType.MIXED_HCV_MAIN,
                             PipelineType.MIXED_HCV_MIDI):
            if self.config.mixed_hcv_pipeline_id is None:
                return None
            if pipeline_type == PipelineType.MIXED_HCV_MAIN:
                input_datasets = dict(fastq1=sample_watcher.fastq_datasets[0],
                                      fastq2=sample_watcher.fastq_datasets[1])
                sample_name = sample_watcher.sample_group.names[0]
            else:
                input_datasets = dict(fastq1=sample_watcher.fastq_datasets[2],
                                      fastq2=sample_watcher.fastq_datasets[3])
                sample_name = sample_watcher.sample_group.names[1]
            return self.find_or_launch_run(
                self.config.mixed_hcv_pipeline_id,
                input_datasets,
                'Mixed HCV on ' + (trim_name(sample_name) if sample_name is not None else 'unknown'),
                folder_watcher.batch)
        if pipeline_type == PipelineType.MAIN:
            group_position = 0
            run_name = 'MiCall main'
            pipeline_id = self.config.micall_main_pipeline_id
        elif pipeline_type == PipelineType.MIDI:
            group_position = 1
            run_name = 'MiCall main'
            pipeline_id = self.config.micall_main_pipeline_id
        elif pipeline_type == PipelineType.DENOVO_MAIN:
            group_position = 0
            run_name = 'MiCall denovo main'
            pipeline_id = self.config.denovo_main_pipeline_id
        else:
            assert pipeline_type == PipelineType.DENOVO_MIDI
            group_position = 1
            run_name = 'MiCall denovo main'
            pipeline_id = self.config.denovo_main_pipeline_id
        if pipeline_id is None:
            return None
        fastq1, fastq2 = sample_watcher.fastq_datasets[
                         group_position*2:(group_position+1)*2]
        sample_name = sample_watcher.sample_group.names[group_position]
        run_name += ' on ' + (trim_name(sample_name) if sample_name is not None else 'unknown')
        sample_info = self.get_sample_info(pipeline_id,
                                           sample_watcher,
                                           folder_watcher,
                                           group_position)
        if folder_watcher.bad_cycles_dataset is None:
            if folder_watcher.filter_quality_run is None:
                raise RuntimeError('Filter quality run not available for bad cycles dataset')
            filter_run_id = folder_watcher.filter_quality_run['id']
            run_datasets: Sequence[RunDataset] = self.kive_retry(
                lambda: self.session.endpoints.containerruns.get(
                    f'{filter_run_id}/dataset_list/'))
            bad_cycles_run_dataset, = [
                run_dataset
                for run_dataset in run_datasets
                if run_dataset['argument_name'] == 'bad_cycles_csv']
            folder_watcher.bad_cycles_dataset = self.kive_retry(
                lambda: self.session.get(bad_cycles_run_dataset['dataset']).json())

        inputs = dict(fastq1=fastq1,
                      fastq2=fastq2,
                      bad_cycles_csv=folder_watcher.bad_cycles_dataset)
        if sample_info is not None:
            inputs['sample_info_csv'] = sample_info
        return self.find_or_launch_run(
            pipeline_id,
            inputs,
            run_name,
            folder_watcher.batch)

    def get_sample_info(self,
                        pipeline_id: int,
                        sample_watcher: SampleWatcher,
                        folder_watcher: FolderWatcher,
                        group_position: int):
        pipeline_args = self.get_kive_arguments(pipeline_id)
        if 'sample_info_csv' not in pipeline_args:
            return None

        # We need a sample info dataset.
        if group_position < len(sample_watcher.sample_info_datasets):
            return sample_watcher.sample_info_datasets[group_position]
        assert group_position == len(sample_watcher.sample_info_datasets)

        fastq_name = sample_watcher.sample_group.names[group_position]
        sample_name = trim_name(fastq_name if fastq_name is not None else 'unknown')
        project_code = sample_watcher.sample_group.project_codes[group_position]
        micall_version = self.get_container_version(pipeline_id)
        info_file = StringIO()
        writer = DictWriter(info_file, ['sample', 'project', 'run_name', 'micall_version'])
        writer.writeheader()
        writer.writerow(dict(sample=sample_name,
                             project=project_code,
                             run_name=folder_watcher.run_name,
                             micall_version=micall_version))
        bytes_file = BytesIO(info_file.getvalue().encode('utf8'))
        info_dataset = self.find_or_upload_dataset(
            bytes_file,
            f'{sample_name}_info.csv')
        sample_watcher.sample_info_datasets.append(info_dataset)

        return info_dataset

    def run_resistance_pipeline(self, sample_watcher: SampleWatcher, folder_watcher: FolderWatcher, input_pipeline_types: Sequence[PipelineType], description: str) -> Optional[Run]:
        pipeline_id = self.config.micall_resistance_pipeline_id
        if pipeline_id is None:
            return None
        main_runs = (
            (pipeline_type, run)
            for pipeline_type in input_pipeline_types
            for run in [sample_watcher.runs.get(pipeline_type)]
            if run is not None)
        input_dataset_urls = {
            (pipeline_type, run_dataset['argument_name']): run_dataset['dataset']
            for pipeline_type, run in main_runs
            for run_dataset in run['datasets']}
        main_type, midi_type = input_pipeline_types
        input_types = [(main_type, 'amino_csv'),
                       (midi_type, 'amino_csv'),
                       (main_type, 'nuc_csv')]
        if input_types[1] not in input_dataset_urls:
            input_types[1] = input_types[0]
        selected_urls = [input_dataset_urls[input_type]
                         for input_type in input_types]
        input_datasets = [self.kive_retry(lambda: self.session.get(url).json())
                          for url in selected_urls]
        inputs_dict = dict(zip(('main_amino_csv', 'midi_amino_csv', 'main_nuc_csv'),
                               input_datasets))
        run = self.find_or_launch_run(
            pipeline_id,
            inputs_dict,
            description + ' on ' + sample_watcher.sample_group.enum,
            folder_watcher.batch)
        return run

    def run_proviral_pipeline(self, sample_watcher: SampleWatcher, folder_watcher: FolderWatcher, description: str) -> Optional[Run]:
        pipeline_id = self.config.proviral_pipeline_id
        if pipeline_id is None:
            return None
        main_run = sample_watcher.runs.get(PipelineType.DENOVO_MAIN)
        if main_run is None:
            return None
        input_dataset_urls = {
            run_dataset['argument_name']: run_dataset['dataset']
            for run_dataset in main_run['datasets']
            if run_dataset['argument_name'] in ('sample_info_csv',
                                                'unstitched_cascade_csv',
                                                'unstitched_conseq_csv',
                                                'unstitched_contigs_csv')}
        input_datasets = {
            argument_name: self.kive_retry(lambda: self.session.get(url).json())
            for argument_name, url in input_dataset_urls.items()}
        input_datasets['cascade_csv'] = input_datasets.pop('unstitched_cascade_csv')
        input_datasets['conseqs_csv'] = input_datasets.pop('unstitched_conseq_csv')
        input_datasets['contigs_csv'] = input_datasets.pop('unstitched_contigs_csv')
        run = self.find_or_launch_run(
            pipeline_id,
            input_datasets,
            description + ' on ' + sample_watcher.sample_group.enum,
            folder_watcher.batch)
        return run

    def find_or_launch_run(self,
                           pipeline_id: int,
                           inputs: Mapping[str, RunDataset],
                           run_name: str,
                           run_batch: Optional[Batch]) -> Optional[Run]:
        """ Look for a matching container run, or start a new one.

        :return: the run dictionary
        """
        app_url = self.get_kive_app(pipeline_id)
        app_args = self.get_kive_arguments(pipeline_id)
        run_name = trim_run_name(run_name)
        filters = ['name', run_name, 'app_id', pipeline_id, 'states', 'NLRSC']
        for arg in inputs.values():
            filters.append('input_id')
            filters.append(arg['id'])

        old_runs = self.session.endpoints.containerruns.filter(*filters)
        untyped_run = self.find_name_and_permissions_match(old_runs, run_name, 'container run')
        run: Optional[Run] = untyped_run  # type: ignore

        if run:
            if run['state'] == 'C':
                run_id = run['id']
                run_datasets: Sequence[RunCreationDataset] = self.session.endpoints.containerruns.get(
                    f'{run_id}/dataset_list/')
                if any(run_dataset.get('dataset_purged')
                       for run_dataset in run_datasets):
                    run = None
        if run is None:
            try:
                run_datasets = [dict(argument=app_arg,
                                     dataset=inputs[name]['url'])
                                for name, app_arg in app_args.items()]
            except KeyError as e:
                raise ValueError(f"Pipeline input error: {repr(e)}."
                                 f" The specified app with id {pipeline_id} appears to expect a different set of inputs."
                                 f" Does the run name {repr(run_name)} make sense for it?")
            if run_batch is None:
                raise ValueError("Run batch is required.")
            run_params = dict(name=run_name,
                              batch=run_batch['url'],
                              groups_allowed=ALLOWED_GROUPS,
                              app=app_url,
                              datasets=run_datasets)
            try:
                run = self.session.endpoints.containerruns.post(json=run_params)
            except Exception as ex:
                raise RuntimeError(
                    'Failed to launch run {}.'.format(run_name)) from ex
        return run

    def kive_retry(self, target: Callable[[], T]) -> T:
        """ Add a single retry to a Kive API call.

        Tries to call the target function, then refreshes the session login if
        the call fails and tries a second time.
        """
        try:
            return target()
        except KiveClientException:
            self.session.login(self.config.kive_user, self.config.kive_password)
            return target()

    def fetch_run_status(self, run: Run, folder_watcher: FolderWatcher, pipeline_type: PipelineType, sample_watchers: Sequence[Optional[SampleWatcher]]) -> Optional[Run]:
        self.check_session()
        new_status: Run = self.kive_retry(lambda: self.session.endpoints.containerruns.get(run['id']))
        is_complete = new_status['state'] == 'C'
        if new_status['state'] == 'X':
            new_run: Optional[Run] = None
            for sample_watcher in sample_watchers:
                if sample_watcher is None:
                    if pipeline_type != PipelineType.FILTER_QUALITY:
                        raise RuntimeError('Filter quality run failed, but no sample watcher provided.')
                    new_run = self.run_filter_quality_pipeline(folder_watcher)
                else:
                    new_run = self.run_pipeline(folder_watcher,
                                                pipeline_type,
                                                sample_watcher)
            return new_run
        if new_status['state'] == 'F':
            raise KiveRunFailedException(f'Run id {new_status["id"]} failed.')
        if is_complete and pipeline_type != PipelineType.FILTER_QUALITY:
            run_datasets = self.kive_retry(
                lambda: self.session.endpoints.containerruns.get(
                    f"{run['id']}/dataset_list/"))
            for sample_watcher in sample_watchers:
                if sample_watcher is None:
                    raise RuntimeError('A non-filter quality run completed, but no sample watcher provided.')

                for other_run in sample_watcher.runs.values():
                    if other_run['id'] == run['id']:
                        other_run['datasets'] = run_datasets
                        break
                sample_name = (sample_watcher.sample_group.names[1]
                               if pipeline_type in (PipelineType.MIDI,
                                                    PipelineType.MIXED_HCV_MIDI,
                                                    PipelineType.DENOVO_MIDI)
                               else sample_watcher.sample_group.names[0])
                results_path = self.get_results_path(folder_watcher)
                pipeline_group = PIPELINE_GROUPS[pipeline_type]
                scratch_path = get_scratch_path(results_path, pipeline_group)
                scratch_path /= trim_name(sample_name if sample_name is not None else 'unknown')
                disk_operations.mkdir_p(scratch_path, parents=True, exist_ok=True)
                for output_name in DOWNLOADED_RESULTS:
                    matches = [run_dataset
                               for run_dataset in run_datasets
                               if (run_dataset['argument_name'] == output_name and
                                   run_dataset['argument_type'] == 'O')]
                    if not matches:
                        continue
                    filename = get_output_filename(output_name)
                    dataset_url, = [match['dataset'] for match in matches]
                    self.kive_retry(
                        lambda: self.download_file(dataset_url + 'download/',
                                                   scratch_path / filename))

        if is_complete:
            return None
        return run

    def download_file(self, dataset_url: str, target_path: Path) -> None:
        with disk_operations.disk_file_operation(target_path, 'wb') as f:
            self.session.download_file(f, dataset_url)

    def get_results_path(self, folder_watcher: FolderWatcher) -> Path:
        version_name = f'version_{self.config.pipeline_version}'
        results_path = folder_watcher.run_folder / "Results" / version_name
        return results_path

    def upload_filter_quality(self, folder_watcher: FolderWatcher) -> None:
        read_sizes = parse_read_sizes(folder_watcher.run_folder / "RunInfo.xml")
        read_lengths = ReadLengths4(
            forward_read=read_sizes.read1,
            index1=read_sizes.index1,
            index2=read_sizes.index2,
            reverse_read=read_sizes.read2
        )
        error_path = folder_watcher.run_folder / "InterOp" / "ErrorMetricsOut.bin"
        quality_csv = StringIO()

        if not error_path.exists():
            logger.warning("ErrorMetricsOut.bin not found in %s",
                           folder_watcher.run_folder)
            disk_operations.write_text(folder_watcher.run_folder / "errorprocessing",
                                     "ErrorMetricsOut.bin not found.\n")
            return

        start_time = None
        for attempt_count in count(1):
            try:
                reader = InterOpReader(folder_watcher.run_folder)
                records = reader.read_error_records()
                write_phix_csv(quality_csv, records, read_lengths)
                break
            except OSError as ex:
                logger.info(f"Caught error while trying to read the phix data: {ex}")

                if start_time is None:
                    start_time = datetime.now()
                wait_for_retry(attempt_count, start_time)
            except Exception:
                logger.error("Finding error metrics in %s",
                         folder_watcher.run_folder,
                         exc_info=True)
                disk_operations.write_text(folder_watcher.run_folder / "errorprocessing",
                                           "Finding error metrics failed.\n")
                return

        quality_csv_bytes = BytesIO()
        quality_csv_bytes.write(quality_csv.getvalue().encode('utf8'))
        quality_csv_bytes.seek(0)
        folder_watcher.quality_dataset = self.find_or_upload_dataset(
            quality_csv_bytes,
            folder_watcher.run_name + '_quality.csv',
            'Error rates for {} run.'.format(folder_watcher.run_name))

    @property
    def session(self) -> KiveAPI:
        if self._session is None:
            new = open_kive(self.config.kive_server)
            new.login(self.config.kive_user, self.config.kive_password)
            self._session = new
        return self._session

    def check_session(self):
        if self._session is None:
            # retrieve external file directory
            directories = self.session.endpoints.externalfiledirectories.get()
            self.external_directory_name = self.external_directory_path = None
            search_path = self.config.raw_data / 'MiSeq' / 'runs'
            for directory in directories:
                directory_path = Path(directory['path'])
                if (directory_path == search_path or
                        directory_path in search_path.parents):
                    self.external_directory_name = directory['name']
                    self.external_directory_path = directory_path
                    break
            else:
                logger.warning('No external directory found to match %r.',
                               str(search_path))

    def find_or_upload_dataset(self,
                               dataset_file: IO[bytes],
                               dataset_name: str,
                               description: str = '') -> Optional[RunDataset]:
        dataset = self.find_kive_dataset(dataset_file, dataset_name)
        if dataset is None:
            dataset_file.seek(0)
            dataset = self.upload_kive_dataset(
                dataset_file,
                dataset_name,
                description)
        return dataset
