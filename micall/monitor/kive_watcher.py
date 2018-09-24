import hashlib
import logging
from datetime import datetime, timedelta
from queue import Full

from io import StringIO, BytesIO

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
FOLDER_SCAN_INTERVAL = timedelta(hours=1)


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


def find_samples(raw_data_folder, sample_queue, wait=True):
    while True:
        is_complete = scan_samples(raw_data_folder, sample_queue)
        if is_complete and not wait:
            break


def scan_samples(raw_data_folder, sample_queue):
    next_scan = now() + FOLDER_SCAN_INTERVAL
    flag_paths = sorted(
        raw_data_folder.glob("MiSeq/runs/*/needsprocessing"),
        reverse=True)
    for flag_path in flag_paths:
        run_path = flag_path.parent
        base_calls_path = run_path / "Data/Intensities/BaseCalls"
        fastq_files = base_calls_path.glob("*_R1_*.fastq.gz")
        sample_sheet_path = run_path / "SampleSheet.csv"
        file_names = [f.name for f in fastq_files]
        for sample_group in find_groups(file_names, sample_sheet_path):
            is_sent = False
            while not is_sent and now() < next_scan:
                try:
                    sample_queue.put((base_calls_path, sample_group), timeout=60)
                    logger.debug('Put %s, %s in queue.', base_calls_path, sample_group)
                    is_sent = True
                except Full:
                    pass
            if now() >= next_scan:
                return False
    return True


class KiveWatcher:
    def __init__(self, config=None):
        self.config = config
        self.session = None
        self.folder_watchers = {}  # {base_calls_folder: FolderWatcher}
        self.pipelines = {}  # {pipeline_id: pipeline}
        self.input_pipeline_ids = dict(
            quality_csv=config.micall_filter_quality_pipeline_id,
            bad_cycles_csv=config.micall_main_pipeline_id,
            main_amino_csv=config.micall_resistance_pipeline_id,
            midi_amino_csv=config.micall_resistance_pipeline_id)

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
        batch_name = folder_watcher.run_name + ' ' + self.config.pipeline_version
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
        dataset = self.session.add_dataset(
            name=dataset_name,
            description=description,
            handle=source_file,
            cdt=cdt,
            groups=ALLOWED_GROUPS)
        # TODO: external data sets
        # if (self.external_directory_name is None or
        #         not filename.startswith(self.external_directory_name)):
        # else:
        #     external_path = os.path.relpath(filename,
        #                                     self.external_directory_path)
        #     dataset = self.kive.add_dataset(
        #         name=dataset_name,
        #         description=description,
        #         handle=None,
        #         externalfiledirectory=self.external_directory_name,
        #         external_path=external_path,
        #         cdt=cdt,
        #         groups=settings.kive_groups_allowed)
        return dataset

    def add_sample_group(self, base_calls, sample_group):
        self.check_session()
        folder_watcher = self.folder_watchers.get(base_calls)
        if folder_watcher is None:
            folder_watcher = FolderWatcher(base_calls, self)
            self.folder_watchers[base_calls] = folder_watcher
            self.create_batch(folder_watcher)
            self.upload_filter_quality(folder_watcher)

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

    def poll_runs(self):
        for folder_watcher in self.folder_watchers.values():
            folder_watcher.poll_runs()

    def run_pipeline(self, folder_watcher, sample_watcher, pipeline_type):
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
            inputs = [run.get_results()['amino_csv']
                      for run in sample_watcher.main_runs]
            if len(inputs) == 1:
                inputs *= 2
            return self.session.run_pipeline(
                resistance_pipeline,
                inputs,
                'MiCall resistance on ' + sample_watcher.sample_group.enum,
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
            folder_watcher.bad_cycles_dataset.cdt = bad_cycles_input.compounddatatype

        run_name = 'MiCall main on ' + '_'.join(sample_name.split('_')[:2])
        return self.session.run_pipeline(
            main_pipeline,
            [fastq1, fastq2, folder_watcher.bad_cycles_dataset],
            run_name,
            runbatch=folder_watcher.batch,
            groups=ALLOWED_GROUPS)

    @staticmethod
    def fetch_run_status(run):
        return run.is_complete()

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

    def find_or_upload_dataset(self, dataset_file, dataset_name, description, compounddatatype):
        quality_dataset = self.find_kive_dataset(dataset_file,
                                                 dataset_name,
                                                 compounddatatype)
        if quality_dataset is None:
            dataset_file.seek(0)
            quality_dataset = self.upload_kive_dataset(
                dataset_file,
                dataset_name,
                compounddatatype,
                description)
        return quality_dataset
