import hashlib
import logging
from datetime import datetime, timedelta
from queue import Full

from io import StringIO, BytesIO

from micall.drivers.run_info import parse_read_sizes
from micall.monitor import error_metrics_parser
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
ALLOWED_GROUPS = ['Everyone']
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

    def is_full(self):
        return False

    def get_kive_pipeline(self, pipeline_id):
        kive_pipeline = self.session.get_pipeline(pipeline_id)
        return kive_pipeline

    def create_batch(self, run_name):
        batch_name = run_name + ' ' + self.config.pipeline_version
        description = 'MiCall batch for folder {}, pipeline version {}.'.format(
            run_name,
            self.config.pipeline_version)
        batch = self.session.create_run_batch(batch_name,
                                              description=description,
                                              users=[],
                                              groups=ALLOWED_GROUPS)
        return batch

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
        if self.session is None:
            self.session = open_kive(self.config.kive_server)
            self.session.login(self.config.kive_user, self.config.kive_password)
        run_folder = (base_calls / "../../..").resolve()
        run_name = '_'.join(run_folder.name.split('_')[:2])
        run_batch = self.create_batch(run_name)
        read_sizes = parse_read_sizes(run_folder / "RunInfo.xml")
        read_lengths = [read_sizes.read1,
                        read_sizes.index1,
                        read_sizes.index2,
                        read_sizes.read2]
        quality_pipeline = self.get_kive_pipeline(
            self.config.micall_filter_quality_pipeline_id)
        quality_input = quality_pipeline.inputs[0]
        error_path = run_folder / "InterOp/ErrorMetricsOut.bin"
        quality_csv = StringIO()
        with error_path.open('rb') as error_file:
            records = error_metrics_parser.read_errors(error_file)
            error_metrics_parser.write_phix_csv(quality_csv,
                                                records,
                                                read_lengths)
        quality_csv_bytes = BytesIO()
        quality_csv_bytes.write(quality_csv.getvalue().encode('utf8'))
        quality_csv_bytes.seek(0)
        quality_file_name = run_name + '_quality.csv'
        quality_dataset = self.find_kive_dataset(quality_csv_bytes,
                                                 quality_file_name,
                                                 quality_input.compounddatatype)
        if quality_dataset is None:
            quality_dataset = self.upload_kive_dataset(
                quality_csv_bytes,
                quality_file_name,
                quality_input.compounddatatype,
                'Error rates for {} run.'.format(run_name))
        self.session.run_pipeline(quality_pipeline,
                                  [quality_dataset],
                                  'MiCall filter quality on ' + run_name,
                                  runbatch=run_batch,
                                  groups=ALLOWED_GROUPS)
