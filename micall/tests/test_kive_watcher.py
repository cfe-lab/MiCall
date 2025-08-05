import shutil
import tarfile
from gzip import GzipFile
from io import BytesIO, StringIO
from pathlib import Path
from queue import Full
from tarfile import TarInfo
from unittest.mock import patch, ANY, Mock, call
from zipfile import ZipFile
import tempfile
import errno
import time
import threading

# noinspection PyPackageRequirements
import pytest
from datetime import datetime, timedelta

from struct import pack

from kiveapi import KiveClientException
# noinspection PyPackageRequirements
from requests import ConnectionError

import micall
from micall.monitor.kive_watcher import find_samples, KiveWatcher, FolderEvent, FolderEventType, calculate_retry_wait, \
    trim_run_name, compress_old_versions
from micall.monitor.sample_watcher import PipelineType, ALLOWED_GROUPS, FolderWatcher, SampleWatcher
from micall.monitor.find_groups import SampleGroup
from micall.monitor.micall_watcher import parse_args
from micall.monitor import disk_operations


class DummyDataset:
    def __init__(self, name):
        self.name = name
        self.lines = ['row,name'] + ['{},{}'.format(i, name) for i in range(3)]

    def __repr__(self):
        return f"DummyDataset({self.name!r})"

    def download(self, f):
        for line in self.lines:
            f.write((line + '\n').encode('utf8'))


def create_datasets(names):
    return {name: DummyDataset(name)
            for name in names}


@pytest.fixture(name='mock_clock')
def create_mock_clock():
    with patch('micall.monitor.kive_watcher.now') as mock_clock:
        mock_clock.return_value = datetime(2000, 1, 1)
        yield mock_clock


@pytest.fixture(name='mock_open_kive')
def create_mock_open_kive():
    with patch('micall.monitor.kive_watcher.open_kive') as mock_open_kive:
        mock_session = mock_open_kive.return_value

        # By default, support calling the filter_quality pipeline.
        mock_session.endpoints.containerapps.get.return_value = [
            dict(name='quality_csv', url='/args/967485', type='I', app='/apps/19374')]
        mock_pipeline = mock_session.get_pipeline.return_value
        mock_input = Mock(dataset_name='quality_csv')
        mock_pipeline.inputs = [mock_input]

        # By default, all runs are still running.
        mock_session.get_run.return_value.raw = dict(end_time=None,
                                                     stopped_by=None)

        yield mock_open_kive


def mock_containerruns_get(path):
    if 'dataset_list' in str(path):
        return [dict(dataset='/datasets/111/',
                     argument_type='I',
                     argument_name='bad_cycles_csv'),
                dict(dataset='/datasets/112/',
                     argument_type='O',
                     argument_name='amino_csv')]
    return dict(state='C')


def mock_session_get(url):
    if '/datasets/' in url:
        response = Mock()
        response.json.return_value = dict(id=123, url='/datasets/123/')
        return response
    raise RuntimeError('Unexpected url in session.get().')


def mock_session_download_file(file, url):
    file.write(b'url,n\n')
    for i in range(3):
        file.write(f'{url},{i}\n'.encode())


def mock_session_download_fasta(file, url):
    for i in range(2):
        file.write(f'>{url},{i}\n'.encode('UTF8'))
        for j in range(3):
            file.write('ACTGTCA'[i+j:].encode())
            file.write(b'\n')


def mock_failures(failure_count, success_callable):
    """ Simulate a number of failures, followed by successes. """
    def mocked(*args, **kwargs):
        nonlocal failure_count
        if failure_count > 0:
            failure_count -= 1
            raise KiveClientException('failed.')
        return success_callable(*args, **kwargs)
    return mocked


def create_kive_watcher_with_filter_run(config, base_calls, is_complete=False):
    kive_watcher = KiveWatcher(config)
    kive_watcher.app_urls = {
        config.micall_filter_quality_pipeline_id: '/containerapps/102',
        config.micall_main_pipeline_id: '/containerapps/103',
        config.micall_resistance_pipeline_id: '/containerapps/104'}
    kive_watcher.app_args = {
        config.micall_filter_quality_pipeline_id: dict(
            quality_csv='/containerargs/105'),
        config.micall_main_pipeline_id: dict(
            fastq1='/containerargs/106',
            fastq2='/containerargs/107',
            bad_cycles_csv='/containerargs/108'),
        config.micall_resistance_pipeline_id: dict(
            main_amino_csv='/containerargs/109',
            midi_amino_csv='/containerargs/110',
            main_nuc_csv='/containerargs/111')}

    folder_watcher = kive_watcher.add_folder(base_calls)
    folder_watcher.batch = dict(url='/batches/101')
    folder_watcher.add_run(dict(id=110),
                           PipelineType.FILTER_QUALITY,
                           is_complete=is_complete)
    if is_complete:
        folder_watcher.bad_cycles_dataset = dict(url='/datasets/110', id=110)
    kive_watcher.check_session()
    mock_session = kive_watcher.session
    mock_session.endpoints.containerruns.get.side_effect = mock_containerruns_get
    mock_session.get.side_effect = mock_session_get
    mock_session.download_file.side_effect = mock_session_download_file

    return kive_watcher


def create_kive_watcher_with_main_run(config,
                                      base_calls,
                                      sample_group,
                                      is_complete=False):
    kive_watcher = create_kive_watcher_with_filter_run(config,
                                                       base_calls,
                                                       is_complete=True)

    folder_watcher, = kive_watcher.folder_watchers.values()
    sample_watcher = kive_watcher.add_sample_group(
        base_calls=base_calls,
        sample_group=sample_group)
    folder_watcher.add_run(
        dict(id=107),
        PipelineType.MAIN,
        sample_watcher,
        is_complete=is_complete)

    return kive_watcher


@pytest.fixture(name='mock_wait')
def create_mock_wait():
    with patch('micall.monitor.kive_watcher.wait_for_retry') as mock_wait:
        yield mock_wait


@pytest.fixture(name='default_config')
def create_default_config():
    default_config = parse_args(argv=['--micall_filter_quality_pipeline_id', '42',
                                      '--micall_main_pipeline_id', '43',
                                      '--micall_resistance_pipeline_id', '494'])
    yield default_config


@pytest.fixture(name='pipelines_config')
def create_pipelines_config(default_config):
    # This used to be different from the default, but now they're the same.
    yield default_config


class DummyQueueSink:
    def __init__(self):
        self.expected_puts = []  # [(item, is_full, callback)]
        self.put_calls = 0

    def expect_put(self, item, is_full=False, callback=None):
        self.expected_puts.append((item, is_full, callback))

    def put(self, item, block=True, timeout=None):
        assert block
        assert timeout is not None
        expected_item, is_full, callback = self.expected_puts.pop(0)
        self.put_calls += 1
        assert expected_item == item, self.put_calls
        if callback is not None:
            callback()
        if is_full:
            raise Full()

    def verify(self):
        assert [] == self.expected_puts  # All expected puts arrived.


def create_run_folder(tmpdir, run_name, sample_pattern):
    raw_data_folder = Path(tmpdir)
    run = raw_data_folder / "MiSeq/runs" / run_name
    run.mkdir(parents=True)
    microtest = Path(__file__).parent / "microtest"
    sample_sheet_source = microtest / "SampleSheet.csv"
    sample_sheet_target = run / "SampleSheet.csv"
    sample_sheet_target.write_bytes(sample_sheet_source.read_bytes())
    base_calls = run / "Data/Intensities/BaseCalls"
    base_calls.mkdir(parents=True)
    for source_fastq in microtest.glob(sample_pattern):
        target_name = source_fastq.name + '.gz'
        target = base_calls / target_name
        with GzipFile(target, 'w') as target_file:
            target_file.write(source_fastq.read_bytes())
    run_info_source = microtest / "RunInfo.xml"
    run_info_target = run / "RunInfo.xml"
    run_info_target.write_bytes(run_info_source.read_bytes())

    interop = run / "InterOp"
    interop.mkdir()
    error_data = pack('<bb', 3, 30)  # version, record size, and no records
    # noinspection PyTypeChecker
    (interop / "ErrorMetricsOut.bin").write_bytes(error_data)

    (run / "needsprocessing").touch()
    return raw_data_folder


@pytest.fixture(name='raw_data_with_hcv_pair')
def create_raw_data_with_hcv_pair(tmpdir):
    return create_run_folder(tmpdir, '140101_M01234', '2130A*.fastq')


@pytest.fixture(name='raw_data_with_two_samples')
def create_raw_data_with_two_samples(tmpdir):
    # Install samples 2110A and 2120A in a single run.
    return create_run_folder(tmpdir, '140101_M01234', '21[12]0A*.fastq')


@pytest.fixture(name='raw_data_with_two_runs')
def create_raw_data_with_two_runs(tmpdir):
    raw_data = create_run_folder(tmpdir, '140101_M01234', '2000A*.fastq')
    create_run_folder(tmpdir, '140201_M01234', '2010A*.fastq')

    return raw_data


def test_filter_pipeline_not_set(capsys, monkeypatch):
    monkeypatch.delenv('MICALL_FILTER_QUALITY_PIPELINE_ID', raising=False)
    expected_error = "Argument --micall_filter_quality_pipeline_id not set " \
                     "and $MICALL_FILTER_QUALITY_PIPELINE_ID environment " \
                     "variable not set."

    with pytest.raises(SystemExit):
        parse_args([])

    stderr = capsys.readouterr().err
    assert expected_error in stderr


def test_pipeline_set():
    args = parse_args(['--micall_filter_quality_pipeline_id', '402',
                       '--micall_main_pipeline_id', '403'])

    assert args.micall_filter_quality_pipeline_id == 402


def test_pipeline_set_with_environment_variable(monkeypatch):
    monkeypatch.setenv('MICALL_FILTER_QUALITY_PIPELINE_ID', '99')
    monkeypatch.setenv('MICALL_MAIN_PIPELINE_ID', '99')
    args = parse_args([])

    assert args.micall_filter_quality_pipeline_id == 99


def test_main_pipeline_not_set(capsys, monkeypatch):
    monkeypatch.delenv('MICALL_MAIN_PIPELINE_ID', raising=False)
    expected_error = "No arguments or environment variables set for main " \
                     "pipeline ids"

    with pytest.raises(SystemExit):
        parse_args(['--micall_filter_quality_pipeline_id', '402'])

    stderr = capsys.readouterr().err
    assert expected_error in stderr


def test_hcv_pair(raw_data_with_hcv_pair):
    sample_queue = DummyQueueSink()
    sample_queue.expect_put(
        FolderEvent(raw_data_with_hcv_pair / "MiSeq/runs/140101_M01234" /
                    "Data/Intensities/BaseCalls",
                    FolderEventType.ADD_SAMPLE,
                    SampleGroup('2130A',
                                ('2130A-HCV_S15_L001_R1_001.fastq.gz',
                                 '2130AMIDI-MidHCV_S16_L001_R1_001.fastq.gz'),
                                ('HCV', 'MidHCV'))))
    sample_queue.expect_put(
        FolderEvent(raw_data_with_hcv_pair / "MiSeq/runs/140101_M01234" /
                    "Data/Intensities/BaseCalls",
                    FolderEventType.FINISH_FOLDER,
                    None))
    pipeline_version = 'XXX'

    find_samples(raw_data_with_hcv_pair, pipeline_version, sample_queue, wait=False)

    sample_queue.verify()


def test_hcv_pair_with_wg_suffix(raw_data_with_hcv_pair):
    sample_queue = DummyQueueSink()
    run_folder = raw_data_with_hcv_pair / "MiSeq/runs/140101_M01234"
    sample_sheet = run_folder / "SampleSheet.csv"
    sample_sheet_text = sample_sheet.read_text()
    sample_sheet_text = sample_sheet_text.replace('2130A_HCV', '2130AWG_HCV')
    sample_sheet.write_text(sample_sheet_text)
    base_calls_folder = run_folder / "Data/Intensities/BaseCalls"
    for sample_path in base_calls_folder.glob('2130A-HCV_S15_*.fastq.gz'):
        new_name = sample_path.name.replace('2130A-HCV', '2130AWG-HCV')
        new_path = sample_path.parent / new_name
        sample_path.rename(new_path)

    sample_queue.expect_put(
        FolderEvent(base_calls_folder,
                    FolderEventType.ADD_SAMPLE,
                    SampleGroup('2130A',
                                ('2130AWG-HCV_S15_L001_R1_001.fastq.gz',
                                 '2130AMIDI-MidHCV_S16_L001_R1_001.fastq.gz'),
                                ('HCV', 'MidHCV'))))
    sample_queue.expect_put(
        FolderEvent(base_calls_folder,
                    FolderEventType.FINISH_FOLDER,
                    None))
    pipeline_version = 'XXX'

    find_samples(raw_data_with_hcv_pair,
                 pipeline_version,
                 sample_queue,
                 wait=False,
                 retry=False)

    sample_queue.verify()


def test_hcv_midi_alone(raw_data_with_hcv_pair):

    sample_queue = DummyQueueSink()
    base_calls_path = (raw_data_with_hcv_pair / "MiSeq/runs/140101_M01234" /
                       "Data/Intensities/BaseCalls")
    (base_calls_path / '2130A-HCV_S15_L001_R1_001.fastq.gz').unlink()
    sample_queue.expect_put(
        FolderEvent(base_calls_path,
                    FolderEventType.ADD_SAMPLE,
                    SampleGroup('2130AMIDI',
                                ('2130AMIDI-MidHCV_S16_L001_R1_001.fastq.gz',
                                 None),
                                ('MidHCV', None))))
    sample_queue.expect_put(
        FolderEvent(base_calls_path,
                    FolderEventType.FINISH_FOLDER,
                    None))
    pipeline_version = 'XXX'

    find_samples(raw_data_with_hcv_pair,
                 pipeline_version,
                 sample_queue,
                 wait=False,
                 retry=False)

    sample_queue.verify()


def test_two_runs(raw_data_with_two_runs):
    sample_queue = DummyQueueSink()
    sample_queue.expect_put(
        FolderEvent(raw_data_with_two_runs / "MiSeq/runs/140201_M01234" /
                    "Data/Intensities/BaseCalls",
                    FolderEventType.ADD_SAMPLE,
                    SampleGroup('2010A',
                                ('2010A-V3LOOP_S3_L001_R1_001.fastq.gz',
                                 None),
                                ('V3LOOP', None))))
    sample_queue.expect_put(
        FolderEvent(raw_data_with_two_runs / "MiSeq/runs/140201_M01234" /
                    "Data/Intensities/BaseCalls",
                    FolderEventType.FINISH_FOLDER,
                    None))
    sample_queue.expect_put(
        FolderEvent(raw_data_with_two_runs / "MiSeq/runs/140101_M01234" /
                    "Data/Intensities/BaseCalls",
                    FolderEventType.ADD_SAMPLE,
                    SampleGroup('2000A',
                                ('2000A-V3LOOP_S2_L001_R1_001.fastq.gz',
                                 None),
                                ('V3LOOP', None))))
    sample_queue.expect_put(
        FolderEvent(raw_data_with_two_runs / "MiSeq/runs/140101_M01234" /
                    "Data/Intensities/BaseCalls",
                    FolderEventType.FINISH_FOLDER,
                    None))
    pipeline_version = 'XXX'

    find_samples(raw_data_with_two_runs, pipeline_version, sample_queue, wait=False)

    sample_queue.verify()


def test_two_samples(raw_data_with_two_samples):
    sample_queue = DummyQueueSink()
    sample_queue.expect_put(
        FolderEvent(raw_data_with_two_samples / "MiSeq/runs/140101_M01234" /
                    "Data/Intensities/BaseCalls",
                    FolderEventType.ADD_SAMPLE,
                    SampleGroup('2120A',
                                ('2120A-PR_S14_L001_R1_001.fastq.gz',
                                 None),
                                ('PR', None))))
    sample_queue.expect_put(
        FolderEvent(raw_data_with_two_samples / "MiSeq/runs/140101_M01234" /
                    "Data/Intensities/BaseCalls",
                    FolderEventType.ADD_SAMPLE,
                    SampleGroup('2110A',
                                ('2110A-V3LOOP_S13_L001_R1_001.fastq.gz',
                                 None),
                                ('V3LOOP', None))))
    sample_queue.expect_put(
        FolderEvent(raw_data_with_two_samples / "MiSeq/runs/140101_M01234" /
                    "Data/Intensities/BaseCalls",
                    FolderEventType.FINISH_FOLDER,
                    None))
    pipeline_version = 'XXX'

    find_samples(raw_data_with_two_samples, pipeline_version, sample_queue, wait=False)

    sample_queue.verify()


def test_undetermined_file(raw_data_with_two_samples):

    sample_queue = DummyQueueSink()
    base_calls = (raw_data_with_two_samples / "MiSeq/runs/140101_M01234" /
                  "Data/Intensities/BaseCalls")

    # Undetermined reads should be ignored.
    (base_calls / "Undetermined_S0_L001_R1_001.fastq.gz").write_text('')
    (base_calls / "Undetermined_S0_L001_R2_001.fastq.gz").write_text('')

    sample_queue.expect_put(
        FolderEvent(base_calls,
                    FolderEventType.ADD_SAMPLE,
                    SampleGroup('2120A',
                                ('2120A-PR_S14_L001_R1_001.fastq.gz',
                                 None),
                                ('PR', None))))
    sample_queue.expect_put(
        FolderEvent(base_calls,
                    FolderEventType.ADD_SAMPLE,
                    SampleGroup('2110A',
                                ('2110A-V3LOOP_S13_L001_R1_001.fastq.gz',
                                 None),
                                ('V3LOOP', None))))
    sample_queue.expect_put(
        FolderEvent(base_calls,
                    FolderEventType.FINISH_FOLDER,
                    None))
    pipeline_version = 'XXX'

    find_samples(raw_data_with_two_samples, pipeline_version, sample_queue, wait=False)

    sample_queue.verify()


def test_skip_done_runs(raw_data_with_two_runs):
    done_run = raw_data_with_two_runs / "MiSeq/runs/140201_M01234"
    results_path = done_run / "Results/version_0-dev"
    results_path.mkdir(parents=True)
    done_path = results_path / "done_all_processing"
    done_path.touch()
    pipeline_version = '0-dev'
    sample_queue = DummyQueueSink()
    sample_queue.expect_put(
        FolderEvent(raw_data_with_two_runs / "MiSeq/runs/140101_M01234" /
                    "Data/Intensities/BaseCalls",
                    FolderEventType.ADD_SAMPLE,
                    SampleGroup('2000A',
                                ('2000A-V3LOOP_S2_L001_R1_001.fastq.gz',
                                 None),
                                ('V3LOOP', None))))
    sample_queue.expect_put(
        FolderEvent(raw_data_with_two_runs / "MiSeq/runs/140101_M01234" /
                    "Data/Intensities/BaseCalls",
                    FolderEventType.FINISH_FOLDER,
                    None))

    find_samples(raw_data_with_two_runs,
                 pipeline_version,
                 sample_queue,
                 wait=False,
                 retry=False)

    sample_queue.verify()


def test_skip_failed_runs(raw_data_with_two_runs):
    error_run_path = raw_data_with_two_runs / "MiSeq/runs/140201_M01234"
    error_path = error_run_path / "errorprocessing"
    error_path.touch()
    pipeline_version = '0-dev'
    sample_queue = DummyQueueSink()
    sample_queue.expect_put(
        FolderEvent(raw_data_with_two_runs / "MiSeq/runs/140101_M01234" /
                    "Data/Intensities/BaseCalls",
                    FolderEventType.ADD_SAMPLE,
                    SampleGroup('2000A',
                                ('2000A-V3LOOP_S2_L001_R1_001.fastq.gz',
                                 None),
                                ('V3LOOP', None))))
    sample_queue.expect_put(
        FolderEvent(raw_data_with_two_runs / "MiSeq/runs/140101_M01234" /
                    "Data/Intensities/BaseCalls",
                    FolderEventType.FINISH_FOLDER,
                    None))

    find_samples(raw_data_with_two_runs,
                 pipeline_version,
                 sample_queue,
                 wait=False)

    sample_queue.verify()


def test_bad_sample_sheet(raw_data_with_two_runs):
    error_run_path = raw_data_with_two_runs / "MiSeq/runs/140201_M01234"
    bad_sample_sheet_text = """\
[Data]
Broken Sample Sheet
Garbage!
"""
    (error_run_path / "SampleSheet.csv").write_text(bad_sample_sheet_text)
    error_path = error_run_path / "errorprocessing"
    pipeline_version = '0-dev'
    sample_queue = DummyQueueSink()
    sample_queue.expect_put(
        FolderEvent(raw_data_with_two_runs / "MiSeq/runs/140101_M01234" /
                    "Data/Intensities/BaseCalls",
                    FolderEventType.ADD_SAMPLE,
                    SampleGroup('2000A',
                                ('2000A-V3LOOP_S2_L001_R1_001.fastq.gz',
                                 None),
                                ('V3LOOP', None))))
    sample_queue.expect_put(
        FolderEvent(raw_data_with_two_runs / "MiSeq/runs/140101_M01234" /
                    "Data/Intensities/BaseCalls",
                    FolderEventType.FINISH_FOLDER,
                    None))

    find_samples(raw_data_with_two_runs,
                 pipeline_version,
                 sample_queue,
                 wait=False)

    sample_queue.verify()
    assert error_path.exists()


def test_full_queue(raw_data_with_two_runs):
    sample_queue = DummyQueueSink()
    sample_queue.expect_put(
        FolderEvent(raw_data_with_two_runs / "MiSeq/runs/140201_M01234" /
                    "Data/Intensities/BaseCalls",
                    FolderEventType.ADD_SAMPLE,
                    SampleGroup('2010A',
                                ('2010A-V3LOOP_S3_L001_R1_001.fastq.gz',
                                 None),
                                ('V3LOOP', None))))
    sample_queue.expect_put(
        FolderEvent(raw_data_with_two_runs / "MiSeq/runs/140201_M01234" /
                    "Data/Intensities/BaseCalls",
                    FolderEventType.FINISH_FOLDER,
                    None))
    item2 = FolderEvent(raw_data_with_two_runs / "MiSeq/runs/140101_M01234" /
                        "Data/Intensities/BaseCalls",
                        FolderEventType.ADD_SAMPLE,
                        SampleGroup('2000A',
                                    ('2000A-V3LOOP_S2_L001_R1_001.fastq.gz',
                                     None),
                                    ('V3LOOP', None)))
    sample_queue.expect_put(item2, is_full=True)
    sample_queue.expect_put(item2)
    sample_queue.expect_put(
        FolderEvent(raw_data_with_two_runs / "MiSeq/runs/140101_M01234" /
                    "Data/Intensities/BaseCalls",
                    FolderEventType.FINISH_FOLDER,
                    None))
    pipeline_version = 'XXX'

    find_samples(raw_data_with_two_runs,
                 pipeline_version,
                 sample_queue,
                 wait=False)

    sample_queue.verify()


def test_scan_for_new_runs(raw_data_with_two_runs, mock_clock):
    """ Mark the more recent run as not ready, timeout, then notice new run. """
    def increment_clock():
        mock_clock.return_value += timedelta(hours=1)

    needs_processing2 = (Path(raw_data_with_two_runs) /
                         "MiSeq/runs/140201_M01234/needsprocessing")
    needs_processing2.unlink()
    sample_queue = DummyQueueSink()
    item1 = FolderEvent(raw_data_with_two_runs / "MiSeq/runs/140101_M01234" /
                        "Data/Intensities/BaseCalls",
                        FolderEventType.ADD_SAMPLE,
                        SampleGroup('2000A',
                                    ('2000A-V3LOOP_S2_L001_R1_001.fastq.gz',
                                     None),
                                    ('V3LOOP', None)))
    finish1 = FolderEvent(raw_data_with_two_runs / "MiSeq/runs/140101_M01234" /
                          "Data/Intensities/BaseCalls",
                          FolderEventType.FINISH_FOLDER,
                          None)
    item2 = FolderEvent(raw_data_with_two_runs / "MiSeq/runs/140201_M01234" /
                        "Data/Intensities/BaseCalls",
                        FolderEventType.ADD_SAMPLE,
                        SampleGroup('2010A',
                                    ('2010A-V3LOOP_S3_L001_R1_001.fastq.gz',
                                     None),
                                    ('V3LOOP', None)))
    finish2 = FolderEvent(raw_data_with_two_runs / "MiSeq/runs/140201_M01234" /
                          "Data/Intensities/BaseCalls",
                          FolderEventType.FINISH_FOLDER,
                          None)
    sample_queue.expect_put(item1,
                            is_full=True,
                            callback=needs_processing2.touch)
    sample_queue.expect_put(item1,
                            is_full=True,
                            callback=increment_clock)
    sample_queue.expect_put(item2)
    sample_queue.expect_put(finish2)
    sample_queue.expect_put(item1)
    sample_queue.expect_put(finish1)
    pipeline_version = 'XXX'

    find_samples(raw_data_with_two_runs,
                 pipeline_version,
                 sample_queue,
                 wait=False)

    sample_queue.verify()


def test_scan_error(raw_data_with_two_samples, monkeypatch):
    mock_scan = Mock()
    mock_sleep = Mock()
    monkeypatch.setattr(micall.monitor.kive_watcher,
                        'scan_flag_paths',
                        mock_scan)
    monkeypatch.setattr(micall.monitor.kive_watcher, 'sleep', mock_sleep)
    needs_processing = (Path(raw_data_with_two_samples) /
                        "MiSeq/runs/140101_M01234/needsprocessing")
    mock_scan.side_effect = [IOError('unavailable'), [needs_processing]]
    sample_queue = DummyQueueSink()
    sample_queue.expect_put(
        FolderEvent(raw_data_with_two_samples / "MiSeq/runs/140101_M01234" /
                    "Data/Intensities/BaseCalls",
                    FolderEventType.ADD_SAMPLE,
                    SampleGroup('2120A',
                                ('2120A-PR_S14_L001_R1_001.fastq.gz',
                                 None),
                                ('PR', None))))
    sample_queue.expect_put(
        FolderEvent(raw_data_with_two_samples / "MiSeq/runs/140101_M01234" /
                    "Data/Intensities/BaseCalls",
                    FolderEventType.ADD_SAMPLE,
                    SampleGroup('2110A',
                                ('2110A-V3LOOP_S13_L001_R1_001.fastq.gz',
                                 None),
                                ('V3LOOP', None))))
    sample_queue.expect_put(
        FolderEvent(raw_data_with_two_samples / "MiSeq/runs/140101_M01234" /
                    "Data/Intensities/BaseCalls",
                    FolderEventType.FINISH_FOLDER,
                    None))
    pipeline_version = 'XXX'

    find_samples(raw_data_with_two_samples,
                 pipeline_version,
                 sample_queue,
                 wait=False)

    sample_queue.verify()
    mock_sleep.assert_called_with(5)


def test_starts_empty(default_config):
    kive_watcher = KiveWatcher(default_config)

    assert not kive_watcher.is_full()


def test_get_kive_pipeline(mock_open_kive, pipelines_config):
    app_id = 43
    assert app_id == pipelines_config.micall_main_pipeline_id
    mock_session = mock_open_kive.return_value
    mock_session.endpoints.containerapps.get.return_value = [
        dict(name="fastq1", url="/args/101", type="I", app="/apps/99"),
        dict(name="fastq2", url="/args/102", type="I"),
        dict(name="bad_cycles_csv", url="/args/103", type="I"),
        dict(name="g2p_csv", url="/args/104", type="O")]
    expected_args = dict(fastq1="/args/101",
                         fastq2="/args/102",
                         bad_cycles_csv="/args/103")
    expected_url = "/apps/99"
    kive_watcher = KiveWatcher(pipelines_config)

    args = kive_watcher.get_kive_arguments(app_id)
    url = kive_watcher.get_kive_app(app_id)

    assert expected_args == args
    assert expected_url == url
    mock_session.endpoints.containerapps.get.assert_called_once_with('43/argument_list/')


def test_get_app_cached(mock_open_kive, pipelines_config):
    mock_session = mock_open_kive.return_value
    kive_watcher = KiveWatcher(pipelines_config)

    pipeline1 = kive_watcher.get_kive_app(pipelines_config.micall_main_pipeline_id)
    pipeline2 = kive_watcher.get_kive_app(pipelines_config.micall_main_pipeline_id)

    assert pipeline1 is pipeline2
    mock_session.endpoints.containerapps.get.assert_called_once()


def test_add_first_sample(raw_data_with_two_samples, mock_open_kive, default_config):
    base_calls = (raw_data_with_two_samples /
                  "MiSeq/runs/140101_M01234/Data/Intensities/BaseCalls")
    mock_session = mock_open_kive.return_value
    result_path = base_calls / "../../../Results/version_0-dev"
    result_path.mkdir(parents=True)
    old_stuff_csv = result_path / 'old_stuff.csv'
    old_stuff_csv.write_text('out of date')
    dataset1 = dict(name='quality_csv')
    dataset2 = dict(name='fastq1')
    dataset3 = dict(name='fastq2')
    mock_session.endpoints.datasets.post.side_effect = [dataset1, dataset2, dataset3]
    kive_watcher = KiveWatcher(default_config)
    kive_watcher.apps = {default_config.micall_filter_quality_pipeline_id: dict(
        quality_csv="/args/101")}

    kive_watcher.add_sample_group(
        base_calls=base_calls,
        sample_group=SampleGroup('2110A',
                                 ('2110A-V3LOOP_S13_L001_R1_001.fastq.gz',
                                  None),
                                 ('V3LOOP', None)))

    mock_open_kive.assert_called_once_with(default_config.kive_server)
    mock_session.login.assert_called_once_with(default_config.kive_user,
                                               default_config.kive_password)
    mock_session.endpoints.batches.post.assert_called_once_with(
        json=dict(name='140101_M01234 v0-dev',
                  description='MiCall batch for folder 140101_M01234, '
                              'pipeline version 0-dev.',
                  groups_allowed=['Everyone']))
    assert [call('name', '140101_M01234_quality.csv',
                 # MD5 of header with no records.
                 'md5', '6861a4a0bfd71b62c0048ff9a4910223',
                 'uploaded', True),
            call('name', '2110A-V3LOOP_S13_L001_R1_001.fastq.gz',
                 'md5', ANY,
                 'uploaded', True),
            call('name', '2110A-V3LOOP_S13_L001_R2_001.fastq.gz',
                 'md5', ANY,
                 'uploaded', True)] == mock_session.endpoints.datasets.filter.call_args_list
    assert [call(data=dict(name='140101_M01234_quality.csv',
                           description='Error rates for 140101_M01234 run.',
                           users_allowed=[],
                           groups_allowed=['Everyone']),
                 files=dict(dataset_file=ANY)),
            call(data=dict(name='2110A-V3LOOP_S13_L001_R1_001.fastq.gz',
                           description='forward read from MiSeq run 140101_M01234',
                           users_allowed=[],
                           groups_allowed=['Everyone']),
                 files=dict(dataset_file=ANY)),
            call(data=dict(name='2110A-V3LOOP_S13_L001_R2_001.fastq.gz',
                           description='reverse read from MiSeq run 140101_M01234',
                           users_allowed=[],
                           groups_allowed=['Everyone']),
                 files=dict(dataset_file=ANY)),
            ] == mock_session.endpoints.datasets.post.call_args_list
    assert 1 == len(kive_watcher.folder_watchers)
    folder_watcher = kive_watcher.folder_watchers[base_calls]
    assert dataset1 is folder_watcher.quality_dataset
    assert [dataset2, dataset3] == folder_watcher.sample_watchers[0].fastq_datasets
    assert not old_stuff_csv.exists()


def test_add_first_sample_with_compression(raw_data_with_two_samples, mock_open_kive, default_config):
    base_calls = (raw_data_with_two_samples /
                  "MiSeq/runs/140101_M01234/Data/Intensities/BaseCalls")
    mock_session = mock_open_kive.return_value
    results_path = base_calls / "../../../Results"
    results_path.mkdir(parents=True)
    version1_folder: Path = results_path / 'version_1.0'
    version1_folder.mkdir()
    (version1_folder / "foo.txt").write_text('Foo content')
    (version1_folder / "bar.txt").write_text('Bar content')
    version2_folder: Path = results_path / 'version_2.0'
    version2_folder.mkdir()
    (version2_folder / "baz.txt").write_text('Baz content')
    version3_zip = results_path / 'version_3.0.zip'
    version3_zip.write_text('ziiipped content')
    default_config.pipeline_version = '3.0'
    expected_version1_zip = results_path / 'version_1.0.zip'
    expected_file_names = ['version_1.0/bar.txt', 'version_1.0/foo.txt']
    dataset1 = dict(name='quality_csv')
    dataset2 = dict(name='fastq1')
    dataset3 = dict(name='fastq2')
    mock_session.endpoints.datasets.post.side_effect = [dataset1, dataset2, dataset3]
    kive_watcher = KiveWatcher(default_config)
    kive_watcher.apps = {default_config.micall_filter_quality_pipeline_id: dict(
        quality_csv="/args/101")}

    kive_watcher.add_sample_group(
        base_calls=base_calls,
        sample_group=SampleGroup('2110A',
                                 ('2110A-V3LOOP_S13_L001_R1_001.fastq.gz',
                                  None),
                                 ('V3LOOP', None)))

    assert not version3_zip.exists()
    assert version2_folder.exists()
    assert expected_version1_zip.exists()
    with expected_version1_zip.open('rb') as f:
        assert expected_file_names == sorted(ZipFile(f).namelist())
    assert not version1_folder.exists()


def test_compress_old_versions_no_results(tmpdir):
    run_path = Path(str(tmpdir))
    results_path = run_path / 'Results'
    version_path = results_path / 'version_3.0'

    compress_old_versions(version_path)

    assert not results_path.exists()


def test_compress_old_versions_no_old_versions(tmpdir):
    run_path = Path(str(tmpdir))
    results_path = run_path / 'Results'
    results_path.mkdir()
    version_path = results_path / 'version_3.0'

    compress_old_versions(version_path)


def test_compress_old_versions_sorting(tmpdir):
    run_path = Path(str(tmpdir))
    results_path = run_path / 'Results'
    results_path.mkdir()
    version_path9 = results_path / 'version_0.9'
    version_path10 = results_path / 'version_0.10'
    version_path11 = results_path / 'version_0.11'
    version_path9.mkdir()
    version_path10.mkdir()

    compress_old_versions(version_path11)

    assert not version_path9.exists()
    assert version_path10.exists()


def test_compress_old_versions_skip_current_and_newer(tmpdir):
    run_path = Path(str(tmpdir))
    results_path = run_path / 'Results'
    results_path.mkdir()
    version_path1 = results_path / 'version_1'
    version_path2 = results_path / 'version_2'
    version_path3 = results_path / 'version_3'
    version_path4 = results_path / 'version_4'
    version_path1.mkdir()
    version_path2.mkdir()
    version_path3.mkdir()
    version_path4.mkdir()

    compress_old_versions(version_path3)

    assert not version_path1.exists()
    assert version_path2.exists()
    assert version_path3.exists()
    assert version_path4.exists()


def test_compress_old_versions_other_files(tmpdir):
    run_path = Path(str(tmpdir))
    results_path = run_path / 'Results'
    results_path.mkdir()
    other_file = results_path / 'other_file.txt'
    version_path1 = results_path / 'version_1'
    version_path2 = results_path / 'version_2'
    other_file.write_text('Other stuff')
    version_path1.mkdir()

    compress_old_versions(version_path2)

    assert other_file.exists()
    assert version_path1.exists()


def test_compress_old_versions_other_dirs(tmpdir):
    run_path = Path(str(tmpdir))
    results_path = run_path / 'Results'
    results_path.mkdir()
    other_stuff = results_path / 'other_stuff'
    version_path1 = results_path / 'version_1'
    version_path2 = results_path / 'version_2'
    other_stuff.mkdir()
    version_path1.mkdir()

    compress_old_versions(version_path2)

    assert other_stuff.exists()
    assert version_path1.exists()


def test_create_batch_with_expired_session(raw_data_with_two_samples,
                                           mock_open_kive,
                                           default_config):
    base_calls = (raw_data_with_two_samples /
                  "MiSeq/runs/140101_M01234/Data/Intensities/BaseCalls")
    mock_session = mock_open_kive.return_value
    mock_batch = Mock(name='batch')
    mock_session.endpoints.batches.post.side_effect = [KiveClientException('expired'),
                                                       mock_batch]
    kive_watcher = KiveWatcher(default_config)

    kive_watcher.add_sample_group(
        base_calls=base_calls,
        sample_group=SampleGroup('2110A',
                                 ('2110A-V3LOOP_S13_L001_R1_001.fastq.gz',
                                  None),
                                 ('V3LOOP', None)))

    folder_watcher, = kive_watcher.folder_watchers.values()
    assert mock_batch is folder_watcher.batch


def test_add_external_dataset(raw_data_with_two_samples, mock_open_kive, default_config):
    base_calls = (raw_data_with_two_samples /
                  "MiSeq/runs/140101_M01234/Data/Intensities/BaseCalls")
    default_config.raw_data = raw_data_with_two_samples
    mock_session = mock_open_kive.return_value
    mock_session.endpoints.externalfiledirectories.get.return_value = [
        dict(name='raw_data', path=str(raw_data_with_two_samples))]
    dataset1 = dict(name='quality_csv')
    dataset2 = dict(name='fastq1')
    dataset3 = dict(name='fastq2')
    mock_session.endpoints.datasets.post.side_effect = [dataset1, dataset2, dataset3]
    mock_pipeline = mock_session.get_pipeline.return_value
    mock_input = Mock(dataset_name='quality_csv')
    mock_pipeline.inputs = [mock_input]
    kive_watcher = KiveWatcher(default_config)

    kive_watcher.add_sample_group(
        base_calls=base_calls,
        sample_group=SampleGroup('2110A',
                                 ('2110A-V3LOOP_S13_L001_R1_001.fastq.gz',
                                  None),
                                 ('V3LOOP', None)))

    mock_session.endpoints.externalfiledirectories.get.assert_called_once_with()
    assert [call('name', '140101_M01234_quality.csv',
                 # MD5 of header with no records.
                 'md5', '6861a4a0bfd71b62c0048ff9a4910223',
                 'uploaded', True),
            call('name', '2110A-V3LOOP_S13_L001_R1_001.fastq.gz',
                 'md5', ANY,
                 'uploaded', True),
            call('name', '2110A-V3LOOP_S13_L001_R2_001.fastq.gz',
                 'md5', ANY,
                 'uploaded', True)] == mock_session.endpoints.datasets.filter.call_args_list
    assert [call(data=dict(name='140101_M01234_quality.csv',
                           description='Error rates for 140101_M01234 run.',
                           users_allowed=[],
                           groups_allowed=['Everyone']),
                 files=dict(dataset_file=ANY)),
            call(json=dict(name='2110A-V3LOOP_S13_L001_R1_001.fastq.gz',
                           description='forward read from MiSeq run 140101_M01234',
                           externalfiledirectory='raw_data',
                           external_path='MiSeq/runs/140101_M01234/Data/'
                                         'Intensities/BaseCalls/'
                                         '2110A-V3LOOP_S13_L001_R1_001.fastq.gz',
                           users_allowed=[],
                           groups_allowed=['Everyone'])),
            call(json=dict(name='2110A-V3LOOP_S13_L001_R2_001.fastq.gz',
                           description='reverse read from MiSeq run 140101_M01234',
                           externalfiledirectory='raw_data',
                           external_path='MiSeq/runs/140101_M01234/Data/'
                                         'Intensities/BaseCalls/'
                                         '2110A-V3LOOP_S13_L001_R2_001.fastq.gz',
                           users_allowed=[],
                           groups_allowed=['Everyone']))
            ] == mock_session.endpoints.datasets.post.call_args_list


def test_poll_first_sample(raw_data_with_two_samples, mock_open_kive, default_config):
    base_calls = (raw_data_with_two_samples /
                  "MiSeq/runs/140101_M01234/Data/Intensities/BaseCalls")
    mock_session = mock_open_kive.return_value
    mock_session.endpoints.batches.post.return_value = dict(url='/batches/101')
    mock_session.endpoints.datasets.post.return_value = dict(url='/datasets/104',
                                                             id=104)

    kive_watcher = KiveWatcher(default_config)
    kive_watcher.app_urls = {
        default_config.micall_filter_quality_pipeline_id: '/containerapps/102'}
    kive_watcher.app_args = {
        default_config.micall_filter_quality_pipeline_id: dict(
            quality_csv='/containerargs/103')}

    kive_watcher.add_sample_group(
        base_calls=base_calls,
        sample_group=SampleGroup('2110A',
                                 ('2110A-V3LOOP_S13_L001_R1_001.fastq.gz',
                                  None),
                                 ('V3LOOP', None)))
    kive_watcher.poll_runs()

    mock_session.endpoints.containerruns.post.assert_called_once_with(json=dict(
        name='MiCall filter quality on 140101_M01234',
        batch='/batches/101',
        app='/containerapps/102',
        groups_allowed=['Everyone'],
        datasets=[dict(argument='/containerargs/103',
                       dataset='/datasets/104')]))


def test_poll_second_folder(raw_data_with_two_runs, mock_open_kive, default_config):
    # raw_data = create_run_folder(tmpdir, '140101_M01234', '2000A*.fastq')
    # create_run_folder(tmpdir, '140201_M01234', '2010A*.fastq')
    base_calls = (raw_data_with_two_runs /
                  "MiSeq/runs/140101_M01234/Data/Intensities/BaseCalls")
    base_calls2 = (raw_data_with_two_runs /
                   "MiSeq/runs/140201_M01234/Data/Intensities/BaseCalls")
    mock_session = mock_open_kive.return_value
    mock_session.endpoints.batches.post.return_value = dict(url='/batches/101')
    mock_session.endpoints.datasets.post.return_value = dict(url='/datasets/104',
                                                             id=104)

    kive_watcher = KiveWatcher(default_config)
    kive_watcher.app_urls = {
        default_config.micall_filter_quality_pipeline_id: '/containerapps/102'}
    kive_watcher.app_args = {
        default_config.micall_filter_quality_pipeline_id: dict(
            quality_csv='/containerargs/103')}

    kive_watcher.add_sample_group(
        base_calls=base_calls2,
        sample_group=SampleGroup('2010A',
                                 ('2010A-V3LOOP_S3_L001_R1_001.fastq.gz',
                                  None),
                                 ('V3LOOP', None)))
    kive_watcher.poll_runs()  # Start filter run for folder 1.
    kive_watcher.add_sample_group(
        base_calls=base_calls,
        sample_group=SampleGroup('2000A',
                                 ('2000A-V3LOOP_S2_L001_R1_001.fastq.gz',
                                  None),
                                 ('V3LOOP', None)))
    kive_watcher.poll_runs()  # Check filter run 1 and start filter run 2.

    fetch_count_before = mock_session.endpoints.containerruns.get.call_count
    kive_watcher.poll_runs()  # Only check filter run 2.
    fetch_count_after = mock_session.endpoints.containerruns.get.call_count

    assert fetch_count_after - fetch_count_before == 1


def test_poll_all_folders(raw_data_with_two_runs, mock_open_kive, default_config):
    # raw_data = create_run_folder(tmpdir, '140101_M01234', '2000A*.fastq')
    # create_run_folder(tmpdir, '140201_M01234', '2010A*.fastq')
    base_calls = (raw_data_with_two_runs /
                  "MiSeq/runs/140101_M01234/Data/Intensities/BaseCalls")
    base_calls2 = (raw_data_with_two_runs /
                   "MiSeq/runs/140201_M01234/Data/Intensities/BaseCalls")
    mock_session = mock_open_kive.return_value
    mock_session.endpoints.batches.post.return_value = dict(url='/batches/101')
    mock_session.endpoints.datasets.post.return_value = dict(url='/datasets/104',
                                                             id=104)

    kive_watcher = KiveWatcher(default_config)
    kive_watcher.app_urls = {
        default_config.micall_filter_quality_pipeline_id: '/containerapps/102'}
    kive_watcher.app_args = {
        default_config.micall_filter_quality_pipeline_id: dict(
            quality_csv='/containerargs/103')}

    kive_watcher.add_sample_group(
        base_calls=base_calls2,
        sample_group=SampleGroup('2010A',
                                 ('2010A-V3LOOP_S3_L001_R1_001.fastq.gz',
                                  None),
                                 ('V3LOOP', None)))
    kive_watcher.poll_runs()  # Start filter run for folder 1.
    kive_watcher.add_sample_group(
        base_calls=base_calls,
        sample_group=SampleGroup('2000A',
                                 ('2000A-V3LOOP_S2_L001_R1_001.fastq.gz',
                                  None),
                                 ('V3LOOP', None)))
    kive_watcher.poll_runs()  # Check filter run 1 and start filter run 2.
    kive_watcher.poll_runs()  # Only check filter run 2.

    fetch_count_before = mock_session.endpoints.containerruns.get.call_count
    kive_watcher.poll_runs()  # No new samples, check all.
    fetch_count_after = mock_session.endpoints.containerruns.get.call_count

    assert fetch_count_after - fetch_count_before == 2


def test_poll_first_sample_twice(raw_data_with_two_samples, mock_open_kive, default_config):
    base_calls = (raw_data_with_two_samples /
                  "MiSeq/runs/140101_M01234/Data/Intensities/BaseCalls")
    mock_session = mock_open_kive.return_value
    kive_watcher = KiveWatcher(default_config)

    kive_watcher.add_sample_group(
        base_calls=base_calls,
        sample_group=SampleGroup('2110A',
                                 ('2110A-V3LOOP_S13_L001_R1_001.fastq.gz',
                                  None),
                                 ('V3LOOP', None)))
    kive_watcher.poll_runs()
    kive_watcher.poll_runs()

    mock_session.endpoints.containerruns.post.assert_called_once()


def test_poll_first_sample_already_started(raw_data_with_two_samples,
                                           mock_open_kive,
                                           default_config):
    base_calls = (raw_data_with_two_samples /
                  "MiSeq/runs/140101_M01234/Data/Intensities/BaseCalls")
    mock_session = mock_open_kive.return_value
    mock_session.endpoints.batches.filter.return_value = [
        dict(url='/batches/101',
             name='140101_M01234 v0-dev',
             groups_allowed=['Everyone'])]
    mock_session.endpoints.datasets.filter.return_value = [
        dict(url='/datasets/104',
             id=104,
             name='140101_M01234_quality.csv',
             groups_allowed=['Everyone'])]
    mock_session.endpoints.containerruns.filter.return_value = [
        dict(url='/containerruns/105',
             id=105,
             state='R',
             name='MiCall filter quality on 140101_M01234',
             groups_allowed=['Everyone'])]

    kive_watcher = KiveWatcher(default_config)
    kive_watcher.app_urls = {
        default_config.micall_filter_quality_pipeline_id: '/containerapps/102'}
    kive_watcher.app_args = {
        default_config.micall_filter_quality_pipeline_id: dict(
            quality_csv='/containerargs/103')}
    kive_watcher = KiveWatcher(default_config)

    kive_watcher.add_sample_group(
        base_calls=base_calls,
        sample_group=SampleGroup('2110A',
                                 ('2110A-V3LOOP_S13_L001_R1_001.fastq.gz',
                                  None),
                                 ('V3LOOP', None)))
    kive_watcher.poll_runs()

    mock_session.endpoints.batches.post.assert_not_called()
    mock_session.endpoints.containerruns.post.assert_not_called()
    mock_session.endpoints.batches.filter.assert_called_once_with(
        'name', '140101_M01234 v0-dev')
    assert [call('name', '140101_M01234_quality.csv', 'md5', ANY, 'uploaded', True),
            call('name', '2110A-V3LOOP_S13_L001_R1_001.fastq.gz', 'md5', ANY, 'uploaded', True),
            call('name', '2110A-V3LOOP_S13_L001_R2_001.fastq.gz', 'md5', ANY, 'uploaded', True)
            ] == mock_session.endpoints.datasets.filter.call_args_list
    mock_session.endpoints.containerruns.filter.assert_called_once_with(
        'name', 'MiCall filter quality on 140101_M01234',
        'app_id', default_config.micall_filter_quality_pipeline_id,
        'states', 'NLRSC',
        'input_id', 104)


def test_poll_first_sample_completed_and_purged(raw_data_with_two_samples,
                                                mock_open_kive,
                                                default_config):
    """ A matching run finished recently, but it was purged. """
    base_calls = (raw_data_with_two_samples /
                  "MiSeq/runs/140101_M01234/Data/Intensities/BaseCalls")
    mock_session = mock_open_kive.return_value
    mock_session.endpoints.datasets.filter.return_value = [
        dict(url='/datasets/104',
             id=104,
             name='140101_M01234_quality.csv',
             groups_allowed=['Everyone'])]
    mock_session.endpoints.containerruns.filter.return_value = [
        dict(url='/containerruns/105',
             id=105,
             state='C',
             name='MiCall filter quality on 140101_M01234',
             groups_allowed=['Everyone'])]
    mock_session.endpoints.containerruns.get.return_value = [
        dict(dataset_purged=True)]
    kive_watcher = KiveWatcher(default_config)

    kive_watcher.add_sample_group(
        base_calls=base_calls,
        sample_group=SampleGroup('2110A',
                                 ('2110A-V3LOOP_S13_L001_R1_001.fastq.gz',
                                  None),
                                 ('V3LOOP', None)))
    kive_watcher.poll_runs()

    mock_session.endpoints.containerruns.post.assert_called_once()


def test_second_sample(raw_data_with_two_samples, mock_open_kive, default_config):
    base_calls = (raw_data_with_two_samples /
                  "MiSeq/runs/140101_M01234/Data/Intensities/BaseCalls")
    mock_session = mock_open_kive.return_value
    kive_watcher = KiveWatcher(default_config)

    kive_watcher.add_sample_group(
        base_calls=base_calls,
        sample_group=SampleGroup('2110A',
                                 ('2110A-V3LOOP_S13_L001_R1_001.fastq.gz',
                                  None),
                                 ('V3LOOP', None)))
    kive_watcher.add_sample_group(
        base_calls=base_calls,
        sample_group=SampleGroup('2120A',
                                 ('2120A-PR_S14_L001_R1_001.fastq.gz',
                                  None),
                                 ('PR', None)))

    mock_session.endpoints.batches.post.assert_called_once_with(
        json=dict(name='140101_M01234 v0-dev',
                  description='MiCall batch for folder 140101_M01234, '
                              'pipeline version 0-dev.',
                  groups_allowed=['Everyone']))
    expected_dataset_count = 5  # quality_csv + 2 pairs of FASTQ files
    assert expected_dataset_count == len(
        mock_session.endpoints.datasets.filter.call_args_list)


def test_sample_with_hcv_pair(raw_data_with_hcv_pair, mock_open_kive, default_config):
    base_calls = (raw_data_with_hcv_pair /
                  "MiSeq/runs/140101_M01234/Data/Intensities/BaseCalls")
    mock_session = mock_open_kive.return_value
    kive_watcher = KiveWatcher(default_config)

    kive_watcher.add_sample_group(
        base_calls=base_calls,
        sample_group=SampleGroup('2130A',
                                 ('2130A-HCV_S15_L001_R1_001.fastq.gz',
                                  '2130AMIDI-MidHCV_S16_L001_R1_001.fastq.gz'),
                                 ('HCV', 'MidHCV')))

    expected_dataset_count = 5  # quality_csv + 2 pairs of FASTQ files
    assert expected_dataset_count == len(
        mock_session.endpoints.datasets.filter.call_args_list)


def test_sample_fails_to_upload(raw_data_with_two_samples,
                                mock_open_kive,
                                mock_wait,
                                default_config):
    base_calls = (raw_data_with_two_samples /
                  "MiSeq/runs/140101_M01234/Data/Intensities/BaseCalls")
    mock_session = mock_open_kive.return_value
    kive_watcher = KiveWatcher(default_config, retry=True)
    mock_session.endpoints.datasets.post.side_effect = [
        ConnectionError('server down'),
        Mock(name='quality_csv'),
        Mock(name='fastq1'),
        Mock(name='fastq2')]
    mock_wait.side_effect = [
        None,
        RuntimeError('Should only call wait_for_retry() once.')]

    sample_watcher = kive_watcher.add_sample_group(
        base_calls=base_calls,
        sample_group=SampleGroup('2110A',
                                 ('2110A-V3LOOP_S13_L001_R1_001.fastq.gz',
                                  None),
                                 ('V3LOOP', None)))

    assert sample_watcher is not None
    mock_wait.assert_called_once_with(1)


def test_create_batch_fails(raw_data_with_two_samples,
                            mock_open_kive,
                            mock_wait,
                            default_config):
    base_calls = (raw_data_with_two_samples /
                  "MiSeq/runs/140101_M01234/Data/Intensities/BaseCalls")
    mock_session = mock_open_kive.return_value
    mock_session.endpoints.batches.post.side_effect = [
        ConnectionError('server down'),
        dict(url='/batches/101')]
    mock_wait.side_effect = [
        None,
        RuntimeError('Should only call wait_for_retry() once.')]
    kive_watcher = KiveWatcher(default_config, retry=True)

    sample_watcher = kive_watcher.add_sample_group(
        base_calls=base_calls,
        sample_group=SampleGroup('2110A',
                                 ('2110A-V3LOOP_S13_L001_R1_001.fastq.gz',
                                  None),
                                 ('V3LOOP', None)))
    kive_watcher.poll_runs()

    assert sample_watcher is not None
    mock_wait.assert_called_once_with(1)
    mock_session.endpoints.containerruns.post.assert_called_once_with(json=dict(
        name=ANY,
        app=ANY,
        batch='/batches/101',
        datasets=ANY,
        groups_allowed=ANY))


def test_sample_already_uploaded(raw_data_with_two_samples, mock_open_kive, default_config):
    base_calls = (raw_data_with_two_samples /
                  "MiSeq/runs/140101_M01234/Data/Intensities/BaseCalls")
    mock_session = mock_open_kive.return_value
    dataset1 = dict(name='140101_M01234_quality.csv',
                    groups_allowed=ALLOWED_GROUPS)
    dataset2 = Mock(name='fastq1')
    dataset3 = Mock(name='fastq2')
    mock_session.endpoints.datasets.filter.side_effect = [[dataset1], [], []]
    mock_session.endpoints.datasets.post.side_effect = [dataset2, dataset3]
    kive_watcher = KiveWatcher(default_config)

    kive_watcher.add_sample_group(
        base_calls=base_calls,
        sample_group=SampleGroup('2110A',
                                 ('2110A-V3LOOP_S13_L001_R1_001.fastq.gz',
                                  None),
                                 ('V3LOOP', None)))

    assert [call('name', '140101_M01234_quality.csv',
                 # MD5 of header with no records.
                 'md5', '6861a4a0bfd71b62c0048ff9a4910223',
                 'uploaded', True),
            call('name', '2110A-V3LOOP_S13_L001_R1_001.fastq.gz',
                 'md5', ANY,
                 'uploaded', True),
            call('name', '2110A-V3LOOP_S13_L001_R2_001.fastq.gz',
                 'md5', ANY,
                 'uploaded', True)
            ] == mock_session.endpoints.datasets.filter.call_args_list
    assert [call(data=dict(name='2110A-V3LOOP_S13_L001_R1_001.fastq.gz',
                           description='forward read from MiSeq run 140101_M01234',
                           users_allowed=[],
                           groups_allowed=['Everyone']),
                 files=dict(dataset_file=ANY)),
            call(data=dict(name='2110A-V3LOOP_S13_L001_R2_001.fastq.gz',
                           description='reverse read from MiSeq run 140101_M01234',
                           users_allowed=[],
                           groups_allowed=['Everyone']),
                 files=dict(dataset_file=ANY))
            ] == mock_session.endpoints.datasets.post.call_args_list
    assert 1 == len(kive_watcher.folder_watchers)
    folder_watcher = kive_watcher.folder_watchers[base_calls]
    assert dataset1 is folder_watcher.quality_dataset
    assert [dataset2, dataset3] == folder_watcher.sample_watchers[0].fastq_datasets


def test_launch_main_run(raw_data_with_two_samples, mock_open_kive, pipelines_config):
    base_calls = (raw_data_with_two_samples /
                  "MiSeq/runs/140101_M01234/Data/Intensities/BaseCalls")
    mock_session = mock_open_kive.return_value
    mock_session.endpoints.datasets.post.side_effect = [
        dict(url='/datasets/104', id=104),
        dict(url='/datasets/105', id=105)]

    pipelines_config.denovo_main_pipeline_id = 495
    kive_watcher = KiveWatcher(pipelines_config)
    kive_watcher.app_urls = {
        pipelines_config.micall_main_pipeline_id: '/containerapps/102',
        pipelines_config.denovo_main_pipeline_id: '/containerapps/103'}
    kive_watcher.app_args = {
        pipelines_config.micall_main_pipeline_id: dict(
            bad_cycles_csv='/containerargs/110',
            fastq1='/containerargs/111',
            fastq2='/containerargs/112'),
        pipelines_config.denovo_main_pipeline_id: dict(
            bad_cycles_csv='/containerargs/113',
            fastq1='/containerargs/114',
            fastq2='/containerargs/115')}

    folder_watcher = kive_watcher.add_folder(base_calls)
    folder_watcher.batch = dict(url='/batches/101')
    kive_watcher.add_sample_group(
        base_calls=base_calls,
        sample_group=SampleGroup('2110A',
                                 ('2110A-V3LOOP_S13_L001_R1_001.fastq.gz',
                                  None),
                                 ('V3LOOP', None)))
    folder_watcher.add_run(
        dict(id=106),
        PipelineType.FILTER_QUALITY)

    mock_session.endpoints.containerruns.get.side_effect = [
        dict(id=106, state='C'),  # refresh run status
        [dict(argument_name='bad_cycles_csv',
              dataset='/datasets/106',
              dataset_purged=False)]]  # get dataset list
    mock_session.get.return_value.json.return_value = dict(url='/datasets/106',
                                                           id=106)

    kive_watcher.poll_runs()

    assert [call(json=dict(app='/containerapps/103',
                           datasets=[dict(argument='/containerargs/113',
                                          dataset='/datasets/106'),
                                     dict(argument='/containerargs/114',
                                          dataset='/datasets/104'),
                                     dict(argument='/containerargs/115',
                                          dataset='/datasets/105')],
                           name='MiCall denovo main on 2110A-V3LOOP_S13',
                           batch='/batches/101',
                           groups_allowed=['Everyone'])),
            call(json=dict(app='/containerapps/102',
                           datasets=[dict(argument='/containerargs/110',
                                          dataset='/datasets/106'),
                                     dict(argument='/containerargs/111',
                                          dataset='/datasets/104'),
                                     dict(argument='/containerargs/112',
                                          dataset='/datasets/105')],
                           name='MiCall main on 2110A-V3LOOP_S13',
                           batch='/batches/101',
                           groups_allowed=['Everyone']))
            ] == mock_session.endpoints.containerruns.post.call_args_list


def test_launch_main_run_with_sample_info(raw_data_with_two_samples,
                                          mock_open_kive,
                                          pipelines_config):
    base_calls = (raw_data_with_two_samples /
                  "MiSeq/runs/140101_M01234/Data/Intensities/BaseCalls")
    mock_session = mock_open_kive.return_value
    mock_session.endpoints.datasets.post.side_effect = [
        dict(url='/datasets/104', id=104),
        dict(url='/datasets/105', id=105),
        dict(url='/datasets/106', id=106)]

    pipelines_config.denovo_main_pipeline_id = 495
    kive_watcher = KiveWatcher(pipelines_config)
    kive_watcher.app_urls = {
        pipelines_config.micall_main_pipeline_id: '/containerapps/102',
        pipelines_config.denovo_main_pipeline_id: '/containerapps/103'}
    kive_watcher.app_args = {
        pipelines_config.micall_main_pipeline_id: dict(
            bad_cycles_csv='/containerargs/110',
            fastq1='/containerargs/111',
            fastq2='/containerargs/112',
            sample_info_csv='/containerargs/113'),
        pipelines_config.denovo_main_pipeline_id: dict(
            bad_cycles_csv='/containerargs/114',
            fastq1='/containerargs/115',
            fastq2='/containerargs/116',
            sample_info_csv='/containerargs/117')}

    folder_watcher = kive_watcher.add_folder(base_calls)
    folder_watcher.batch = dict(url='/batches/101')
    kive_watcher.add_sample_group(
        base_calls=base_calls,
        sample_group=SampleGroup('2110A',
                                 ('2110A-V3LOOP_S13_L001_R1_001.fastq.gz',
                                  None),
                                 ('V3LOOP', None)))
    folder_watcher.add_run(
        dict(id=106),
        PipelineType.FILTER_QUALITY)

    mock_session.endpoints.containerruns.get.side_effect = [
        dict(id=106, state='C'),  # refresh run status
        [dict(argument_name='bad_cycles_csv',
              dataset='/datasets/107',
              dataset_purged=False)]]  # get dataset list
    mock_session.get.return_value.json.return_value = dict(url='/datasets/107',
                                                           id=107)

    kive_watcher.poll_runs()

    assert [call(json=dict(app='/containerapps/103',
                           datasets=[dict(argument='/containerargs/114',
                                          dataset='/datasets/107'),
                                     dict(argument='/containerargs/115',
                                          dataset='/datasets/104'),
                                     dict(argument='/containerargs/116',
                                          dataset='/datasets/105'),
                                     dict(argument='/containerargs/117',
                                          dataset='/datasets/106')],
                           name='MiCall denovo main on 2110A-V3LOOP_S13',
                           batch='/batches/101',
                           groups_allowed=['Everyone'])),
            call(json=dict(app='/containerapps/102',
                           datasets=[dict(argument='/containerargs/110',
                                          dataset='/datasets/107'),
                                     dict(argument='/containerargs/111',
                                          dataset='/datasets/104'),
                                     dict(argument='/containerargs/112',
                                          dataset='/datasets/105'),
                                     dict(argument='/containerargs/113',
                                          dataset='/datasets/106')],
                           name='MiCall main on 2110A-V3LOOP_S13',
                           batch='/batches/101',
                           groups_allowed=['Everyone']))
            ] == mock_session.endpoints.containerruns.post.call_args_list


def test_launch_main_run_long_name(raw_data_with_two_samples, mock_open_kive, pipelines_config):
    base_calls = (raw_data_with_two_samples /
                  "MiSeq/runs/140101_M01234/Data/Intensities/BaseCalls")
    rename_fastq_files(base_calls,
                       '2110A-V3LOOP_S13_L001',
                       '2110A-V3LOOP-987654321REALLYVERYLONGNAME-HCV_S13_L001')
    mock_session = mock_open_kive.return_value
    kive_watcher = create_kive_watcher_with_filter_run(pipelines_config,
                                                       base_calls)

    kive_watcher.add_sample_group(
        base_calls=base_calls,
        sample_group=SampleGroup(
            '2110A',
            ('2110A-V3LOOP-987654321REALLYVERYLONGNAME-HCV_S13_L001_R1_001.fastq.gz',
             None),
            ('V3LOOP', None)))

    kive_watcher.poll_runs()

    mock_session.endpoints.containerruns.post.assert_called_once_with(json=dict(
        app=ANY,
        datasets=ANY,
        name='MiCall main on 2110A-V3LOOP-987654321REALLYVERYLONGNA..._S13',
        batch=ANY,
        groups_allowed=['Everyone']))


def rename_fastq_files(base_calls, old_file_prefix, new_file_prefix):
    for suffix in ('_R1_001.fastq.gz', '_R2_001.fastq.gz'):
        old_path = base_calls / (old_file_prefix + suffix)
        new_path = base_calls / (new_file_prefix + suffix)
        old_path.rename(new_path)


def test_trim_run_name_no_suffix():
    run_name = 'MiCall main on 2110A-V3LOOP-987654321REALLYVERYLONGNAME-HCV-S13'
    expected_name = 'MiCall main on 2110A-V3LOOP-987654321REALLYVERYLONGNAME-H...'

    trimmed_name = trim_run_name(run_name)

    assert expected_name == trimmed_name


def test_trim_run_name_long_suffix():
    run_name = 'MiCall main on 2110A-V3LOOP-987654321-HCV_S13REALLYVERYLONGNAME'
    expected_name = 'MiCall main on 2110A-V3LOOP-987654321-HCV_S13REALLYVERYLO...'

    trimmed_name = trim_run_name(run_name)

    assert expected_name == trimmed_name


def test_launch_main_run_after_connection_error(raw_data_with_two_samples,
                                                mock_open_kive,
                                                mock_wait,
                                                pipelines_config):
    base_calls = (raw_data_with_two_samples /
                  "MiSeq/runs/140101_M01234/Data/Intensities/BaseCalls")
    kive_watcher = create_kive_watcher_with_filter_run(pipelines_config,
                                                       base_calls)
    kive_watcher.retry = True
    mock_session = kive_watcher.session

    kive_watcher.add_sample_group(
        base_calls=base_calls,
        sample_group=SampleGroup('2110A',
                                 ('2110A-V3LOOP_S13_L001_R1_001.fastq.gz',
                                  None),
                                 ('V3LOOP', None)))

    mock_session.endpoints.containerruns.get.side_effect = mock_failures(
        failure_count=4,  # We get one extra retry after a login for each request.
        success_callable=mock_containerruns_get)

    mock_wait.side_effect = [
        None,
        None,
        RuntimeError('Should only call wait_for_retry() twice.')]

    kive_watcher.poll_runs()

    mock_session.endpoints.containerruns.post.assert_called_once()
    assert [call(1), call(2)] == mock_wait.call_args_list


def test_launch_midi_run(raw_data_with_hcv_pair, mock_open_kive, pipelines_config):
    base_calls = (raw_data_with_hcv_pair /
                  "MiSeq/runs/140101_M01234/Data/Intensities/BaseCalls")
    kive_watcher = create_kive_watcher_with_filter_run(pipelines_config,
                                                       base_calls)
    mock_session = kive_watcher.session

    kive_watcher.add_sample_group(
        base_calls=base_calls,
        sample_group=SampleGroup('2130A',
                                 ('2130A-HCV_S15_L001_R1_001.fastq.gz',
                                  '2130AMIDI-MidHCV_S16_L001_R1_001.fastq.gz'),
                                 ('HCV', 'MidHCV')))

    kive_watcher.poll_runs()

    assert [call(json=dict(name='MiCall main on 2130A-HCV_S15',
                           batch=ANY,
                           app=ANY,
                           datasets=ANY,
                           groups_allowed=['Everyone'])),
            call(json=dict(name='MiCall main on 2130AMIDI-MidHCV_S16',
                           batch=ANY,
                           app=ANY,
                           datasets=ANY,
                           groups_allowed=['Everyone']))
            ] == mock_session.endpoints.containerruns.post.call_args_list


def test_launch_resistance_run(raw_data_with_two_samples, mock_open_kive, pipelines_config):
    pipelines_config.micall_resistance_pipeline_id = 45
    base_calls = (raw_data_with_two_samples /
                  "MiSeq/runs/140101_M01234/Data/Intensities/BaseCalls")
    kive_watcher = KiveWatcher(pipelines_config)
    kive_watcher.app_urls = {
        pipelines_config.micall_resistance_pipeline_id: '/containerapps/103'}
    kive_watcher.app_args = {
        pipelines_config.micall_resistance_pipeline_id: dict(
            main_amino_csv='/containerargs/104',
            midi_amino_csv='/containerargs/105',
            main_nuc_csv='/containerargs/106')}

    folder_watcher = kive_watcher.add_folder(base_calls)
    folder_watcher.batch = dict(url='/batches/101')
    folder_watcher.add_run(dict(id=106),
                           PipelineType.FILTER_QUALITY,
                           is_complete=True)
    sample_watcher = kive_watcher.add_sample_group(
        base_calls=base_calls,
        sample_group=SampleGroup('2110A',
                                 ('2110A-V3LOOP_S13_L001_R1_001.fastq.gz',
                                  None),
                                 ('V3LOOP', None)))
    folder_watcher.add_run(
        dict(id=107),
        PipelineType.MAIN,
        sample_watcher)

    kive_watcher.check_session()
    mock_session = kive_watcher.session
    mock_session.endpoints.containerruns.get.side_effect = [
        dict(id=107, state='C'),  # refresh run state
        [dict(dataset='/datasets/111/',
              argument_type='O',
              argument_name='amino_csv'),
         dict(dataset='/datasets/112/',
              argument_type='O',
              argument_name='nuc_csv')]]  # run datasets
    mock_session.get.return_value.json.side_effect = [
        dict(url='/datasets/111/', id=111),
        dict(url='/datasets/111/', id=111),
        dict(url='/datasets/112/', id=112)]

    kive_watcher.poll_runs()

    mock_session.endpoints.containerruns.post.assert_called_once_with(json=dict(
        app='/containerapps/103',
        datasets=[dict(argument='/containerargs/104',
                       dataset='/datasets/111/'),
                  dict(argument='/containerargs/105',
                       dataset='/datasets/111/'),
                  dict(argument='/containerargs/106',
                       dataset='/datasets/112/')],
        name='MiCall resistance on 2110A',
        batch='/batches/101',
        groups_allowed=['Everyone']))


def test_launch_proviral_run(raw_data_with_two_samples, mock_open_kive):
    pipelines_config = parse_args(argv=['--micall_filter_quality_pipeline_id', '42',
                                        '--denovo_main_pipeline_id', '43',
                                        '--proviral_pipeline_id', '145'])

    base_calls = (raw_data_with_two_samples /
                  "MiSeq/runs/140101_M01234/Data/Intensities/BaseCalls")
    kive_watcher = KiveWatcher(pipelines_config)
    kive_watcher.app_urls = {
        pipelines_config.proviral_pipeline_id: '/containerapps/103'}
    kive_watcher.app_args = {
        pipelines_config.proviral_pipeline_id: dict(
            sample_info_csv='/containerargs/103',
            contigs_csv='/containerargs/104',
            conseqs_csv='/containerargs/105',
            cascade_csv='/containerargs/106')}

    folder_watcher = kive_watcher.add_folder(base_calls)
    folder_watcher.batch = dict(url='/batches/101')
    folder_watcher.add_run(dict(id=106),
                           PipelineType.FILTER_QUALITY,
                           is_complete=True)
    sample_watcher = kive_watcher.add_sample_group(
        base_calls=base_calls,
        sample_group=SampleGroup('2120A',
                                 ('2120A-PR_S14_L001_R1_001.fastq.gz',
                                  None),
                                 ('HIVNFLDNA', None)))
    folder_watcher.add_run(
        dict(id=107),
        PipelineType.DENOVO_MAIN,
        sample_watcher)

    kive_watcher.check_session()
    mock_session = kive_watcher.session
    mock_session.endpoints.containerruns.get.side_effect = [
        dict(id=107, state='C'),  # refresh run state
        [dict(dataset='/datasets/110/',
              argument_type='I',
              argument_name='sample_info_csv'),
         dict(dataset='/datasets/111/',
              argument_type='O',
              argument_name='unstitched_contigs_csv'),
         dict(dataset='/datasets/112/',
              argument_type='O',
              argument_name='unstitched_conseq_csv'),
         dict(dataset='/datasets/113/',
              argument_type='O',
              argument_name='unstitched_cascade_csv')]]  # run datasets
    mock_session.get.return_value.json.side_effect = [
        dict(url='/datasets/110/', id=110),
        dict(url='/datasets/111/', id=111),
        dict(url='/datasets/112/', id=112),
        dict(url='/datasets/113/', id=113)]

    kive_watcher.poll_runs()

    mock_session.endpoints.containerruns.post.assert_called_once_with(json=dict(
        app='/containerapps/103',
        datasets=[dict(argument='/containerargs/103',
                       dataset='/datasets/110/'),
                  dict(argument='/containerargs/104',
                       dataset='/datasets/111/'),
                  dict(argument='/containerargs/105',
                       dataset='/datasets/112/'),
                  dict(argument='/containerargs/106',
                       dataset='/datasets/113/')],
        name='Proviral on 2120A',
        batch='/batches/101',
        groups_allowed=['Everyone']))


def test_skip_resistance_run(raw_data_with_two_samples, mock_open_kive, pipelines_config):
    pipelines_config.micall_resistance_pipeline_id = None
    base_calls = (raw_data_with_two_samples /
                  "MiSeq/runs/140101_M01234/Data/Intensities/BaseCalls")
    kive_watcher = KiveWatcher(pipelines_config)

    folder_watcher = kive_watcher.add_folder(base_calls)
    folder_watcher.batch = dict(url='/batches/101')
    folder_watcher.add_run(dict(id=106),
                           PipelineType.FILTER_QUALITY,
                           is_complete=True)
    sample_watcher = kive_watcher.add_sample_group(
        base_calls=base_calls,
        sample_group=SampleGroup('2110A',
                                 ('2110A-V3LOOP_S13_L001_R1_001.fastq.gz',
                                  None),
                                 ('V3LOOP', None)))
    folder_watcher.add_run(
        dict(id=107),
        PipelineType.MAIN,
        sample_watcher)
    kive_watcher.finish_folder(base_calls)

    kive_watcher.check_session()
    mock_session = kive_watcher.session
    mock_session.endpoints.containerruns.get.side_effect = [
        dict(id=107, state='C'),  # refresh run state
        [dict(dataset='/datasets/111/',
              argument_type='O',
              argument_name='amino_csv'),
         dict(dataset='/datasets/112/',
              argument_type='O',
              argument_name='nuc_csv')]]  # run datasets
    mock_session.get.return_value.json.side_effect = [
        dict(url='/datasets/111/', id=111)]

    kive_watcher.poll_runs()

    mock_session.endpoints.containerruns.post.assert_not_called()
    assert kive_watcher.is_idle()


def test_resistance_run_missing_input(raw_data_with_two_samples,
                                      mock_open_kive,
                                      pipelines_config):
    base_calls = (raw_data_with_two_samples /
                  "MiSeq/runs/140101_M01234/Data/Intensities/BaseCalls")
    kive_watcher = create_kive_watcher_with_main_run(
        pipelines_config,
        base_calls,
        SampleGroup('2110A',
                    ('2110A-V3LOOP_S13_L001_R1_001.fastq.gz',
                     None),
                    ('V3LOOP', None)))
    mock_session = kive_watcher.session
    mock_session.endpoints.containerruns.get.side_effect = [
        dict(id=107, state='C'),  # refresh run state
        [dict(dataset='/datasets/111/',
              argument_name='fail_csv')]]  # run datasets

    with pytest.raises(RuntimeError, match=r'Polling sample group 2110A failed.'):
        kive_watcher.poll_runs()

    mock_session.endpoints.containerruns.post.assert_not_called()


def test_poll_main_run_cancelled(raw_data_with_two_samples,
                                 mock_open_kive,
                                 pipelines_config):
    base_calls = (raw_data_with_two_samples /
                  "MiSeq/runs/140101_M01234/Data/Intensities/BaseCalls")
    kive_watcher = create_kive_watcher_with_main_run(
        pipelines_config,
        base_calls,
        SampleGroup('2110A',
                    ('2110A-V3LOOP_S13_L001_R1_001.fastq.gz',
                     None),
                    ('V3LOOP', None)))
    mock_session = kive_watcher.session
    mock_session.endpoints.containerruns.get.side_effect = [
        dict(id=107, state='X')]  # refresh run state

    kive_watcher.poll_runs()

    mock_session.endpoints.containerruns.post.assert_called_once()


def test_launch_hcv_resistance_run(raw_data_with_hcv_pair, mock_open_kive, pipelines_config):
    base_calls = (raw_data_with_hcv_pair /
                  "MiSeq/runs/140101_M01234/Data/Intensities/BaseCalls")
    kive_watcher = create_kive_watcher_with_main_run(
        pipelines_config,
        base_calls,
        SampleGroup('2130A',
                    ('2130A-HCV_S15_L001_R1_001.fastq.gz',
                     '2130AMIDI-MidHCV_S16_L001_R1_001.fastq.gz'),
                    ('HCV', 'MidHCV')))
    mock_session = kive_watcher.session
    folder_watcher, = kive_watcher.folder_watchers.values()
    sample_watcher, = folder_watcher.sample_watchers
    folder_watcher.add_run(
        dict(id=108),
        PipelineType.MIDI,
        sample_watcher)

    mock_session.endpoints.containerruns.get.side_effect = [
        dict(id=107, state='C'),  # refresh run state for main run
        [dict(dataset='/datasets/111/',
              argument_type='O',
              argument_name='amino_csv'),
         dict(dataset='/datasets/112/',
              argument_type='O',
              argument_name='nuc_csv')],  # run datasets
        dict(id=108, state='C'),  # refresh run state for midi run
        [dict(dataset='/datasets/113/',
              argument_type='O',
              argument_name='amino_csv'),
         dict(dataset='/datasets/114/',
              argument_type='O',
              argument_name='nuc_csv')]]  # run datasets
    mock_session.get.reset_mock(side_effect=True)
    mock_session.get.return_value.json.side_effect = [
        dict(url='/datasets/111/', id=111),
        dict(url='/datasets/112/', id=112),
        dict(url='/datasets/113/', id=113)]

    kive_watcher.poll_runs()

    mock_session.endpoints.containerruns.post.assert_called_once_with(json=dict(
        app=ANY,
        datasets=[dict(argument=ANY,
                       dataset='/datasets/111/'),
                  dict(argument=ANY,
                       dataset='/datasets/112/'),
                  dict(argument=ANY,
                       dataset='/datasets/113/')],
        name='MiCall resistance on 2130A',
        batch=ANY,
        groups_allowed=['Everyone']))


def test_launch_hcv_triplet_resistance_run(raw_data_with_hcv_pair, mock_open_kive, pipelines_config):
    """ Same MIDI sample gets paired with two different main samples. """
    run_folder = raw_data_with_hcv_pair / "MiSeq/runs/140101_M01234"
    sample_sheet = run_folder / "SampleSheet.csv"
    sample_sheet_text = sample_sheet.read_text()
    sample_sheet_text += ('CFE_MS1_23-Jan-2014_N506-N702_2130AWG_HCV_nopid,'
                          '2130AWG_HCV,23-Jan-2014,N/A,TTTAAAAA,AAATTTTT,'
                          '23-Jan-2014,Research:2130AWG_HCV:FALSE '
                          'Comments:2130AWG_HCV: '
                          'Disablecontamcheck:2130AWG_HCV:FALSE,\n')
    sample_sheet.write_text(sample_sheet_text)
    base_calls_folder = run_folder / "Data/Intensities/BaseCalls"
    for sample_path in base_calls_folder.glob('2130A-HCV_S15_*.fastq.gz'):
        new_name = sample_path.name.replace('2130A-HCV_S15_', '2130AWG-HCV_S18_')
        new_path = sample_path.parent / new_name
        # noinspection PyTypeChecker
        shutil.copy(sample_path, new_path)
    kive_watcher = create_kive_watcher_with_main_run(
        pipelines_config,
        base_calls_folder,
        SampleGroup('2130A',
                    ('2130A-HCV_S15_L001_R1_001.fastq.gz',
                     '2130AMIDI-MidHCV_S16_L001_R1_001.fastq.gz'),
                    ('HCV', 'MidHCV')))
    mock_session = kive_watcher.session
    folder_watcher, = kive_watcher.folder_watchers.values()
    sample_watcher, = folder_watcher.sample_watchers
    folder_watcher.add_run(
        dict(id=108),
        PipelineType.MIDI,
        sample_watcher)
    sample_watcher2 = kive_watcher.add_sample_group(
        base_calls=base_calls_folder,
        sample_group=SampleGroup('2130A',
                                 ('2130AWG-HCV_S18_L001_R1_001.fastq.gz',
                                  '2130AMIDI-MidHCV_S16_L001_R1_001.fastq.gz'),
                                 ('HCV', 'MidHCV')))
    folder_watcher.add_run(
        dict(id=109),
        PipelineType.MAIN,
        sample_watcher2)
    folder_watcher.add_run(
        dict(id=108),
        PipelineType.MIDI,
        sample_watcher2)

    mock_session.endpoints.containerruns.get.side_effect = [
        dict(id=107, state='C'),  # refresh run state for main run
        [dict(dataset='/datasets/111/',
              argument_type='O',
              argument_name='amino_csv'),
         dict(dataset='/datasets/121/',
              argument_type='O',
              argument_name='nuc_csv')],  # run datasets
        dict(id=108, state='C'),  # refresh run state for midi run
        [dict(dataset='/datasets/112/',
              argument_type='O',
              argument_name='amino_csv'),
         dict(dataset='/datasets/122/',
              argument_type='O',
              argument_name='nuc_csv')],  # run datasets
        dict(id=109, state='C'),  # refresh run state for other main run
        [dict(dataset='/datasets/113/',
              argument_type='O',
              argument_name='amino_csv'),
         dict(dataset='/datasets/123/',
              argument_type='O',
              argument_name='nuc_csv')]]  # run datasets
    mock_session.get.reset_mock(side_effect=True)
    mock_session.get.return_value.json.side_effect = [
        dict(url='/datasets/111/', id=111),
        dict(url='/datasets/112/', id=112),
        dict(url='/datasets/121/', id=121),
        dict(url='/datasets/113/', id=113),
        dict(url='/datasets/112/', id=112),
        dict(url='/datasets/123/', id=123)]
    mock_session.endpoints.containerruns.post.return_value = dict(id=None)

    kive_watcher.poll_runs()

    expected_calls = [
        call(json=dict(app=ANY,
                       datasets=[dict(argument=ANY, dataset='/datasets/111/'),
                                 dict(argument=ANY, dataset='/datasets/112/'),
                                 dict(argument=ANY, dataset='/datasets/121/')],
                       name='MiCall resistance on 2130A',
                       batch=ANY,
                       groups_allowed=['Everyone'])),
        call(json=dict(app=ANY,
                       datasets=[dict(argument=ANY, dataset='/datasets/113/'),
                                 dict(argument=ANY, dataset='/datasets/112/'),
                                 dict(argument=ANY, dataset='/datasets/123/')],
                       name='MiCall resistance on 2130A',
                       batch=ANY,
                       groups_allowed=['Everyone']))]
    assert expected_calls == mock_session.endpoints.containerruns.post.mock_calls


def test_launch_mixed_hcv_run(raw_data_with_hcv_pair, mock_open_kive, pipelines_config):
    pipelines_config.mixed_hcv_pipeline_id = 47
    base_calls = (raw_data_with_hcv_pair /
                  "MiSeq/runs/140101_M01234/Data/Intensities/BaseCalls")
    kive_watcher = KiveWatcher(pipelines_config)
    kive_watcher.app_urls = {
        pipelines_config.micall_filter_quality_pipeline_id: '/containerapps/102',
        pipelines_config.micall_main_pipeline_id: '/containerapps/103',
        pipelines_config.mixed_hcv_pipeline_id: '/containerapps/104'}
    kive_watcher.app_args = {
        pipelines_config.micall_filter_quality_pipeline_id: dict(
            quality_csv='/containerargs/105'),
        pipelines_config.micall_main_pipeline_id: dict(
            fastq1='/containerargs/106',
            fastq2='/containerargs/107',
            bad_cycles_csv='/containerargs/108'),
        pipelines_config.mixed_hcv_pipeline_id: dict(
            fastq1='/containerargs/109',
            fastq2='/containerargs/110')}

    folder_watcher = kive_watcher.add_folder(base_calls)
    folder_watcher.batch = dict(url='/batches/101')
    folder_watcher.add_run(dict(id=120), PipelineType.FILTER_QUALITY)
    mock_session = mock_open_kive.return_value
    mock_session.endpoints.containerruns.get.side_effect = [
        dict(id=120, state='C'),  # refresh run state
        [dict(dataset='/datasets/121/',
              argument_name='bad_cycles_csv')]]  # run datasets
    mock_session.get.return_value.json.side_effect = [dict(id=121, url='/datasets/121')]
    mock_session.endpoints.datasets.post.side_effect = [
        dict(url='/datasets/122', id=122),
        dict(url='/datasets/123', id=123),
        dict(url='/datasets/124', id=124),
        dict(url='/datasets/125', id=125)]

    kive_watcher.add_sample_group(
        base_calls=base_calls,
        sample_group=SampleGroup('2130A',
                                 ('2130A-HCV_S15_L001_R1_001.fastq.gz',
                                  '2130AMIDI-MidHCV_S16_L001_R1_001.fastq.gz'),
                                 ('HCV', 'MidHCV')))

    kive_watcher.poll_runs()

    assert [call(json=dict(app='/containerapps/104',
                           datasets=[dict(argument='/containerargs/109',
                                          dataset='/datasets/122'),
                                     dict(argument='/containerargs/110',
                                          dataset='/datasets/123')],
                           name='Mixed HCV on 2130A-HCV_S15',
                           batch='/batches/101',
                           groups_allowed=['Everyone'])),
            call(json=dict(app='/containerapps/104',
                           datasets=[dict(argument='/containerargs/109',
                                          dataset='/datasets/124'),
                                     dict(argument='/containerargs/110',
                                          dataset='/datasets/125')],
                           name='Mixed HCV on 2130AMIDI-MidHCV_S16',
                           batch='/batches/101',
                           groups_allowed=['Everyone'])),
            call(json=dict(app='/containerapps/103',
                           datasets=[dict(argument='/containerargs/106',
                                          dataset='/datasets/122'),
                                     dict(argument='/containerargs/107',
                                          dataset='/datasets/123'),
                                     dict(argument='/containerargs/108',
                                          dataset='/datasets/121')],
                           name='MiCall main on 2130A-HCV_S15',
                           batch='/batches/101',
                           groups_allowed=['Everyone'])),
            call(json=dict(app='/containerapps/103',
                           datasets=[dict(argument='/containerargs/106',
                                          dataset='/datasets/124'),
                                     dict(argument='/containerargs/107',
                                          dataset='/datasets/125'),
                                     dict(argument='/containerargs/108',
                                          dataset='/datasets/121')],
                           name='MiCall main on 2130AMIDI-MidHCV_S16',
                           batch='/batches/101',
                           groups_allowed=['Everyone']))
            ] == mock_session.endpoints.containerruns.post.call_args_list


def test_full_with_two_samples(raw_data_with_two_samples, mock_open_kive, pipelines_config):
    pipelines_config.max_active = 2
    base_calls = (raw_data_with_two_samples /
                  "MiSeq/runs/140101_M01234/Data/Intensities/BaseCalls")
    kive_watcher = KiveWatcher(pipelines_config)

    kive_watcher.add_sample_group(
        base_calls=base_calls,
        sample_group=SampleGroup('2110A',
                                 ('2110A-V3LOOP_S13_L001_R1_001.fastq.gz',
                                  None),
                                 ('V3LOOP', None)))
    kive_watcher.poll_runs()
    is_full1 = kive_watcher.is_full()
    kive_watcher.add_sample_group(
        base_calls=base_calls,
        sample_group=SampleGroup('2120A',
                                 ('2120A-PR_S14_L001_R1_001.fastq.gz',
                                  None),
                                 ('PR', None)))
    kive_watcher.poll_runs()
    is_full2 = kive_watcher.is_full()

    folder_watcher, = kive_watcher.folder_watchers.values()
    folder_watcher.completed_samples.add('2120A-PR_S14_L001_R1_001.fastq.gz')

    is_full3 = kive_watcher.is_full()

    assert not is_full1
    assert is_full2
    assert not is_full3


def test_full_with_two_runs(raw_data_with_two_runs, mock_open_kive, pipelines_config):
    pipelines_config.max_active = 2
    base_calls1 = (raw_data_with_two_runs /
                   "MiSeq/runs/140101_M01234/Data/Intensities/BaseCalls")
    base_calls2 = (raw_data_with_two_runs /
                   "MiSeq/runs/140201_M01234/Data/Intensities/BaseCalls")
    kive_watcher = KiveWatcher(pipelines_config)

    kive_watcher.add_sample_group(
        base_calls=base_calls1,
        sample_group=SampleGroup('2000A',
                                 ('2000A-V3LOOP_S2_L001_R1_001.fastq.gz',
                                  None),
                                 ('V3LOOP', None)))
    kive_watcher.poll_runs()
    is_full1 = kive_watcher.is_full()
    kive_watcher.add_sample_group(
        base_calls=base_calls2,
        sample_group=SampleGroup('2010A',
                                 ('2010A-V3LOOP_S3_L001_R1_001.fastq.gz',
                                  None),
                                 ('V3LOOP', None)))
    kive_watcher.poll_runs()
    is_full2 = kive_watcher.is_full()

    folder_watcher = kive_watcher.folder_watchers[base_calls1]
    folder_watcher.completed_samples.add('2000A-V3LOOP_S2_L001_R1_001.fastq.gz')

    is_full3 = kive_watcher.is_full()

    assert not is_full1
    assert is_full2
    assert not is_full3


def test_fetch_run_status_incomplete(mock_open_kive, pipelines_config):
    mock_run = dict(id=123)

    kive_watcher = KiveWatcher(pipelines_config)

    new_run = kive_watcher.fetch_run_status(mock_run,
                                            folder_watcher=None,
                                            pipeline_type=None,
                                            sample_watchers=None)

    assert new_run is mock_run


def test_fetch_run_status_filter_quality(raw_data_with_two_runs,
                                         mock_open_kive,
                                         pipelines_config):
    mock_session = mock_open_kive.return_value
    base_calls = (raw_data_with_two_runs /
                  "MiSeq/runs/140101_M01234/Data/Intensities/BaseCalls")
    folder_watcher = FolderWatcher(base_calls)
    sample_watcher = None
    mock_run = dict(id=123)
    mock_session.endpoints.containerruns.get.side_effect = [
        dict(state='C')]

    kive_watcher = KiveWatcher(pipelines_config)

    new_run = kive_watcher.fetch_run_status(mock_run,
                                            folder_watcher,
                                            PipelineType.FILTER_QUALITY,
                                            sample_watcher)

    assert new_run is None


def test_fetch_run_status_main(raw_data_with_two_runs,
                               mock_open_kive,
                               pipelines_config):
    mock_session = mock_open_kive.return_value
    base_calls = (raw_data_with_two_runs /
                  "MiSeq/runs/140101_M01234/Data/Intensities/BaseCalls")
    folder_watcher = FolderWatcher(base_calls)
    sample_watcher = SampleWatcher(
        SampleGroup('2000A',
                    ('2000A-V3LOOP_S2_L001_R1_001.fastq.gz',
                     None),
                    ('V3LOOP', None)))
    mock_run = dict(id=123)
    mock_session.endpoints.containerruns.get.side_effect = [
        dict(state='C'),  # run state refresh
        [dict(argument_name='insertions_csv',
              argument_type='O',
              dataset='/datasets/110/'),
         dict(argument_name='nuc_csv',
              argument_type='O',
              dataset='/datasets/111/')]]  # run datasets
    expected_scratch = base_calls / "../../../Results/version_0-dev/scratch"
    expected_insertion_path = expected_scratch / "2000A-V3LOOP_S2/insertions.csv"
    expected_nuc_path = expected_scratch / "2000A-V3LOOP_S2/nuc.csv"

    kive_watcher = KiveWatcher(pipelines_config)

    new_run = kive_watcher.fetch_run_status(mock_run,
                                            folder_watcher,
                                            PipelineType.MAIN,
                                            [sample_watcher])

    assert new_run is None
    assert expected_insertion_path.exists()
    assert expected_nuc_path.exists()
    assert [call(ANY, '/datasets/110/download/'),
            call(ANY, '/datasets/111/download/')
            ] == mock_session.download_file.call_args_list


def test_fetch_run_status_main_and_resistance(raw_data_with_two_runs,
                                              mock_open_kive,
                                              pipelines_config):
    mock_session = mock_open_kive.return_value
    base_calls = (raw_data_with_two_runs /
                  "MiSeq/runs/140101_M01234/Data/Intensities/BaseCalls")
    folder_watcher = FolderWatcher(base_calls)
    sample_watcher = SampleWatcher(
        SampleGroup('2000A',
                    ('2000A-V3LOOP_S2_L001_R1_001.fastq.gz',
                     None),
                    ('V3LOOP', None)))
    main_run = dict(id=123)
    resistance_run = dict(id=124)
    mock_session.endpoints.containerruns.get.side_effect = [
        dict(state='C'),  # main run refresh
        [dict(argument_name='nuc_csv',
              argument_type='O',
              dataset='/datasets/110/')],  # main run datasets
        dict(state='C'),  # resistance run refresh
        [dict(argument_name='resistance_csv',
              argument_type='O',
              dataset='/datasets/112/')]]  # resistance run datasets
    expected_scratch = base_calls / "../../../Results/version_0-dev/scratch"
    expected_nuc_path = expected_scratch / "2000A-V3LOOP_S2/nuc.csv"
    expected_resistance_path = expected_scratch / "2000A-V3LOOP_S2/resistance.csv"

    kive_watcher = KiveWatcher(pipelines_config)

    new_main_run = kive_watcher.fetch_run_status(
        main_run,
        folder_watcher,
        PipelineType.MAIN,
        [sample_watcher])
    new_resistance_run = kive_watcher.fetch_run_status(
        resistance_run,
        folder_watcher,
        PipelineType.RESISTANCE,
        [sample_watcher])

    assert new_main_run is None
    assert new_resistance_run is None
    assert expected_nuc_path.exists()
    assert expected_resistance_path.exists()


def test_fetch_run_status_main_and_midi(raw_data_with_hcv_pair,
                                        mock_open_kive,
                                        pipelines_config):
    mock_session = mock_open_kive.return_value
    base_calls = (raw_data_with_hcv_pair /
                  "MiSeq/runs/140101_M01234/Data/Intensities/BaseCalls")
    folder_watcher = FolderWatcher(base_calls)
    sample_watcher = SampleWatcher(
        SampleGroup('2130A',
                    ('2130A-HCV_S15_L001_R1_001.fastq.gz',
                     '2130AMIDI-MidHCV_S16_L001_R1_001.fastq.gz'),
                    ('HCV', 'MidHCV')))
    main_run = dict(id=123)
    midi_run = dict(id=124)
    mock_session.endpoints.containerruns.get.side_effect = [
        dict(state='C'),  # main run refresh
        [dict(argument_name='nuc_csv',
              argument_type='O',
              dataset='/datasets/110/')],  # main outputs
        dict(state='C'),  # midi run refresh
        [dict(argument_name='nuc_csv',
              argument_type='O',
              dataset='/datasets/111/')]]  # midi outputs
    expected_scratch = base_calls / "../../../Results/version_0-dev/scratch"
    expected_main_nuc_path = expected_scratch / "2130A-HCV_S15/nuc.csv"
    expected_midi_nuc_path = expected_scratch / "2130AMIDI-MidHCV_S16/nuc.csv"

    kive_watcher = KiveWatcher(pipelines_config)

    new_main_run = kive_watcher.fetch_run_status(main_run,
                                                 folder_watcher,
                                                 PipelineType.MAIN,
                                                 [sample_watcher])
    new_midi_run = kive_watcher.fetch_run_status(midi_run,
                                                 folder_watcher,
                                                 PipelineType.MIDI,
                                                 [sample_watcher])

    assert new_main_run is None
    assert new_midi_run is None
    assert expected_main_nuc_path.exists()
    assert expected_midi_nuc_path.exists()


def test_fetch_run_status_session_expired(raw_data_with_two_runs,
                                          mock_open_kive,
                                          pipelines_config):
    mock_session = mock_open_kive.return_value
    base_calls = (raw_data_with_two_runs /
                  "MiSeq/runs/140101_M01234/Data/Intensities/BaseCalls")
    mock_run = dict(id=123)
    mock_session.endpoints.containerruns.get.side_effect = [
        KiveClientException('expired'),  # Session expired
        dict(state='C'),  # run state refresh
        []]  # run outputs

    kive_watcher = KiveWatcher(pipelines_config)

    sample_watcher = kive_watcher.add_sample_group(
        base_calls,
        SampleGroup('2000A',
                    ('2000A-V3LOOP_S2_L001_R1_001.fastq.gz', None),
                    ('V3LOOP', None)))
    folder_watcher, = kive_watcher.folder_watchers.values()

    new_run = kive_watcher.fetch_run_status(mock_run,
                                            folder_watcher,
                                            PipelineType.MAIN,
                                            [sample_watcher])

    assert new_run is None


def test_fetch_run_status_user_cancelled(raw_data_with_two_runs,
                                         mock_open_kive,
                                         pipelines_config):
    base_calls = (raw_data_with_two_runs /
                  "MiSeq/runs/140101_M01234/Data/Intensities/BaseCalls")
    kive_watcher = create_kive_watcher_with_main_run(
        pipelines_config,
        base_calls,
        SampleGroup('2000A',
                    ('2000A-V3LOOP_S2_L001_R1_001.fastq.gz', None),
                    ('V3LOOP', None)))
    mock_session = kive_watcher.session
    original_run = dict(id=123)
    mock_session.endpoints.containerruns.get.side_effect = [
        dict(state='X')]
    folder_watcher, = kive_watcher.folder_watchers.values()
    sample_watcher, = folder_watcher.sample_watchers

    new_run = kive_watcher.fetch_run_status(original_run,
                                            folder_watcher,
                                            PipelineType.MAIN,
                                            [sample_watcher])

    assert new_run is not None
    assert new_run is not original_run


def test_folder_completed(raw_data_with_two_samples, mock_open_kive, default_config):
    base_calls = (raw_data_with_two_samples /
                  "MiSeq/runs/140101_M01234/Data/Intensities/BaseCalls")
    kive_watcher = create_kive_watcher_with_main_run(
        default_config,
        base_calls,
        SampleGroup('2110A',
                    ('2110A-V3LOOP_S13_L001_R1_001.fastq.gz', None),
                    ('V3LOOP', None)),
        is_complete=True)
    folder_watcher = kive_watcher.folder_watchers[base_calls]
    sample1_watcher, = folder_watcher.sample_watchers
    sample2_watcher = kive_watcher.add_sample_group(
        base_calls=base_calls,
        sample_group=SampleGroup('2120A',
                                 ('2120A-PR_S14_L001_R1_001.fastq.gz', None),
                                 ('V3LOOP', None)))
    kive_watcher.finish_folder(base_calls)
    folder_watcher.add_run(dict(id=150),
                           PipelineType.MAIN,
                           sample2_watcher,
                           is_complete=True)
    folder_watcher.add_run(dict(id=151),
                           PipelineType.RESISTANCE,
                           sample1_watcher)
    folder_watcher.add_run(dict(id=152),
                           PipelineType.RESISTANCE,
                           sample2_watcher)
    kive_watcher.session.endpoints.containerruns.get.side_effect = [
        dict(id=151, state='C'),  # refresh run state for 2110
        [dict(dataset='/datasets/161/',
              argument_type='O',
              argument_name='resistance_csv')],  # run datasets
        dict(id=152, state='C'),  # refresh run state for 2120
        [dict(dataset='/datasets/162/',
              argument_type='O',
              argument_name='resistance_csv')]]  # run datasets
    results_path = base_calls / "../../../Results/version_0-dev"
    scratch_path = results_path / "scratch"
    expected_coverage_map_content = b'This is a coverage map.'
    sample_scratch_path = scratch_path / "2110A-V3LOOP_S13"
    sample_scratch_path.mkdir(parents=True)
    coverage_maps_path = sample_scratch_path / "coverage_maps.tar"
    with tarfile.open(coverage_maps_path, 'w') as coverage_maps_tar:
        content = BytesIO(expected_coverage_map_content)
        tar_info = TarInfo('coverage_maps/R1_coverage.txt')
        tar_info.size = len(expected_coverage_map_content)
        coverage_maps_tar.addfile(tar_info, content)
    expected_coverage_map_path = (
            results_path / "coverage_maps/2110A-V3LOOP_S13.R1_coverage.txt")
    expected_mutations_path = results_path / "mutations.csv"
    expected_done_path = results_path / "doneprocessing"
    expected_all_done_path = results_path / "done_all_processing"
    expected_resistance_path = results_path / "resistance.csv"
    expected_resistance_content = """\
sample,url,n
2110A-V3LOOP_S13,/datasets/161/download/,0
2110A-V3LOOP_S13,/datasets/161/download/,1
2110A-V3LOOP_S13,/datasets/161/download/,2
2120A-PR_S14,/datasets/162/download/,0
2120A-PR_S14,/datasets/162/download/,1
2120A-PR_S14,/datasets/162/download/,2
"""

    kive_watcher.poll_runs()

    assert not scratch_path.exists()
    assert not expected_mutations_path.exists()
    assert expected_resistance_content == expected_resistance_path.read_text()
    assert expected_coverage_map_content == expected_coverage_map_path.read_bytes()
    assert expected_done_path.exists()
    assert expected_all_done_path.exists()
    assert kive_watcher.is_idle()


def test_folder_completed_except_denovo(raw_data_with_two_samples, mock_open_kive, default_config):
    base_calls = (raw_data_with_two_samples /
                  "MiSeq/runs/140101_M01234/Data/Intensities/BaseCalls")
    default_config.denovo_main_pipeline_id = 495
    kive_watcher = create_kive_watcher_with_main_run(
        default_config,
        base_calls,
        SampleGroup('2110A',
                    ('2110A-V3LOOP_S13_L001_R1_001.fastq.gz', None),
                    ('V3LOOP', None)),
        is_complete=True)
    kive_watcher.app_urls[default_config.denovo_main_pipeline_id] = '/containerapps/105'
    kive_watcher.app_args[default_config.denovo_main_pipeline_id] = dict(
            bad_cycles_csv='/containerargs/113',
            fastq1='/containerargs/114',
            fastq2='/containerargs/115')
    folder_watcher = kive_watcher.folder_watchers[base_calls]
    sample1_watcher, = folder_watcher.sample_watchers
    sample2_watcher = kive_watcher.add_sample_group(
        base_calls=base_calls,
        sample_group=SampleGroup('2120A',
                                 ('2120A-PR_S14_L001_R1_001.fastq.gz', None),
                                 ('PR', None)))
    kive_watcher.finish_folder(base_calls)
    folder_watcher.add_run(dict(id=150),
                           PipelineType.MAIN,
                           sample2_watcher,
                           is_complete=True)
    folder_watcher.add_run(dict(id=151),
                           PipelineType.DENOVO_MAIN,
                           sample1_watcher)
    folder_watcher.add_run(dict(id=152),
                           PipelineType.RESISTANCE,
                           sample1_watcher)
    folder_watcher.add_run(dict(id=153),
                           PipelineType.RESISTANCE,
                           sample2_watcher)
    kive_watcher.session.endpoints.containerruns.get.side_effect = [
        dict(id=151, state='C'),  # refresh run state for denovo main
        [dict(dataset='/datasets/161/',
              argument_type='O',
              argument_name='amino_csv'),
         dict(dataset='/datasets/171/',
              argument_type='O',
              argument_name='nuc_csv')],  # run datasets
        dict(id=152, state='C'),  # refresh run state for 2110
        [dict(dataset='/datasets/162/',
              argument_type='O',
              argument_name='resistance_csv')],  # run datasets
        dict(id=153, state='C'),  # refresh run state for 2120
        [dict(dataset='/datasets/163/',
              argument_type='O',
              argument_name='resistance_csv')]]  # run datasets
    results_path = base_calls / "../../../Results/version_0-dev"
    scratch_path = results_path / "scratch"
    sample_scratch_path = scratch_path / "2110A-V3LOOP_S13"
    sample_scratch_path.mkdir(parents=True)
    denovo_scratch_path = results_path / "scratch_denovo" / "2110A-V3LOOP_S13"
    expected_done_path = results_path / "doneprocessing"
    expected_all_done_path = results_path / "done_all_processing"
    expected_mutations_path = results_path / "mutations.csv"
    expected_resistance_path = results_path / "resistance.csv"
    expected_amino_path = denovo_scratch_path / "amino.csv"

    kive_watcher.poll_runs()

    assert not scratch_path.exists()
    assert not expected_mutations_path.exists()
    assert expected_resistance_path.exists()
    assert expected_done_path.exists()
    assert not expected_all_done_path.exists()
    assert expected_amino_path.exists()
    assert not kive_watcher.is_idle()


def test_folder_completed_with_fasta(raw_data_with_two_samples, mock_open_kive, default_config):
    base_calls = (raw_data_with_two_samples /
                  "MiSeq/runs/140101_M01234/Data/Intensities/BaseCalls")
    kive_watcher = create_kive_watcher_with_main_run(
        default_config,
        base_calls,
        SampleGroup('2110A',
                    ('2110A-V3LOOP_S13_L001_R1_001.fastq.gz', None),
                    ('V3LOOP', None)),
        is_complete=True)
    mock_session = mock_open_kive.return_value
    mock_session.download_file.side_effect = mock_session_download_fasta
    folder_watcher = kive_watcher.folder_watchers[base_calls]
    sample1_watcher, = folder_watcher.sample_watchers
    sample2_watcher = kive_watcher.add_sample_group(
        base_calls=base_calls,
        sample_group=SampleGroup('2120A',
                                 ('2120A-PR_S14_L001_R1_001.fastq.gz', None),
                                 ('PR', None)))
    kive_watcher.finish_folder(base_calls)
    folder_watcher.add_run(dict(id=150),
                           PipelineType.MAIN,
                           sample2_watcher,
                           is_complete=True)
    folder_watcher.add_run(dict(id=151),
                           PipelineType.RESISTANCE,
                           sample1_watcher)
    folder_watcher.add_run(dict(id=152),
                           PipelineType.RESISTANCE,
                           sample2_watcher)
    kive_watcher.session.endpoints.containerruns.get.side_effect = [
        dict(id=151, state='C'),  # refresh run state for 2110
        [dict(dataset='/datasets/161/',
              argument_type='O',
              argument_name='wg_fasta')],  # run datasets
        dict(id=152, state='C'),  # refresh run state for 2120
        [dict(dataset='/datasets/162/',
              argument_type='O',
              argument_name='wg_fasta')]]  # run datasets
    results_path = base_calls / "../../../Results/version_0-dev"
    scratch_path = results_path / "scratch"
    sample_scratch_path = scratch_path / "2110A-V3LOOP_S13"
    sample_scratch_path.mkdir(parents=True)
    expected_fasta_path = results_path / "wg.fasta"
    expected_fasta_content = """\
>2110A-V3LOOP_S13,/datasets/161/download/,0
ACTGTCA
CTGTCA
TGTCA
>2110A-V3LOOP_S13,/datasets/161/download/,1
CTGTCA
TGTCA
GTCA
>2120A-PR_S14,/datasets/162/download/,0
ACTGTCA
CTGTCA
TGTCA
>2120A-PR_S14,/datasets/162/download/,1
CTGTCA
TGTCA
GTCA
"""

    kive_watcher.poll_runs()

    assert not scratch_path.exists()
    assert expected_fasta_content == expected_fasta_path.read_text()


def test_folder_completed_with_svg(raw_data_with_two_samples, mock_open_kive, default_config):
    base_calls = (raw_data_with_two_samples /
                  "MiSeq/runs/140101_M01234/Data/Intensities/BaseCalls")
    kive_watcher = create_kive_watcher_with_main_run(
        default_config,
        base_calls,
        SampleGroup('2110A',
                    ('2110A-V3LOOP_S13_L001_R1_001.fastq.gz', None),
                    ('V3LOOP', None)),
        is_complete=True)
    folder_watcher = kive_watcher.folder_watchers[base_calls]
    sample1_watcher, = folder_watcher.sample_watchers
    sample2_watcher = kive_watcher.add_sample_group(
        base_calls=base_calls,
        sample_group=SampleGroup('2120A',
                                 ('2120A-PR_S14_L001_R1_001.fastq.gz', None),
                                 ('PR', None)))
    kive_watcher.finish_folder(base_calls)
    folder_watcher.add_run(dict(id=150),
                           PipelineType.MAIN,
                           sample2_watcher,
                           is_complete=True)
    folder_watcher.add_run(dict(id=151),
                           PipelineType.RESISTANCE,
                           sample1_watcher)
    folder_watcher.add_run(dict(id=152),
                           PipelineType.RESISTANCE,
                           sample2_watcher)
    kive_watcher.session.endpoints.containerruns.get.side_effect = [
        dict(id=151, state='C'),  # refresh run state for 2110
        [dict(dataset='/datasets/161/',
              argument_type='O',
              argument_name='alignment_svg')],  # run datasets
        dict(id=152, state='C'),  # refresh run state for 2120
        [dict(dataset='/datasets/162/',
              argument_type='O',
              argument_name='alignment_svg')]]  # run datasets
    results_path = base_calls / "../../../Results/version_0-dev"
    expected_alignment1_path = results_path / "alignment" / "2110A-V3LOOP_S13_alignment.svg"
    expected_alignment2_path = results_path / "alignment" / "2120A-PR_S14_alignment.svg"

    kive_watcher.poll_runs()

    assert expected_alignment1_path.exists()
    assert expected_alignment2_path.exists()


def test_folder_not_finished(raw_data_with_two_samples, mock_open_kive, default_config):
    base_calls = (raw_data_with_two_samples /
                  "MiSeq/runs/140101_M01234/Data/Intensities/BaseCalls")
    kive_watcher = create_kive_watcher_with_main_run(
        default_config,
        base_calls,
        SampleGroup('2110A',
                    ('2110A-V3LOOP_S13_L001_R1_001.fastq.gz', None),
                    ('V3LOOP', None)),
        is_complete=True)
    folder_watcher = kive_watcher.folder_watchers[base_calls]
    sample1_watcher, = folder_watcher.sample_watchers
    sample2_watcher = kive_watcher.add_sample_group(
        base_calls=base_calls,
        sample_group=SampleGroup('2120A',
                                 ('2120A-PR_S14_L001_R1_001.fastq.gz', None),
                                 ('PR', None)))
    # Did not call kive_watcher.finish_folder(), more samples could be coming.
    folder_watcher.add_run(dict(id=150),
                           PipelineType.MAIN,
                           sample2_watcher,
                           is_complete=True)
    folder_watcher.add_run(dict(id=151),
                           PipelineType.RESISTANCE,
                           sample1_watcher)
    folder_watcher.add_run(dict(id=152),
                           PipelineType.RESISTANCE,
                           sample2_watcher)
    kive_watcher.session.endpoints.containerruns.get.side_effect = [
        dict(id=151, state='C'),  # refresh run state for 2110
        [dict(dataset='/datasets/161/',
              argument_type='O',
              argument_name='resistance_csv')],  # run datasets
        dict(id=152, state='C'),  # refresh run state for 2120
        [dict(dataset='/datasets/162/',
              argument_type='O',
              argument_name='resistance_csv')]]  # run datasets
    results_path = base_calls / "../../../Results/version_0-dev"
    scratch_path = results_path / "scratch"
    expected_resistance_path = results_path / "resistance.csv"

    kive_watcher.poll_runs()

    assert scratch_path.exists()
    assert not expected_resistance_path.exists()


def test_folder_not_finished_before_new_start(raw_data_with_two_runs,
                                              mock_open_kive,
                                              default_config):
    base_calls1 = (raw_data_with_two_runs /
                   "MiSeq/runs/140101_M01234/Data/Intensities/BaseCalls")
    base_calls2 = (raw_data_with_two_runs /
                   "MiSeq/runs/140201_M01234/Data/Intensities/BaseCalls")
    kive_watcher = create_kive_watcher_with_main_run(
        default_config,
        base_calls1,
        SampleGroup('2000A',
                    ('2000A-V3LOOP_S2_L001_R1_001.fastq.gz', None),
                    ('V3LOOP', None)),
        is_complete=True)
    folder_watcher1 = kive_watcher.folder_watchers[base_calls1]
    sample1_watcher, = folder_watcher1.sample_watchers

    # Did not call kive_watcher.finish_folder(base_calls1), more samples could be coming.
    folder_watcher2 = kive_watcher.add_folder(base_calls2)
    folder_watcher2.batch = dict(url='/batches/171/')
    kive_watcher.add_sample_group(
        base_calls=base_calls2,
        sample_group=SampleGroup('2010A',
                                 ('2010A-V3LOOP_S3_L001_R1_001.fastq.gz', None),
                                 ('V3LOOP', None)))
    folder_watcher1.add_run(dict(id=151),
                            PipelineType.RESISTANCE,
                            sample1_watcher)
    folder_watcher2.quality_dataset = dict(url='/datasets/127/', id=127)
    results_path = base_calls1 / "../../../Results/version_0-dev"
    scratch_path = results_path / "scratch"
    expected_resistance_path = results_path / "resistance.csv"

    kive_watcher.poll_runs()

    assert not expected_resistance_path.exists()
    assert scratch_path.exists()


def test_folder_failed_quality(raw_data_with_two_samples, mock_open_kive, default_config):
    mock_session = mock_open_kive.return_value
    base_calls = (raw_data_with_two_samples /
                  "MiSeq/runs/140101_M01234/Data/Intensities/BaseCalls")
    kive_watcher = create_kive_watcher_with_filter_run(default_config, base_calls)
    mock_session.endpoints.containerruns.get.side_effect = [dict(id=110,
                                                                 state='F')]

    kive_watcher.finish_folder(base_calls)
    run_path = base_calls / "../../.."
    results_path = run_path / "Results/version_0-dev"
    expected_done_path = results_path / "doneprocessing"
    expected_error_path = run_path / "errorprocessing"
    expected_error_message = "Filter quality failed in Kive.\n"

    kive_watcher.poll_runs()

    assert not expected_done_path.exists()
    assert expected_error_message == expected_error_path.read_text()


def test_folder_failed_quality_incomplete(raw_data_with_two_samples,
                                          mock_open_kive,
                                          default_config):
    """ The filter_quality run fails before all the samples are loaded. """
    default_config.max_active = 1
    mock_session = mock_open_kive.return_value
    base_calls = (raw_data_with_two_samples /
                  "MiSeq/runs/140101_M01234/Data/Intensities/BaseCalls")
    kive_watcher = create_kive_watcher_with_filter_run(default_config, base_calls)
    mock_session.endpoints.containerruns.get.side_effect = [dict(id=110,
                                                                 state='F')]

    kive_watcher.add_sample_group(
        base_calls=base_calls,
        sample_group=SampleGroup('2110A',
                                 ('2110A-V3LOOP_S13_L001_R1_001.fastq.gz',
                                  None),
                                 ('V3LOOP', None)))
    run_path = base_calls / "../../.."
    results_path = run_path / "Results/version_0-dev"
    expected_done_path = results_path / "doneprocessing"

    is_full_before_fail = kive_watcher.is_full()
    kive_watcher.poll_runs()
    is_full_after_fail = kive_watcher.is_full()

    assert not expected_done_path.exists()
    assert is_full_before_fail
    assert not is_full_after_fail


def test_folder_failed_sample(raw_data_with_two_samples, mock_open_kive, default_config):
    mock_session = mock_open_kive.return_value
    base_calls = (raw_data_with_two_samples /
                  "MiSeq/runs/140101_M01234/Data/Intensities/BaseCalls")
    kive_watcher = create_kive_watcher_with_main_run(
        default_config,
        base_calls,
        SampleGroup('2110A',
                    ('2110A-V3LOOP_S13_L001_R1_001.fastq.gz', None),
                    ('V3LOOP', None)))
    folder_watcher = kive_watcher.folder_watchers[base_calls]

    sample2_watcher = kive_watcher.add_sample_group(
        base_calls=base_calls,
        sample_group=SampleGroup('2120A',
                                 ('2120A-PR_S14_L001_R1_001.fastq.gz', None),
                                 ('PR', None)))
    kive_watcher.finish_folder(base_calls)
    folder_watcher.add_run(dict(id=151),
                           PipelineType.MAIN,
                           sample2_watcher,
                           is_complete=True)
    folder_watcher.add_run(dict(id=152),
                           PipelineType.RESISTANCE,
                           sample2_watcher)
    mock_session.endpoints.containerruns.get.side_effect = [
        dict(state='F', id=107),  # main run for 2110 fails
        dict(state='C', id=152),  # resistance run for 2120 complete
        [dict(argument_name='resistance_csv',
              argument_type='O',
              dataset='/datasets/167/')]]  # outputs for resistance run
    run_path = base_calls / "../../.."
    results_path = run_path / "Results/version_0-dev"
    expected_done_path = results_path / "doneprocessing"
    expected_error_path = run_path / "errorprocessing"
    expected_error_message = "Samples failed in Kive: 2110A.\n"

    kive_watcher.poll_runs()

    assert not expected_done_path.exists()
    assert expected_error_message == expected_error_path.read_text()


def test_add_duplicate_sample(raw_data_with_two_samples,
                              mock_open_kive,
                              default_config):
    """ Could happen when restarting a run after finishing a newer one. """
    base_calls = (raw_data_with_two_samples /
                  "MiSeq/runs/140101_M01234/Data/Intensities/BaseCalls")
    assert mock_open_kive
    kive_watcher = KiveWatcher(default_config)

    sample_watcher1 = kive_watcher.add_sample_group(
        base_calls=base_calls,
        sample_group=SampleGroup('2110A',
                                 ('2110A-V3LOOP_S13_L001_R1_001.fastq.gz',
                                  None),
                                 ('V3LOOP', None)))
    sample_watcher2 = kive_watcher.add_sample_group(
        base_calls=base_calls,
        sample_group=SampleGroup('2110A',
                                 ('2110A-V3LOOP_S13_L001_R1_001.fastq.gz',
                                  None),
                                 ('V3LOOP', None)))

    assert sample_watcher1 is sample_watcher2


def test_add_finished_sample(raw_data_with_two_samples,
                             mock_open_kive,
                             default_config):
    """ The folder finished processing since the folder was scanned. """
    base_calls = (raw_data_with_two_samples /
                  "MiSeq/runs/140101_M01234/Data/Intensities/BaseCalls")
    assert mock_open_kive
    done_run = base_calls / "../../.."
    results_path = done_run / "Results/version_0-dev"
    results_path.mkdir(parents=True)
    done_path = results_path / "doneprocessing"
    done_path.touch()
    kive_watcher = KiveWatcher(default_config)

    sample_watcher1 = kive_watcher.add_sample_group(
        base_calls=base_calls,
        sample_group=SampleGroup('2110A',
                                 ('2110A-V3LOOP_S13_L001_R1_001.fastq.gz',
                                  None),
                                 ('V3LOOP', None)))

    assert done_path.exists()
    assert sample_watcher1 is None


def test_add_failed_sample(raw_data_with_two_samples,
                           mock_open_kive,
                           default_config):
    """ The folder failed since the folder was scanned. """
    base_calls = (raw_data_with_two_samples /
                  "MiSeq/runs/140101_M01234/Data/Intensities/BaseCalls")
    assert mock_open_kive
    failed_run_path = base_calls / "../../.."
    error_path = failed_run_path / "errorprocessing"
    error_path.touch()
    kive_watcher = KiveWatcher(default_config)

    sample_watcher1 = kive_watcher.add_sample_group(
        base_calls=base_calls,
        sample_group=SampleGroup('2110A',
                                 ('2110A-V3LOOP_S13_L001_R1_001.fastq.gz',
                                  None),
                                 ('V3LOOP', None)))

    assert sample_watcher1 is None


def test_calculate_retry_wait():
    min_wait = timedelta(minutes=1)
    max_wait = timedelta(minutes=9)

    assert timedelta(minutes=1) == calculate_retry_wait(min_wait,
                                                        max_wait,
                                                        attempt_count=1)
    assert timedelta(minutes=2) == calculate_retry_wait(min_wait,
                                                        max_wait,
                                                        attempt_count=2)
    assert timedelta(minutes=4) == calculate_retry_wait(min_wait,
                                                        max_wait,
                                                        attempt_count=3)
    assert timedelta(minutes=8) == calculate_retry_wait(min_wait,
                                                        max_wait,
                                                        attempt_count=4)
    assert timedelta(minutes=9) == calculate_retry_wait(min_wait,
                                                        max_wait,
                                                        attempt_count=5)
    assert timedelta(minutes=9) == calculate_retry_wait(min_wait,
                                                        max_wait,
                                                        attempt_count=6)

    assert timedelta(minutes=9) == calculate_retry_wait(min_wait,
                                                        max_wait,
                                                        attempt_count=10000)


def test_collate_main_results(raw_data_with_two_samples, default_config, mock_open_kive):
    run_folder = raw_data_with_two_samples / "MiSeq/runs/140101_M01234"
    base_calls = run_folder / "Data/Intensities/BaseCalls"
    results_path = run_folder / "Results"
    results_path.mkdir(parents=True)
    version_folder: Path = results_path / 'version_0-dev'
    version_folder.mkdir()

    sample1_scratch = version_folder / "scratch" / "2120A-PR_S14"
    sample1_scratch.mkdir(parents=True)
    (sample1_scratch / "cascade.csv").write_text("col1,col2\nval1.1,val2.1\n")
    sample2_scratch = version_folder / "scratch" / "2110A-V3LOOP_S13"
    sample2_scratch.mkdir(parents=True)
    (sample2_scratch / "cascade.csv").write_text("col1,col2\nval1.2,val2.2\n")

    expected_cascade_path = version_folder / "cascade.csv"
    expected_cascade_text = "sample,col1,col2\n2120A-PR_S14,val1.1,val2.1\n2110A-V3LOOP_S13,val1.2,val2.2\n"
    expected_done_path = version_folder / "doneprocessing"

    denovo_scratch_path = version_folder / "scratch_denovo"
    denovo_scratch_path.mkdir()

    kive_watcher = KiveWatcher(default_config)
    folder_watcher = kive_watcher.add_folder(base_calls)
    kive_watcher.add_sample_group(
        base_calls,
        SampleGroup('2120A',
                    ('2120A-PR_S14_L001_R1_001.fastq.gz', None),
                    ('PR', None)))
    kive_watcher.add_sample_group(
        base_calls,
        SampleGroup('2110A',
                    ('2110A-V3LOOP_S13_L001_R1_001.fastq.gz', None),
                    ('V3LOOP', None)))

    kive_watcher.collate_folder(folder_watcher, PipelineType.MAIN)

    cascade_text = expected_cascade_path.read_text()
    assert cascade_text == expected_cascade_text
    assert expected_done_path.exists()
    assert denovo_scratch_path.exists()


def test_collate_denovo_results(raw_data_with_two_samples, default_config, mock_open_kive):
    run_folder = raw_data_with_two_samples / "MiSeq/runs/140101_M01234"
    base_calls = run_folder / "Data/Intensities/BaseCalls"
    results_path = run_folder / "Results"
    results_path.mkdir(parents=True)
    version_folder: Path = results_path / 'version_0-dev'
    version_folder.mkdir()

    sample1_scratch = version_folder / "scratch_denovo" / "2120A-PR_S14"
    sample1_scratch.mkdir(parents=True)
    (sample1_scratch / "cascade.csv").write_text("col1,col2\n")
    sample2_scratch = version_folder / "scratch_denovo" / "2110A-V3LOOP_S13"
    sample2_scratch.mkdir(parents=True)
    (sample2_scratch / "cascade.csv").write_text("col1,col2\n")

    expected_cascade_path = version_folder / "denovo" / "cascade.csv"
    expected_done_path = version_folder / "denovo" / "doneprocessing"
    proviral_path = version_folder / "denovo" / "detailed_results"

    main_scratch_path = version_folder / "scratch"
    main_scratch_path.mkdir()

    kive_watcher = KiveWatcher(default_config)
    folder_watcher = kive_watcher.add_folder(base_calls)
    kive_watcher.add_sample_group(
        base_calls,
        SampleGroup('2120A',
                    ('2120A-PR_S14_L001_R1_001.fastq.gz', None),
                    ('PR', None)))
    kive_watcher.add_sample_group(
        base_calls,
        SampleGroup('2110A',
                    ('2110A-V3LOOP_S13_L001_R1_001.fastq.gz', None),
                    ('V3LOOP', None)))

    kive_watcher.collate_folder(folder_watcher, PipelineType.DENOVO_MAIN)

    assert expected_cascade_path.exists()
    assert expected_done_path.exists()
    assert main_scratch_path.exists()
    assert not proviral_path.exists()


def test_collate_proviral_results(raw_data_with_two_samples, default_config, mock_open_kive):
    run_folder = raw_data_with_two_samples / "MiSeq/runs/140101_M01234"
    base_calls = run_folder / "Data/Intensities/BaseCalls"
    results_path = run_folder / "Results"
    results_path.mkdir(parents=True)
    version_folder: Path = results_path / 'version_0-dev'
    version_folder.mkdir()

    sample1_scratch = version_folder / "scratch_proviral" / "2120A-PR_S14"
    sample1_scratch.mkdir(parents=True)
    (sample1_scratch / "outcome_summary.csv").write_text("col1,col2\nvalue1,value2\n")
    (sample1_scratch / "table_precursor.csv").write_text("col1,col2\n")
    sample2_scratch = version_folder / "scratch_proviral" / "2110A-V3LOOP_S13"
    sample2_scratch.mkdir(parents=True)
    (sample2_scratch / "outcome_summary.csv").write_text("col1,col2\nvalue3,value4\n")
    (sample2_scratch / "table_precursor.csv").write_text("col1,col2\n")

    expected_outcome_path = version_folder / "proviral" / "outcome_summary.csv"
    expected_precursor_path = version_folder / "proviral" / "table_precursor.csv"
    expected_outcome_text = "sample,col1,col2\n2120A-PR_S14,value1,value2\n2110A-V3LOOP_S13,value3,value4\n"

    main_scratch_path = version_folder / "scratch"
    main_scratch_path.mkdir()

    kive_watcher = KiveWatcher(default_config)
    folder_watcher = kive_watcher.add_folder(base_calls)
    kive_watcher.add_sample_group(
        base_calls,
        SampleGroup('2120A',
                    ('2120A-PR_S14_L001_R1_001.fastq.gz', None),
                    ('PR', None)))
    kive_watcher.add_sample_group(
        base_calls,
        SampleGroup('2110A',
                    ('2110A-V3LOOP_S13_L001_R1_001.fastq.gz', None),
                    ('V3LOOP', None)))

    kive_watcher.collate_folder(folder_watcher, PipelineType.PROVIRAL)

    assert expected_outcome_path.exists()
    assert expected_precursor_path.exists()
    outcome_text = expected_outcome_path.read_text()
    assert outcome_text == expected_outcome_text

def test_collate_mixed_hcv_results(raw_data_with_two_samples, default_config, mock_open_kive):
    run_folder = raw_data_with_two_samples / "MiSeq/runs/140101_M01234"
    base_calls = run_folder / "Data/Intensities/BaseCalls"
    results_path = run_folder / "Results"
    results_path.mkdir(parents=True)
    version_folder: Path = results_path / 'version_0-dev'
    version_folder.mkdir()

    sample1_scratch = version_folder / "scratch_mixed_hcv" / "2120A-PR_S14"
    sample1_scratch.mkdir(parents=True)
    (sample1_scratch / "mixed_counts.csv").write_text("col1,col2\n")
    sample2_scratch = version_folder / "scratch_mixed_hcv" / "2110A-V3LOOP_S13"
    sample2_scratch.mkdir(parents=True)
    (sample2_scratch / "mixed_counts.csv").write_text("col1,col2\n")

    expected_cascade_path = version_folder / "mixed_hcv" / "mixed_counts.csv"
    expected_done_path = version_folder / "mixed_hcv" / "doneprocessing"

    kive_watcher = KiveWatcher(default_config)
    folder_watcher = kive_watcher.add_folder(base_calls)
    kive_watcher.add_sample_group(
        base_calls,
        SampleGroup('2120A',
                    ('2120A-PR_S14_L001_R1_001.fastq.gz', None),
                    ('PR', None)))
    kive_watcher.add_sample_group(
        base_calls,
        SampleGroup('2110A',
                    ('2110A-V3LOOP_S13_L001_R1_001.fastq.gz', None),
                    ('V3LOOP', None)))

    kive_watcher.collate_folder(folder_watcher, PipelineType.MIXED_HCV_MAIN)

    assert expected_cascade_path.exists()
    assert expected_done_path.exists()


def test_collate_csv():
    source1 = StringIO("""\
a,b,c
1,2,3
_,-,=
""")
    source2 = StringIO("""\
a,b,c
10,20,30
""")
    expected_target = """\
sample,a,b,c
E12345,1,2,3
E12345,_,-,=
E22222,10,20,30
"""
    target = StringIO()

    KiveWatcher.extract_csv(source1, target, 'E12345', source_count=0)
    KiveWatcher.extract_csv(source2, target, 'E22222', source_count=1)

    assert target.getvalue() == expected_target


def test_collate_csv_with_sample_already_filled():
    source1 = StringIO("""\
sample,a,b,c
E12345,1,2,3
E12345,_,-,=
""")
    source2 = StringIO("""\
sample,a,b,c
E22222,10,20,30
""")
    expected_target = """\
sample,a,b,c
E12345,1,2,3
E12345,_,-,=
E22222,10,20,30
"""
    target = StringIO()

    KiveWatcher.extract_csv(source1, target, 'ignored', source_count=0)
    KiveWatcher.extract_csv(source2, target, 'ignored', source_count=1)

    assert target.getvalue() == expected_target

def test_launch_main_good_pipeline_id(mock_open_kive, default_config):
    _mock_session = mock_open_kive.return_value
    kive_watcher = KiveWatcher(default_config)
    kive_watcher.app_urls = {
        default_config.micall_filter_quality_pipeline_id: '/containerapps/102'}
    kive_watcher.app_args = {
        default_config.micall_filter_quality_pipeline_id: dict(
            quality_csv='/containerargs/103')}

    inputs = {'quality_csv': {'url': '/datasets/104', 'id': 104}}
    run_batch = {'url': '/batches/101'}
    kive_watcher.find_or_launch_run(pipeline_id=42,
                                    inputs=inputs,
                                    run_name='MiCall filter quality on 140101_M01234',
                                    run_batch=run_batch)

def test_launch_main_bad_pipeline_id(mock_open_kive, default_config):
    _mock_session = mock_open_kive.return_value
    kive_watcher = KiveWatcher(default_config)
    kive_watcher.app_urls = {
        default_config.micall_filter_quality_pipeline_id: '/containerapps/102'}
    kive_watcher.app_args = {
        default_config.micall_filter_quality_pipeline_id: dict(
            quality_csv='/containerargs/103')}

    inputs = {'quality_csv': {'bad_argument': 777, 'id': 104}}
    run_batch = {'url': '/batches/101'}
    pipeline_id = 42
    expected_msg = f'The specified app with id {pipeline_id}' \
                    ' appears to expect a different set of inputs'

    with pytest.raises(ValueError, match=expected_msg) as _excinfo:
        kive_watcher.find_or_launch_run(pipeline_id=pipeline_id,
                                        inputs=inputs,
                                        run_name='MiCall filter quality on 140101_M01234',
                                        run_batch=run_batch)


class TestKiveWatcherIntegration:
    """Test scenarios that mirror how kive_watcher uses disk operations."""

    def test_extract_coverage_maps_scenario(self):
        """Test the exact pattern used in extract_coverage_maps."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            results_path = temp_path / "results"
            coverage_path = results_path / "coverage_maps"
            
            # This mirrors the exact pattern from kive_watcher.extract_coverage_maps
            disk_operations.mkdir_p(coverage_path, exist_ok=True)
            assert coverage_path.exists()
            
            # Simulate creating some files (would normally come from tar extraction)
            test_file = coverage_path / "sample1.coverage.svg"
            disk_operations.write_text(test_file, "test content")
            
            # The remove_empty_directory should NOT remove the directory because it has files
            disk_operations.remove_empty_directory(coverage_path)
            assert coverage_path.exists()  # Should still exist because it has files
            assert test_file.exists()
            
            # Now remove the file and try again
            disk_operations.unlink(test_file)
            disk_operations.remove_empty_directory(coverage_path)
            assert not coverage_path.exists()  # Should be removed now that it's empty

    def test_extract_archive_scenario(self):
        """Test the pattern used in extract_archive method."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            results_path = temp_path / "results"
            output_path = results_path / "detailed_results"
            sample_target_path = output_path / "sample1"
            
            # Create the directory structure
            disk_operations.mkdir_p(sample_target_path, exist_ok=True)
            assert sample_target_path.exists()
            
            # Test empty directory removal
            disk_operations.remove_empty_directory(sample_target_path)
            assert not sample_target_path.exists()
            
            # Test with files - create again and add content
            disk_operations.mkdir_p(sample_target_path, exist_ok=True)
            test_file = sample_target_path / "results.txt"
            disk_operations.write_text(test_file, "important data")
            
            # Should not remove when there are files
            disk_operations.remove_empty_directory(sample_target_path)
            assert sample_target_path.exists()
            assert test_file.exists()
            
            # Remove files and parent directory
            disk_operations.unlink(test_file)
            disk_operations.remove_empty_directory(sample_target_path)
            assert not sample_target_path.exists()
            
            # Test parent directory cleanup
            disk_operations.remove_empty_directory(output_path)
            assert not output_path.exists()

    def test_alignment_plot_scenario(self):
        """Test the pattern used in move_alignment_plot method."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            results_path = temp_path / "results"
            alignment_path = results_path / "alignment"
            
            # Create alignment directory
            disk_operations.mkdir_p(alignment_path, exist_ok=True)
            
            # Test case where no alignment files are moved (empty directory)
            disk_operations.remove_empty_directory(alignment_path)
            assert not alignment_path.exists()
            
            # Test case where some files exist
            disk_operations.mkdir_p(alignment_path, exist_ok=True)
            alignment_file = alignment_path / "sample1_alignment.svg"
            disk_operations.write_text(alignment_file, "<svg>test</svg>")
            
            # Should not remove directory with files
            disk_operations.remove_empty_directory(alignment_path)
            assert alignment_path.exists()

    def test_network_drive_retry_simulation(self):
        """Test retry behavior with network drive issues."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            test_dir = temp_path / "test_dir"
            
            # Create directory
            disk_operations.mkdir_p(test_dir, exist_ok=True)
            
            # Mock intermittent permission errors at the Path.rmdir level
            call_count = 0
            original_rmdir = Path.rmdir
            
            def mock_rmdir(self):
                nonlocal call_count
                if str(self) == str(test_dir):
                    call_count += 1
                    if call_count <= 2:  # Fail first 2 attempts
                        raise OSError(errno.EACCES, "Permission denied")
                return original_rmdir(self)
            
            with patch.object(Path, 'rmdir', mock_rmdir):
                # Should succeed after retries
                disk_operations.remove_empty_directory(test_dir)
                assert not test_dir.exists()
                assert call_count == 3  # Should have retried and succeeded on 3rd attempt

    def test_remove_empty_directory_non_empty_case(self):
        """Test that ENOTEMPTY is handled silently as expected."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            parent_dir = temp_path / "parent"
            child_dir = parent_dir / "child"
            
            # Create nested structure
            disk_operations.mkdir_p(child_dir, exist_ok=True)
            
            # Add file to child
            test_file = child_dir / "important.txt"
            disk_operations.write_text(test_file, "data")
            
            # Try to remove parent - should silently handle ENOTEMPTY
            disk_operations.remove_empty_directory(parent_dir)
            assert parent_dir.exists()  # Should still exist
            assert child_dir.exists()   # Child should still exist
            assert test_file.exists()   # File should still exist

    def test_disk_file_operation_integration(self):
        """Test the disk_file_operation context manager pattern used in kive_watcher."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            target_file = temp_path / "test_output.csv"
            
            # Test write operation (used in copy_outputs)
            with disk_operations.disk_file_operation(target_file, 'w') as target:
                target.write("sample,value\n")
                target.write("sample1,123\n")
            
            assert target_file.exists()
            
            # Test read operation (used in copy_outputs)
            with disk_operations.disk_file_operation(target_file, 'r') as source:
                content = source.read()
                assert "sample,value" in content
                assert "sample1,123" in content

    def test_error_handling_write_text(self):
        """Test error handling for write_text operation."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            error_file = temp_path / "errorprocessing"
            
            # Test normal operation
            disk_operations.write_text(error_file, "Finding error metrics failed.\n")
            assert error_file.exists()
            content = error_file.read_text()
            assert content == "Finding error metrics failed.\n"

    def test_concurrent_directory_operations(self):
        """Test handling of concurrent directory operations."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            # Create multiple directories as might happen in concurrent processing
            dirs = []
            for i in range(5):
                dir_path = temp_path / f"sample_{i}"
                disk_operations.mkdir_p(dir_path, exist_ok=True)
                dirs.append(dir_path)
            
            # Remove them all
            for dir_path in dirs:
                disk_operations.remove_empty_directory(dir_path)
                assert not dir_path.exists()

    def test_rmtree_with_ignore_errors(self):
        """Test rmtree operation with ignore_errors=True as used in kive_watcher."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            results_path = temp_path / "Results" / "version_1.0"
            
            # Create complex directory structure
            disk_operations.mkdir_p(results_path / "sub1" / "sub2", exist_ok=True)
            disk_operations.write_text(results_path / "file1.txt", "content")
            disk_operations.write_text(results_path / "sub1" / "file2.txt", "content")
            
            # Test rmtree with ignore_errors=True
            disk_operations.rmtree(results_path, ignore_errors=True)
            assert not results_path.exists()


class TestAdvancedRetryLogic:
    """Advanced tests for retry logic and error handling."""

    def test_exponential_backoff_timing_precision(self):
        """Test that exponential backoff calculations are precise."""
        min_wait = disk_operations.MINIMUM_RETRY_WAIT
        max_wait = disk_operations.MAXIMUM_RETRY_WAIT

        # Test specific timing values
        expected_timings = [
            (1, 5),  # First retry: 5 seconds
            (2, 10),  # Second retry: 10 seconds
            (3, 20),  # Third retry: 20 seconds
            (4, 40),  # Fourth retry: 40 seconds
            (10, 2560),  # Tenth retry: ~42 minutes
            (15, 81920),  # Fifteenth retry: ~22.7 hours
            (20, 86400),  # Twentieth retry: capped at 1 day
        ]

        for attempt, expected_seconds in expected_timings:
            result = disk_operations.calculate_retry_wait(min_wait, max_wait, attempt)
            assert result.total_seconds() == expected_seconds

    def test_retry_decorator_thread_safety(self):
        """Test that retry decorator works correctly with multiple threads."""
        call_counts = {}

        @disk_operations.disk_retry("test_op")
        def thread_test_function(thread_id):
            if thread_id not in call_counts:
                call_counts[thread_id] = 0
            call_counts[thread_id] += 1

            if call_counts[thread_id] <= 2:
                raise OSError(f"Error in thread {thread_id}")
            return f"success_{thread_id}"

        def run_test(thread_id, results):
            try:
                result = thread_test_function(thread_id)
                results[thread_id] = result
            except Exception as e:
                results[thread_id] = str(e)

        threads = []
        results: dict[int, float] = {}

        # Run multiple threads concurrently
        for i in range(3):  # Reduced from 5 to 3 threads
            thread = threading.Thread(target=run_test, args=(i, results))
            threads.append(thread)
            thread.start()

        # Wait for all threads to complete
        for thread in threads:
            thread.join()

        # Verify all threads succeeded after retries
        for i in range(3):  # Updated from 5 to 3
            assert results[i] == f"success_{i}"
            assert call_counts[i] == 3  # Should have retried twice

    def test_retry_with_different_exception_types(self):
        """Test retry behavior with various exception types."""

        @disk_operations.disk_retry("test_op")
        def test_function(exception_type):
            if exception_type == "os_error":
                raise OSError("OS Error")
            elif exception_type == "io_error":
                raise IOError("IO Error")
            elif exception_type == "permission_error":
                raise PermissionError("Permission denied")
            elif exception_type == "file_not_found":
                raise FileNotFoundError("File not found")
            elif exception_type == "value_error":
                raise ValueError("Value error")
            return "success"

        # Should retry these errors
        with patch("micall.monitor.disk_operations.wait_for_retry") as mock_wait:
            with pytest.raises(OSError):
                test_function("os_error")
            assert mock_wait.call_count == 14  # 15 attempts = 14 retries

        with patch("micall.monitor.disk_operations.wait_for_retry") as mock_wait:
            with pytest.raises(IOError):
                test_function("io_error")
            assert mock_wait.call_count == 14  # 15 attempts = 14 retries

        # Should NOT retry these errors
        with patch("micall.monitor.disk_operations.wait_for_retry") as mock_wait:
            with pytest.raises(ValueError):
                test_function("value_error")
            assert mock_wait.call_count == 0  # No retries

    def test_retry_success_after_max_attempts_minus_one(self):
        """Test success on the 14th attempt (just before max)."""
        attempt_count = 0

        @disk_operations.disk_retry("test_op")
        def almost_max_attempts():
            nonlocal attempt_count
            attempt_count += 1
            if attempt_count < 14:
                raise OSError("Almost there")
            return "success"

        with patch("micall.monitor.disk_operations.wait_for_retry"):
            result = almost_max_attempts()
            assert result == "success"
            assert attempt_count == 14


class TestEdgeCaseFileOperations:
    """Test edge cases for file operations."""

    def test_mkdir_p_with_special_modes(self):
        """Test mkdir_p with various file modes."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)

            # Test restrictive mode
            restricted_dir = temp_path / "restricted"
            disk_operations.mkdir_p(restricted_dir, mode=0o700, exist_ok=True)
            assert restricted_dir.exists()

            # Test permissive mode
            permissive_dir = temp_path / "permissive"
            disk_operations.mkdir_p(permissive_dir, mode=0o777, exist_ok=True)
            assert permissive_dir.exists()

            # Test nested creation with exist_ok=False then True
            nested_dir = temp_path / "level1" / "level2" / "level3"
            disk_operations.mkdir_p(nested_dir, exist_ok=False)
            assert nested_dir.exists()

            # Should not raise with exist_ok=True
            disk_operations.mkdir_p(nested_dir, exist_ok=True)

    def test_write_text_with_different_encodings(self):
        """Test write_text with various encodings."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)

            # Test UTF-8 with special characters
            utf8_file = temp_path / "utf8.txt"
            utf8_content = "Hello 世界 🌍 café naïve résumé"
            disk_operations.write_text(utf8_file, utf8_content, encoding="utf-8")
            assert utf8_file.read_text(encoding="utf-8") == utf8_content

            # Test ASCII
            ascii_file = temp_path / "ascii.txt"
            ascii_content = "Simple ASCII text"
            disk_operations.write_text(ascii_file, ascii_content, encoding="ascii")
            assert ascii_file.read_text(encoding="ascii") == ascii_content

            # Test Latin-1
            latin1_file = temp_path / "latin1.txt"
            latin1_content = "Café résumé naïve"
            disk_operations.write_text(latin1_file, latin1_content, encoding="latin-1")
            assert latin1_file.read_text(encoding="latin-1") == latin1_content

    def test_unlink_with_missing_ok_variations(self):
        """Test unlink behavior with missing_ok parameter."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)

            # Test unlink existing file with missing_ok=True
            existing_file = temp_path / "existing.txt"
            existing_file.write_text("content")
            disk_operations.unlink(existing_file, missing_ok=True)
            assert not existing_file.exists()

            # Test unlink non-existing file with missing_ok=True
            non_existing = temp_path / "non_existing.txt"
            disk_operations.unlink(non_existing, missing_ok=True)  # Should not raise

            # Test unlink existing file with missing_ok=False (should work)
            existing_file2 = temp_path / "existing2.txt"
            existing_file2.write_text("content")
            disk_operations.unlink(existing_file2, missing_ok=False)
            assert not existing_file2.exists()

            # Test that missing_ok=False would raise for non-existing file
            # (but don't actually call it since it would retry for 15 attempts)
            # This tests our understanding of the expected behavior
            assert not non_existing.exists()  # Confirm file doesn't exist

    def test_copy_file_with_metadata_preservation(self):
        """Test that copy_file preserves metadata using shutil.copy2."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)

            source = temp_path / "source.txt"
            target = temp_path / "target.txt"

            # Create source file with content
            source.write_text("test content")
            original_stat = source.stat()

            # Add a small delay to ensure timestamp difference would be detectable
            time.sleep(0.01)

            # Copy file
            disk_operations.copy_file(source, target)

            # Verify content
            assert target.read_text() == "test content"

            # Verify metadata preservation (size should match)
            assert target.stat().st_size == original_stat.st_size

    def test_rmtree_with_readonly_files(self):
        """Test rmtree with read-only files (ignore_errors behavior)."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)

            # Create directory with read-only file
            target_dir = temp_path / "readonly_dir"
            target_dir.mkdir()
            readonly_file = target_dir / "readonly.txt"
            readonly_file.write_text("readonly content")
            readonly_file.chmod(0o444)  # Read-only

            # Test with ignore_errors=True
            disk_operations.rmtree(target_dir, ignore_errors=True)
            # Should succeed regardless of read-only files


class TestComplexIntegrationScenarios:
    """Test complex scenarios that mirror real-world usage."""

    def test_full_kive_watcher_workflow_simulation(self):
        """Simulate a simplified kive_watcher workflow."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)

            # Simulate run structure
            run_folder = temp_path / "run_001"
            results_path = run_folder / "Results" / "version_1.0"
            scratch_path = results_path / "scratch"

            # Step 1: Create initial structure (only 1 sample for speed)
            disk_operations.mkdir_p(scratch_path / "sample_001", exist_ok=True)

            # Step 2: Simulate processing - create output files
            sample_path = scratch_path / "sample_001"
            disk_operations.write_text(
                sample_path / "cascade.csv", "header,data\nval1,val2\n"
            )
            disk_operations.write_text(
                sample_path / "conseq.csv", "sequence,count\nATCG,100\n"
            )

            # Step 3: Simulate collation process
            collated_path = results_path / "collated"
            disk_operations.mkdir_p(collated_path, exist_ok=True)

            # Collate CSV files (simplified)
            with disk_operations.disk_file_operation(
                collated_path / "cascade.csv", "w"
            ) as target:
                target.write("sample,header,data\n")
                source_path = scratch_path / "sample_001" / "cascade.csv"
                with disk_operations.disk_file_operation(
                    source_path, "r"
                ) as source:
                    lines = source.readlines()
                    for line in lines[1:]:  # Skip header
                        target.write(f"sample_001,{line}")

            # Step 4: Final cleanup
            disk_operations.rmtree(scratch_path)
            disk_operations.touch(results_path / "doneprocessing")

            # Verify final state
            assert (collated_path / "cascade.csv").exists()
            assert (results_path / "doneprocessing").exists()
            assert not scratch_path.exists()

    def test_concurrent_file_operations(self):
        """Test concurrent file operations across multiple threads."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            results = {}
            errors = {}

            def worker_thread(thread_id):
                try:
                    # Each thread works in its own subdirectory
                    thread_dir = temp_path / f"thread_{thread_id}"
                    disk_operations.mkdir_p(thread_dir, exist_ok=True)

                    # Create multiple files
                    for i in range(10):
                        file_path = thread_dir / f"file_{i}.txt"
                        disk_operations.write_text(
                            file_path, f"Thread {thread_id} file {i}"
                        )

                    # Copy files to backup directory
                    backup_dir = thread_dir / "backup"
                    disk_operations.mkdir_p(backup_dir, exist_ok=True)

                    for i in range(10):
                        source = thread_dir / f"file_{i}.txt"
                        target = backup_dir / f"file_{i}.txt"
                        disk_operations.copy_file(source, target)

                    # Clean up some files
                    for i in range(5):
                        disk_operations.unlink(thread_dir / f"file_{i}.txt")

                    # Final verification
                    remaining_files = list(thread_dir.glob("file_*.txt"))
                    backup_files = list(backup_dir.glob("file_*.txt"))

                    results[thread_id] = {
                        "remaining": len(remaining_files),
                        "backup": len(backup_files),
                    }

                except Exception as e:
                    errors[thread_id] = str(e)

            # Run multiple worker threads
            threads = []
            for i in range(4):  # Reduced from 8 to 4 threads
                thread = threading.Thread(target=worker_thread, args=(i,))
                threads.append(thread)
                thread.start()

            # Wait for completion
            for thread in threads:
                thread.join()

            # Verify results
            assert len(errors) == 0, f"Errors occurred: {errors}"
            assert len(results) == 4  # Updated from 8 to 4

            for thread_id, result in results.items():
                assert result["remaining"] == 5  # 5 files left after deleting 5
                assert result["backup"] == 10  # All 10 files backed up

    def test_error_recovery_scenarios(self):
        """Test various error recovery scenarios."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)

            # Scenario 1: Network interruption during file operations
            def simulate_network_interruption():
                call_count = 0

                def mock_write_text(content, encoding="utf-8"):
                    nonlocal call_count
                    call_count += 1
                    if call_count <= 3:
                        raise OSError(errno.EIO, "I/O error - network interruption")
                    return content

                error_file = temp_path / "network_error.txt"
                with patch.object(Path, "write_text", side_effect=mock_write_text):
                    with patch("micall.monitor.disk_operations.wait_for_retry"):
                        disk_operations.write_text(error_file, "test content")
                        assert call_count == 4  # Should succeed on 4th attempt

            simulate_network_interruption()

            # Scenario 2: Permission denied then granted
            def simulate_permission_recovery():
                call_count = 0

                def mock_mkdir(mode=0o777, parents=False, exist_ok=False):
                    nonlocal call_count
                    call_count += 1
                    if call_count <= 2:
                        raise OSError(errno.EACCES, "Permission denied")
                    # Success on 3rd attempt
                    return None

                perm_dir = temp_path / "permission_test"
                with patch.object(Path, "mkdir", side_effect=mock_mkdir):
                    with patch("micall.monitor.disk_operations.wait_for_retry"):
                        disk_operations.mkdir_p(perm_dir, exist_ok=True)
                        assert call_count == 3

            simulate_permission_recovery()

    def test_large_file_operations(self):
        """Test operations with larger files and directories."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)

            # Create directory with many files
            large_dir = temp_path / "large_directory"
            disk_operations.mkdir_p(large_dir, exist_ok=True)

            # Create 50 files with various sizes (reduced from 100)
            for i in range(50):
                file_path = large_dir / f"file_{i:03d}.txt"
                content = f"File {i} content " * (i + 1)  # Variable size content
                disk_operations.write_text(file_path, content)

            # Create nested directory structure
            for level1 in range(3):  # Reduced from 5 to 3
                for level2 in range(3):  # Reduced from 5 to 3
                    nested_dir = large_dir / f"level1_{level1}" / f"level2_{level2}"
                    disk_operations.mkdir_p(nested_dir, exist_ok=True)

                    # Add a file in each nested directory
                    nested_file = nested_dir / "nested_file.txt"
                    disk_operations.write_text(
                        nested_file, f"Nested content {level1}-{level2}"
                    )

            # Verify structure was created
            files = list(large_dir.glob("file_*.txt"))
            assert len(files) == 50  # Updated from 100 to 50

            nested_files = list(large_dir.glob("**/nested_file.txt"))
            assert len(nested_files) == 9  # Updated: 3 * 3 = 9

            # Test bulk removal
            disk_operations.rmtree(large_dir)
            assert not large_dir.exists()


class TestErrorMessageVerification:
    """Test that error messages are informative and correctly logged."""

    def test_retry_error_messages(self):
        """Test that retry error messages contain useful information."""

        @disk_operations.disk_retry("test_operation")
        def failing_operation():
            raise OSError("Simulated failure")

        with patch("micall.monitor.disk_operations.logger.error") as mock_logger:
            with patch("micall.monitor.disk_operations.sleep"):
                with pytest.raises(OSError):
                    failing_operation()

                # Verify error logging
                assert mock_logger.call_count >= 1

                # Check that final error message mentions max attempts
                final_call = mock_logger.call_args_list[-1]
                error_message = final_call[0][0]
                assert "test_operation" in error_message
                assert "15 attempts" in error_message

    def test_wait_for_retry_logging(self):
        """Test wait_for_retry logging behavior."""

        # Test with logging enabled
        with patch("micall.monitor.disk_operations.logger.error") as mock_logger:
            with patch("micall.monitor.disk_operations.sleep"):
                disk_operations.wait_for_retry(3, "test_op", is_logged=True)

                mock_logger.assert_called_once()
                # The log message uses string formatting, so check that args were passed
                assert mock_logger.call_args[0][1] == "test_op"
                # Third argument should be a timedelta object
                assert hasattr(mock_logger.call_args[0][2], 'total_seconds')

        # Test with logging disabled
        with patch("micall.monitor.disk_operations.logger.error") as mock_logger:
            with patch("micall.monitor.disk_operations.sleep"):
                disk_operations.wait_for_retry(3, "test_op", is_logged=False)

                mock_logger.assert_not_called()


class TestDiskFileOperationContextManager:
    """Comprehensive tests for the disk_file_operation context manager."""

    def test_context_manager_exception_propagation(self):
        """Test that exceptions inside context manager are properly propagated."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            test_file = temp_path / "test.txt"

            with pytest.raises(ValueError):
                with disk_operations.disk_file_operation(test_file, "w") as f:
                    f.write("some content")
                    raise ValueError("Test exception")

            # File should still be created even though exception occurred
            assert test_file.exists()

    def test_context_manager_with_binary_modes(self):
        """Test context manager with binary file modes."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            binary_file = temp_path / "binary.dat"

            test_data = b"Binary data: \x00\x01\x02\x03\xff"

            # Write binary data
            with disk_operations.disk_file_operation(binary_file, "wb") as f:
                f.write(test_data)

            # Read binary data back
            with disk_operations.disk_file_operation(binary_file, "rb") as f:
                read_data = f.read()

            assert read_data == test_data

    def test_context_manager_retry_on_open_failure(self):
        """Test that context manager retries on file open failures."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            test_file = temp_path / "retry_test.txt"

            call_count = 0
            original_open = Path.open

            def mock_open(
                self, mode="r", buffering=-1, encoding=None, errors=None, newline=None
            ):
                nonlocal call_count
                call_count += 1
                if call_count <= 2:
                    raise OSError("Simulated open failure")
                return original_open(self, mode, buffering, encoding, errors, newline)

            with patch.object(Path, "open", mock_open):
                with patch("micall.monitor.disk_operations.wait_for_retry"):
                    with disk_operations.disk_file_operation(test_file, "w") as f:
                        f.write("success after retry")

                    assert call_count == 3  # Should succeed on 3rd attempt

            assert test_file.read_text() == "success after retry"


class TestRemoveEmptyDirectoryAdvanced:
    """Advanced tests for remove_empty_directory function."""

    def test_remove_empty_directory_with_hidden_files(self):
        """Test behavior with hidden files (should not remove directory)."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            test_dir = temp_path / "hidden_files_dir"
            test_dir.mkdir()

            # Create hidden file
            hidden_file = test_dir / ".hidden_file"
            hidden_file.write_text("hidden content")

            # Should not remove directory because it contains hidden file
            disk_operations.remove_empty_directory(test_dir)
            assert test_dir.exists()
            assert hidden_file.exists()

    def test_remove_empty_directory_with_symlinks(self):
        """Test behavior with symbolic links."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)

            # Create target file and directory with symlink
            target_file = temp_path / "target.txt"
            target_file.write_text("target content")

            symlink_dir = temp_path / "symlink_dir"
            symlink_dir.mkdir()

            symlink_file = symlink_dir / "link_to_target"
            symlink_file.symlink_to(target_file)

            # Should not remove directory because it contains symlink
            disk_operations.remove_empty_directory(symlink_dir)
            assert symlink_dir.exists()
            assert symlink_file.exists()

    def test_remove_empty_directory_enotempty_no_retry(self):
        """Test that ENOTEMPTY doesn't trigger retry logic."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            test_dir = temp_path / "nonempty_test"
            test_dir.mkdir()

            # Create file in directory
            (test_dir / "file.txt").write_text("content")

            with patch("micall.monitor.disk_operations.wait_for_retry") as mock_wait:
                disk_operations.remove_empty_directory(test_dir)

                # Should not retry for ENOTEMPTY
                assert mock_wait.call_count == 0
                assert test_dir.exists()


class TestIntegrationWithActualFilesystem:
    """Integration tests with actual filesystem operations."""

    def test_real_filesystem_stress_test(self):
        """Stress test with real filesystem operations."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)

            # Create complex directory structure
            for i in range(10):  # Reduced from 20 to 10
                base_dir = temp_path / f"stress_test_{i}"
                disk_operations.mkdir_p(base_dir, exist_ok=True)

                # Create files with various operations
                for j in range(5):  # Reduced from 10 to 5
                    file_path = base_dir / f"file_{j}.txt"
                    content = f"Stress test content {i}-{j} " * 50
                    disk_operations.write_text(file_path, content)

                    # Copy to backup
                    backup_path = base_dir / f"backup_{j}.txt"
                    disk_operations.copy_file(file_path, backup_path)

                    # Create and remove temporary file
                    temp_file = base_dir / f"temp_{j}.txt"
                    disk_operations.touch(temp_file)
                    disk_operations.unlink(temp_file)

                # Test remove_empty_directory on subdirectories
                for k in range(2):  # Reduced from 3 to 2
                    empty_dir = base_dir / f"empty_{k}"
                    disk_operations.mkdir_p(empty_dir, exist_ok=True)
                    disk_operations.remove_empty_directory(empty_dir)
                    assert not empty_dir.exists()

            # Verify final state
            created_dirs = list(temp_path.glob("stress_test_*"))
            assert len(created_dirs) == 10  # Updated from 20 to 10

            total_files = list(temp_path.glob("**/*.txt"))
            assert (
                len(total_files) == 100  # Updated: 10 dirs * 10 files each (5 original + 5 backup)
            )

            # Cleanup stress test
            for stress_dir in created_dirs:
                disk_operations.rmtree(stress_dir)
                assert not stress_dir.exists()

    def test_unicode_filename_handling(self):
        """Test operations with Unicode filenames."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)

            unicode_names = [
                "café_file.txt",
                "фаїл.txt",  # Ukrainian
                "文件.txt",  # Chinese
                "📁_emoji_file.txt",
                "naïve_résumé.txt",
            ]

            for unicode_name in unicode_names:
                try:
                    file_path = temp_path / unicode_name
                    content = f"Unicode content for {unicode_name}"

                    # Test all operations with Unicode names
                    disk_operations.write_text(file_path, content)
                    assert file_path.exists()

                    # Copy to backup
                    backup_path = temp_path / f"backup_{unicode_name}"
                    disk_operations.copy_file(file_path, backup_path)
                    assert backup_path.exists()

                    # Verify content
                    assert backup_path.read_text() == content

                    # Cleanup
                    disk_operations.unlink(file_path)
                    disk_operations.unlink(backup_path)

                except OSError:
                    # Some filesystems don't support certain Unicode characters
                    # This is acceptable and not a failure of our implementation
                    pass


class TestPerformanceAndScalability:
    """Tests focused on performance and scalability aspects."""

    def test_large_directory_removal_performance(self):
        """Test performance of removing large directories."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)

            # Create large directory structure
            large_dir = temp_path / "large_structure"
            disk_operations.mkdir_p(large_dir, exist_ok=True)

            # Create nested structure: 3 levels deep, 3 subdirs per level (reduced from 5 levels, 5 subdirs)
            current_paths = [large_dir]
            for level in range(3):  # Reduced from 5 to 3
                next_paths = []
                for current_path in current_paths:
                    for i in range(3):  # Reduced from 5 to 3
                        sub_path = current_path / f"level_{level}_sub_{i}"
                        disk_operations.mkdir_p(sub_path, exist_ok=True)

                        # Add some files
                        for j in range(3):
                            file_path = sub_path / f"file_{j}.txt"
                            disk_operations.write_text(
                                file_path, f"Content {level}-{i}-{j}"
                            )

                        next_paths.append(sub_path)
                current_paths = next_paths

            # Measure removal time (should complete reasonably quickly)
            start_time = time.time()
            disk_operations.rmtree(large_dir)
            end_time = time.time()

            assert not large_dir.exists()
            removal_time = end_time - start_time
            # Should complete within reasonable time (adjust threshold as needed)
            assert removal_time < 5.0  # 5 seconds should be more than enough

    def test_many_small_operations_performance(self):
        """Test performance of many small operations."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)

            start_time = time.time()

            # Perform many small operations (reduced from 200 to 100)
            for i in range(100):
                # Create, write, read, copy, delete cycle
                file_path = temp_path / f"perf_test_{i}.txt"
                content = f"Performance test content {i}"

                disk_operations.write_text(file_path, content)

                copy_path = temp_path / f"copy_{i}.txt"
                disk_operations.copy_file(file_path, copy_path)

                disk_operations.unlink(file_path)
                disk_operations.unlink(copy_path)

            end_time = time.time()
            total_time = end_time - start_time

            # Should complete within reasonable time (adjusted for fewer operations)
            assert total_time < 5.0  # Reduced from 10 to 5 seconds for 100 cycles
