import tarfile
from gzip import GzipFile
from io import BytesIO
from pathlib import Path
from queue import Full
from tarfile import TarInfo
from unittest.mock import patch, ANY, Mock, call, MagicMock

import pytest
from datetime import datetime, timedelta

from struct import pack

from kiveapi import KiveClientException, KiveRunFailedException
from requests import ConnectionError

from micall.monitor.kive_watcher import find_samples, KiveWatcher, FolderEvent, FolderEventType, calculate_retry_wait
from micall.monitor.sample_watcher import PipelineType, ALLOWED_GROUPS, FolderWatcher, SampleWatcher
from micall.monitor.find_groups import SampleGroup
from micall_watcher import parse_args


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
        mock_pipeline = mock_session.get_pipeline.return_value
        mock_input = Mock(dataset_name='quality_csv')
        mock_pipeline.inputs = [mock_input]

        # By default, all runs are still running.
        mock_session.get_run.return_value.raw = dict(end_time=None,
                                                     stopped_by=None)

        yield mock_open_kive


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
        assert not self.expected_puts  # All expected puts arrived.


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


def test_find_default_pipelines(mock_open_kive):
    mock_session = mock_open_kive.return_value
    mock_session.get.side_effect = [
        Mock(name='get external file directories',
             **{'json.return_value': []}),
        Mock(name='get pipeline families',
             **{'json.return_value': [
                 dict(id=42,
                      name='Other MiCall',
                      members=[dict(id=420)]),
                 dict(id=43,
                      name='MiCall Filter Quality',
                      members=[dict(id=435),
                               dict(id=432)]),
                 dict(id=44,
                      name='MiCall Main',
                      members=[dict(id=440)]),
                 dict(id=45,
                      name='MiCall Resistance',
                      members=[dict(id=450)])]})]
    expected_filter_quality_pipeline_id = 435
    expected_main_pipeline_id = 440
    expected_resistance_pipeline_id = 450
    args = parse_args([])  # No parameters: all defaults.

    kive_watcher = KiveWatcher(args)

    assert expected_filter_quality_pipeline_id == \
        kive_watcher.config.micall_filter_quality_pipeline_id
    assert expected_main_pipeline_id == \
        kive_watcher.config.micall_main_pipeline_id
    assert expected_resistance_pipeline_id == \
        kive_watcher.config.micall_resistance_pipeline_id


def test_default_pipeline_not_found(mock_open_kive):
    mock_session = mock_open_kive.return_value
    mock_session.get.side_effect = [
        Mock(name='get external file directories',
             **{'json.return_value': []}),
        Mock(name='get pipeline families',
             **{'json.return_value': [
                 dict(id=43,
                      name='MiCall Filter Quality',
                      members=[dict(id=435)]),
                 dict(id=44,
                      name='MiCrawl Main',  # <== Typo
                      members=[dict(id=440)]),
                 dict(id=45,
                      name='MiCall Resistance',
                      members=[dict(id=450)])]})]
    args = parse_args([])  # No parameters: all defaults.

    with pytest.raises(
            RuntimeError,
            match=r"Argument micall_main_pipeline_id not set, and no "
                  r"pipeline found named 'micall main'\."):
        KiveWatcher(args)


def test_hcv_pair(raw_data_with_hcv_pair):
    sample_queue = DummyQueueSink()
    sample_queue.expect_put(
        FolderEvent(raw_data_with_hcv_pair / "MiSeq/runs/140101_M01234" /
                    "Data/Intensities/BaseCalls",
                    FolderEventType.ADD_SAMPLE,
                    SampleGroup('2130A',
                                ('2130A-HCV_S15_L001_R1_001.fastq.gz',
                                 '2130AMIDI-MidHCV_S16_L001_R1_001.fastq.gz'))))
    sample_queue.expect_put(
        FolderEvent(raw_data_with_hcv_pair / "MiSeq/runs/140101_M01234" /
                    "Data/Intensities/BaseCalls",
                    FolderEventType.FINISH_FOLDER,
                    None))
    pipeline_version = 'XXX'

    find_samples(raw_data_with_hcv_pair, pipeline_version, sample_queue, wait=False)

    sample_queue.verify()


def test_two_runs(raw_data_with_two_runs):
    sample_queue = DummyQueueSink()
    sample_queue.expect_put(
        FolderEvent(raw_data_with_two_runs / "MiSeq/runs/140201_M01234" /
                    "Data/Intensities/BaseCalls",
                    FolderEventType.ADD_SAMPLE,
                    SampleGroup('2010A',
                                ('2010A-V3LOOP_S3_L001_R1_001.fastq.gz',
                                 None))))
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
                                 None))))
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
                                 None))))
    sample_queue.expect_put(
        FolderEvent(raw_data_with_two_samples / "MiSeq/runs/140101_M01234" /
                    "Data/Intensities/BaseCalls",
                    FolderEventType.ADD_SAMPLE,
                    SampleGroup('2110A',
                                ('2110A-V3LOOP_S13_L001_R1_001.fastq.gz',
                                 None))))
    sample_queue.expect_put(
        FolderEvent(raw_data_with_two_samples / "MiSeq/runs/140101_M01234" /
                    "Data/Intensities/BaseCalls",
                    FolderEventType.FINISH_FOLDER,
                    None))
    pipeline_version = 'XXX'

    find_samples(raw_data_with_two_samples, pipeline_version, sample_queue, wait=False)

    sample_queue.verify()


def test_skip_done_runs(raw_data_with_two_runs):
    done_run = raw_data_with_two_runs / "MiSeq/runs/140201_M01234"
    results_path = done_run / "Results/version_0-dev"
    results_path.mkdir(parents=True)
    done_path = results_path / "doneprocessing"
    done_path.touch()
    pipeline_version = '0-dev'
    sample_queue = DummyQueueSink()
    sample_queue.expect_put(
        FolderEvent(raw_data_with_two_runs / "MiSeq/runs/140101_M01234" /
                    "Data/Intensities/BaseCalls",
                    FolderEventType.ADD_SAMPLE,
                    SampleGroup('2000A',
                                ('2000A-V3LOOP_S2_L001_R1_001.fastq.gz',
                                 None))))
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
                                 None))))
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


def test_full_queue(raw_data_with_two_runs):
    sample_queue = DummyQueueSink()
    sample_queue.expect_put(
        FolderEvent(raw_data_with_two_runs / "MiSeq/runs/140201_M01234" /
                    "Data/Intensities/BaseCalls",
                    FolderEventType.ADD_SAMPLE,
                    SampleGroup('2010A',
                                ('2010A-V3LOOP_S3_L001_R1_001.fastq.gz',
                                 None))))
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
                                     None)))
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
                                     None)))
    finish1 = FolderEvent(raw_data_with_two_runs / "MiSeq/runs/140101_M01234" /
                          "Data/Intensities/BaseCalls",
                          FolderEventType.FINISH_FOLDER,
                          None)
    item2 = FolderEvent(raw_data_with_two_runs / "MiSeq/runs/140201_M01234" /
                        "Data/Intensities/BaseCalls",
                        FolderEventType.ADD_SAMPLE,
                        SampleGroup('2010A',
                                    ('2010A-V3LOOP_S3_L001_R1_001.fastq.gz',
                                     None)))
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


def test_starts_empty(default_config):
    kive_watcher = KiveWatcher(default_config)

    assert not kive_watcher.is_full()


def test_get_kive_pipeline(mock_open_kive, pipelines_config):
    mock_session = mock_open_kive.return_value
    expected_pipeline = mock_session.get_pipeline.return_value
    kive_watcher = KiveWatcher(pipelines_config)

    pipeline1 = kive_watcher.get_kive_pipeline(
        pipelines_config.micall_main_pipeline_id)

    assert expected_pipeline is pipeline1
    mock_session.get_pipeline.assert_called_once_with(
        pipelines_config.micall_main_pipeline_id)


def test_get_pipeline_cached(mock_open_kive, pipelines_config):
    mock_session = mock_open_kive.return_value
    kive_watcher = KiveWatcher(pipelines_config)

    pipeline1 = kive_watcher.get_kive_pipeline(pipelines_config.micall_main_pipeline_id)
    pipeline2 = kive_watcher.get_kive_pipeline(pipelines_config.micall_main_pipeline_id)

    assert pipeline1 is pipeline2
    mock_session.get_pipeline.assert_called_once_with(
        pipelines_config.micall_main_pipeline_id)


def test_get_kive_input(mock_open_kive, pipelines_config):
    kive_watcher = KiveWatcher(pipelines_config)
    mock_pipeline = Mock()
    kive_watcher.pipelines[
        pipelines_config.micall_main_pipeline_id] = mock_pipeline
    expected_input = Mock(dataset_name='bad_cycles_csv')
    mock_pipeline.inputs = [Mock(), expected_input]

    kive_input = kive_watcher.get_kive_input('bad_cycles_csv')

    assert expected_input is kive_input
    mock_open_kive.assert_called_once()


def test_get_kive_input_wrong_pipeline(mock_open_kive, pipelines_config):
    pipelines_config.micall_resistance_pipeline_id = 44
    kive_watcher = KiveWatcher(pipelines_config)
    mock_pipeline = Mock(name='resistance_pipeline')
    kive_watcher.pipelines[
        pipelines_config.micall_resistance_pipeline_id] = mock_pipeline
    mock_pipeline.inputs = [Mock(), Mock(dataset_name='bad_cycles_csv')]

    with pytest.raises(
            ValueError,
            match=r'Input main_amino_csv not found on pipeline id 44\.'):
        kive_watcher.get_kive_input('main_amino_csv')

    mock_open_kive.assert_called_once()


def test_add_first_sample(raw_data_with_two_samples, mock_open_kive, default_config):
    base_calls = (raw_data_with_two_samples /
                  "MiSeq/runs/140101_M01234/Data/Intensities/BaseCalls")
    mock_session = mock_open_kive.return_value
    result_path = base_calls / "../../../Results/version_0-dev"
    result_path.mkdir(parents=True)
    old_stuff_csv = result_path / 'old_stuff.csv'
    old_stuff_csv.write_text('out of date')
    dataset1 = Mock(name='quality_csv')
    dataset2 = Mock(name='fastq1')
    dataset3 = Mock(name='fastq2')
    mock_session.add_dataset.side_effect = [dataset1, dataset2, dataset3]
    mock_pipeline = mock_session.get_pipeline.return_value
    mock_input = Mock(dataset_name='quality_csv')
    mock_pipeline.inputs = [mock_input]
    kive_watcher = KiveWatcher(default_config)

    kive_watcher.add_sample_group(
        base_calls=base_calls,
        sample_group=SampleGroup('2110A',
                                 ('2110A-V3LOOP_S13_L001_R1_001.fastq.gz',
                                  None)))

    mock_open_kive.assert_called_once_with(default_config.kive_server)
    mock_session.login.assert_called_once_with(default_config.kive_user,
                                               default_config.kive_password)
    mock_session.create_run_batch.assert_called_once_with(
        '140101_M01234 v0-dev',
        description='MiCall batch for folder 140101_M01234, pipeline version 0-dev.',
        users=[],
        groups=['Everyone'])
    mock_session.get_pipeline.assert_called_once_with(
        default_config.micall_filter_quality_pipeline_id)
    assert [call(cdt=mock_input.compounddatatype,
                 name='140101_M01234_quality.csv',
                 uploaded=True,
                 # MD5 of header with no records.
                 md5='6861a4a0bfd71b62c0048ff9a4910223'),
            call(cdt=None,
                 name='2110A-V3LOOP_S13_L001_R1_001.fastq.gz',
                 uploaded=True,
                 md5=ANY),
            call(cdt=None,
                 name='2110A-V3LOOP_S13_L001_R2_001.fastq.gz',
                 uploaded=True,
                 md5=ANY)] == mock_session.find_datasets.call_args_list
    assert [call(name='140101_M01234_quality.csv',
                 description='Error rates for 140101_M01234 run.',
                 handle=ANY,
                 cdt=mock_input.compounddatatype,
                 groups=['Everyone']),
            call(name='2110A-V3LOOP_S13_L001_R1_001.fastq.gz',
                 description='forward read from MiSeq run 140101_M01234',
                 handle=ANY,
                 cdt=None,
                 groups=['Everyone']),
            call(name='2110A-V3LOOP_S13_L001_R2_001.fastq.gz',
                 description='reverse read from MiSeq run 140101_M01234',
                 handle=ANY,
                 cdt=None,
                 groups=['Everyone'])] == mock_session.add_dataset.call_args_list
    assert 1 == len(kive_watcher.folder_watchers)
    folder_watcher = kive_watcher.folder_watchers[base_calls]
    assert dataset1 is folder_watcher.quality_dataset
    assert [dataset2, dataset3] == folder_watcher.sample_watchers[0].fastq_datasets
    assert not old_stuff_csv.exists()


def test_create_batch_with_expired_session(raw_data_with_two_samples,
                                           mock_open_kive,
                                           default_config):
    base_calls = (raw_data_with_two_samples /
                  "MiSeq/runs/140101_M01234/Data/Intensities/BaseCalls")
    mock_session = mock_open_kive.return_value
    mock_batch = Mock(name='batch')
    mock_session.create_run_batch.side_effect = [KiveClientException('expired'),
                                                 mock_batch]
    kive_watcher = KiveWatcher(default_config)

    kive_watcher.add_sample_group(
        base_calls=base_calls,
        sample_group=SampleGroup('2110A',
                                 ('2110A-V3LOOP_S13_L001_R1_001.fastq.gz',
                                  None)))

    folder_watcher, = kive_watcher.folder_watchers.values()
    assert mock_batch is folder_watcher.batch


def test_add_external_dataset(raw_data_with_two_samples, mock_open_kive, default_config):
    base_calls = (raw_data_with_two_samples /
                  "MiSeq/runs/140101_M01234/Data/Intensities/BaseCalls")
    default_config.raw_data = raw_data_with_two_samples
    mock_session = mock_open_kive.return_value
    mock_session.get.return_value.json.return_value = [
        dict(name='raw_data', path=str(raw_data_with_two_samples))]
    dataset1 = Mock(name='quality_csv')
    dataset2 = Mock(name='fastq1')
    dataset3 = Mock(name='fastq2')
    mock_session.add_dataset.side_effect = [dataset1, dataset2, dataset3]
    mock_pipeline = mock_session.get_pipeline.return_value
    mock_input = Mock(dataset_name='quality_csv')
    mock_pipeline.inputs = [mock_input]
    kive_watcher = KiveWatcher(default_config)

    kive_watcher.add_sample_group(
        base_calls=base_calls,
        sample_group=SampleGroup('2110A',
                                 ('2110A-V3LOOP_S13_L001_R1_001.fastq.gz',
                                  None)))

    mock_session.get.assert_called_once_with('/api/externalfiledirectories',
                                             is_json=True)
    assert [call(cdt=mock_input.compounddatatype,
                 name='140101_M01234_quality.csv',
                 uploaded=True,
                 # MD5 of header with no records.
                 md5='6861a4a0bfd71b62c0048ff9a4910223'),
            call(cdt=None,
                 name='2110A-V3LOOP_S13_L001_R1_001.fastq.gz',
                 uploaded=True,
                 md5=ANY),
            call(cdt=None,
                 name='2110A-V3LOOP_S13_L001_R2_001.fastq.gz',
                 uploaded=True,
                 md5=ANY)] == mock_session.find_datasets.call_args_list
    assert [call(name='140101_M01234_quality.csv',
                 description='Error rates for 140101_M01234 run.',
                 handle=ANY,
                 cdt=mock_input.compounddatatype,
                 groups=['Everyone']),
            call(name='2110A-V3LOOP_S13_L001_R1_001.fastq.gz',
                 description='forward read from MiSeq run 140101_M01234',
                 handle=None,
                 externalfiledirectory='raw_data',
                 external_path='MiSeq/runs/140101_M01234/Data/Intensities/'
                               'BaseCalls/2110A-V3LOOP_S13_L001_R1_001.fastq.gz',
                 cdt=None,
                 groups=['Everyone']),
            call(name='2110A-V3LOOP_S13_L001_R2_001.fastq.gz',
                 description='reverse read from MiSeq run 140101_M01234',
                 handle=None,
                 externalfiledirectory='raw_data',
                 external_path='MiSeq/runs/140101_M01234/Data/Intensities/'
                               'BaseCalls/2110A-V3LOOP_S13_L001_R2_001.fastq.gz',
                 cdt=None,
                 groups=['Everyone'])] == mock_session.add_dataset.call_args_list


def test_poll_first_sample(raw_data_with_two_samples, mock_open_kive, default_config):
    base_calls = (raw_data_with_two_samples /
                  "MiSeq/runs/140101_M01234/Data/Intensities/BaseCalls")
    mock_session = mock_open_kive.return_value
    kive_watcher = KiveWatcher(default_config)

    kive_watcher.add_sample_group(
        base_calls=base_calls,
        sample_group=SampleGroup('2110A',
                                 ('2110A-V3LOOP_S13_L001_R1_001.fastq.gz',
                                  None)))
    kive_watcher.poll_runs()

    mock_session.run_pipeline.assert_called_once_with(
        mock_session.get_pipeline.return_value,
        [mock_session.add_dataset.return_value],
        'MiCall filter quality on 140101_M01234',
        runbatch=mock_session.create_run_batch.return_value,
        groups=['Everyone'])


def test_poll_first_sample_twice(raw_data_with_two_samples, mock_open_kive, default_config):
    base_calls = (raw_data_with_two_samples /
                  "MiSeq/runs/140101_M01234/Data/Intensities/BaseCalls")
    mock_session = mock_open_kive.return_value
    mock_run = mock_session.run_pipeline.return_value
    mock_run.is_complete.return_value = False
    kive_watcher = KiveWatcher(default_config)

    kive_watcher.add_sample_group(
        base_calls=base_calls,
        sample_group=SampleGroup('2110A',
                                 ('2110A-V3LOOP_S13_L001_R1_001.fastq.gz',
                                  None)))
    kive_watcher.poll_runs()
    kive_watcher.poll_runs()

    mock_session.run_pipeline.assert_called_once()


def test_poll_first_sample_already_running(raw_data_with_two_samples,
                                           mock_open_kive,
                                           default_config):
    base_calls = (raw_data_with_two_samples /
                  "MiSeq/runs/140101_M01234/Data/Intensities/BaseCalls")
    mock_session = mock_open_kive.return_value
    quality_dataset_id = 100
    dataset1 = Mock(name='quality_csv',
                    groups_allowed=ALLOWED_GROUPS,
                    dataset_id=quality_dataset_id)
    dataset1.name = '140101_M01234_quality.csv'
    unfinished_dataset = Mock(name='unfinished_csv', dataset_id=None)
    mock_session.find_datasets.side_effect = [[dataset1], [], []]
    filter_run = MagicMock(
        name='filter_run',
        pipeline_id=default_config.micall_filter_quality_pipeline_id,
        raw=dict(inputs=[dict(index=1, dataset=quality_dataset_id)]),
        **{'get_results.return_value': dict(purged_csv=unfinished_dataset),
           'is_complete.return_value': False})
    mock_session.find_runs.return_value = [filter_run]
    kive_watcher = KiveWatcher(default_config)

    kive_watcher.add_sample_group(
        base_calls=base_calls,
        sample_group=SampleGroup('2110A',
                                 ('2110A-V3LOOP_S13_L001_R1_001.fastq.gz',
                                  None)))
    kive_watcher.poll_runs()

    mock_session.run_pipeline.assert_not_called()
    assert not kive_watcher.other_runs  # Remove the run once it is used.


def test_poll_first_sample_with_other_running(raw_data_with_two_samples,
                                              mock_open_kive,
                                              default_config):
    base_calls = (raw_data_with_two_samples /
                  "MiSeq/runs/140101_M01234/Data/Intensities/BaseCalls")
    mock_session = mock_open_kive.return_value
    quality_dataset_id = 100
    dataset1 = Mock(name='quality_csv',
                    groups_allowed=ALLOWED_GROUPS,
                    dataset_id=quality_dataset_id)
    dataset1.name = '140101_M01234_quality.csv'
    mock_session.find_datasets.side_effect = [[dataset1], [], []]
    other_run = MagicMock(
        name='other_run',
        pipeline_id=default_config.micall_main_pipeline_id,
        raw=dict(inputs=[dict(index=1, dataset=quality_dataset_id)]))
    mock_session.find_runs.return_value = [other_run]
    kive_watcher = KiveWatcher(default_config)

    kive_watcher.add_sample_group(
        base_calls=base_calls,
        sample_group=SampleGroup('2110A',
                                 ('2110A-V3LOOP_S13_L001_R1_001.fastq.gz',
                                  None)))
    kive_watcher.poll_runs()

    mock_session.run_pipeline.assert_called_once()


def test_poll_first_sample_completed_and_purged(raw_data_with_two_samples,
                                                mock_open_kive,
                                                default_config):
    """ A matching run finished recently, but it was purged. """
    base_calls = (raw_data_with_two_samples /
                  "MiSeq/runs/140101_M01234/Data/Intensities/BaseCalls")
    mock_session = mock_open_kive.return_value
    quality_dataset_id = 100
    dataset1 = Mock(name='quality_csv',
                    groups_allowed=ALLOWED_GROUPS,
                    dataset_id=quality_dataset_id)
    dataset1.name = '140101_M01234_quality.csv'
    purged_dataset = Mock(name='purged_csv',
                          dataset_id=None)
    mock_session.find_datasets.side_effect = [[dataset1], [], []]
    filter_run = MagicMock(
        name='filter_run',
        pipeline_id=default_config.micall_filter_quality_pipeline_id,
        raw=dict(inputs=[dict(index=1, dataset=quality_dataset_id)]),
        **{'get_results.return_value': dict(purged_csv=purged_dataset),
           'is_complete.return_value': True})
    mock_session.find_runs.return_value = [filter_run]
    kive_watcher = KiveWatcher(default_config)

    kive_watcher.add_sample_group(
        base_calls=base_calls,
        sample_group=SampleGroup('2110A',
                                 ('2110A-V3LOOP_S13_L001_R1_001.fastq.gz',
                                  None)))
    kive_watcher.poll_runs()

    mock_session.run_pipeline.assert_called_once()


def test_second_sample(raw_data_with_two_samples, mock_open_kive, default_config):
    base_calls = (raw_data_with_two_samples /
                  "MiSeq/runs/140101_M01234/Data/Intensities/BaseCalls")
    mock_session = mock_open_kive.return_value
    kive_watcher = KiveWatcher(default_config)

    kive_watcher.add_sample_group(
        base_calls=base_calls,
        sample_group=SampleGroup('2110A',
                                 ('2110A-V3LOOP_S13_L001_R1_001.fastq.gz',
                                  None)))
    kive_watcher.add_sample_group(
        base_calls=base_calls,
        sample_group=SampleGroup('2120A',
                                 ('2120A-PR_S14_L001_R1_001.fastq.gz',
                                  None)))

    mock_session.create_run_batch.assert_called_once_with(
        '140101_M01234 v0-dev',
        description='MiCall batch for folder 140101_M01234, pipeline version 0-dev.',
        users=[],
        groups=['Everyone'])
    expected_dataset_count = 5  # quality_csv + 2 pairs of FASTQ files
    assert expected_dataset_count == len(
        mock_session.find_datasets.call_args_list)


def test_sample_with_hcv_pair(raw_data_with_hcv_pair, mock_open_kive, default_config):
    base_calls = (raw_data_with_hcv_pair /
                  "MiSeq/runs/140101_M01234/Data/Intensities/BaseCalls")
    mock_session = mock_open_kive.return_value
    kive_watcher = KiveWatcher(default_config)

    kive_watcher.add_sample_group(
        base_calls=base_calls,
        sample_group=SampleGroup('2130A',
                                 ('2130A-HCV_S15_L001_R1_001.fastq.gz',
                                  '2130AMIDI-MidHCV_S16_L001_R1_001.fastq.gz')))

    expected_dataset_count = 5  # quality_csv + 2 pairs of FASTQ files
    assert expected_dataset_count == len(
        mock_session.find_datasets.call_args_list)


def test_sample_fails_to_upload(raw_data_with_two_samples,
                                mock_open_kive,
                                mock_wait,
                                default_config):
    base_calls = (raw_data_with_two_samples /
                  "MiSeq/runs/140101_M01234/Data/Intensities/BaseCalls")
    mock_session = mock_open_kive.return_value
    kive_watcher = KiveWatcher(default_config, retry=True)
    mock_session.add_dataset.side_effect = [ConnectionError('server down'),
                                            Mock(name='quality_csv'),
                                            Mock(name='fastq1'),
                                            Mock(name='fastq2')]

    sample_watcher = kive_watcher.add_sample_group(
        base_calls=base_calls,
        sample_group=SampleGroup('2110A',
                                 ('2110A-V3LOOP_S13_L001_R1_001.fastq.gz',
                                  None)))

    assert sample_watcher is not None
    mock_wait.assert_called_once_with(1)


def test_create_batch_fails(raw_data_with_two_samples,
                            mock_open_kive,
                            mock_wait,
                            default_config):
    base_calls = (raw_data_with_two_samples /
                  "MiSeq/runs/140101_M01234/Data/Intensities/BaseCalls")
    mock_session = mock_open_kive.return_value
    mock_batch = Mock(name='batch')
    mock_session.create_run_batch.side_effect = [ConnectionError('server down'),
                                                 mock_batch]
    kive_watcher = KiveWatcher(default_config, retry=True)

    sample_watcher = kive_watcher.add_sample_group(
        base_calls=base_calls,
        sample_group=SampleGroup('2110A',
                                 ('2110A-V3LOOP_S13_L001_R1_001.fastq.gz',
                                  None)))
    kive_watcher.poll_runs()

    assert sample_watcher is not None
    mock_wait.assert_called_once_with(1)
    mock_session.run_pipeline.assert_called_once_with(
        ANY,
        ANY,
        ANY,
        runbatch=mock_batch,
        groups=ANY)


def test_sample_already_uploaded(raw_data_with_two_samples, mock_open_kive, default_config):
    base_calls = (raw_data_with_two_samples /
                  "MiSeq/runs/140101_M01234/Data/Intensities/BaseCalls")
    mock_session = mock_open_kive.return_value
    dataset1 = Mock(name='quality_csv', groups_allowed=ALLOWED_GROUPS)
    dataset1.name = '140101_M01234_quality.csv'
    dataset2 = Mock(name='fastq1')
    dataset3 = Mock(name='fastq2')
    mock_session.find_datasets.side_effect = [[dataset1], [], []]
    mock_session.add_dataset.side_effect = [dataset2, dataset3]
    mock_pipeline = mock_session.get_pipeline.return_value
    mock_input = Mock(dataset_name='quality_csv')
    mock_pipeline.inputs = [mock_input]
    kive_watcher = KiveWatcher(default_config)

    kive_watcher.add_sample_group(
        base_calls=base_calls,
        sample_group=SampleGroup('2110A',
                                 ('2110A-V3LOOP_S13_L001_R1_001.fastq.gz',
                                  None)))

    assert [call(cdt=mock_input.compounddatatype,
                 name='140101_M01234_quality.csv',
                 uploaded=True,
                 # MD5 of header with no records.
                 md5='6861a4a0bfd71b62c0048ff9a4910223'),
            call(cdt=None,
                 name='2110A-V3LOOP_S13_L001_R1_001.fastq.gz',
                 uploaded=True,
                 md5=ANY),
            call(cdt=None,
                 name='2110A-V3LOOP_S13_L001_R2_001.fastq.gz',
                 uploaded=True,
                 md5=ANY)] == mock_session.find_datasets.call_args_list
    assert [call(name='2110A-V3LOOP_S13_L001_R1_001.fastq.gz',
                 description='forward read from MiSeq run 140101_M01234',
                 handle=ANY,
                 cdt=None,
                 groups=['Everyone']),
            call(name='2110A-V3LOOP_S13_L001_R2_001.fastq.gz',
                 description='reverse read from MiSeq run 140101_M01234',
                 handle=ANY,
                 cdt=None,
                 groups=['Everyone'])] == mock_session.add_dataset.call_args_list
    assert 1 == len(kive_watcher.folder_watchers)
    folder_watcher = kive_watcher.folder_watchers[base_calls]
    assert dataset1 is folder_watcher.quality_dataset
    assert [dataset2, dataset3] == folder_watcher.sample_watchers[0].fastq_datasets


def test_launch_main_run(raw_data_with_two_samples, mock_open_kive, pipelines_config):
    base_calls = (raw_data_with_two_samples /
                  "MiSeq/runs/140101_M01234/Data/Intensities/BaseCalls")
    fastq1 = Mock(name='fastq1')
    fastq2 = Mock(name='fastq2')
    bad_cycles_csv = Mock(name='bad_cycles_csv')
    mock_session = mock_open_kive.return_value
    mock_main_pipeline = Mock(name='main_pipeline')
    mock_session.get_pipeline.return_value = mock_main_pipeline
    mock_main_pipeline.inputs = [Mock(dataset_name='fastq1'),
                                 Mock(dataset_name='fastq1'),
                                 Mock(dataset_name='bad_cycles_csv')]
    mock_session.add_dataset.side_effect = [fastq1, fastq2]
    kive_watcher = KiveWatcher(pipelines_config)

    folder_watcher = kive_watcher.add_folder(base_calls)
    folder_watcher.batch = Mock('batch')
    kive_watcher.add_sample_group(
        base_calls=base_calls,
        sample_group=SampleGroup('2110A',
                                 ('2110A-V3LOOP_S13_L001_R1_001.fastq.gz',
                                  None)))
    folder_watcher.add_run(
        Mock(name='quality_run',
             **{'is_complete.return_value': True,
                'get_results.return_value': dict(bad_cycles_csv=bad_cycles_csv)}),
        PipelineType.FILTER_QUALITY)

    kive_watcher.poll_runs()

    assert [call(pipelines_config.micall_main_pipeline_id)
            ] == mock_session.get_pipeline.call_args_list
    mock_session.run_pipeline.assert_called_once_with(
        mock_main_pipeline,
        [fastq1,
         fastq2,
         bad_cycles_csv],
        'MiCall main on 2110A-V3LOOP_S13',
        runbatch=folder_watcher.batch,
        groups=['Everyone'])


def test_launch_main_run_after_connection_error(raw_data_with_two_samples,
                                                mock_open_kive,
                                                mock_wait,
                                                pipelines_config):
    base_calls = (raw_data_with_two_samples /
                  "MiSeq/runs/140101_M01234/Data/Intensities/BaseCalls")
    fastq1 = Mock(name='fastq1')
    fastq2 = Mock(name='fastq2')
    bad_cycles_csv = Mock(name='bad_cycles_csv')
    mock_session = mock_open_kive.return_value
    mock_main_pipeline = Mock(name='main_pipeline')
    mock_session.get_pipeline.return_value = mock_main_pipeline
    mock_main_pipeline.inputs = [Mock(dataset_name='fastq1'),
                                 Mock(dataset_name='fastq1'),
                                 Mock(dataset_name='bad_cycles_csv')]
    mock_session.add_dataset.side_effect = [fastq1, fastq2]
    kive_watcher = KiveWatcher(pipelines_config, retry=True)

    folder_watcher = kive_watcher.add_folder(base_calls)
    folder_watcher.batch = Mock('batch')
    kive_watcher.add_sample_group(
        base_calls=base_calls,
        sample_group=SampleGroup('2110A',
                                 ('2110A-V3LOOP_S13_L001_R1_001.fastq.gz',
                                  None)))
    folder_watcher.add_run(
        Mock(name='quality_run',
             **{'is_complete.side_effect': [ConnectionError('server down'),
                                            ConnectionError('server down'),
                                            True],
                'get_results.return_value': dict(bad_cycles_csv=bad_cycles_csv)}),
        PipelineType.FILTER_QUALITY)

    kive_watcher.poll_runs()

    mock_session.run_pipeline.assert_called_once_with(
        mock_main_pipeline,
        [fastq1,
         fastq2,
         bad_cycles_csv],
        'MiCall main on 2110A-V3LOOP_S13',
        runbatch=folder_watcher.batch,
        groups=['Everyone'])
    assert [call(1), call(2)] == mock_wait.call_args_list


def test_launch_midi_run(raw_data_with_hcv_pair, mock_open_kive, pipelines_config):
    base_calls = (raw_data_with_hcv_pair /
                  "MiSeq/runs/140101_M01234/Data/Intensities/BaseCalls")
    main_fastq1 = Mock(name='main_fastq1')
    main_fastq2 = Mock(name='main_fastq2')
    midi_fastq1 = Mock(name='midi_fastq1')
    midi_fastq2 = Mock(name='midi_fastq2')
    bad_cycles_csv = Mock(name='bad_cycles_csv')
    mock_session = mock_open_kive.return_value
    mock_main_pipeline = Mock(name='main_pipeline')
    mock_session.get_pipeline.return_value = mock_main_pipeline
    mock_main_pipeline.inputs = [Mock(dataset_name='fastq1'),
                                 Mock(dataset_name='fastq1'),
                                 Mock(dataset_name='bad_cycles_csv')]
    mock_session.add_dataset.side_effect = [main_fastq1,
                                            main_fastq2,
                                            midi_fastq1,
                                            midi_fastq2]
    kive_watcher = KiveWatcher(pipelines_config)

    folder_watcher = kive_watcher.add_folder(base_calls)
    folder_watcher.batch = Mock('batch')
    kive_watcher.add_sample_group(
        base_calls=base_calls,
        sample_group=SampleGroup('2130A',
                                 ('2130A-HCV_S15_L001_R1_001.fastq.gz',
                                  '2130AMIDI-MidHCV_S16_L001_R1_001.fastq.gz')))
    folder_watcher.add_run(
        Mock(name='quality_run',
             **{'is_complete.return_value': True,
                'get_results.return_value': dict(bad_cycles_csv=bad_cycles_csv)}),
        PipelineType.FILTER_QUALITY)

    kive_watcher.poll_runs()

    assert [call(pipelines_config.micall_main_pipeline_id)
            ] == mock_session.get_pipeline.call_args_list
    assert [call(mock_main_pipeline,
                 [main_fastq1, main_fastq2, bad_cycles_csv],
                 'MiCall main on 2130A-HCV_S15',
                 runbatch=folder_watcher.batch,
                 groups=['Everyone']),
            call(mock_main_pipeline,
                 [midi_fastq1, midi_fastq2, bad_cycles_csv],
                 'MiCall main on 2130AMIDI-MidHCV_S16',
                 runbatch=folder_watcher.batch,
                 groups=['Everyone'])
            ] == mock_session.run_pipeline.call_args_list


def test_launch_resistance_run(raw_data_with_two_samples, mock_open_kive, pipelines_config):
    pipelines_config.micall_resistance_pipeline_id = 45
    base_calls = (raw_data_with_two_samples /
                  "MiSeq/runs/140101_M01234/Data/Intensities/BaseCalls")
    amino_csv = Mock(name='amino_csv')
    mock_session = mock_open_kive.return_value
    mock_resistance_pipeline = Mock(name='resistance_pipeline')
    mock_session.get_pipeline.return_value = mock_resistance_pipeline
    mock_resistance_pipeline.inputs = [Mock(dataset_name='main_amino_csv'),
                                       Mock(dataset_name='midi_amino_csv')]
    kive_watcher = KiveWatcher(pipelines_config)

    folder_watcher = kive_watcher.add_folder(base_calls)
    folder_watcher.batch = Mock('batch')
    sample_watcher = kive_watcher.add_sample_group(
        base_calls=base_calls,
        sample_group=SampleGroup('2110A',
                                 ('2110A-V3LOOP_S13_L001_R1_001.fastq.gz',
                                  None)))
    folder_watcher.add_run(Mock(name='filter_quality_run'),
                           PipelineType.FILTER_QUALITY,
                           is_complete=True)
    folder_watcher.add_run(
        Mock(name='main_run',
             **{'is_complete.return_value': True,
                'get_results.return_value': dict(amino_csv=amino_csv)}),
        PipelineType.MAIN,
        sample_watcher)

    kive_watcher.poll_runs()

    assert [call(pipelines_config.micall_resistance_pipeline_id)
            ] == mock_session.get_pipeline.call_args_list
    mock_session.run_pipeline.assert_called_once_with(
        mock_resistance_pipeline,
        [amino_csv, amino_csv],
        'MiCall resistance on 2110A',
        runbatch=folder_watcher.batch,
        groups=['Everyone'])


def test_resistance_run_missing_input(raw_data_with_two_samples,
                                      mock_open_kive,
                                      pipelines_config):
    pipelines_config.micall_resistance_pipeline_id = 45
    base_calls = (raw_data_with_two_samples /
                  "MiSeq/runs/140101_M01234/Data/Intensities/BaseCalls")
    fail_csv = Mock(name='fail_csv')
    mock_session = mock_open_kive.return_value
    mock_resistance_pipeline = Mock(name='resistance_pipeline')
    mock_session.get_pipeline.return_value = mock_resistance_pipeline
    mock_resistance_pipeline.inputs = [Mock(dataset_name='main_amino_csv'),
                                       Mock(dataset_name='midi_amino_csv')]
    kive_watcher = KiveWatcher(pipelines_config)

    folder_watcher = kive_watcher.add_folder(base_calls)
    folder_watcher.batch = Mock('batch')
    sample_watcher = kive_watcher.add_sample_group(
        base_calls=base_calls,
        sample_group=SampleGroup('2110A',
                                 ('2110A-V3LOOP_S13_L001_R1_001.fastq.gz',
                                  None)))
    folder_watcher.add_run(Mock(name='filter_quality_run'),
                           PipelineType.FILTER_QUALITY,
                           is_complete=True)
    folder_watcher.add_run(
        Mock(name='main_run',
             **{'is_complete.return_value': True,
                'get_results.return_value': dict(fail_csv=fail_csv)}),
        PipelineType.MAIN,
        sample_watcher)

    with pytest.raises(RuntimeError, match=r'Polling sample group 2110A failed.'):
        kive_watcher.poll_runs()

    assert [] == mock_session.get_pipeline.call_args_list
    mock_session.run_pipeline.assert_not_called()


def test_poll_main_run_cancelled(raw_data_with_two_samples,
                                 mock_open_kive,
                                 pipelines_config):
    pipelines_config.micall_resistance_pipeline_id = 45
    base_calls = (raw_data_with_two_samples /
                  "MiSeq/runs/140101_M01234/Data/Intensities/BaseCalls")
    mock_session = mock_open_kive.return_value
    mock_main_pipeline = Mock(name='resistance_pipeline')
    mock_main_pipeline.inputs = [Mock(dataset_name='fastq1'),
                                 Mock(dataset_name='fastq2'),
                                 Mock(dataset_name='bad_cycles_csv')]
    mock_session.get_pipeline.return_value = mock_main_pipeline
    original_run = Mock(
        name='original_run',
        **{'is_complete.side_effect': KiveRunFailedException("Run 9 cancelled")})
    new_run = Mock(name='new_run')
    mock_session.run_pipeline.return_value = new_run
    kive_watcher = KiveWatcher(pipelines_config)

    folder_watcher = kive_watcher.add_folder(base_calls)
    folder_watcher.batch = Mock('batch')
    sample_watcher = kive_watcher.add_sample_group(
        base_calls=base_calls,
        sample_group=SampleGroup('2110A',
                                 ('2110A-V3LOOP_S13_L001_R1_001.fastq.gz',
                                  None)))
    folder_watcher.add_run(
        Mock(name='filter_quality_run',
             **{'get_results.return_value': dict(
                 bad_cycles_csv=Mock(name='bad_cycles_csv'))}),
        PipelineType.FILTER_QUALITY,
        is_complete=True)
    folder_watcher.add_run(
        original_run,
        PipelineType.MAIN,
        sample_watcher)

    kive_watcher.poll_runs()

    mock_session.run_pipeline.assert_called_once()
    assert new_run in folder_watcher.active_runs
    assert original_run not in folder_watcher.active_runs


def test_launch_hcv_resistance_run(raw_data_with_hcv_pair, mock_open_kive, pipelines_config):
    pipelines_config.micall_resistance_pipeline_id = 45
    base_calls = (raw_data_with_hcv_pair /
                  "MiSeq/runs/140101_M01234/Data/Intensities/BaseCalls")
    main_amino_csv = Mock(name='main_amino_csv')
    midi_amino_csv = Mock(name='midi_amino_csv')
    mock_session = mock_open_kive.return_value
    mock_resistance_pipeline = Mock(name='main_pipeline')
    mock_session.get_pipeline.return_value = mock_resistance_pipeline
    mock_resistance_pipeline.inputs = [Mock(dataset_name='main_amino_csv'),
                                       Mock(dataset_name='midi_amino_csv')]
    kive_watcher = KiveWatcher(pipelines_config)

    folder_watcher = kive_watcher.add_folder(base_calls)
    folder_watcher.batch = Mock('batch')
    sample_watcher = kive_watcher.add_sample_group(
        base_calls=base_calls,
        sample_group=SampleGroup('2130A',
                                 ('2130A-HCV_S15_L001_R1_001.fastq.gz',
                                  '2130AMIDI-MidHCV_S16_L001_R1_001.fastq.gz')))
    folder_watcher.add_run(Mock(name='filter_quality_run'),
                           PipelineType.FILTER_QUALITY,
                           is_complete=True)
    folder_watcher.add_run(
        Mock(name='main_run',
             **{'is_complete.return_value': True,
                'get_results.return_value': dict(amino_csv=main_amino_csv)}),
        PipelineType.MAIN,
        sample_watcher)
    folder_watcher.add_run(
        Mock(name='midi_run',
             **{'is_complete.return_value': True,
                'get_results.return_value': dict(amino_csv=midi_amino_csv)}),
        PipelineType.MIDI,
        sample_watcher)

    kive_watcher.poll_runs()

    assert [call(pipelines_config.micall_resistance_pipeline_id)
            ] == mock_session.get_pipeline.call_args_list
    mock_session.run_pipeline.assert_called_once_with(
        mock_resistance_pipeline,
        [main_amino_csv, midi_amino_csv],
        'MiCall resistance on 2130A',
        runbatch=folder_watcher.batch,
        groups=['Everyone'])


def test_launch_mixed_hcv_run(raw_data_with_hcv_pair, mock_open_kive, pipelines_config):
    pipelines_config.mixed_hcv_pipeline_id = 47
    base_calls = (raw_data_with_hcv_pair /
                  "MiSeq/runs/140101_M01234/Data/Intensities/BaseCalls")
    main_fastq1 = Mock(name='main_fastq1')
    main_fastq2 = Mock(name='main_fastq2')
    midi_fastq1 = Mock(name='midi_fastq1')
    midi_fastq2 = Mock(name='midi_fastq2')
    bad_cycles_csv = Mock(name='bad_cycles_csv')
    mock_session = mock_open_kive.return_value
    mock_mixed_hcv_pipeline = Mock(name='mixed_hcv_pipeline')
    mock_main_pipeline = Mock(name='main_pipeline')
    mock_session.get_pipeline.side_effect = [mock_mixed_hcv_pipeline,
                                             mock_main_pipeline]
    mock_session.add_dataset.side_effect = [main_fastq1,
                                            main_fastq2,
                                            midi_fastq1,
                                            midi_fastq2]
    mock_main_pipeline.inputs = [Mock(dataset_name='fastq1'),
                                 Mock(dataset_name='fastq1'),
                                 Mock(dataset_name='bad_cycles_csv')]
    mock_mixed_hcv_pipeline.inputs = [Mock(dataset_name='FASTQ1'),
                                      Mock(dataset_name='FASTQ2')]
    kive_watcher = KiveWatcher(pipelines_config)

    folder_watcher = kive_watcher.add_folder(base_calls)
    folder_watcher.batch = Mock('batch')
    kive_watcher.add_sample_group(
        base_calls=base_calls,
        sample_group=SampleGroup('2130A',
                                 ('2130A-HCV_S15_L001_R1_001.fastq.gz',
                                  '2130AMIDI-MidHCV_S16_L001_R1_001.fastq.gz')))
    folder_watcher.add_run(
        Mock(name='quality_run',
             **{'is_complete.return_value': True,
                'get_results.return_value': dict(bad_cycles_csv=bad_cycles_csv)}),
        PipelineType.FILTER_QUALITY)

    kive_watcher.poll_runs()

    assert [call(pipelines_config.mixed_hcv_pipeline_id),
            call(pipelines_config.micall_main_pipeline_id)
            ] == mock_session.get_pipeline.call_args_list
    assert [call(mock_mixed_hcv_pipeline,
                 [main_fastq1, main_fastq2],
                 'Mixed HCV on 2130A-HCV_S15',
                 runbatch=folder_watcher.batch,
                 groups=['Everyone']),
            call(mock_mixed_hcv_pipeline,
                 [midi_fastq1, midi_fastq2],
                 'Mixed HCV on 2130AMIDI-MidHCV_S16',
                 runbatch=folder_watcher.batch,
                 groups=['Everyone']),
            call(mock_main_pipeline,
                 [main_fastq1, main_fastq2, bad_cycles_csv],
                 'MiCall main on 2130A-HCV_S15',
                 runbatch=folder_watcher.batch,
                 groups=['Everyone']),
            call(mock_main_pipeline,
                 [midi_fastq1, midi_fastq2, bad_cycles_csv],
                 'MiCall main on 2130AMIDI-MidHCV_S16',
                 runbatch=folder_watcher.batch,
                 groups=['Everyone'])
            ] == mock_session.run_pipeline.call_args_list


def test_full_with_two_samples(raw_data_with_two_samples, mock_open_kive, pipelines_config):
    assert mock_open_kive
    pipelines_config.max_active = 2
    base_calls = (raw_data_with_two_samples /
                  "MiSeq/runs/140101_M01234/Data/Intensities/BaseCalls")
    kive_watcher = KiveWatcher(pipelines_config)

    kive_watcher.add_sample_group(
        base_calls=base_calls,
        sample_group=SampleGroup('2110A',
                                 ('2110A-V3LOOP_S13_L001_R1_001.fastq.gz',
                                  None)))
    is_full1 = kive_watcher.is_full()
    kive_watcher.add_sample_group(
        base_calls=base_calls,
        sample_group=SampleGroup('2120A',
                                 ('2120A-PR_S14_L001_R1_001.fastq.gz',
                                  None)))
    is_full2 = kive_watcher.is_full()

    folder_watcher, = kive_watcher.folder_watchers.values()
    folder_watcher.completed_samples.add('2120A-PR_S14_L001_R1_001.fastq.gz')

    is_full3 = kive_watcher.is_full()

    assert not is_full1
    assert is_full2
    assert not is_full3


def test_full_with_two_runs(raw_data_with_two_runs, mock_open_kive, pipelines_config):
    assert mock_open_kive
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
                                  None)))
    is_full1 = kive_watcher.is_full()
    kive_watcher.add_sample_group(
        base_calls=base_calls2,
        sample_group=SampleGroup('2010A',
                                 ('2010A-V3LOOP_S3_L001_R1_001.fastq.gz',
                                  None)))
    is_full2 = kive_watcher.is_full()

    folder_watcher = kive_watcher.folder_watchers[base_calls1]
    folder_watcher.completed_samples.add('2000A-V3LOOP_S2_L001_R1_001.fastq.gz')

    is_full3 = kive_watcher.is_full()

    assert not is_full1
    assert is_full2
    assert not is_full3


def test_fetch_run_status_incomplete(mock_open_kive, pipelines_config):
    assert mock_open_kive
    mock_run = Mock(name='run',
                    **{'is_complete.return_value': False})

    kive_watcher = KiveWatcher(pipelines_config)

    new_run = kive_watcher.fetch_run_status(mock_run,
                                            folder_watcher=None,
                                            pipeline_type=None,
                                            sample_watcher=None)

    assert new_run is mock_run


def test_fetch_run_status_filter_quality(raw_data_with_two_runs,
                                         mock_open_kive,
                                         pipelines_config):
    assert mock_open_kive
    base_calls = (raw_data_with_two_runs /
                  "MiSeq/runs/140101_M01234/Data/Intensities/BaseCalls")
    folder_watcher = FolderWatcher(base_calls)
    sample_watcher = None
    mock_run = Mock(name='run',
                    **{'is_complete.return_value': True})

    kive_watcher = KiveWatcher(pipelines_config)

    new_run = kive_watcher.fetch_run_status(mock_run,
                                            folder_watcher,
                                            PipelineType.FILTER_QUALITY,
                                            sample_watcher)

    assert new_run is None


def test_fetch_run_status_main(raw_data_with_two_runs,
                               mock_open_kive,
                               pipelines_config):
    assert mock_open_kive
    base_calls = (raw_data_with_two_runs /
                  "MiSeq/runs/140101_M01234/Data/Intensities/BaseCalls")
    folder_watcher = FolderWatcher(base_calls)
    sample_watcher = SampleWatcher(
        SampleGroup('2000A',
                    ('2000A-V3LOOP_S2_L001_R1_001.fastq.gz',
                     None)))
    mock_run = Mock(**{
        'is_complete.return_value': True,
        'get_results.return_value': create_datasets(['coord_ins_csv',
                                                     'nuc_csv'])})
    expected_scratch = base_calls / "../../../Results/version_0-dev/scratch"
    expected_coord_ins_path = expected_scratch / "2000A-V3LOOP_S2/coord_ins.csv"
    expected_nuc_path = expected_scratch / "2000A-V3LOOP_S2/nuc.csv"
    expected_nuc_content = """\
row,name
0,nuc_csv
1,nuc_csv
2,nuc_csv
"""

    kive_watcher = KiveWatcher(pipelines_config)

    new_run = kive_watcher.fetch_run_status(mock_run,
                                            folder_watcher,
                                            PipelineType.MAIN,
                                            sample_watcher)

    assert new_run is None
    assert expected_coord_ins_path.exists()
    assert expected_nuc_content == expected_nuc_path.read_text()


def test_fetch_run_status_main_and_resistance(raw_data_with_two_runs,
                                              mock_open_kive,
                                              pipelines_config):
    assert mock_open_kive
    base_calls = (raw_data_with_two_runs /
                  "MiSeq/runs/140101_M01234/Data/Intensities/BaseCalls")
    folder_watcher = FolderWatcher(base_calls)
    sample_watcher = SampleWatcher(
        SampleGroup('2000A',
                    ('2000A-V3LOOP_S2_L001_R1_001.fastq.gz',
                     None)))
    mock_run = Mock(**{
        'is_complete.return_value': True,
        'get_results.side_effect': [create_datasets(['nuc_csv']),
                                    create_datasets(['resistance_csv'])]})
    expected_scratch = base_calls / "../../../Results/version_0-dev/scratch"
    expected_nuc_path = expected_scratch / "2000A-V3LOOP_S2/nuc.csv"
    expected_resistance_path = expected_scratch / "2000A-V3LOOP_S2/resistance.csv"

    kive_watcher = KiveWatcher(pipelines_config)

    new_main_run = kive_watcher.fetch_run_status(
        mock_run,
        folder_watcher,
        PipelineType.MAIN,
        sample_watcher)
    new_resistance_run = kive_watcher.fetch_run_status(
        mock_run,
        folder_watcher,
        PipelineType.RESISTANCE,
        sample_watcher)

    assert new_main_run is None
    assert new_resistance_run is None
    assert expected_nuc_path.exists()
    assert expected_resistance_path.exists()


def test_fetch_run_status_main_and_midi(raw_data_with_hcv_pair,
                                        mock_open_kive,
                                        pipelines_config):
    assert mock_open_kive
    base_calls = (raw_data_with_hcv_pair /
                  "MiSeq/runs/140101_M01234/Data/Intensities/BaseCalls")
    folder_watcher = FolderWatcher(base_calls)
    sample_watcher = SampleWatcher(
        SampleGroup('2130A', ('2130A-HCV_S15_L001_R1_001.fastq.gz',
                              '2130AMIDI-MidHCV_S16_L001_R1_001.fastq.gz')))
    mock_run = Mock(**{'is_complete.return_value': True,
                       'get_results.return_value': create_datasets(['nuc_csv'])})
    expected_scratch = base_calls / "../../../Results/version_0-dev/scratch"
    expected_main_nuc_path = expected_scratch / "2130A-HCV_S15/nuc.csv"
    expected_midi_nuc_path = expected_scratch / "2130AMIDI-MidHCV_S16/nuc.csv"

    kive_watcher = KiveWatcher(pipelines_config)

    new_main_run = kive_watcher.fetch_run_status(mock_run,
                                                 folder_watcher,
                                                 PipelineType.MAIN,
                                                 sample_watcher)
    new_midi_run = kive_watcher.fetch_run_status(mock_run,
                                                 folder_watcher,
                                                 PipelineType.MIDI,
                                                 sample_watcher)

    assert new_main_run is None
    assert new_midi_run is None
    assert expected_main_nuc_path.exists()
    assert expected_midi_nuc_path.exists()


def test_fetch_run_status_session_expired(raw_data_with_two_runs,
                                          mock_open_kive,
                                          pipelines_config):
    assert mock_open_kive
    base_calls = (raw_data_with_two_runs /
                  "MiSeq/runs/140101_M01234/Data/Intensities/BaseCalls")
    mock_run = Mock(**{'is_complete.side_effect': [KiveClientException('expired'),
                                                   True],
                       'get_results.return_value': {}})

    kive_watcher = KiveWatcher(pipelines_config)

    sample_watcher = kive_watcher.add_sample_group(
        base_calls,
        SampleGroup('2000A', ('2000A-V3LOOP_S2_L001_R1_001.fastq.gz', None)))
    folder_watcher, = kive_watcher.folder_watchers.values()

    new_run = kive_watcher.fetch_run_status(mock_run,
                                            folder_watcher,
                                            PipelineType.MAIN,
                                            sample_watcher)

    assert new_run is None


def test_fetch_run_status_user_cancelled(raw_data_with_two_runs,
                                         mock_open_kive,
                                         pipelines_config):
    assert mock_open_kive
    base_calls = (raw_data_with_two_runs /
                  "MiSeq/runs/140101_M01234/Data/Intensities/BaseCalls")
    original_run = Mock(
        name='original_run',
        **{'is_complete.side_effect': KiveRunFailedException("Run 9 cancelled")})

    kive_watcher = KiveWatcher(pipelines_config)

    sample_watcher = kive_watcher.add_sample_group(
        base_calls,
        SampleGroup('2000A', ('2000A-V3LOOP_S2_L001_R1_001.fastq.gz', None)))
    folder_watcher, = kive_watcher.folder_watchers.values()
    folder_watcher.bad_cycles_dataset = Mock(name='bad_cycles_csv')

    new_run = kive_watcher.fetch_run_status(original_run,
                                            folder_watcher,
                                            PipelineType.MAIN,
                                            sample_watcher)

    assert new_run is not None
    assert new_run is not original_run


def test_folder_completed(raw_data_with_two_samples, mock_open_kive, default_config):
    base_calls = (raw_data_with_two_samples /
                  "MiSeq/runs/140101_M01234/Data/Intensities/BaseCalls")
    assert mock_open_kive
    resistance_run1 = Mock(name='resistance_run1',
                           **{'is_complete.return_value': True,
                              'get_results.return_value': create_datasets(['resistance_csv'])})
    resistance_run2 = Mock(name='resistance_run2',
                           **{'is_complete.return_value': True,
                              'get_results.return_value': create_datasets(['resistance_csv'])})
    kive_watcher = KiveWatcher(default_config)

    folder_watcher = kive_watcher.add_folder(base_calls)
    sample1_watcher = kive_watcher.add_sample_group(
        base_calls=base_calls,
        sample_group=SampleGroup('2110A',
                                 ('2110A-V3LOOP_S13_L001_R1_001.fastq.gz',
                                  None)))
    sample2_watcher = kive_watcher.add_sample_group(
        base_calls=base_calls,
        sample_group=SampleGroup('2120A',
                                 ('2120A-PR_S14_L001_R1_001.fastq.gz',
                                  None)))
    kive_watcher.finish_folder(base_calls)
    folder_watcher.add_run(Mock(name='filter_quality_run'),
                           PipelineType.FILTER_QUALITY,
                           is_complete=True)
    folder_watcher.add_run(Mock(name='main_run1'),
                           PipelineType.MAIN,
                           sample1_watcher,
                           is_complete=True)
    folder_watcher.add_run(Mock(name='main_run2'),
                           PipelineType.MAIN,
                           sample2_watcher,
                           is_complete=True)
    folder_watcher.add_run(resistance_run1,
                           PipelineType.RESISTANCE,
                           sample1_watcher)
    folder_watcher.add_run(resistance_run2,
                           PipelineType.RESISTANCE,
                           sample2_watcher)
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
    expected_done_path = results_path / "doneprocessing"
    expected_mutations_path = results_path / "mutations.csv"
    expected_resistance_path = results_path / "resistance.csv"
    expected_resistance_content = """\
sample,row,name
2110A-V3LOOP_S13,0,resistance_csv
2110A-V3LOOP_S13,1,resistance_csv
2110A-V3LOOP_S13,2,resistance_csv
2120A-PR_S14,0,resistance_csv
2120A-PR_S14,1,resistance_csv
2120A-PR_S14,2,resistance_csv
"""

    kive_watcher.poll_runs()

    assert not scratch_path.exists()
    assert not expected_mutations_path.exists()
    assert expected_resistance_content == expected_resistance_path.read_text()
    assert expected_coverage_map_content == expected_coverage_map_path.read_bytes()
    assert expected_done_path.exists()


def test_folder_not_finished(raw_data_with_two_samples, mock_open_kive, default_config):
    assert mock_open_kive
    base_calls = (raw_data_with_two_samples /
                  "MiSeq/runs/140101_M01234/Data/Intensities/BaseCalls")
    resistance_run1 = Mock(name='resistance_run1',
                           **{'is_complete.return_value': True,
                              'get_results.return_value': create_datasets(['resistance_csv'])})
    resistance_run2 = Mock(name='resistance_run2',
                           **{'is_complete.return_value': True,
                              'get_results.return_value': create_datasets(['resistance_csv'])})
    kive_watcher = KiveWatcher(default_config)

    folder_watcher = kive_watcher.add_folder(base_calls)
    sample1_watcher = kive_watcher.add_sample_group(
        base_calls=base_calls,
        sample_group=SampleGroup('2110A',
                                 ('2110A-V3LOOP_S13_L001_R1_001.fastq.gz',
                                  None)))
    sample2_watcher = kive_watcher.add_sample_group(
        base_calls=base_calls,
        sample_group=SampleGroup('2120A',
                                 ('2120A-PR_S14_L001_R1_001.fastq.gz',
                                  None)))
    # Did not call kive_watcher.finish_folder(), more samples could be coming.
    folder_watcher.add_run(Mock(name='filter_quality_run'),
                           PipelineType.FILTER_QUALITY,
                           is_complete=True)
    folder_watcher.add_run(Mock(name='main_run1'),
                           PipelineType.MAIN,
                           sample1_watcher,
                           is_complete=True)
    folder_watcher.add_run(Mock(name='main_run2'),
                           PipelineType.MAIN,
                           sample2_watcher,
                           is_complete=True)
    folder_watcher.add_run(resistance_run1,
                           PipelineType.RESISTANCE,
                           sample1_watcher)
    folder_watcher.add_run(resistance_run2,
                           PipelineType.RESISTANCE,
                           sample2_watcher)
    results_path = base_calls / "../../../Results/version_0-dev"
    scratch_path = results_path / "scratch"
    expected_resistance_path = results_path / "resistance.csv"

    kive_watcher.poll_runs()

    assert scratch_path.exists()
    assert not expected_resistance_path.exists()


def test_folder_not_finished_before_new_start(raw_data_with_two_runs,
                                              mock_open_kive,
                                              default_config):
    assert mock_open_kive
    base_calls1 = (raw_data_with_two_runs /
                   "MiSeq/runs/140101_M01234/Data/Intensities/BaseCalls")
    base_calls2 = (raw_data_with_two_runs /
                   "MiSeq/runs/140201_M01234/Data/Intensities/BaseCalls")
    resistance_run = Mock(name='resistance_run',
                          **{'is_complete.return_value': True,
                             'get_results.return_value': create_datasets(['resistance_csv'])})
    kive_watcher = KiveWatcher(default_config)

    folder_watcher1 = kive_watcher.add_folder(base_calls1)
    sample1_watcher = kive_watcher.add_sample_group(
        base_calls=base_calls1,
        sample_group=SampleGroup('2000A',
                                 ('2000A-V3LOOP_S2_L001_R1_001.fastq.gz',
                                  None)))
    # Did not call kive_watcher1.finish_folder(), more samples could be coming.
    folder_watcher2 = kive_watcher.add_folder(base_calls2)
    kive_watcher.add_sample_group(
        base_calls=base_calls2,
        sample_group=SampleGroup('2010A',
                                 ('2010A-V3LOOP_S3_L001_R1_001.fastq.gz',
                                  None)))
    folder_watcher1.add_run(Mock(name='filter_quality_run'),
                            PipelineType.FILTER_QUALITY,
                            is_complete=True)
    folder_watcher1.add_run(Mock(name='main_run'),
                            PipelineType.MAIN,
                            sample1_watcher,
                            is_complete=True)
    folder_watcher1.add_run(resistance_run,
                            PipelineType.RESISTANCE,
                            sample1_watcher)
    folder_watcher2.quality_dataset = Mock(name='quality2_csv')
    results_path = base_calls1 / "../../../Results/version_0-dev"
    scratch_path = results_path / "scratch"
    expected_resistance_path = results_path / "resistance.csv"

    kive_watcher.poll_runs()

    assert not expected_resistance_path.exists()
    assert scratch_path.exists()


def test_folder_failed(raw_data_with_two_samples, mock_open_kive, default_config):
    assert mock_open_kive
    base_calls = (raw_data_with_two_samples /
                  "MiSeq/runs/140101_M01234/Data/Intensities/BaseCalls")
    main_run1 = Mock(
        name='main_run1',
        **{'is_complete.side_effect': KiveRunFailedException('failed')})
    resistance_run2 = Mock(
        name='resistance_run2',
        **{'is_complete.return_value': True,
           'get_results.return_value': create_datasets(['resistance_csv'])})
    kive_watcher = KiveWatcher(default_config)

    folder_watcher = kive_watcher.add_folder(base_calls)
    sample1_watcher = kive_watcher.add_sample_group(
        base_calls=base_calls,
        sample_group=SampleGroup('2110A',
                                 ('2110A-V3LOOP_S13_L001_R1_001.fastq.gz',
                                  None)))
    sample2_watcher = kive_watcher.add_sample_group(
        base_calls=base_calls,
        sample_group=SampleGroup('2120A',
                                 ('2120A-PR_S14_L001_R1_001.fastq.gz',
                                  None)))
    kive_watcher.finish_folder(base_calls)
    folder_watcher.add_run(Mock(name='filter_quality_run'),
                           PipelineType.FILTER_QUALITY,
                           is_complete=True)
    folder_watcher.add_run(main_run1,
                           PipelineType.MAIN,
                           sample1_watcher)
    folder_watcher.add_run(Mock(name='main_run2'),
                           PipelineType.MAIN,
                           sample2_watcher,
                           is_complete=True)
    folder_watcher.add_run(resistance_run2,
                           PipelineType.RESISTANCE,
                           sample2_watcher)
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
                                  None)))
    sample_watcher2 = kive_watcher.add_sample_group(
        base_calls=base_calls,
        sample_group=SampleGroup('2110A',
                                 ('2110A-V3LOOP_S13_L001_R1_001.fastq.gz',
                                  None)))

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
                                  None)))

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
                                  None)))

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
