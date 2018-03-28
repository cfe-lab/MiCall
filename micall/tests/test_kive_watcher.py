from gzip import GzipFile
from pathlib import Path
from queue import Full
from unittest.mock import patch, ANY, Mock, call

import pytest
from datetime import datetime, timedelta

from struct import pack

from micall.monitor.kive_watcher import find_samples, KiveWatcher
from micall.resistance.resistance import SampleGroup
from micall_watcher import parse_args


@pytest.fixture(name='mock_clock')
def create_mock_clock():
    with patch('micall.monitor.kive_watcher.now') as mock_clock:
        mock_clock.return_value = datetime(2000, 1, 1)
        yield mock_clock


@pytest.fixture(name='mock_open_kive')
def create_mock_open_kive():
    with patch('micall.monitor.kive_watcher.open_kive') as mock_open_kive:
        yield mock_open_kive


@pytest.fixture(name='default_config')
def create_default_config():
    default_config = parse_args(argv=['--micall_filter_quality_pipeline_id', '42'])
    yield default_config


@pytest.fixture(name='pipelines_config')
def create_pipelines_config(default_config):
    default_config.micall_main_pipeline_id = 43
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


def create_run(tmpdir, run_name, sample_pattern):
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
    (run / "qc_uploaded").touch()
    return raw_data_folder


@pytest.fixture(name='raw_data_with_hcv_pair')
def create_raw_data_with_hcv_pair(tmpdir):
    return create_run(tmpdir, '140101_M01234', '2130A*.fastq')


@pytest.fixture(name='raw_data_with_two_samples')
def create_raw_data_with_two_samples(tmpdir):
    # Install samples 2110A and 2120A in a single run.
    return create_run(tmpdir, '140101_M01234', '21[12]0A*.fastq')


@pytest.fixture(name='raw_data_with_two_runs')
def create_raw_data_with_two_runs(tmpdir):
    raw_data = create_run(tmpdir, '140101_M01234', '2000A*.fastq')
    create_run(tmpdir, '140201_M01234', '2010A*.fastq')

    return raw_data


def test_hcv_pair(raw_data_with_hcv_pair):
    sample_queue = DummyQueueSink()
    sample_queue.expect_put(
        (raw_data_with_hcv_pair / "MiSeq/runs/140101_M01234" /
         "Data/Intensities/BaseCalls",
         SampleGroup('2130A', ('2130A-HCV_S15_L001_R1_001.fastq.gz',
                               '2130AMIDI-MidHCV_S16_L001_R1_001.fastq.gz'))))

    find_samples(raw_data_with_hcv_pair, sample_queue, wait=False)

    sample_queue.verify()


def test_two_runs(raw_data_with_two_runs):
    sample_queue = DummyQueueSink()
    sample_queue.expect_put(
        (raw_data_with_two_runs / "MiSeq/runs/140201_M01234" /
         "Data/Intensities/BaseCalls",
         SampleGroup('2010A', ('2010A-V3LOOP_S3_L001_R1_001.fastq.gz',
                               None))))
    sample_queue.expect_put(
        (raw_data_with_two_runs / "MiSeq/runs/140101_M01234" /
         "Data/Intensities/BaseCalls",
         SampleGroup('2000A', ('2000A-V3LOOP_S2_L001_R1_001.fastq.gz',
                               None))))

    find_samples(raw_data_with_two_runs, sample_queue, wait=False)

    sample_queue.verify()


def test_full_queue(raw_data_with_two_runs):
    sample_queue = DummyQueueSink()
    sample_queue.expect_put(
        (raw_data_with_two_runs / "MiSeq/runs/140201_M01234" /
         "Data/Intensities/BaseCalls",
         SampleGroup('2010A', ('2010A-V3LOOP_S3_L001_R1_001.fastq.gz',
                               None))))
    item2 = (raw_data_with_two_runs / "MiSeq/runs/140101_M01234" /
             "Data/Intensities/BaseCalls",
             SampleGroup('2000A', ('2000A-V3LOOP_S2_L001_R1_001.fastq.gz',
                                   None)))
    sample_queue.expect_put(item2, is_full=True)
    sample_queue.expect_put(item2)

    find_samples(raw_data_with_two_runs, sample_queue, wait=False)

    sample_queue.verify()


def test_scan_for_new_runs(raw_data_with_two_runs, mock_clock):
    """ Mark the more recent run as not ready, timeout, then notice new run. """
    def increment_clock():
        mock_clock.return_value += timedelta(hours=1)

    needs_processing2 = (Path(raw_data_with_two_runs) /
                         "MiSeq/runs/140201_M01234/needsprocessing")
    needs_processing2.unlink()
    sample_queue = DummyQueueSink()
    item1 = (raw_data_with_two_runs / "MiSeq/runs/140101_M01234" /
             "Data/Intensities/BaseCalls",
             SampleGroup('2000A', ('2000A-V3LOOP_S2_L001_R1_001.fastq.gz',
                                   None)))
    item2 = (raw_data_with_two_runs / "MiSeq/runs/140201_M01234" /
             "Data/Intensities/BaseCalls",
             SampleGroup('2010A', ('2010A-V3LOOP_S3_L001_R1_001.fastq.gz',
                                   None)))
    sample_queue.expect_put(item1,
                            is_full=True,
                            callback=needs_processing2.touch)
    sample_queue.expect_put(item1,
                            is_full=True,
                            callback=increment_clock)
    sample_queue.expect_put(item2)
    sample_queue.expect_put(item1)

    find_samples(raw_data_with_two_runs, sample_queue, wait=False)

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


def test_get_kive_input(pipelines_config):
    kive_watcher = KiveWatcher(pipelines_config)
    mock_pipeline = Mock()
    kive_watcher.pipelines[
        pipelines_config.micall_main_pipeline_id] = mock_pipeline
    expected_input = Mock(dataset_name='bad_cycles_csv')
    mock_pipeline.inputs = [Mock(), expected_input]

    kive_input = kive_watcher.get_kive_input('bad_cycles_csv')

    assert expected_input is kive_input


def test_get_kive_input_wrong_pipeline(pipelines_config):
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


def test_first_sample(raw_data_with_two_samples, mock_open_kive, default_config):
    base_calls = (raw_data_with_two_samples /
                  "MiSeq/runs/140101_M01234/Data/Intensities/BaseCalls")
    mock_session = mock_open_kive.return_value
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
        description='MiCall batch for folder 140101_M01234, pipeline version v0-dev.',
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
    mock_session.run_pipeline.assert_called_once_with(
        mock_session.get_pipeline.return_value,
        [mock_session.add_dataset.return_value],
        'MiCall filter quality on 140101_M01234',
        runbatch=mock_session.create_run_batch.return_value,
        groups=['Everyone'])


def test_second_sample(raw_data_with_two_samples, mock_open_kive, default_config):
    base_calls = (raw_data_with_two_samples /
                  "MiSeq/runs/140101_M01234/Data/Intensities/BaseCalls")
    mock_session = mock_open_kive.return_value
    mock_session.get_pipeline.return_value.inputs = [
        Mock(dataset_name='quality_csv')]
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
        description='MiCall batch for folder 140101_M01234, pipeline version v0-dev.',
        users=[],
        groups=['Everyone'])
    mock_session.run_pipeline.assert_called_once_with(
        mock_session.get_pipeline.return_value,
        [mock_session.add_dataset.return_value],
        'MiCall filter quality on 140101_M01234',
        runbatch=mock_session.create_run_batch.return_value,
        groups=['Everyone'])
    expected_dataset_count = 5  # quality_csv + 2 pairs of FASTQ files
    assert expected_dataset_count == len(
        mock_session.find_datasets.call_args_list)
