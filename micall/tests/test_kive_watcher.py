from gzip import GzipFile
from pathlib import Path
from queue import Full
from unittest.mock import patch, ANY, Mock, call

import pytest
from datetime import datetime, timedelta

from struct import pack

from micall.monitor.kive_watcher import find_samples, KiveWatcher
from micall.monitor.sample_watcher import PipelineType
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
    (run / "qc_uploaded").touch()
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
    assert 1 == len(kive_watcher.folder_watchers)
    folder_watcher = kive_watcher.folder_watchers[base_calls]
    assert dataset1 is folder_watcher.quality_dataset
    assert [dataset2, dataset3] == folder_watcher.sample_watchers[0].fastq_datasets


def test_poll_first_sample(raw_data_with_two_samples, mock_open_kive, default_config):
    base_calls = (raw_data_with_two_samples /
                  "MiSeq/runs/140101_M01234/Data/Intensities/BaseCalls")
    mock_session = mock_open_kive.return_value
    mock_pipeline = mock_session.get_pipeline.return_value
    mock_input = Mock(dataset_name='quality_csv')
    mock_pipeline.inputs = [mock_input]
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
    mock_run.is_complete.assert_called_once_with()


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
    expected_dataset_count = 5  # quality_csv + 2 pairs of FASTQ files
    assert expected_dataset_count == len(
        mock_session.find_datasets.call_args_list)


def test_sample_with_hcv_pair(raw_data_with_hcv_pair, mock_open_kive, default_config):
    base_calls = (raw_data_with_hcv_pair /
                  "MiSeq/runs/140101_M01234/Data/Intensities/BaseCalls")
    mock_session = mock_open_kive.return_value
    mock_session.get_pipeline.return_value.inputs = [
        Mock(dataset_name='quality_csv')]
    kive_watcher = KiveWatcher(default_config)

    kive_watcher.add_sample_group(
        base_calls=base_calls,
        sample_group=SampleGroup('2130A',
                                 ('2130A-HCV_S15_L001_R1_001.fastq.gz',
                                  '2130AMIDI-MidHCV_S16_L001_R1_001.fastq.gz')))

    expected_dataset_count = 5  # quality_csv + 2 pairs of FASTQ files
    assert expected_dataset_count == len(
        mock_session.find_datasets.call_args_list)


def test_launch_main_run(raw_data_with_two_samples, mock_open_kive, pipelines_config):
    base_calls = (raw_data_with_two_samples /
                  "MiSeq/runs/140101_M01234/Data/Intensities/BaseCalls")
    quality_csv = Mock(name='quality_csv')
    fastq1 = Mock(name='fastq1')
    fastq2 = Mock(name='fastq2')
    bad_cycles_csv = Mock(name='bad_cycles_csv')
    mock_session = mock_open_kive.return_value
    mock_quality_pipeline = Mock(name='quality_pipeline')
    mock_main_pipeline = Mock(name='main_pipeline')
    mock_session.get_pipeline.side_effect = [mock_quality_pipeline, mock_main_pipeline]
    mock_input = Mock(dataset_name='quality_csv')
    mock_quality_pipeline.inputs = [mock_input]
    mock_main_pipeline.inputs = [Mock(dataset_name='fastq1'),
                                 Mock(dataset_name='fastq1'),
                                 Mock(dataset_name='bad_cycles_csv')]
    mock_session.add_dataset.side_effect = [quality_csv, fastq1, fastq2]
    kive_watcher = KiveWatcher(pipelines_config)

    kive_watcher.add_sample_group(
        base_calls=base_calls,
        sample_group=SampleGroup('2110A',
                                 ('2110A-V3LOOP_S13_L001_R1_001.fastq.gz',
                                  None)))
    folder_watcher, = kive_watcher.folder_watchers.values()
    filter_quality_run = Mock(
        name='quality_run',
        **{'get_results.return_value': dict(bad_cycles_csv=bad_cycles_csv)})
    folder_watcher.filter_quality_run = filter_quality_run
    sample_watcher, = folder_watcher.sample_watchers

    run = kive_watcher.run_pipeline(folder_watcher, sample_watcher, PipelineType.MAIN)

    assert [call(pipelines_config.micall_filter_quality_pipeline_id),
            call(pipelines_config.micall_main_pipeline_id)
            ] == mock_session.get_pipeline.call_args_list
    mock_session.run_pipeline.assert_called_once_with(
        mock_main_pipeline,
        [fastq1,
         fastq2,
         bad_cycles_csv],
        'MiCall main on 2110A-V3LOOP_S13',
        runbatch=mock_session.create_run_batch.return_value,
        groups=['Everyone'])
    assert mock_session.run_pipeline.return_value is run
    assert bad_cycles_csv is folder_watcher.bad_cycles_dataset


def test_launch_midi_run(raw_data_with_hcv_pair, mock_open_kive, pipelines_config):
    base_calls = (raw_data_with_hcv_pair /
                  "MiSeq/runs/140101_M01234/Data/Intensities/BaseCalls")
    quality_csv = Mock(name='quality_csv')
    main_fastq1 = Mock(name='main_fastq1')
    main_fastq2 = Mock(name='main_fastq2')
    midi_fastq1 = Mock(name='midi_fastq1')
    midi_fastq2 = Mock(name='midi_fastq2')
    bad_cycles_csv = Mock(name='bad_cycles_csv')
    mock_session = mock_open_kive.return_value
    mock_quality_pipeline = Mock(name='quality_pipeline')
    mock_main_pipeline = Mock(name='main_pipeline')
    mock_session.get_pipeline.side_effect = [mock_quality_pipeline, mock_main_pipeline]
    mock_input = Mock(dataset_name='quality_csv')
    mock_quality_pipeline.inputs = [mock_input]
    mock_session.add_dataset.side_effect = [quality_csv,
                                            main_fastq1,
                                            main_fastq2,
                                            midi_fastq1,
                                            midi_fastq2]
    kive_watcher = KiveWatcher(pipelines_config)

    kive_watcher.add_sample_group(
        base_calls=base_calls,
        sample_group=SampleGroup('2130A',
                                 ('2130A-HCV_S15_L001_R1_001.fastq.gz',
                                  '2130AMIDI-MidHCV_S16_L001_R1_001.fastq.gz')))
    folder_watcher, = kive_watcher.folder_watchers.values()
    folder_watcher.bad_cycles_dataset = bad_cycles_csv
    sample_watcher, = folder_watcher.sample_watchers

    run = kive_watcher.run_pipeline(folder_watcher, sample_watcher, PipelineType.MIDI)

    assert [call(pipelines_config.micall_filter_quality_pipeline_id),
            call(pipelines_config.micall_main_pipeline_id)
            ] == mock_session.get_pipeline.call_args_list
    mock_session.run_pipeline.assert_called_once_with(
        mock_main_pipeline,
        [midi_fastq1,
         midi_fastq2,
         bad_cycles_csv],
        'MiCall main on 2130AMIDI-MidHCV_S16',
        runbatch=mock_session.create_run_batch.return_value,
        groups=['Everyone'])
    assert mock_session.run_pipeline.return_value is run


def test_launch_resistance_run(raw_data_with_two_samples, mock_open_kive, pipelines_config):
    pipelines_config.micall_resistance_pipeline_id = 45
    base_calls = (raw_data_with_two_samples /
                  "MiSeq/runs/140101_M01234/Data/Intensities/BaseCalls")
    amino_csv = Mock(name='amino_csv')
    mock_session = mock_open_kive.return_value
    mock_quality_pipeline = Mock(name='quality_pipeline')
    mock_resistance_pipeline = Mock(name='resistance_pipeline')
    mock_session.get_pipeline.side_effect = [mock_quality_pipeline, mock_resistance_pipeline]
    mock_input = Mock(dataset_name='quality_csv')
    mock_quality_pipeline.inputs = [mock_input]
    mock_resistance_pipeline.inputs = [Mock(dataset_name='main_amino_csv'),
                                       Mock(dataset_name='midi_amino_csv')]
    kive_watcher = KiveWatcher(pipelines_config)

    kive_watcher.add_sample_group(
        base_calls=base_calls,
        sample_group=SampleGroup('2110A',
                                 ('2110A-V3LOOP_S13_L001_R1_001.fastq.gz',
                                  None)))
    folder_watcher, = kive_watcher.folder_watchers.values()
    sample_watcher, = folder_watcher.sample_watchers
    main_run = Mock(
        name='main_run',
        **{'get_results.return_value': dict(amino_csv=amino_csv)})
    sample_watcher.main_runs.append(main_run)

    run = kive_watcher.run_pipeline(folder_watcher,
                                    sample_watcher,
                                    PipelineType.RESISTANCE)

    assert [call(pipelines_config.micall_filter_quality_pipeline_id),
            call(pipelines_config.micall_resistance_pipeline_id)
            ] == mock_session.get_pipeline.call_args_list
    mock_session.run_pipeline.assert_called_once_with(
        mock_resistance_pipeline,
        [amino_csv, amino_csv],
        'MiCall resistance on 2110A',
        runbatch=mock_session.create_run_batch.return_value,
        groups=['Everyone'])
    assert mock_session.run_pipeline.return_value is run


def test_launch_hcv_resistance_run(raw_data_with_hcv_pair, mock_open_kive, pipelines_config):
    pipelines_config.micall_resistance_pipeline_id = 45
    base_calls = (raw_data_with_hcv_pair /
                  "MiSeq/runs/140101_M01234/Data/Intensities/BaseCalls")
    main_amino_csv = Mock(name='main_amino_csv')
    midi_amino_csv = Mock(name='midi_amino_csv')
    mock_session = mock_open_kive.return_value
    mock_quality_pipeline = Mock(name='quality_pipeline')
    mock_resistance_pipeline = Mock(name='main_pipeline')
    mock_session.get_pipeline.side_effect = [mock_quality_pipeline,
                                             mock_resistance_pipeline]
    mock_input = Mock(dataset_name='quality_csv')
    mock_quality_pipeline.inputs = [mock_input]
    mock_resistance_pipeline.inputs = [Mock(dataset_name='main_amino_csv'),
                                       Mock(dataset_name='midi_amino_csv')]
    kive_watcher = KiveWatcher(pipelines_config)

    kive_watcher.add_sample_group(
        base_calls=base_calls,
        sample_group=SampleGroup('2130A',
                                 ('2130A-HCV_S15_L001_R1_001.fastq.gz',
                                  '2130AMIDI-MidHCV_S16_L001_R1_001.fastq.gz')))
    folder_watcher, = kive_watcher.folder_watchers.values()
    sample_watcher, = folder_watcher.sample_watchers
    main_run = Mock(
        name='main_run',
        **{'get_results.return_value': dict(amino_csv=main_amino_csv)})
    midi_run = Mock(
        name='midi_run',
        **{'get_results.return_value': dict(amino_csv=midi_amino_csv)})
    sample_watcher.main_runs.extend((main_run, midi_run))

    run = kive_watcher.run_pipeline(folder_watcher,
                                    sample_watcher,
                                    PipelineType.RESISTANCE)

    assert [call(pipelines_config.micall_filter_quality_pipeline_id),
            call(pipelines_config.micall_resistance_pipeline_id)
            ] == mock_session.get_pipeline.call_args_list
    mock_session.run_pipeline.assert_called_once_with(
        mock_resistance_pipeline,
        [main_amino_csv, midi_amino_csv],
        'MiCall resistance on 2130A',
        runbatch=mock_session.create_run_batch.return_value,
        groups=['Everyone'])
    assert mock_session.run_pipeline.return_value is run


def test_full_with_two_samples(raw_data_with_two_samples, mock_open_kive, pipelines_config):
    pipelines_config.max_active = 2
    base_calls = (raw_data_with_two_samples /
                  "MiSeq/runs/140101_M01234/Data/Intensities/BaseCalls")
    mock_session = mock_open_kive.return_value
    mock_quality_pipeline = Mock(name='quality_pipeline')
    mock_resistance_pipeline = Mock(name='resistance_pipeline')
    mock_session.get_pipeline.side_effect = [mock_quality_pipeline, mock_resistance_pipeline]
    mock_input = Mock(dataset_name='quality_csv')
    mock_quality_pipeline.inputs = [mock_input]
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
    pipelines_config.max_active = 2
    base_calls1 = (raw_data_with_two_runs /
                   "MiSeq/runs/140101_M01234/Data/Intensities/BaseCalls")
    base_calls2 = (raw_data_with_two_runs /
                   "MiSeq/runs/140201_M01234/Data/Intensities/BaseCalls")
    mock_session = mock_open_kive.return_value
    mock_quality_pipeline = Mock(name='quality_pipeline')
    mock_resistance_pipeline = Mock(name='resistance_pipeline')
    mock_session.get_pipeline.side_effect = [mock_quality_pipeline, mock_resistance_pipeline]
    mock_input = Mock(dataset_name='quality_csv')
    mock_quality_pipeline.inputs = [mock_input]
    kive_watcher = KiveWatcher(pipelines_config)

    # raw_data = create_run_folder(tmpdir, '140101_M01234', '2000A*.fastq')
    # create_run_folder(tmpdir, '140201_M01234', '2010A*.fastq')
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
