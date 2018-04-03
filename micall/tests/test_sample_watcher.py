from pathlib import Path

from micall.monitor.sample_watcher import FolderWatcher, SampleWatcher, PipelineType
from micall.resistance.resistance import SampleGroup


class DummySession:
    def __init__(self):
        self.active_runs = []

    def run_pipeline(self, *args):
        self.active_runs.append(args)
        return args

    def fetch_run_status(self, run):
        return run not in self.active_runs


def test_folder_watcher_repr():
    base_calls_folder = '/path/Data/Intensities/BaseCalls'
    expected_repr = "FolderWatcher('/path/Data/Intensities/BaseCalls')"
    watcher = FolderWatcher(base_calls_folder)

    assert expected_repr == repr(watcher)


def test_folder_watcher_repr_with_pathlib():
    base_calls_folder = Path('/path/Data/Intensities/BaseCalls')
    expected_repr = "FolderWatcher('/path/Data/Intensities/BaseCalls')"
    watcher = FolderWatcher(base_calls_folder)

    assert expected_repr == repr(watcher)


def test_folder_watcher_run_details():
    base_calls_folder = '/path/140101_M01234_JUNK/Data/Intensities/BaseCalls'
    expected_run_folder = Path('/path/140101_M01234_JUNK')
    expected_run_name = '140101_M01234'
    watcher = FolderWatcher(base_calls_folder)

    assert expected_run_folder == watcher.run_folder
    assert expected_run_name == watcher.run_name


def test_sample_watcher_repr_pair():
    expected_repr = "SampleWatcher(SampleGroup('1234A', ('...', '...')))"
    watcher = SampleWatcher(SampleGroup('1234A', ('foo', 'bar')))

    assert expected_repr == repr(watcher)


def test_sample_watcher_repr_single():
    expected_repr = "SampleWatcher(SampleGroup('1234A', ('...', None)))"
    watcher = SampleWatcher(SampleGroup('1234A', ('foo', None)))

    assert expected_repr == repr(watcher)


def test_launch_filter_quality():
    base_calls_folder = '/path/Data/Intensities/BaseCalls'
    session = DummySession()
    folder_watcher = FolderWatcher(base_calls_folder, runner=session)
    sample_watcher = SampleWatcher(SampleGroup('1234A', ('1234A-V3LOOP_R1_001.fastq.gz', None)))
    folder_watcher.sample_watchers.append(sample_watcher)

    folder_watcher.poll_runs()  # Start filter_quality

    assert [(folder_watcher, None, PipelineType.FILTER_QUALITY)] == session.active_runs
    assert {'1234A-V3LOOP_R1_001.fastq.gz'} == folder_watcher.active_samples
    assert not folder_watcher.is_complete


def test_filter_quality_running():
    base_calls_folder = '/path/Data/Intensities/BaseCalls'
    session = DummySession()
    folder_watcher = FolderWatcher(base_calls_folder, runner=session)
    sample_watcher = SampleWatcher(SampleGroup('1234A', ('1234A-V3LOOP_R1_001.fastq.gz', None)))
    folder_watcher.sample_watchers.append(sample_watcher)

    folder_watcher.poll_runs()  # Start filter_quality
    folder_watcher.poll_runs()  # Start filter_quality

    assert [(folder_watcher, None, PipelineType.FILTER_QUALITY)] == session.active_runs


def test_filter_quality_finished():
    base_calls_folder = '/path/Data/Intensities/BaseCalls'
    session = DummySession()
    folder_watcher = FolderWatcher(base_calls_folder, runner=session)
    sample_watcher = SampleWatcher(SampleGroup('1234A', ('1234A-V3LOOP_R1_001.fastq.gz', None)))
    folder_watcher.sample_watchers.append(sample_watcher)

    folder_watcher.poll_runs()   # Start filter_quality
    session.active_runs.clear()  # Finish filter_quality

    folder_watcher.poll_runs()   # start main

    assert [(folder_watcher, sample_watcher, PipelineType.MAIN)] == session.active_runs


def test_main_running():
    base_calls_folder = '/path/Data/Intensities/BaseCalls'
    session = DummySession()
    folder_watcher = FolderWatcher(base_calls_folder, runner=session)
    sample_watcher = SampleWatcher(SampleGroup('1234A', ('1234A-V3LOOP_R1_001.fastq.gz', None)))
    folder_watcher.sample_watchers.append(sample_watcher)

    folder_watcher.poll_runs()   # Start filter_quality
    session.active_runs.clear()  # Finish filter_quality

    folder_watcher.poll_runs()   # start main
    folder_watcher.poll_runs()   # main still running

    assert [(folder_watcher, sample_watcher, PipelineType.MAIN)] == session.active_runs


def test_main_finished():
    base_calls_folder = '/path/Data/Intensities/BaseCalls'
    session = DummySession()
    folder_watcher = FolderWatcher(base_calls_folder, runner=session)
    sample_watcher = SampleWatcher(SampleGroup('1234A', ('1234A-V3LOOP_R1_001.fastq.gz', None)))
    folder_watcher.sample_watchers.append(sample_watcher)

    folder_watcher.poll_runs()   # Start filter_quality
    session.active_runs.clear()  # Finish filter_quality
    folder_watcher.poll_runs()   # Start main
    session.active_runs.clear()  # Finish main
    folder_watcher.poll_runs()   # Start resistance

    assert [(folder_watcher, sample_watcher, PipelineType.RESISTANCE)] == session.active_runs


def test_resistance_running():
    base_calls_folder = '/path/Data/Intensities/BaseCalls'
    session = DummySession()
    folder_watcher = FolderWatcher(base_calls_folder, runner=session)
    sample_watcher = SampleWatcher(SampleGroup('1234A', ('1234A-V3LOOP_R1_001.fastq.gz', None)))
    folder_watcher.sample_watchers.append(sample_watcher)

    folder_watcher.poll_runs()   # Start filter_quality
    session.active_runs.clear()  # Finish filter_quality
    folder_watcher.poll_runs()   # Start main
    session.active_runs.clear()  # Finish main
    folder_watcher.poll_runs()   # Start resistance
    folder_watcher.poll_runs()   # resistance still running

    assert [(folder_watcher, sample_watcher, PipelineType.RESISTANCE)] == session.active_runs
    assert {'1234A-V3LOOP_R1_001.fastq.gz'} == folder_watcher.active_samples
    assert not folder_watcher.is_complete


def test_resistance_finished():
    base_calls_folder = '/path/Data/Intensities/BaseCalls'
    session = DummySession()
    folder_watcher = FolderWatcher(base_calls_folder, runner=session)
    sample_watcher = SampleWatcher(SampleGroup('1234A', ('1234A-V3LOOP_R1_001.fastq.gz', None)))
    folder_watcher.sample_watchers.append(sample_watcher)

    folder_watcher.poll_runs()   # Start filter_quality
    session.active_runs.clear()  # Finish filter_quality
    folder_watcher.poll_runs()   # Start main
    session.active_runs.clear()  # Finish main
    folder_watcher.poll_runs()   # Start resistance
    session.active_runs.clear()  # Finish resistance
    folder_watcher.poll_runs()   # Finish sample

    assert [] == session.active_runs
    assert set() == folder_watcher.active_samples
    assert folder_watcher.is_complete


def test_hcv_filter_quality_finished():
    base_calls_folder = '/path/Data/Intensities/BaseCalls'
    session = DummySession()
    folder_watcher = FolderWatcher(base_calls_folder, runner=session)
    sample_watcher = SampleWatcher(SampleGroup(
        '2130A',
        ('2130A-HCV_S15_L001_R1_001.fastq.gz',
         '2130AMIDI-MidHCV_S16_L001_R1_001.fastq.gz')))
    folder_watcher.sample_watchers.append(sample_watcher)

    folder_watcher.poll_runs()   # Start filter_quality
    session.active_runs.clear()  # Finish filter_quality

    folder_watcher.poll_runs()   # start main

    assert [(folder_watcher, sample_watcher, PipelineType.MAIN),
            (folder_watcher, sample_watcher, PipelineType.MIDI)] == session.active_runs
