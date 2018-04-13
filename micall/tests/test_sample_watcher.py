from pathlib import Path

from micall.monitor.sample_watcher import FolderWatcher, SampleWatcher, PipelineType
from micall.resistance.resistance import SampleGroup


class DummySession:
    def __init__(self, skipped_types=frozenset()):
        """ Initialize.

        :param set skipped_types: {PipelineType} types that won't run
        """
        self.skipped_types = skipped_types
        self.active_runs = []  # [(folder_watcher, sample_watcher, pipeline_type)]

    def run_pipeline(self, folder_watcher, pipeline_type, sample_watcher):
        if pipeline_type in self.skipped_types:
            return None
        run = (folder_watcher, sample_watcher, pipeline_type)
        self.active_runs.append(run)
        return run

    def fetch_run_status(self,
                         run,
                         _folder_watcher,
                         _pipeline_type,
                         _sample_watcher):
        if run not in self.active_runs:
            return None
        return run

    def finish_run(self, run):
        self.active_runs.remove(run)

    def finish_all_runs(self):
        self.active_runs.clear()


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
    folder_watcher.poll_runs()  # filter_quality still running

    assert [(folder_watcher, None, PipelineType.FILTER_QUALITY)] == session.active_runs


def test_filter_quality_finished():
    base_calls_folder = '/path/Data/Intensities/BaseCalls'
    session = DummySession()
    folder_watcher = FolderWatcher(base_calls_folder, runner=session)
    sample_watcher = SampleWatcher(SampleGroup('1234A', ('1234A-V3LOOP_R1_001.fastq.gz', None)))
    folder_watcher.sample_watchers.append(sample_watcher)

    folder_watcher.poll_runs()   # Start filter_quality
    session.finish_all_runs()  # Finish filter_quality

    folder_watcher.poll_runs()   # start main

    assert [(folder_watcher, sample_watcher, PipelineType.MAIN)] == session.active_runs


def test_main_running():
    base_calls_folder = '/path/Data/Intensities/BaseCalls'
    session = DummySession()
    folder_watcher = FolderWatcher(base_calls_folder, runner=session)
    sample_watcher = SampleWatcher(SampleGroup('1234A', ('1234A-V3LOOP_R1_001.fastq.gz', None)))
    folder_watcher.sample_watchers.append(sample_watcher)

    folder_watcher.poll_runs()   # Start filter_quality
    session.finish_all_runs()  # Finish filter_quality

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
    session.finish_all_runs()  # Finish filter_quality
    folder_watcher.poll_runs()   # Start main
    session.finish_all_runs()  # Finish main
    folder_watcher.poll_runs()   # Start resistance

    assert [(folder_watcher, sample_watcher, PipelineType.RESISTANCE)] == session.active_runs
    assert 1 == len(folder_watcher.active_runs)


def test_resistance_running():
    base_calls_folder = '/path/Data/Intensities/BaseCalls'
    session = DummySession()
    folder_watcher = FolderWatcher(base_calls_folder, runner=session)
    sample_watcher = SampleWatcher(SampleGroup('1234A', ('1234A-V3LOOP_R1_001.fastq.gz', None)))
    folder_watcher.sample_watchers.append(sample_watcher)

    folder_watcher.poll_runs()   # Start filter_quality
    session.finish_all_runs()  # Finish filter_quality
    folder_watcher.poll_runs()   # Start main
    session.finish_all_runs()  # Finish main
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
    session.finish_all_runs()  # Finish filter_quality
    folder_watcher.poll_runs()   # Start main
    session.finish_all_runs()  # Finish main
    folder_watcher.poll_runs()   # Start resistance
    session.finish_all_runs()  # Finish resistance
    folder_watcher.poll_runs()   # Finish sample

    assert not session.active_runs
    assert not folder_watcher.active_samples
    assert not folder_watcher.active_runs
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
    session.finish_all_runs()  # Finish filter_quality

    folder_watcher.poll_runs()   # start main, midi, and mixed HCV

    assert [(folder_watcher, sample_watcher, PipelineType.MIXED_HCV_MAIN),
            (folder_watcher, sample_watcher, PipelineType.MIXED_HCV_MIDI),
            (folder_watcher, sample_watcher, PipelineType.MAIN),
            (folder_watcher, sample_watcher, PipelineType.MIDI)
            ] == session.active_runs
    assert 4 == len(folder_watcher.active_runs)


def test_hcv_mixed_hcv_running():
    base_calls_folder = '/path/Data/Intensities/BaseCalls'
    session = DummySession()
    folder_watcher = FolderWatcher(base_calls_folder, runner=session)
    sample_watcher = SampleWatcher(SampleGroup(
        '2130A',
        ('2130A-HCV_S15_L001_R1_001.fastq.gz',
         '2130AMIDI-MidHCV_S16_L001_R1_001.fastq.gz')))
    folder_watcher.sample_watchers.append(sample_watcher)

    folder_watcher.poll_runs()   # Start filter_quality
    session.finish_all_runs()  # Finish filter_quality

    folder_watcher.poll_runs()   # start main, midi, and mixed HCV
    folder_watcher.poll_runs()   # main, midi, and mixed HCV still running

    assert [(folder_watcher, sample_watcher, PipelineType.MIXED_HCV_MAIN),
            (folder_watcher, sample_watcher, PipelineType.MIXED_HCV_MIDI),
            (folder_watcher, sample_watcher, PipelineType.MAIN),
            (folder_watcher, sample_watcher, PipelineType.MIDI)
            ] == session.active_runs
    assert 4 == len(folder_watcher.active_runs)


def test_hcv_mixed_hcv_finished():
    base_calls_folder = '/path/Data/Intensities/BaseCalls'
    session = DummySession()
    folder_watcher = FolderWatcher(base_calls_folder, runner=session)
    sample_watcher = SampleWatcher(SampleGroup(
        '2130A',
        ('2130A-HCV_S15_L001_R1_001.fastq.gz',
         '2130AMIDI-MidHCV_S16_L001_R1_001.fastq.gz')))
    folder_watcher.sample_watchers.append(sample_watcher)

    folder_watcher.poll_runs()   # Start filter_quality
    session.finish_all_runs()  # Finish filter_quality

    folder_watcher.poll_runs()   # start main, midi, and mixed HCV
    # Finish mixed HCV
    session.finish_run(
        (folder_watcher, sample_watcher, PipelineType.MIXED_HCV_MAIN))
    session.finish_run(
        (folder_watcher, sample_watcher, PipelineType.MIXED_HCV_MIDI))
    folder_watcher.poll_runs()   # main and midi still running

    assert [(folder_watcher, sample_watcher, PipelineType.MAIN),
            (folder_watcher, sample_watcher, PipelineType.MIDI)
            ] == session.active_runs
    assert 2 == len(folder_watcher.active_runs)


def test_hcv_mixed_hcv_not_finished():
    base_calls_folder = '/path/Data/Intensities/BaseCalls'
    session = DummySession()
    folder_watcher = FolderWatcher(base_calls_folder, runner=session)
    sample_watcher = SampleWatcher(SampleGroup(
        '2130A',
        ('2130A-HCV_S15_L001_R1_001.fastq.gz',
         '2130AMIDI-MidHCV_S16_L001_R1_001.fastq.gz')))
    folder_watcher.sample_watchers.append(sample_watcher)

    folder_watcher.poll_runs()   # Start filter_quality
    session.finish_all_runs()  # Finish filter_quality

    folder_watcher.poll_runs()   # start main, midi, and mixed HCV
    session.finish_run(
        (folder_watcher, sample_watcher, PipelineType.MAIN))  # Finish main
    session.finish_run(
        (folder_watcher, sample_watcher, PipelineType.MIDI))  # Finish midi
    folder_watcher.poll_runs()   # mixed HCV still running, resistance started
    session.finish_run(
        (folder_watcher, sample_watcher, PipelineType.RESISTANCE))  # Finish res
    folder_watcher.poll_runs()   # mixed HCV still running, resistance finished

    assert [(folder_watcher, sample_watcher, PipelineType.MIXED_HCV_MAIN),
            (folder_watcher, sample_watcher, PipelineType.MIXED_HCV_MIDI)
            ] == session.active_runs
    assert 2 == len(folder_watcher.active_runs)
    assert not folder_watcher.is_complete


def test_mixed_hcv_skipped():
    base_calls_folder = '/path/Data/Intensities/BaseCalls'
    session = DummySession(skipped_types={PipelineType.MIXED_HCV_MAIN,
                                          PipelineType.MIXED_HCV_MIDI})
    folder_watcher = FolderWatcher(base_calls_folder, runner=session)
    sample_watcher = SampleWatcher(SampleGroup(
        '2130A',
        ('2130A-HCV_S15_L001_R1_001.fastq.gz',
         '2130AMIDI-MidHCV_S16_L001_R1_001.fastq.gz')))
    folder_watcher.sample_watchers.append(sample_watcher)

    folder_watcher.poll_runs()   # Start filter_quality
    session.finish_all_runs()  # Finish filter_quality

    folder_watcher.poll_runs()   # start main and midi

    assert [(folder_watcher, sample_watcher, PipelineType.MAIN),
            (folder_watcher, sample_watcher, PipelineType.MIDI)
            ] == session.active_runs
    assert 2 == len(folder_watcher.active_runs)


def test_mixed_hcv_skipped_and_complete():
    base_calls_folder = '/path/Data/Intensities/BaseCalls'
    session = DummySession(skipped_types={PipelineType.MIXED_HCV_MAIN,
                                          PipelineType.MIXED_HCV_MIDI})
    folder_watcher = FolderWatcher(base_calls_folder, runner=session)
    sample_watcher = SampleWatcher(SampleGroup(
        '2130A',
        ('2130A-HCV_S15_L001_R1_001.fastq.gz',
         '2130AMIDI-MidHCV_S16_L001_R1_001.fastq.gz')))
    folder_watcher.sample_watchers.append(sample_watcher)

    folder_watcher.poll_runs()   # Start filter_quality
    session.finish_all_runs()  # Finish filter_quality
    folder_watcher.poll_runs()   # start main and midi
    session.finish_all_runs()  # Finish main and midi
    folder_watcher.poll_runs()   # start resistance
    session.finish_all_runs()  # Finish resistance
    folder_watcher.poll_runs()   # done

    assert not session.active_runs
    assert not folder_watcher.active_runs
    assert folder_watcher.is_complete
