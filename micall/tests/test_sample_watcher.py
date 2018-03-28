from pathlib import Path

from micall.monitor.sample_watcher import FolderWatcher


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
