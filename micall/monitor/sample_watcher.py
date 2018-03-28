from pathlib import Path

ALLOWED_GROUPS = ['Everyone']


class FolderWatcher:
    def __init__(self,
                 base_calls_folder,
                 kive_session=None,
                 pipeline_version='unknown'):
        self.base_calls_folder = Path(base_calls_folder)
        self.session = kive_session
        self.pipeline_version = pipeline_version
        self.run_folder = (self.base_calls_folder / '../../..').resolve()
        self.run_name = '_'.join(self.run_folder.name.split('_')[:2])
        self.sample_groups = []
        self.samples = {}  # {file_name: SampleWatcher}
        self.batch = None
        self.filter_quality_run = None
        self.filter_quality_finished = False

    def __repr__(self):
        return f'FolderWatcher({str(self.base_calls_folder)!r})'


class SampleWatcher:
    pass
