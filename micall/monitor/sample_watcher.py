from enum import Enum
from pathlib import Path

ALLOWED_GROUPS = ['Everyone']
# noinspection PyArgumentList
PipelineType = Enum('PipelineType', 'FILTER_QUALITY MAIN MIDI RESISTANCE')


class FolderWatcher:
    def __init__(self,
                 base_calls_folder,
                 runner=None):
        """ Set up an instance.

        :param base_calls_folder: path to the BaseCalls folder under a MiSeq
            run folder
        :param runner: an object for running Kive pipelines. Must have these
            methods:
            run_pipeline(folder_watcher, sample_watcher, pipeline_type) => run
            fetch_run_status(run) => True if successfully finished, raise if
                run failed
        """
        self.base_calls_folder = Path(base_calls_folder)
        self.runner = runner
        self.run_folder = (self.base_calls_folder / '../../..').resolve()
        self.run_name = '_'.join(self.run_folder.name.split('_')[:2])
        self.sample_watchers = []
        self.batch = None
        self.quality_dataset = None
        self.filter_quality_run = None
        self.bad_cycles_dataset = None
        self.active_runs = set()
        self.completed_samples = set()  # {fastq1_name}

    def __repr__(self):
        return f'FolderWatcher({str(self.base_calls_folder)!r})'

    @property
    def active_samples(self):
        started_samples = {name
                           for sample_watcher in self.sample_watchers
                           for name in sample_watcher.sample_group.names
                           if name is not None}
        return started_samples - self.completed_samples

    @property
    def is_complete(self):
        return not self.active_samples

    def poll_runs(self):
        if self.filter_quality_run is None:
            self.filter_quality_run = self.runner.run_pipeline(
                self,
                None,
                PipelineType.FILTER_QUALITY)
            self.active_runs.add(self.filter_quality_run)
            # Launched filter run, nothing more to check.
            return
        if self.filter_quality_run in self.active_runs:
            if not self.runner.fetch_run_status(self.filter_quality_run):
                # Still running, nothing more to check.
                return
            self.active_runs.remove(self.filter_quality_run)
        for sample_watcher in self.sample_watchers:
            is_finished = self.poll_sample_runs(sample_watcher)
            if is_finished:
                self.completed_samples.update(
                    name
                    for name in sample_watcher.sample_group.names
                    if name is not None)

    def poll_sample_runs(self, sample_watcher):
        if not sample_watcher.main_runs:
            sample_watcher.main_runs.append(self.runner.run_pipeline(
                self,
                sample_watcher,
                PipelineType.MAIN))
            if sample_watcher.sample_group.names[1] is not None:
                sample_watcher.main_runs.append(self.runner.run_pipeline(
                    self,
                    sample_watcher,
                    PipelineType.MIDI))
            self.active_runs.update(sample_watcher.main_runs)
            # Launched main and midi runs, nothing more to check on sample.
            return False
        for run in sample_watcher.main_runs:
            if run in self.active_runs:
                if not self.runner.fetch_run_status(run):
                    # Still running, nothing more to check on sample
                    return False
                self.active_runs.remove(run)
        if sample_watcher.resistance_run is None:
            sample_watcher.resistance_run = self.runner.run_pipeline(
                self,
                sample_watcher,
                PipelineType.RESISTANCE)
            # Launched resistance run, nothing more to check on sample.
            return False
        return self.runner.fetch_run_status(sample_watcher.resistance_run)


class SampleWatcher:
    def __init__(self, sample_group):
        self.sample_group = sample_group
        self.fastq_datasets = []
        self.main_runs = []
        self.resistance_run = None

    def __repr__(self):
        enum = self.sample_group.enum
        midi_name = self.sample_group.names[1] and '...'
        return f"SampleWatcher(SampleGroup({enum!r}, ('...', {midi_name!r})))"
