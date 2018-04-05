from enum import Enum
from pathlib import Path

ALLOWED_GROUPS = ['Everyone']
# noinspection PyArgumentList
PipelineType = Enum(
    'PipelineType',
    'FILTER_QUALITY MAIN MIDI RESISTANCE MIXED_HCV_MAIN MIXED_HCV_MIDI')


class FolderWatcher:
    def __init__(self,
                 base_calls_folder,
                 runner=None):
        """ Set up an instance.

        :param base_calls_folder: path to the BaseCalls folder under a MiSeq
            run folder
        :param runner: an object for running Kive pipelines. Must have these
            methods:
            run_pipeline(folder_watcher, sample_watcher, pipeline_type)
                returns run, or None if that pipeline_type is not configured.
            fetch_run_status(
                run,
                folder_watcher,
                sample_watcher,
                pipeline_type) => True if successfully finished, raise if
                run failed, also saves the outputs to temporary files in the
                results folder when the run is finished
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
        self.active_runs = {}  # {run: (sample_watcher, pipeline_type)}
        self.completed_samples = set()  # {fastq1_name}

    def __repr__(self):
        return f'FolderWatcher({str(self.base_calls_folder)!r})'

    @property
    def all_samples(self):
        for sample_watcher in self.sample_watchers:
            for name in sample_watcher.sample_group.names:
                if name is not None:
                    yield name

    @property
    def active_samples(self):
        all_samples = set(self.all_samples)
        return all_samples - self.completed_samples

    @property
    def is_complete(self):
        return not self.active_samples

    def poll_runs(self):
        if self.filter_quality_run is None:
            self.filter_quality_run = self.run_pipeline(
                None,
                PipelineType.FILTER_QUALITY)
            # Launched filter run, nothing more to check.
            return
        if self.filter_quality_run in self.active_runs:
            if not self.fetch_run_status(self.filter_quality_run):
                # Still running, nothing more to check.
                return
        for sample_watcher in self.sample_watchers:
            is_finished = self.poll_sample_runs(sample_watcher)
            if is_finished:
                self.completed_samples.update(
                    name
                    for name in sample_watcher.sample_group.names
                    if name is not None)

    def poll_sample_runs(self, sample_watcher):
        is_mixed_hcv_complete = False
        if not sample_watcher.mixed_hcv_runs:
            if 'HCV' not in sample_watcher.sample_group.names[0]:
                is_mixed_hcv_complete = True
            else:
                mixed_hcv_run = self.run_pipeline(
                    sample_watcher,
                    PipelineType.MIXED_HCV_MAIN)
                if mixed_hcv_run is None:
                    is_mixed_hcv_complete = True
                else:
                    sample_watcher.mixed_hcv_runs.append(mixed_hcv_run)
                    if sample_watcher.sample_group.names[1] is not None:
                        mixed_hcv_run = self.run_pipeline(
                            sample_watcher,
                            PipelineType.MIXED_HCV_MIDI)
                        sample_watcher.mixed_hcv_runs.append(mixed_hcv_run)
        else:
            active_mixed_runs = [
                run
                for run in sample_watcher.mixed_hcv_runs
                if run in self.active_runs and not self.fetch_run_status(run)]
            is_mixed_hcv_complete = not active_mixed_runs

        if not sample_watcher.main_runs:
            sample_watcher.main_runs.append(self.run_pipeline(
                sample_watcher,
                PipelineType.MAIN))
            if sample_watcher.sample_group.names[1] is not None:
                sample_watcher.main_runs.append(self.run_pipeline(
                    sample_watcher,
                    PipelineType.MIDI))
            # Launched main and midi runs, nothing more to check on sample.
            return False
        active_main_runs = [
            run
            for run in sample_watcher.main_runs
            if run in self.active_runs and not self.fetch_run_status(run)]
        if active_main_runs:
            # Still running, nothing more to check on sample
            return False
        if sample_watcher.resistance_run is None:
            sample_watcher.resistance_run = self.run_pipeline(
                sample_watcher,
                PipelineType.RESISTANCE)
            # Launched resistance run, nothing more to check on sample.
            return False
        if sample_watcher.resistance_run in self.active_runs:
            if not self.fetch_run_status(sample_watcher.resistance_run):
                # Still running, nothing more to check on sample.
                return False
        return is_mixed_hcv_complete

    def run_pipeline(self, sample_watcher, pipeline_type):
        run = self.runner.run_pipeline(self, sample_watcher, pipeline_type)
        if run is not None:
            self.active_runs[run] = (sample_watcher, pipeline_type)
        return run

    def fetch_run_status(self, run):
        sample_watcher, pipeline_type = self.active_runs[run]
        is_complete = self.runner.fetch_run_status(run,
                                                   self,
                                                   sample_watcher,
                                                   pipeline_type)
        if is_complete:
            del self.active_runs[run]
        return is_complete


class SampleWatcher:
    def __init__(self, sample_group):
        self.sample_group = sample_group
        self.fastq_datasets = []
        self.mixed_hcv_runs = []
        self.main_runs = []
        self.resistance_run = None

    def __repr__(self):
        enum = self.sample_group.enum
        midi_name = self.sample_group.names[1] and '...'
        return f"SampleWatcher(SampleGroup({enum!r}, ('...', {midi_name!r})))"
