from enum import Enum
from pathlib import Path

from kiveapi import KiveRunFailedException

ALLOWED_GROUPS = ['Everyone']
# noinspection PyArgumentList
PipelineType = Enum(
    'PipelineType',
    'FILTER_QUALITY MAIN MIDI RESISTANCE MIXED_HCV_MAIN MIXED_HCV_MIDI ' +
    'DENOVO_MAIN DENOVO_MIDI DENOVO_RESISTANCE PROVIRAL')

PIPELINE_GROUPS = {
    PipelineType.FILTER_QUALITY: PipelineType.FILTER_QUALITY,
    PipelineType.MAIN: PipelineType.MAIN,
    PipelineType.MIDI: PipelineType.MAIN,
    PipelineType.RESISTANCE: PipelineType.MAIN,
    PipelineType.MIXED_HCV_MAIN: PipelineType.MIXED_HCV_MAIN,
    PipelineType.MIXED_HCV_MIDI: PipelineType.MIXED_HCV_MAIN,
    PipelineType.DENOVO_MAIN: PipelineType.DENOVO_MAIN,
    PipelineType.DENOVO_MIDI: PipelineType.DENOVO_MAIN,
    PipelineType.DENOVO_RESISTANCE: PipelineType.DENOVO_MAIN,
    PipelineType.PROVIRAL: PipelineType.PROVIRAL
}


class FolderWatcher:
    def __init__(self, base_calls_folder, runner=None):
        """ Set up an instance.

        :param base_calls_folder: path to the BaseCalls folder under a MiSeq
            run folder
        :param runner: an object for running Kive pipelines. Must have these
            methods (this is usually a KiveWatcher instance):
            run_pipeline(folder_watcher, pipeline_type, sample_watcher)
                returns run, or None if that pipeline_type is not configured.
            fetch_run_status(
                old_run,
                folder_watcher,
                pipeline_type,
                sample_watcher) => None if successfully finished, raise if
                run failed, new_run if user cancelled the old one, or old_run
                if it's still running. Also saves the outputs to temporary
                files in the results folder when the run is finished
        """
        self.base_calls_folder = Path(base_calls_folder)
        self.runner = runner
        self.run_folder = (self.base_calls_folder / '../../..').resolve()
        self.run_name = '_'.join(self.run_folder.name.split('_')[:2])
        self.sample_watchers = []
        self.is_folder_failed = False
        self.batch = None
        self.quality_dataset = None
        self.filter_quality_run = None
        self.bad_cycles_dataset = None
        self.active_pipeline_groups = set()  # {pipeline_group}
        self.active_runs = {}  # {run_id: ([sample_watcher], pipeline_type)}
        self.new_runs = set()  # {run_id}
        self.completed_samples = set()  # {fastq1_name}
        self.poll_only_new_runs = False

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
        if self.is_folder_failed:
            return set()
        if self.filter_quality_run['id'] in self.active_runs:
            # Individual runs are waiting for filter quality. Return all.
            all_samples = set(self.all_samples)
            return all_samples - self.completed_samples

        active_samples = set()
        for run, (sample_watchers, pipeline_type) in self.active_runs.items():
            sample_watcher = sample_watchers[0]
            if pipeline_type in (PipelineType.MIDI,
                                 PipelineType.MIXED_HCV_MIDI,
                                 PipelineType.DENOVO_MIDI):
                active_samples.add(sample_watcher.sample_group.names[1])
            else:
                active_samples.add(sample_watcher.sample_group.names[0])
        return active_samples

    @property
    def active_run_count(self):
        if self.is_folder_failed:
            return 0
        if self.filter_quality_run['id'] in self.active_runs:
            # Individual runs are waiting for filter quality.
            # Return n * number of samples, because each can launch n runs.
            n = sum(
                pipeline_id is not None
                for pipeline_id in (self.runner.config.micall_main_pipeline_id,
                                    self.runner.config.denovo_main_pipeline_id,
                                    self.runner.config.mixed_hcv_pipeline_id))

            all_samples = set(self.all_samples)
            return n * len(all_samples - self.completed_samples)

        return len(self.active_runs)

    @property
    def is_complete(self):
        return self.is_folder_failed or not self.active_samples

    def is_pipeline_group_finished(self, pipeline_group):
        return not any(
            PIPELINE_GROUPS[pipeline_type] == pipeline_group
            for sample_watchers, pipeline_type in self.active_runs.values())

    def poll_runs(self):
        if self.filter_quality_run is None:
            self.filter_quality_run = self.run_pipeline(
                PipelineType.FILTER_QUALITY)
            # Launched filter run, nothing more to check.
            return
        if self.filter_quality_run['id'] in self.active_runs:
            if not self.fetch_run_status(self.filter_quality_run):
                # Still running, nothing more to check.
                return
        if self.is_folder_failed:
            return
        for sample_watcher in self.sample_watchers:
            try:
                is_finished = self.poll_sample_runs(sample_watcher)
            except Exception as ex:
                raise RuntimeError(
                    f'Polling sample group {sample_watcher.sample_group.enum}'
                    f' failed.') from ex
            if is_finished:
                self.completed_samples.update(
                    name for name in sample_watcher.sample_group.names
                    if name is not None)

    def poll_sample_runs(self, sample_watcher):
        """ Poll all active runs for a sample.

        :param sample_watcher: details about the sample to poll
        :return: True if the sample is complete (success or failure), otherwise
            False
        """
        if sample_watcher.is_failed:
            # Just wait for remaining runs to finish.
            for run in sample_watcher.runs.values():
                if run['id'] in self.active_runs:
                    self.fetch_run_status(run)
            return True
        is_mixed_hcv_complete = False
        mixed_hcv_run = sample_watcher.runs.get(PipelineType.MIXED_HCV_MAIN)
        if mixed_hcv_run is None:
            if 'HCV' not in sample_watcher.sample_group.names[0]:
                is_mixed_hcv_complete = True
            else:
                mixed_hcv_run = self.run_pipeline(PipelineType.MIXED_HCV_MAIN,
                                                  sample_watcher)
                if mixed_hcv_run is None:
                    is_mixed_hcv_complete = True
                else:
                    if sample_watcher.sample_group.names[1] is not None:
                        self.run_pipeline(PipelineType.MIXED_HCV_MIDI,
                                          sample_watcher)
        else:
            mixed_hcv_midi_run = sample_watcher.runs.get(
                PipelineType.MIXED_HCV_MIDI)
            active_mixed_runs = [
                run for run in (mixed_hcv_run, mixed_hcv_midi_run)
                if (run and run['id'] in self.active_runs
                    and not self.fetch_run_status(run))
            ]
            is_mixed_hcv_complete = not active_mixed_runs

        is_denovo_main_complete = False
        denovo_main_run = sample_watcher.runs.get(PipelineType.DENOVO_MAIN)
        if denovo_main_run is None:
            self.run_pipeline(PipelineType.DENOVO_MAIN, sample_watcher)
            if sample_watcher.sample_group.names[1] is not None:
                self.run_pipeline(PipelineType.DENOVO_MIDI, sample_watcher)
        else:
            denovo_midi_run = sample_watcher.runs.get(PipelineType.DENOVO_MIDI)
            active_main_runs = [
                run for run in (denovo_main_run, denovo_midi_run)
                if (run and run['id'] in self.active_runs
                    and not self.fetch_run_status(run))
            ]
            if not active_main_runs:
                denovo_resistance_run = sample_watcher.runs.get(
                    PipelineType.DENOVO_RESISTANCE)
                if denovo_resistance_run is None:
                    self.run_pipeline(PipelineType.DENOVO_RESISTANCE,
                                      sample_watcher)
                else:
                    is_denovo_main_complete = (
                        denovo_resistance_run['id'] not in self.active_runs
                        or self.fetch_run_status(denovo_resistance_run))
                # TODO: Put proviral run here

        main_run = sample_watcher.runs.get(PipelineType.MAIN)
        if main_run is None:
            self.run_pipeline(PipelineType.MAIN, sample_watcher)
            if sample_watcher.sample_group.names[1] is not None:
                self.run_pipeline(PipelineType.MIDI, sample_watcher)
            # Launched main and midi runs, nothing more to check on sample.
            return sample_watcher.is_failed
        midi_run = sample_watcher.runs.get(PipelineType.MIDI)
        active_main_runs = [
            run for run in (main_run, midi_run)
            if run and run['id'] in self.active_runs
            and not self.fetch_run_status(run)
        ]
        if active_main_runs:
            # Still running, nothing more to check on sample
            return sample_watcher.is_failed
        resistance_run = sample_watcher.runs.get(PipelineType.RESISTANCE)
        if resistance_run is None:
            self.run_pipeline(PipelineType.RESISTANCE, sample_watcher)
            # Launched resistance run, nothing more to check on sample.
            return sample_watcher.is_failed
        if resistance_run['id'] in self.active_runs:
            if not self.fetch_run_status(resistance_run):
                # Still running, nothing more to check on sample.
                return sample_watcher.is_failed

        # TODO Take a closer look at this and resistance
        # add is_resistance_complete, and do not return right away
        # Or move proviral up above resistance
        proviral_run = sample_watcher.runs.get(PipelineType.PROVIRAL)
        if proviral_run is None:
            self.run_pipeline(PipelineType.PROVIRAL, sample_watcher)
            # Launched resistance run, nothing more to check on sample.
            return sample_watcher.is_failed
        if proviral_run['id'] in self.active_runs:
            if not self.fetch_run_status(proviral_run):
                # Still running, nothing more to check on sample.
                return sample_watcher.is_failed
        return ((is_denovo_main_complete and is_mixed_hcv_complete)
                or sample_watcher.is_failed)

    @property
    def has_new_runs(self):
        return bool(self.new_runs)

    def run_pipeline(self, pipeline_type, sample_watcher=None):
        if sample_watcher and sample_watcher.is_failed:
            # Don't start runs when the sample has already failed.
            return None
        run = self.runner.run_pipeline(self, pipeline_type, sample_watcher)
        if run is not None:
            self.add_run(run, pipeline_type, sample_watcher)
        return run

    def add_run(self,
                run,
                pipeline_type,
                sample_watcher=None,
                is_complete=False):
        pipeline_group = PIPELINE_GROUPS[pipeline_type]
        self.active_pipeline_groups.add(pipeline_group)
        if not is_complete:
            sample_watchers, _ = self.active_runs.setdefault(
                run['id'], ([], pipeline_type))
            sample_watchers.append(sample_watcher)
            self.new_runs.add(run['id'])
        if pipeline_type == PipelineType.FILTER_QUALITY:
            self.filter_quality_run = run
        else:
            sample_watcher.runs[pipeline_type] = run

    def fetch_run_status(self, run):
        sample_watchers, pipeline_type = self.active_runs[run['id']]
        try:
            self.new_runs.remove(run['id'])
            is_new = True
        except KeyError:
            is_new = False
        is_complete = False
        if self.poll_only_new_runs and not is_new:
            return is_complete
        try:
            new_run = self.runner.fetch_run_status(run, self, pipeline_type,
                                                   sample_watchers)
            if new_run is None:
                is_complete = True
            elif new_run is not run:
                del self.active_runs[run['id']]
                for sample_watcher in sample_watchers:
                    self.add_run(new_run, pipeline_type, sample_watcher)

        except KiveRunFailedException:
            for sample_watcher in sample_watchers:
                if sample_watcher is None:
                    self.is_folder_failed = True
                else:
                    sample_watcher.is_failed = True
            is_complete = True
        if is_complete:
            del self.active_runs[run['id']]
        return is_complete


class SampleWatcher:
    def __init__(self, sample_group):
        self.sample_group = sample_group
        self.fastq_datasets = []
        self.name_datasets = []
        self.runs = {}  # {pipeline_type: run}
        self.is_failed = False

    def __repr__(self):
        enum = self.sample_group.enum
        midi_name = self.sample_group.names[1] and '...'
        return f"SampleWatcher(SampleGroup({enum!r}, ('...', {midi_name!r})))"
