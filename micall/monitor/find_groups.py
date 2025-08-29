from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator, Optional, Sequence

from micall.utils.sample_sheet_parser import read_sample_sheet_and_overrides

@dataclass(frozen=True)
class SampleGroup():
    """Protocol for SampleGroup objects used in the monitor module.

    SampleGroup represents a group of related samples (main and optional MIDI).
    """
    enum: str
    names: Sequence[Optional[str]]
    project_codes: Sequence[Optional[str]]


def find_groups(file_names: Iterable[str], sample_sheet_path: str, included_projects: Optional[set[str]] = None) -> Iterator[SampleGroup]:
    """ Group HCV samples with their MIDI partners.

    :param list[str] file_names: a list of FASTQ file names without paths
    :param sample_sheet_path: path to the SampleSheet.csv file
    :param included_projects: project codes to include, or None to include
        all
    """
    run_info = read_sample_sheet_and_overrides(Path(sample_sheet_path))
    data_split = run_info.get('DataSplit')
    if data_split is None:
        raise RuntimeError(f"Missing 'DataSplit' section in {sample_sheet_path}")

    assert hasattr(data_split, '__iter__')
    if not hasattr(data_split, '__iter__'):
        raise RuntimeError(f"Invalid 'DataSplit' section in {sample_sheet_path}")

    midi_hcv_code = 'MidHCV'
    midi_files: dict[str, str] \
        = {row['sample']: row['filename']
           for row in data_split
           if row['project'] == midi_hcv_code
           and isinstance(row['filename'], str)
           and isinstance(row['sample'], str)}
    wide_names = {row['filename']: (row['sample'], row['project'])
                  for row in data_split
                  if (row['project'] != midi_hcv_code and
                      (included_projects is None or
                       row['project'] in included_projects))}
    trimmed_names = {'_'.join(file_name.split('_')[:2]): file_name
                     for file_name in file_names}
    unused_names = set(trimmed_names.values())
    for trimmed_name, file_name in sorted(trimmed_names.items()):
        sample_entry = wide_names.get(trimmed_name)
        if sample_entry is None:
            # Project was not included.
            continue
        sample_name, project_code = sample_entry
        midi_trimmed = midi_files.get(sample_name + 'MIDI')
        if midi_trimmed is None and sample_name.upper().endswith('WG'):
            sample_name = sample_name[:-2]
            midi_trimmed = midi_files.get(sample_name + 'MIDI')
        if midi_trimmed is not None:
            midi_name = trimmed_names.get(midi_trimmed)
        else:
            midi_name = None
        unused_names.discard(file_name)
        if midi_name is not None:
            unused_names.discard(midi_name)
        midi_project = midi_name and midi_hcv_code
        yield SampleGroup(sample_name,
                          (file_name, midi_name),
                          (project_code, midi_project))

    if unused_names:
        sample_names = {file_name: sample_name
                        for sample_name, file_name in midi_files.items()}
        for trimmed_name, file_name in sorted(trimmed_names.items()):
            if file_name in unused_names:
                unused_names.discard(file_name)
                sample_name = sample_names.get(trimmed_name)
                if sample_name is not None:
                    yield SampleGroup(sample_name,
                                      (file_name, None),
                                      (midi_hcv_code, None))
