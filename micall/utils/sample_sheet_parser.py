#! /usr/bin/env python3

import argparse
from csv import DictReader
from dataclasses import dataclass
from pathlib import Path
from typing import TextIO, Dict
import csv
import multicsv
from io import StringIO
import json

from micall.utils.sample_sheet_v1_parser import sample_sheet_v1_parser
from micall.utils.sample_sheet_v2_parser import sample_sheet_v2_parser, try_parse_sample_project


def _determine_version(file: multicsv.MultiCSVFile) -> int:
    data_section = file.get('Data')
    if data_section is None:
        raise ValueError("Missing 'Data' section in the sample sheet.")

    for row in csv.DictReader(data_section):
        for value in row.values():
            if try_parse_sample_project(value):
                return 2

    return 1


def _sample_sheet_parser(handle: TextIO) -> Dict[str, object]:
    """
    Parse the contents of SampleSheet.csv, convert contents into a
    Python dictionary object.
    Samples are tracked by sample name and sample number (e.g., S9).
    This is to distinguish replicates of the same sample.
    """

    handle = StringIO(handle.read())
    with multicsv.wrap(handle) as csvfile:
        if _determine_version(csvfile) == 1:
            return sample_sheet_v1_parser(handle)
        else:
            return sample_sheet_v2_parser(csvfile)


@dataclass(frozen=True)
class UnknownSamplesInOverrides(ValueError):
    samples: tuple[str, ...]


def _read_sample_sheet_overrides(override_file, run_info):
    reader = DictReader(override_file)
    project_overrides = {row['sample']: row['project']
                         for row in reader}

    # Ensure that overrides is a subset of the samples in the sample sheet.
    sample_names = {row['filename'] for row in run_info['DataSplit']}
    unknown_samples = []
    for sample_name in project_overrides:
        if sample_name not in sample_names:
            unknown_samples.append(sample_name)
    if unknown_samples:
        raise UnknownSamplesInOverrides(tuple(unknown_samples))

    for row in run_info['DataSplit']:
        sample_name = row['filename']
        old_project = row['project']
        new_project = project_overrides.get(sample_name, old_project)
        row['project'] = new_project


def read_sample_sheet_and_overrides(sample_sheet_path: Path) -> dict[str, object]:
    with sample_sheet_path.open() as sample_sheet_file:
        run_info = _sample_sheet_parser(sample_sheet_file)
    overrides_path = sample_sheet_path.parent / 'SampleSheetOverrides.csv'
    if overrides_path.exists():
        with overrides_path.open() as overrides_file:
            _read_sample_sheet_overrides(overrides_file, run_info)
    return run_info


def main():
    parser = argparse.ArgumentParser(description="Read a sample sheet")
    parser.add_argument("samplesheet", type=Path, help="Path to SampleSheet.csv")
    args = parser.parse_args()

    ss = read_sample_sheet_and_overrides(args.samplesheet)
    print(json.dumps(ss, indent='\t'))


if __name__ == "__main__": main() # noqa
