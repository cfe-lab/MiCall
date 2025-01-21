#! /usr/bin/env python3.6

from csv import DictReader
from pathlib import Path
from typing import TextIO, Dict
import csv
import multicsv
from io import StringIO
import json

from micall.utils.sample_sheet_v1_parser import sample_sheet_v1_parser
from micall.utils.sample_sheet_v2_parser import sample_sheet_v2_parser, try_parse_sample_project


def determine_version(file: multicsv.MultiCSVFile) -> int:
    data_section = file.get('Data')
    if data_section is None:
        raise ValueError("Missing 'Data' section in the sample sheet.")

    for row in csv.DictReader(data_section):
        for value in row.values():
            if try_parse_sample_project(value):
                return 2

    return 1


def sample_sheet_parser(handle: TextIO) -> Dict[str, object]:
    """
    Parse the contents of SampleSheet.csv, convert contents into a
    Python dictionary object.
    Samples are tracked by sample name and sample number (e.g., S9).
    This is to distinguish replicates of the same sample.
    """

    handle = StringIO(handle.read())
    with multicsv.wrap(handle) as csvfile:
        if determine_version(csvfile) == 1:
            return sample_sheet_v1_parser(handle)
        else:
            return sample_sheet_v2_parser(csvfile)


def read_sample_sheet_overrides(override_file, run_info):
    reader = DictReader(override_file)
    project_overrides = {row['sample']: row['project']
                         for row in reader}
    for row in run_info['DataSplit']:
        sample_name = row['filename']
        old_project = row['project']
        new_project = project_overrides.get(sample_name, old_project)
        row['project'] = new_project


def read_sample_sheet_and_overrides(sample_sheet_path: Path) -> dict:
    with sample_sheet_path.open() as sample_sheet_file:
        run_info = sample_sheet_parser(sample_sheet_file)
    overrides_path = sample_sheet_path.parent / 'SampleSheetOverrides.csv'
    if overrides_path.exists():
        with overrides_path.open() as overrides_file:
            read_sample_sheet_overrides(overrides_file, run_info)
    return run_info


def main():
    import argparse

    parser = argparse.ArgumentParser(description="Read a sample sheet")
    parser.add_argument("samplesheet")
    args = parser.parse_args()

    with open(args.samplesheet) as f:
        ss = sample_sheet_parser(f)
        print(json.dumps(ss, indent='\t'))


if __name__ == "__main__":
    main()
