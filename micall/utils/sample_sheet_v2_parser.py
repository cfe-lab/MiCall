
from dataclasses import dataclass
from multicsv import MultiCSVFile
from typing import Dict, Any, Sequence, Optional
import csv


SAMPLE_SHEET_MAGIC = 'bccfe'
SAMPLE_SHEET_SEPARATOR = '_'


@dataclass(frozen=True)
class SampleProject:
    magic: str
    major_version: int
    sample_tag: str
    project_codes: Sequence[str]


def try_parse_sample_project(value: Optional[str | list[str]]) -> Optional[SampleProject]:
    if isinstance(value, list):
        # This is a quirk of csv parsing where empty fields are parsed as empty lists.
        assert len(value) == 1, "Cannot parse sample project from list with multiple elements."
        assert value[0] == '', "Cannot parse sample project from list with non-empty element."
        return None

    if value is None or not value.startswith(SAMPLE_SHEET_MAGIC):
        return None

    sample_project_parts = value.split(SAMPLE_SHEET_SEPARATOR)
    if len(sample_project_parts) < 3:
        return None

    (magic, major_version, sample_tag, *project_codes) = sample_project_parts
    return SampleProject(magic=magic,
                         major_version=int(major_version),
                         sample_tag=sample_tag,
                         project_codes=project_codes,
                         )


def parse_sample_project(sample_project: Optional[str]) -> SampleProject:
    ret = try_parse_sample_project(sample_project)
    assert ret is not None
    return ret


def try_parse_sample_name(value: Optional[str]) -> Optional[Sequence[str]]:
    if value is None:
        return None

    sample_name_parts = value.split(SAMPLE_SHEET_SEPARATOR)
    if len(sample_name_parts) < 1:
        return None

    return sample_name_parts


def parse_sample_name(sample_name: Optional[str]) -> Sequence[str]:
    ret = try_parse_sample_name(sample_name)
    assert ret is not None
    return ret


def sample_sheet_v2_parser(file: MultiCSVFile) -> Dict[str, object]:
    ret: Dict[str, Any] = {}

    sample_sheet_v2_verifier(file)

    for row_d in csv.DictReader(file['Data']):
        proj = parse_sample_project(row_d['Sample_Project'])
        ret['sample_sheet_version'] = proj.major_version

    for row in csv.reader(file['Header']):
        if row:
            ret[row[0]] = row[1]

    ret['Reads'] = []
    for row in csv.reader(file['Reads']):
        if row:
            ret['Reads'].append(int(row[0]))

    ret['Data'] = {}
    ret['DataSplit'] = []
    for i, row_d in enumerate(csv.DictReader(file['Data'])):
        number = i + 1
        sample_name = row_d['Sample_Name']
        sample_project = row_d['Sample_Project']
        sample_number = f"S{number}"
        filename = f"{sample_name}_{sample_number}"
        parsed_sample_project = parse_sample_project(sample_project)
        parsed_sample_name = parse_sample_name(sample_name)

        ret['Data'][filename] = {
            'Sample_ID': row_d['Sample_ID'],
            'sample_number': sample_number,
            'orig_sample_name': sample_name,
            'index1': row_d['index'],
            'index2': row_d['index2'],
            'chemistry': ret['Chemistry'],
            'tags': parsed_sample_project.sample_tag,
        }

        for sample, project in zip(parsed_sample_name,
                                   parsed_sample_project.project_codes):
            datasplit = ret['Data'][filename].copy()
            datasplit['sample'] = sample
            datasplit['project'] = project
            datasplit['filename'] = filename
            ret['DataSplit'].append(datasplit)

    return ret


def sample_sheet_v2_verifier(file: MultiCSVFile) -> None:
    if 'Header' not in file.keys():
        raise ValueError("Missing 'Header' section in the sample sheet.")

    # Verify Header section
    for line_number, row in enumerate(csv.reader(file['Header']), start=1):
        if len(row) == 0:
            continue
        if len(row) != 2:
            raise ValueError(
                f"Line {line_number} in [Header] section:"
                " Expected 2 elements (key-value pair), "
                f"but got {len(row)}: {row!r}."
            )

    # Verify Reads section
    if 'Reads' in file.keys():
        for line_number, row in enumerate(csv.reader(file['Reads']), start=1):
            if len(row) == 0:
                continue
            if len(row) != 1:
                raise ValueError(f"Line {line_number} in [Reads] section:"
                                 " Expected 1 element, "
                                 f"but got {len(row)}: {row!r}.")
            try:
                int(row[0])
            except BaseException:
                raise ValueError(f"Line {line_number} in [Reads] section:"
                                 " Expected an integer, "
                                 f"but got {row[0]}: {row!r}.")
    else:
        raise ValueError("Missing 'Reads' section in the sample sheet.")

    required_data_fields = ['Sample_ID', 'Sample_Name', 'index', 'index2', 'Sample_Project']
    data_field_errors = []
    if 'Data' in file.keys():
        for line_number, row in enumerate(csv.reader(file['Data']), start=1):
            if line_number == 1:
                headers = row
                for field in required_data_fields:
                    if field not in headers:
                        data_field_errors.append(f"Line {line_number} in [Data] section: "
                                                 f"Expected field {field!r} not found.")
            else:
                if len(row) != len(headers):
                    data_field_errors.append(f"Line {line_number} in [Data] section: Row length {len(row)} "
                                             f"does not match header length {len(headers)}.")

        reference_parsed_project: Optional[SampleProject] = None
        for row_d in csv.DictReader(file['Data']):
            sample_project = row_d['Sample_Project']
            if reference_parsed_project is not None:
                parsed_project = try_parse_sample_project(sample_project)
                if parsed_project is None:
                    data_field_errors.append(f"Cannot parse sample project field: {sample_project!r}.")
                    continue

                if parsed_project.major_version != reference_parsed_project.major_version:
                    left = parsed_project.major_version
                    right = reference_parsed_project.major_version
                    data_field_errors.append(f"Got two different versions: {left} and {right}.")
            else:
                reference_parsed_project = try_parse_sample_project(sample_project)
                if reference_parsed_project is None:
                    data_field_errors.append(f"Cannot parse sample project field: {sample_project!r}.")
    else:
        raise ValueError("Missing 'Data' section in the sample sheet.")

    if data_field_errors:
        raise ValueError("\n".join(data_field_errors))
