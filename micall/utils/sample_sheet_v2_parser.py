
from multicsv import MultiCSVFile
from typing import Dict, Any
import csv


def sample_sheet_v2_parser(file: MultiCSVFile) -> Dict[str, object]:
    ret: Dict[str, Any] = {}

    sample_sheet_v2_verifier(file)

    for row in csv.reader(file['BCCFE_Settings']):
        ret[row[0]] = row[1]

    for row in csv.reader(file['Header']):
        ret[row[0]] = row[1]

    ret['Reads'] = []
    for row in csv.reader(file['Reads']):
        ret['Reads'].append(int(row[0]))

    ret['Data'] = {}
    for i, row_d in enumerate(csv.DictReader(file['Data'])):
        number = i + 1
        sample_name = row_d['Sample_Name']
        sample_number = f"S{number}"
        filename = f"{sample_name}_{sample_number}"

        ret['Data'][filename] = {
            'Sample_ID': row_d['Sample_ID'],
            'sample_number': sample_number,
            'orig_sample_name': sample_name,
            'index1': row_d['index'],
            'index2': row_d['index2'],
            'chemistry': ret['Chemistry'],
        }

    ret['DataSplit'] = []
    for i, row_d in enumerate(csv.DictReader(file['BCCFE_Data'])):
        ret['DataSplit'].append({
            'Sample_ID': row_d['Sample_ID'],
            'sample': row_d['Enum'],
            'tags': row_d['Tag'],
            'project': row_d['Project'],
        })

    # Combine Data and DataSplit.
    for datasplit in ret['DataSplit']:
        for filename, data in ret['Data'].items():
            if datasplit['Sample_ID'] == data['Sample_ID']:
                data['tags'] = datasplit['tags']
                datasplit['sample_number'] = datasplit
                datasplit['filename'] = filename
                datasplit.update(data)
                break

    ret['sample_sheet_version'] = ret['SampleSheetVersion']

    return ret


def sample_sheet_v2_verifier(file: MultiCSVFile) -> None:
    if 'Header' not in file.keys():
        raise ValueError("Missing 'Header' section in the sample sheet.")

    # Verify BCCFE_Settings section
    found_version = False
    if 'BCCFE_Settings' in file.keys():
        for line_number, row in enumerate(csv.reader(file['BCCFE_Settings']), start=1):
            if len(row) != 2:
                raise ValueError(
                    f"Line {line_number} in [BCCFE_Settings] section: "
                    f"Expected 2 elements (key-value pair), but got {len(row)}."
                )

            if row[0] == 'SampleSheetVersion':
                found_version = True
    else:
        raise ValueError("Missing 'BCCFE_Settings' section in the sample sheet.")

    if not found_version:
        raise ValueError("Did not found version number in the sample sheet.")

    # Verify Header section
    for line_number, row in enumerate(csv.reader(file['Header']), start=1):
        if len(row) != 2:
            raise ValueError(
                f"Line {line_number} in [Header] section: Expected 2 elements (key-value pair), but got {len(row)}."
            )

    # Verify Reads section
    if 'Reads' in file.keys():
        for line_number, row in enumerate(csv.reader(file['Reads']), start=1):
            if len(row) != 1:
                raise ValueError(f"Line {line_number} in [Reads] section: Expected 1 element, but got {len(row)}.")
            try:
                int(row[0])
            except ValueError:
                raise ValueError(f"Line {line_number} in [Reads] section: Expected an integer but got {row[0]}.")
    else:
        raise ValueError("Missing 'Reads' section in the sample sheet.")

    required_data_fields = ['Sample_ID', 'Sample_Name', 'index', 'index2']
    data_field_errors = []
    if 'Data' in file.keys():
        for line_number, row in enumerate(csv.reader(file['Data']), start=1):
            if line_number == 1:
                headers = row
                for field in required_data_fields:
                    if field not in headers:
                        data_field_errors.append(f"Line {line_number} in [Data] section: "
                                                 f"Expected field '{field}' not found.")
            else:
                if len(row) != len(headers):
                    data_field_errors.append(f"Line {line_number} in [Data] section: Row length {len(row)} "
                                             f"does not match header length {len(headers)}.")
    else:
        raise ValueError("Missing 'Data' section in the sample sheet.")

    if data_field_errors:
        raise ValueError("\n".join(data_field_errors))

    bccfe_data_errors = []
    if 'BCCFE_Data' in file.keys():
        for line_number, row in enumerate(csv.reader(file['BCCFE_Data']), start=1):
            if line_number == 1:
                headers = row
                if 'Sample_ID' not in headers:
                    bccfe_data_errors.append(f"Line {line_number} in [BCCFE_Data] section: "
                                             "Expected field 'Sample_ID' not found.")
            else:
                if len(row) != len(headers):
                    bccfe_data_errors.append(f"Line {line_number} in [BCCFE_Data] section: "
                                             f"Row length {len(row)} does not match header length {len(headers)}.")
    else:
        raise ValueError("Missing 'BCCFE_Data' section in the sample sheet.")

    if bccfe_data_errors:
        raise ValueError("\n".join(bccfe_data_errors))
