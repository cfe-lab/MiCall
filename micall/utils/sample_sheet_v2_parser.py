
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

    return ret


def sample_sheet_v2_verifier(file: MultiCSVFile) -> None:
    assert 'Header' in file.keys(), "Missing 'Header' section in the sample sheet."

    for row in csv.reader(file['BCCFE_Settings']):
        assert len(row) == 2, "Expected key-value pairs in [BCCFE_Settings] section."

    for row in csv.reader(file['Header']):
        assert len(row) == 2, "Expected key-value pairs in [Header] section."

    for row in csv.reader(file['Reads']):
        assert len(row) == 1, "Expected a list in [Reads] section."
        try:
            int(row[0])
        except ValueError:
            raise ValueError("Expected an integer in [Reads] section.")

    required_data_fields = ['Sample_ID', 'Sample_Name', 'index']
    for i, row in enumerate(csv.reader(file['Data'])):
        if i == 0:
            # Check that the required header fields are present
            headers = row
            for field in required_data_fields:
                assert field in headers, f"Expected '{field}' in [Data] section header."
        else:
            # Check the data integrity for required fields
            assert len(row) == len(headers), f"Row length mismatch in [Data] section at line {i+1}"

    for i, row in enumerate(csv.reader(file['BCCFE_Data'])):
        if i == 0:
            # Check that the required header fields are present
            headers = row
            assert 'Sample_ID' in headers, "Expected 'Sample_ID' in [BCCFE_Data] section header."
        else:
            # Check the data integrity for each row
            assert len(row) == len(headers), f"Row length mismatch in [BCCFE_Data] section at line {i+1}"
