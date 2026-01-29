import unittest
from io import StringIO
from pathlib import Path

from micall.utils.sample_sheet_parser import _sample_sheet_parser, _read_sample_sheet_overrides

"""
Test suite for the sample sheet parser.

Important cases to test:
 - a version 1 sample sheet
 - a version 2 sample sheet
 - a version 2024 sample sheet
 - entries (wells) with multiple samples
 - defective file whose Data portion does not include a Sample_Name column

Later: do we need to check
 - "is_T_primer"?
"""


# noinspection DuplicatedCode
class VersionOneTest(unittest.TestCase):
    stub_sample_sheet = """
[Header]
IEMFileVersion,3
Investigator Name,RL
Project Name,10-Jul-2014_v1test
Experiment Name,10-Jul-2014_v1test
Date,07/10/2014
Workflow,GenerateFASTQ
Assay,Nextera
Description,Nextera
Chemistry,Amplicon
[Reads]
251
251
[Settings]
[Data]
Sample_ID,Sample_Name,Sample_Plate,Sample_Well,index,index2,Sample_Project,Description,GenomeFolder
CFE_SomeId_10-Jul-2014_N501-N701_Sample1_Proj1,Sample1_Proj1,10-Jul-2014_testing,N/A,ACGTACGT,TGCATGCA,\
10-Jul-2014_testing,Research:Sample1_Proj1:TRUE Comments:Sample1_Proj1:thisiscommentone \
Disablecontamcheck:Sample1_Proj1:FALSE,
CFE_SomeId_10-Jul-2014_N501-N702_Sample2_Proj2,Sample2_Proj2,10-Jul-2014_testing,N/A,AAAAGGGG,CCCCTTTT,\
10-Jul-2014_testing,Research:Sample2_Proj2:FALSE Comments:Sample2_Proj2:thisiscommenttwo \
Chemistry:Sample2_Proj2:BreakingBad Disablecontamcheck:Sample2_Proj2:TRUE,
"""
    clean_filenames = ["Sample1-Proj1_S1", "Sample2-Proj2_S2"]

    def setUp(self):
        self.maxDiff = None
        self.ss = _sample_sheet_parser(StringIO(self.stub_sample_sheet))

    def test_keys_correct(self):
        """
        Test that all of the keys in the resulting dictionary are correct.

        There should be one key per [Header] field, as well as "Data", "DataSplit",
        and sample_sheet_version.
        """
        self.assertSetEqual(
            set(self.ss.keys()),
            {"IEMFileVersion", "Investigator Name", "Project Name", "Experiment Name", "Date", "Workflow",
             "Assay", "Description", "Chemistry", "Data", "DataSplit", "sample_sheet_version",
             "Reads"})

    def test_preamble_correct(self):
        """
        Test that all of the header stuff, as well as the sample sheet version, was set correctly.
        """
        self.assertEqual(self.ss["IEMFileVersion"], "3")
        self.assertEqual(self.ss["Investigator Name"], "RL")
        self.assertEqual(self.ss["Project Name"], "10-Jul-2014_v1test")
        self.assertEqual(self.ss["Experiment Name"], "10-Jul-2014_v1test")
        self.assertEqual(self.ss["Date"], "07/10/2014")
        self.assertEqual(self.ss["Workflow"], "GenerateFASTQ")
        self.assertEqual(self.ss["Assay"], "Nextera")
        self.assertEqual(self.ss["Description"], "Nextera")
        self.assertEqual(self.ss["Chemistry"], "Amplicon")
        self.assertEqual(self.ss["sample_sheet_version"], 1)
        self.assertEqual(self.ss["Reads"], [251, 251])

    def test_data(self):
        """
        Check each entry in the "Data" dictionary.
        """
        self.maxDiff = None
        self.assertDictEqual(
            self.ss["Data"],
            {
                self.clean_filenames[0]:
                    {
                        "index1": "ACGTACGT",
                        "index2": "TGCATGCA",
                        "tags": "N501-N701",
                        "comments": "thisiscommentone",
                        "disable_contamination_check": False,
                        "research": True,
                        "chemistry": "Nextera",
                        "orig_sample_name": "Sample1_Proj1",
                        "is_T_primer": False
                    },
                self.clean_filenames[1]:
                    {
                        "index1": "AAAAGGGG",
                        "index2": "CCCCTTTT",
                        "tags": "N501-N702",
                        "comments": "thisiscommenttwo",
                        "disable_contamination_check": True,
                        "research": False,
                        "chemistry": "BreakingBad",
                        "orig_sample_name": "Sample2_Proj2",
                        "is_T_primer": False
                    }
            }
        )

    def test_datasplit(self):
        """
        Check each entry in the "DataSplit" list.
        """
        self.assertListEqual(
            self.ss["DataSplit"],
            [
                {
                    "sample": "Sample1",
                    "project": "Proj1",
                    "filename": "Sample1-Proj1_S1",
                    "tags": "N501-N701",
                    "index1": "ACGTACGT",
                    "index2": "TGCATGCA",
                    "sample_number": "S1",
                    "chemistry": "Nextera",
                    "disable_contamination_check": False,
                    "research": True,
                    "comments": "thisiscommentone",
                    "orig_sample_name": "Sample1_Proj1",
                    "is_T_primer": False,
                },
                {
                    "sample": "Sample2",
                    "project": "Proj2",
                    "filename": "Sample2-Proj2_S2",
                    "tags": "N501-N702",
                    "index1": "AAAAGGGG",
                    "index2": "CCCCTTTT",
                    "sample_number": "S2",
                    "chemistry": "BreakingBad",
                    "disable_contamination_check": True,
                    "research": False,
                    "comments": "thisiscommenttwo",
                    "orig_sample_name": "Sample2_Proj2",
                    "is_T_primer": False
                }
            ]
        )


# noinspection DuplicatedCode
class MultipleSamplesV1Test(unittest.TestCase):
    stub_sample_sheet = """
[Header]
IEMFileVersion,3
Investigator Name,RL
Project Name,11-Jul-2014_v1multisamplestest
Experiment Name,11-Jul-2014_v1multisamplestest
Date,07/11/2014
Workflow,GenerateFASTQ
Assay,Nextera
Description,Nextera
Chemistry,Amplicon
[Reads]
251
251
[Settings]
[Data]
Sample_ID,Sample_Name,Sample_Plate,Sample_Well,index,index2,Sample_Project,Description,GenomeFolder
CFE_SomeId_11-Jul-2014_N501-N701_Sample1_Proj1;CFE_SomeId_11-Jul-2014_N501-N701_Sample2_Proj2,\
Sample1_Proj1;Sample2_Proj2,11-Jul-2014_testing,N/A,ACGTACGT,TGCATGCA,11-Jul-2014_testing,\
Research:Sample1_Proj1:TRUE;Sample2_Proj2:FALSE \
Comments:Sample1_Proj1:thisiscommentone;Sample2_Proj2:thisiscommenttwo \
Disablecontamcheck:Sample1_Proj1:FALSE;Sample2_Proj2:TRUE,
CFE_SomeId_11-Jul-2014_N501-N702_Sample3_Proj3;CFE_SomeId_11-Jul-2014_N501-N702_Sample4_Proj4,\
Sample3_Proj3;Sample4_Proj4,11-Jul-2014_testing,N/A,AAAAGGGG,CCCCTTTT,11-Jul-2014_testing,\
Research:Sample3_Proj3:FALSE;Sample4_Proj4:FALSE \
Comments:Sample3_Proj3:thisiscommentthree;Sample4_Proj4:thisiscommentfour \
Chemistry:Sample3_Proj3:BreakingBad;Sample4_Proj4:MrWizard \
Disablecontamcheck:Sample3_Proj3:TRUE;Sample4_Proj4:TRUE,
"""
    clean_filenames = ["Sample1-Proj1-Sample2-Proj2_S1", "Sample3-Proj3-Sample4-Proj4_S2"]

    def setUp(self):
        self.ss = _sample_sheet_parser(StringIO(self.stub_sample_sheet))

    def test_keys_correct(self):
        """
        Test that all of the keys in the resulting dictionary are correct.

        There should be one key per [Header] field, as well as "Data", "DataSplit",
        and sample_sheet_version.
        """
        self.assertSetEqual(
            set(self.ss.keys()),
            {"IEMFileVersion", "Investigator Name", "Project Name", "Experiment Name", "Date", "Workflow",
             "Assay", "Description", "Chemistry", "Data", "DataSplit", "sample_sheet_version",
             "Reads"})

    def test_preamble_correct(self):
        """
        Test that all of the header stuff, as well as the sample sheet version, was set correctly.
        """
        self.assertEqual(self.ss["IEMFileVersion"], "3")
        self.assertEqual(self.ss["Investigator Name"], "RL")
        self.assertEqual(self.ss["Project Name"], "11-Jul-2014_v1multisamplestest")
        self.assertEqual(self.ss["Experiment Name"], "11-Jul-2014_v1multisamplestest")
        self.assertEqual(self.ss["Date"], "07/11/2014")
        self.assertEqual(self.ss["Workflow"], "GenerateFASTQ")
        self.assertEqual(self.ss["Assay"], "Nextera")
        self.assertEqual(self.ss["Description"], "Nextera")
        self.assertEqual(self.ss["Chemistry"], "Amplicon")
        self.assertEqual(self.ss["sample_sheet_version"], 1)

    def test_data(self):
        """
        Check each entry in the "Data" dictionary.
        """
        self.assertDictEqual(
            self.ss["Data"],
            {
                self.clean_filenames[0]:
                    {
                        "index1": "ACGTACGT",
                        "index2": "TGCATGCA",
                        "tags": "N501-N701",
                        "comments": "thisiscommentone",
                        "disable_contamination_check": False,
                        "research": True,
                        "chemistry": "Nextera",
                        "orig_sample_name": "Sample1_Proj1;Sample2_Proj2",
                        "is_T_primer": False
                    },
                self.clean_filenames[1]:
                    {
                        "index1": "AAAAGGGG",
                        "index2": "CCCCTTTT",
                        "tags": "N501-N702",
                        "comments": "thisiscommentthree",
                        "disable_contamination_check": True,
                        "research": False,
                        "chemistry": "BreakingBad",
                        "orig_sample_name": "Sample3_Proj3;Sample4_Proj4",
                        "is_T_primer": False
                    }
            }
        )

    def test_datasplit(self):
        """
        Check each entry in the "DataSplit" list.
        """
        self.assertListEqual(
            self.ss["DataSplit"],
            [
                {
                    "sample": "Sample1",
                    "project": "Proj1",
                    "filename": "Sample1-Proj1-Sample2-Proj2_S1",
                    "index1": "ACGTACGT",
                    "index2": "TGCATGCA",
                    "sample_number": "S1",
                    "chemistry": "Nextera",
                    "disable_contamination_check": False,
                    "research": True,
                    "comments": "thisiscommentone",
                    "orig_sample_name": "Sample1_Proj1;Sample2_Proj2",
                    "tags": "N501-N701",
                    "is_T_primer": False,
                },
                {
                    "sample": "Sample2",
                    "project": "Proj2",
                    "filename": "Sample1-Proj1-Sample2-Proj2_S1",
                    "index1": "ACGTACGT",
                    "index2": "TGCATGCA",
                    "sample_number": "S1",
                    "chemistry": "Nextera",
                    "disable_contamination_check": True,
                    "research": False,
                    "comments": "thisiscommenttwo",
                    "orig_sample_name": "Sample1_Proj1;Sample2_Proj2",
                    "tags": "N501-N701",
                    "is_T_primer": False,
                },
                {
                    "sample": "Sample3",
                    "project": "Proj3",
                    "filename": "Sample3-Proj3-Sample4-Proj4_S2",
                    "index1": "AAAAGGGG",
                    "index2": "CCCCTTTT",
                    "sample_number": "S2",
                    "chemistry": "BreakingBad",
                    "disable_contamination_check": True,
                    "research": False,
                    "comments": "thisiscommentthree",
                    "orig_sample_name": "Sample3_Proj3;Sample4_Proj4",
                    "tags": "N501-N702",
                    "is_T_primer": False
                },
                {
                    "sample": "Sample4",
                    "project": "Proj4",
                    "filename": "Sample3-Proj3-Sample4-Proj4_S2",
                    "index1": "AAAAGGGG",
                    "index2": "CCCCTTTT",
                    "sample_number": "S2",
                    "chemistry": "MrWizard",
                    "disable_contamination_check": True,
                    "research": False,
                    "comments": "thisiscommentfour",
                    "orig_sample_name": "Sample3_Proj3;Sample4_Proj4",
                    "tags": "N501-N702",
                    "is_T_primer": False
                }
            ]
        )


class OtherTest(unittest.TestCase):
    def setUp(self):
        self.addTypeEqualityFunc(str, self.assertMultiLineEqual)

    def test_no_sample_name(self):
        """
        Throws an exception if the Data portion has no Sample_Name column.
        """

        stub_sample_sheet = """
[Header]
IEMFileVersion,3
Investigator Name,RL
Project Name,11-Jul-2014_nosamplenametest
Experiment Name,11-Jul-2014_nosamplenametest
Date,07/11/2014
Workflow,GenerateFASTQ
Assay,Nextera
Description,Nextera
Chemistry,Amplicon
[Reads]
251
251
[Settings]
[Data]
There,Is,No,Sample,Name
A,B,C,D,E
"""

        with self.assertRaises(ValueError) as assertion:
            _sample_sheet_parser(StringIO(stub_sample_sheet))
        self.assertEqual("sample sheet data header does not include Sample_Name",
                         assertion.exception.args[0])

    def test_no_index2(self):
        """
        Throws an exception if the Data portion has no Sample_Name column.
        """

        stub_sample_sheet = """
[Header]
IEMFileVersion,3
Investigator Name,RL
Project Name,11-Jul-2014_nosamplenametest
Experiment Name,11-Jul-2014_nosamplenametest
Date,07/11/2014
Workflow,GenerateFASTQ
Assay,Nextera
Description,Nextera
Chemistry,Amplicon
[Reads]
251
251
[Settings]
[Data]
Sample_ID,Sample_Name,Sample_Plate,Sample_Well,index,Sample_Project,Description,GenomeFolder
CFE_SomeId_10-Jul-2014_N501_Sample1_Proj1,Sample1_Proj1,10-Jul-2014_testing,N/A,ACGTACGT,\
10-Jul-2014_testing,Research:Sample1_Proj1:TRUE Comments:Sample1_Proj1:thisiscommentone \
Disablecontamcheck:Sample1_Proj1:FALSE,
CFE_SomeId_10-Jul-2014_N501_Sample2_Proj2,Sample2_Proj2,10-Jul-2014_testing,N/A,AAAAGGGG,\
10-Jul-2014_testing,Research:Sample2_Proj2:FALSE Comments:Sample2_Proj2:thisiscommenttwo \
Chemistry:Sample2_Proj2:BreakingBad Disablecontamcheck:Sample2_Proj2:TRUE,
"""

        ss = _sample_sheet_parser(StringIO(stub_sample_sheet))
        sample = ss['Data']['Sample1-Proj1_S1']
        self.assertEqual('ACGTACGT', sample['index1'])
        self.assertEqual('X', sample['index2'])
        self.assertEqual('N501-X', sample['tags'])

    def test_extra_commas(self):
        """
        Throws an exception if the Data portion has no Sample_Name column.
        """

        stub_sample_sheet = """
[Header],,,,,,,
IEMFileVersion,3,,,,,,,
Investigator Name,RL,,,,,,,
Project Name,10-Jul-2014,,,,,,,
Experiment Name,10-Jul-2014,,,,,,,
Date,07/10/2014,,,,,,,
Workflow,GenerateFASTQ,,,,,,,
Assay,Nextera,,,,,,,
Description,Nextera,,,,,,,
Chemistry,Amplicon,,,,,,,
[Reads],,,,,,,
251,,,,,,,
251,,,,,,,
[Settings],,,,,,,
[Data],,,,,,,
Sample_ID,Sample_Name,Sample_Plate,Sample_Well,index,index2,Sample_Project,Description,GenomeFolder
CFE_SomeId_10-Jul-2014_N501-N701_Sample1_Proj1,Sample1_Proj1,10-Jul-2014_testing,N/A,ACGTACGT,TGCATGCA,\
10-Jul-2014_testing,Research:Sample1_Proj1:TRUE Comments:Sample1_Proj1:thisiscommentone \
Disablecontamcheck:Sample1_Proj1:FALSE,
CFE_SomeId_10-Jul-2014_N501-N702_Sample2_Proj2,Sample2_Proj2,10-Jul-2014_testing,N/A,AAAAGGGG,CCCCTTTT,\
10-Jul-2014_testing,Research:Sample2_Proj2:FALSE Comments:Sample2_Proj2:thisiscommenttwo \
Chemistry:Sample2_Proj2:BreakingBad Disablecontamcheck:Sample2_Proj2:TRUE,
"""

        ss = _sample_sheet_parser(StringIO(stub_sample_sheet))
        self.assertEqual(ss["Experiment Name"], "10-Jul-2014")

    def test_underscores_in_sample_name(self):
        """
        Extracts the correct project code and sample name in presence of underscores.
        """

        stub_sample_sheet = """
[Header]
IEMFileVersion,3
Investigator Name,RL
Project Name,10-Jul-2014_v1test
Experiment Name,10-Jul-2014_v1test
Date,07/10/2014
Workflow,GenerateFASTQ
Assay,Nextera
Description,Nextera
Chemistry,Amplicon
[Reads]
251
251
[Settings]
[Data]
Sample_ID,Sample_Name,Sample_Plate,Sample_Well,index,index2,Sample_Project,Description,GenomeFolder
CFE_SomeId_10-Jul-2014_N501-N701_Sample1_Proj1,Sample1_Proj1,10-Jul-2014_testing,N/A,ACGTACGT,TGCATGCA,\
10-Jul-2014_testing,Research:Sample1_Proj1:TRUE Comments:Sample1_Proj1:thisiscommentone \
Disablecontamcheck:Sample1_Proj1:FALSE,
CFE_SomeId_10-Jul-2014_N501-N702_Sample2_Proj2,Sample2_Proj2,10-Jul-2014_testing,N/A,AAAAGGGG,CCCCTTTT,\
10-Jul-2014_testing,Research:Sample2_Foo_Proj2:FALSE Comments:Sample2_Foo_Proj2:thisiscommenttwo \
Chemistry:Sample2_Foo_Proj2:BreakingBad Disablecontamcheck:Sample2_Foo_Proj2:TRUE,
"""

        ss = _sample_sheet_parser(StringIO(stub_sample_sheet))
        split_rows = ss['DataSplit']
        assert len(split_rows) == 2

        assert split_rows[0]['filename'] == 'Sample1-Proj1_S1'
        assert split_rows[1]['filename'] == 'Sample2-Proj2_S2'

        assert split_rows[0]['project'] == 'Proj1'
        assert split_rows[1]['project'] == 'Proj2'

        assert split_rows[0]['sample'] == 'Sample1'
        assert split_rows[1]['sample'] == 'Sample2'

        assert split_rows[0]['sample_number'] == 'S1'
        assert split_rows[1]['sample_number'] == 'S2'


def test_read_sample_sheet_overrides(tmpdir):
    sample_sheet_path = Path(str(tmpdir)) / 'SampleSheet.csv'
    overrides_path = sample_sheet_path.parent / 'SampleSheetOverrides.csv'

    sample_sheet_path.write_text("""\
[Header]
IEMFileVersion,3
Investigator Name,RL
Project Name,10-Jul-2014_v1test
Experiment Name,10-Jul-2014_v1test
Date,07/10/2014
Workflow,GenerateFASTQ
Assay,Nextera
Description,Nextera
Chemistry,Amplicon
[Reads]
251
251
[Settings]
[Data]
Sample_ID,Sample_Name,Sample_Plate,Sample_Well,index,index2,Sample_Project,Description,GenomeFolder
CFE_SomeId_10-Jul-2014_N501-N701_Sample1_Proj1,Sample1_Proj1,10-Jul-2014_testing,N/A,ACGTACGT,TGCATGCA,\
10-Jul-2014_testing,Research:Sample1_Proj1:TRUE Comments:Sample1_Proj1:thisiscommentone \
Disablecontamcheck:Sample1_Proj1:FALSE,
CFE_SomeId_10-Jul-2014_N501-N702_Sample2_Proj2,Sample2_Proj2,10-Jul-2014_testing,N/A,AAAAGGGG,CCCCTTTT,\
10-Jul-2014_testing,Research:Sample2_Proj2:FALSE Comments:Sample2_Proj2:thisiscommenttwo \
Chemistry:Sample2_Proj2:BreakingBad Disablecontamcheck:Sample2_Proj2:TRUE,
""")
    overrides_path.write_text("""\
sample,project
Sample2-Proj2_S2,AltB
""")

    with sample_sheet_path.open() as f:
        run_info = _sample_sheet_parser(f)

    with overrides_path.open() as f:
        _read_sample_sheet_overrides(f, run_info)

    sample_map = run_info['Data']
    assert len(sample_map) == 2

    split_rows = run_info['DataSplit']
    assert len(split_rows) == 2
    assert split_rows[0]['project'] == 'Proj1'
    assert split_rows[1]['project'] == 'AltB'


# noinspection DuplicatedCode
class VersionTwoTest(unittest.TestCase):
    # Stub is taken from QAI sample sheet writer test case.
    stub_sample_sheet = """\
[Header]
IEMFileVersion,5
Investigator Name,JN
Project Name,20-Jul-2017.M04401
Experiment Name,20-Jul-2017.M04401
Date,07/31/2024
Workflow,Workflow
Assay,Chemistry
Description,Chemistry
Chemistry,Chemistry
[Reads]
251
251
[Settings]
[Data]
Sample_ID,Sample_Name,Sample_Plate,Sample_Well,index,index2,GenomeFolder,Sample_Project
1,1234_1234,20-Jul-2017.M04401,N/A,GTGTTGCT,CATGGTCT,"",bccfe_2_N501-N701_DRT_INT
2,4321,20-Jul-2017.M04401,N/A,GTGTTGCA,CATTGTCA,"",bccfe_2_N501-N702_DRT
"""

    clean_filenames = ["1234_1234_S1", "4321_S2"]

    def setUp(self):
        self.maxDiff = None
        self.ss = _sample_sheet_parser(StringIO(self.stub_sample_sheet))

    def test_preamble_correct(self):
        """
        Test that all of the header stuff, as well as the sample sheet version, was set correctly.
        """

        self.assertEqual(self.ss["IEMFileVersion"], "5")
        self.assertEqual(self.ss["Investigator Name"], "JN")
        self.assertEqual(self.ss["Project Name"], "20-Jul-2017.M04401")
        self.assertEqual(self.ss["Experiment Name"], "20-Jul-2017.M04401")
        self.assertEqual(self.ss["Date"], "07/31/2024")
        self.assertEqual(self.ss["Workflow"], "Workflow")
        self.assertEqual(self.ss["Assay"], "Chemistry")
        self.assertEqual(self.ss["Description"], "Chemistry")
        self.assertEqual(self.ss["Chemistry"], "Chemistry")
        self.assertEqual(self.ss["Reads"], [251, 251])

    def test_data(self):
        """
        Check each entry in the "Data" dictionary.
        """

        self.assertDictEqual(
            self.ss["Data"],
            {
                self.clean_filenames[0]:
                    {
                        "Sample_ID": "1",
                        "sample_number": "S1",
                        "index1": "GTGTTGCT",
                        "index2": "CATGGTCT",
                        "tags": "N501-N701",
                        "chemistry": "Chemistry",
                        "orig_sample_name": "1234_1234",
                    },
                self.clean_filenames[1]:
                    {
                        "Sample_ID": "2",
                        "sample_number": "S2",
                        "index1": "GTGTTGCA",
                        "index2": "CATTGTCA",
                        "tags": "N501-N702",
                        "chemistry": "Chemistry",
                        "orig_sample_name": "4321",
                    },
            }
        )

    def test_datasplit(self):
        """
        Check each entry in the "DataSplit" list.
        """

        self.assertListEqual(
            self.ss["DataSplit"],
            [{'sample': '1234',
              'project': 'DRT',
              'filename': '1234_1234_S1',
              'tags': 'N501-N701',
              'index1': 'GTGTTGCT',
              'index2': 'CATGGTCT',
              'sample_number': 'S1',
              'chemistry': 'Chemistry',
              'orig_sample_name': '1234_1234',
              'Sample_ID': '1',
              },
             {'sample': '1234',
              'project': 'INT',
              'filename': '1234_1234_S1',
              'tags': 'N501-N701',
              'index1': 'GTGTTGCT',
              'index2': 'CATGGTCT',
              'sample_number': 'S1',
              'chemistry': 'Chemistry',
              'orig_sample_name': '1234_1234',
              'Sample_ID': '1',
              },
             {'sample': '4321',
              'project': 'DRT',
              'filename': '4321_S2',
              'tags': 'N501-N702',
              'index1': 'GTGTTGCA',
              'index2': 'CATTGTCA',
              'sample_number': 'S2',
              'chemistry': 'Chemistry',
              'orig_sample_name': '4321',
              'Sample_ID': '2',
              },
             ]
        )


class SampleSheetV2VerifierTests(unittest.TestCase):
    def test_valid_sample_sheet(self):
        valid_sample_sheet = """
[Header]
IEMFileVersion,5
Investigator Name,JN
Project Name,20-Jul-2017.M04401
Experiment Name,20-Jul-2017.M04401
Date,07/31/2024
Workflow,Workflow
Assay,Chemistry
Description,Chemistry
Chemistry,Chemistry
[Reads]
251
251
[Settings]
[Data]
Sample_ID,Sample_Name,Sample_Plate,Sample_Well,index,index2,GenomeFolder,Sample_Project
1,1234_1234,20-Jul-2017.M04401,N/A,GTGTTGCT,CATGGTCT,"",bccfe_2_N501-N701_DRT_INT
2,4321,20-Jul-2017.M04401,N/A,GTGTTGCA,CATTGTCA,"",bccfe_2_N501-N702_DRT
"""
        _sample_sheet_parser(StringIO(valid_sample_sheet))

    def test_missing_header(self):
        invalid_sample_sheet = """
[Reads]
251
251
[Data]
Sample_ID,Sample_Name,Sample_Plate,Sample_Well,index,index2,GenomeFolder,Sample_Project
1,1234_1234,20-Jul-2017.M04401,N/A,GTGTTGCT,CATGGTCT,"",bccfe_2_N501-N701_DRT_INT
2,4321,20-Jul-2017.M04401,N/A,GTGTTGCA,CATTGTCA,"",bccfe_2_N501-N702_DRT
"""
        with self.assertRaises(ValueError) as context:
            _sample_sheet_parser(StringIO(invalid_sample_sheet))
        self.assertIn("Missing 'Header' section in the sample sheet.", str(context.exception))

    def test_missing_required_fields_in_data(self):
        invalid_sample_sheet = """
[Header]
IEMFileVersion,5
Investigator Name,JN
Project Name,TestProject
Experiment Name,TestExperiment
Date,01/01/2021
Workflow,GenerateFASTQ
Assay,Nextera
Description,TestDescription
Chemistry,Amplicon
[Reads]
251
251
[Data]
Sample_Name,Sample_Plate,Sample_Well,index,index2,GenomeFolder,Sample_Project
1234_1234,20-Jul-2017.M04401,N/A,GTGTTGCT,CATGGTCT,"",bccfe_2_N501-N701_DRT_INT
4321,20-Jul-2017.M04401,N/A,GTGTTGCA,CATTGTCA,"",bccfe_2_N501-N702_DRT
"""
        with self.assertRaises(ValueError) as context:
            _sample_sheet_parser(StringIO(invalid_sample_sheet))
        self.assertIn("Expected field 'Sample_ID' not found", str(context.exception))

    def test_row_length_mismatch_in_data(self):
        invalid_sample_sheet = """
[Header]
IEMFileVersion,5
Investigator Name,JN
Project Name,TestProject
Experiment Name,TestExperiment
Date,01/01/2021
Workflow,GenerateFASTQ
Assay,Nextera
Description,TestDescription
Chemistry,Amplicon
[Reads]
251
251
[Data]
Sample_Name,index,index2,Sample_Project
ACGT,TGCA,bccfe_2_N501-N701_INT
CGAT,ATGC,bccfe_2_N501-N701_INT
"""
        with self.assertRaises(ValueError) as context:
            _sample_sheet_parser(StringIO(invalid_sample_sheet))
        self.assertIn("Row length 3 does not match header length 4", str(context.exception))

    def test_row_length_mismatch_in_bccfe_data(self):
        invalid_sample_sheet = """
[Header]
IEMFileVersion,5
Investigator Name,JN
Project Name,TestProject
Experiment Name,TestExperiment
Date,01/01/2021
Workflow,GenerateFASTQ
Assay,Nextera
Description,TestDescription
Chemistry,Amplicon
[Reads]
251
251
[Data]
Sample_ID,Sample_Name,index,index2,Sample_Project
Sample1,ACGT,TGCA,bccfe_2_N501-N701_INT
Sample2,CGAT,ATGC,bccfe_2_N501-N701_INT
"""
        with self.assertRaises(ValueError) as context:
            _sample_sheet_parser(StringIO(invalid_sample_sheet))
        self.assertIn("Row length 4 does not match header length 5", str(context.exception))

    def test_invalid_read_length(self):
        sample_sheet = """
[Header]
IEMFileVersion,5
Investigator Name,JN
Project Name,20-Jul-2017.M04401
Experiment Name,20-Jul-2017.M04401
Date,07/31/2024
Workflow,Workflow
Assay,Chemistry
Description,Chemistry
Chemistry,Chemistry
[Reads]
251
abc
[Settings]
[Data]
Sample_ID,Sample_Name,Sample_Plate,Sample_Well,index,index2,GenomeFolder,Sample_Project
1,1234_1234,20-Jul-2017.M04401,N/A,GTGTTGCT,CATGGTCT,"",bccfe_2_N501-N701_DRT_INT
2,4321,20-Jul-2017.M04401,N/A,GTGTTGCA,CATTGTCA,"",bccfe_2_N501-N702_DRT
"""
        with self.assertRaises(ValueError) as context:
            _sample_sheet_parser(StringIO(sample_sheet))
        self.assertIn("Expected an integer, but got", str(context.exception))

    def test_missing_reads_section(self):
        sample_sheet = """
[Header]
IEMFileVersion,5
Investigator Name,JN
Project Name,TestProject
Experiment Name,TestExperiment
Date,01/01/2021
Workflow,GenerateFASTQ
Assay,Nextera
Description,TestDescription
Chemistry,Amplicon
[Data]
Sample_ID,Sample_Name,Sample_Plate,Sample_Well,index,index2,GenomeFolder,Sample_Project
1,1234_1234,20-Jul-2017.M04401,N/A,GTGTTGCT,CATGGTCT,"",bccfe_2_N501-N701_DRT_INT
2,4321,20-Jul-2017.M04401,N/A,GTGTTGCA,CATTGTCA,"",bccfe_2_N501-N702_DRT
"""
        with self.assertRaises(ValueError) as context:
            _sample_sheet_parser(StringIO(sample_sheet))
        self.assertIn("Missing 'Reads' section in the sample sheet.", str(context.exception))

    def test_not_allowed_missing_index2(self):
        sample_sheet = """
[Header]
IEMFileVersion,5
Investigator Name,JN
Project Name,TestProject
Experiment Name,TestExperiment
Date,01/01/2021
Workflow,GenerateFASTQ
Assay,Nextera
Description,TestDescription
Chemistry,Amplicon
[Reads]
251
251
[Data]
Sample_ID,Sample_Name,Sample_Plate,Sample_Well,index,GenomeFolder,Sample_Project
1,1234_1234,20-Jul-2017.M04401,N/A,GTGTTGCT,"",bccfe_2_N501-N701_DRT_INT
2,4321,20-Jul-2017.M04401,N/A,GTGTTGCA,"",bccfe_2_N501-N702_DRT
"""

        with self.assertRaises(ValueError) as context:
            _sample_sheet_parser(StringIO(sample_sheet))
        self.assertIn("Expected field 'index2' not found", str(context.exception))


class EdgeVersionTwoTests(unittest.TestCase):

    def test_with_optional_fields(self):
        # Valid sample sheet with optional fields like index3 present.
        sample_sheet = """
[Header]
IEMFileVersion,5
Investigator Name,JN
Project Name,CompleteBCCFETest
Experiment Name,CompleteBCCFEExperiment
Date,01/01/2024
Workflow,GenerateFASTQ
Assay,Nextera
Description,CompleteBCCFEDescription
Chemistry,Amplicon
[Reads]
251
251
[Data]
Sample_ID,Sample_Name,Sample_Plate,Sample_Well,index,index2,index3,GenomeFolder,Sample_Project
1,Sample1,20-Jul-2017.M04401,N/A,ACTG,TGCA,ACCC,"",bccfe_2_Tag1_Neurology
2,Sample2,20-Jul-2017.M04401,N/A,CGAT,ATGC,TCCC,"",bccfe_2_Tag2_Oncology
"""
        ss = _sample_sheet_parser(StringIO(sample_sheet))

        self.assertIn("Sample1_S1", ss["Data"])
        self.assertEqual(ss["Data"]["Sample1_S1"]["index2"], "TGCA")

        self.assertIn("Sample2_S2", ss["Data"])
        self.assertEqual(ss["Data"]["Sample2_S2"]["index2"], "ATGC")

    def test_complete_bccfe_data_entries(self):
        # Valid sample sheet with all elements in BCCFE_Data present and correctly parsed.
        sample_sheet = """
[Header]
IEMFileVersion,5
Investigator Name,JN
Project Name,CompleteBCCFETest
Experiment Name,CompleteBCCFEExperiment
Date,01/01/2024
Workflow,GenerateFASTQ
Assay,Nextera
Description,CompleteBCCFEDescription
Chemistry,Amplicon
[Reads]
251
251
[Data]
Sample_ID,Sample_Name,Sample_Plate,Sample_Well,index,index2,GenomeFolder,Sample_Project
1,Sample1,20-Jul-2017.M04401,N/A,ACTG,TGCA,"",bccfe_2_Tag1_Neurology
2,Sample2,20-Jul-2017.M04401,N/A,CGAT,ATGC,"",bccfe_2_Tag2_Oncology
"""
        ss = _sample_sheet_parser(StringIO(sample_sheet))

        self.assertEqual(ss["Data"]["Sample1_S1"]["tags"], "Tag1")
        self.assertEqual(ss["Data"]["Sample2_S2"]["tags"], "Tag2")

        self.assertEqual(len(ss["DataSplit"]), 2)
        self.assertEqual(ss["DataSplit"][0]["project"], "Neurology")
        self.assertEqual(ss["DataSplit"][1]["project"], "Oncology")

    def test_complete_bccfe_data_entries_with_spaces(self):
        # Sample that contains empty lines
        sample_sheet = """
[Header]
IEMFileVersion,5
Investigator Name,JN
Project Name,CompleteBCCFETest
Experiment Name,CompleteBCCFEExperiment
Date,01/01/2024
Workflow,GenerateFASTQ
Assay,Nextera
Description,CompleteBCCFEDescription
Chemistry,Amplicon

[Reads]
251
251

[Settings]

[Data]
Sample_ID,Sample_Name,Sample_Plate,Sample_Well,index,index2,GenomeFolder,Sample_Project
1,Sample1,20-Jul-2017.M04401,N/A,ACTG,TGCA,"",bccfe_2_Tag1_Neurology
2,Sample2,20-Jul-2017.M04401,N/A,CGAT,ATGC,"",bccfe_2_Tag2_Oncology
"""
        ss = _sample_sheet_parser(StringIO(sample_sheet))

        self.assertEqual(ss["Data"]["Sample1_S1"]["tags"], "Tag1")
        self.assertEqual(ss["Data"]["Sample2_S2"]["tags"], "Tag2")

        self.assertEqual(len(ss["DataSplit"]), 2)
        self.assertEqual(ss["DataSplit"][0]["project"], "Neurology")
        self.assertEqual(ss["DataSplit"][1]["project"], "Oncology")
