import unittest
import StringIO
from sample_sheet_parser import sample_sheet_parser

"""
Test suite for the sample sheet parser.

Important cases to test:
 - a version 1 sample sheet
 - a version 2 sample sheet
 - entries (wells) with multiple samples
 - defective file whose Data portion does not include a Sample_Name column

Later: do we need to check
 - "is_T_primer"?
"""

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
CFEidentifier-Sample1_Proj1,Sample1_Proj1,10-Jul-2014_testing,N/A,ACGTACGT,TGCATGCA,10-Jul-2014_testing,Research:Sample1_Proj1:TRUE Comments:Sample1_Proj1:thisiscommentone Disablecontamcheck:Sample1_Proj1:FALSE,
CFEidentifier_Sample2_Proj2,Sample2_Proj2,10-Jul-2014_testing,N/A,AAAAGGGG,CCCCTTTT,10-Jul-2014_testing,Research:Sample2_Proj2:FALSE Comments:Sample2_Proj2:thisiscommenttwo Chemistry:Sample2_Proj2:BreakingBad Disablecontamcheck:Sample2_Proj2:TRUE,
"""
    clean_filenames = ["Sample1-Proj1_S1", "Sample2-Proj2_S2"]

    def setUp(self):
        self.ss = sample_sheet_parser(StringIO.StringIO(self.stub_sample_sheet))

    def test_keys_correct(self):
        """
        Test that all of the keys in the resulting dictionary are correct.

        There should be one key per [Header] field, as well as "Data", "DataSplit",
        and sample_sheet_version.
        """
        self.assertSetEqual(
            set(self.ss.keys()),
            {"IEMFileVersion", "Investigator Name", "Project Name", "Experiment Name", "Date", "Workflow",
             "Assay", "Description", "Chemistry", "Data", "DataSplit", "sample_sheet_version"})

    def test_preamble_correct(self):
        """
        Test that all of the header stuff, as well as the sample sheet version, was set correctly.
        """
        self.assertEquals(self.ss["IEMFileVersion"], "3")
        self.assertEquals(self.ss["Investigator Name"], "RL")
        self.assertEquals(self.ss["Project Name"], "10-Jul-2014_v1test")
        self.assertEquals(self.ss["Experiment Name"], "10-Jul-2014_v1test")
        self.assertEquals(self.ss["Date"], "07/10/2014")
        self.assertEquals(self.ss["Workflow"], "GenerateFASTQ")
        self.assertEquals(self.ss["Assay"], "Nextera")
        self.assertEquals(self.ss["Description"], "Nextera")
        self.assertEquals(self.ss["Chemistry"], "Amplicon")
        self.assertEquals(self.ss["sample_sheet_version"], 1)

    def test_data(self):
        """
        Check each entry in the "Data" dictionary.
        """
        self.assertSetEqual(set(self.ss["Data"].keys()), set(self.clean_filenames))

        self.assertDictEqual(
            self.ss["Data"][self.clean_filenames[0]],
            {
                "index1": "ACGTACGT",
                "index2": "TGCATGCA",
                "comments": "thisiscommentone",
                "disable_contamination_check": False,
                "research": True,
                "chemistry": "Nextera",
                "orig_sample_name": "Sample1_Proj1",
                "is_T_primer": False
            }
        )

        self.assertDictEqual(
            self.ss["Data"][self.clean_filenames[1]],
            {
                "index1": "AAAAGGGG",
                "index2": "CCCCTTTT",
                "comments": "thisiscommenttwo",
                "disable_contamination_check": True,
                "research": False,
                "chemistry": "BreakingBad",
                "orig_sample_name": "Sample2_Proj2",
                "is_T_primer": False
            }
        )

    def test_datasplit(self):
        """
        Check each entry in the "DataSplit" list.
        """
        self.assertEquals(len(self.ss["DataSplit"]), 2)

        self.assertDictEqual(
            self.ss["DataSplit"][0],
            {
                "sample": "Sample1",
                "project": "Proj1",
                "filename": "Sample1-Proj1_S1",
                "index1": "ACGTACGT",
                "index2": "TGCATGCA",
                "sample_number": "S1",
                "chemistry": "Nextera",
                "disable_contamination_check": False,
                "research": True,
                "comments": "thisiscommentone",
                "orig_sample_name": "Sample1_Proj1",
                "is_T_primer": False,
            }
        )

        self.assertDictEqual(
            self.ss["DataSplit"][1],
            {
                "sample": "Sample2",
                "project": "Proj2",
                "filename": "Sample2-Proj2_S2",
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
        )