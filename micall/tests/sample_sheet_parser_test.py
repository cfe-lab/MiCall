import unittest
import StringIO
from micall.utils.sample_sheet_parser import sample_sheet_parser

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
CFE_SomeId_10-Jul-2014_N501-N701_Sample1_Proj1,Sample1_Proj1,10-Jul-2014_testing,N/A,ACGTACGT,TGCATGCA,10-Jul-2014_testing,Research:Sample1_Proj1:TRUE Comments:Sample1_Proj1:thisiscommentone Disablecontamcheck:Sample1_Proj1:FALSE,
CFE_SomeId_10-Jul-2014_N501-N702_Sample2_Proj2,Sample2_Proj2,10-Jul-2014_testing,N/A,AAAAGGGG,CCCCTTTT,10-Jul-2014_testing,Research:Sample2_Proj2:FALSE Comments:Sample2_Proj2:thisiscommenttwo Chemistry:Sample2_Proj2:BreakingBad Disablecontamcheck:Sample2_Proj2:TRUE,
"""
    clean_filenames = ["Sample1-Proj1_S1", "Sample2-Proj2_S2"]

    def setUp(self):
        self.maxDiff=None
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
             "Assay", "Description", "Chemistry", "Data", "DataSplit", "sample_sheet_version",
             "Reads"})

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
        self.assertEquals(self.ss["Reads"], [251, 251])

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


class VersionTwoTest(unittest.TestCase):
    stub_sample_sheet = """
[Header]
IEMFileVersion,3
Investigator Name,RL
Project Name,11-Jul-2014_v2test
Experiment Name,11-Jul-2014_v2test
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
CFE_SomeId_11-Jul-2014_N501-N703_Sample3--Proj3,Sample3--Proj3,11-Jul-2014_testing,N/A,AAAAAAAA,CCCCCCCC,11-Jul-2014_testing,Research:Sample3--Proj3--TRUE Comments:Sample3--Proj3--thisiscommentthree Disablecontamcheck:Sample3--Proj3--FALSE,
CFE_SomeId_11-Jul-2014_N501-N704_Sample4--Proj4,Sample4--Proj4,11-Jul-2014_testing,N/A,GGGGGGGG,TTTTTTTT,11-Jul-2014_testing,Research:Sample4--Proj4--FALSE Comments:Sample4--Proj4--thisiscommentfour Chemistry:Sample4--Proj4--BreakingBad Disablecontamcheck:Sample4--Proj4--TRUE,
"""
    clean_filenames = ["Sample3--Proj3_S1", "Sample4--Proj4_S2"]

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
             "Assay", "Description", "Chemistry", "Data", "DataSplit", "sample_sheet_version",
             "Reads"})

    def test_preamble_correct(self):
        """
        Test that all of the header stuff, as well as the sample sheet version, was set correctly.
        """
        self.assertEquals(self.ss["IEMFileVersion"], "3")
        self.assertEquals(self.ss["Investigator Name"], "RL")
        self.assertEquals(self.ss["Project Name"], "11-Jul-2014_v2test")
        self.assertEquals(self.ss["Experiment Name"], "11-Jul-2014_v2test")
        self.assertEquals(self.ss["Date"], "07/11/2014")
        self.assertEquals(self.ss["Workflow"], "GenerateFASTQ")
        self.assertEquals(self.ss["Assay"], "Nextera")
        self.assertEquals(self.ss["Description"], "Nextera")
        self.assertEquals(self.ss["Chemistry"], "Amplicon")
        self.assertEquals(self.ss["sample_sheet_version"], 2)

    def test_data(self):
        """
        Check each entry in the "Data" dictionary.
        """
        self.assertDictEqual(
            self.ss["Data"],
            {
                self.clean_filenames[0]:
                    {
                        "index1": "AAAAAAAA",
                        "index2": "CCCCCCCC",
                        "tags": "N501-N703",
                        "comments": "thisiscommentthree",
                        "disable_contamination_check": False,
                        "research": True,
                        "chemistry": "Nextera",
                        "orig_sample_name": "Sample3--Proj3",
                        "is_T_primer": False
                    },
                self.clean_filenames[1]:
                    {
                        "index1": "GGGGGGGG",
                        "index2": "TTTTTTTT",
                        "tags": "N501-N704",
                        "comments": "thisiscommentfour",
                        "disable_contamination_check": True,
                        "research": False,
                        "chemistry": "BreakingBad",
                        "orig_sample_name": "Sample4--Proj4",
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
                    "sample": "Sample3",
                    "project": "Proj3",
                    "filename": "Sample3--Proj3_S1",
                    "index1": "AAAAAAAA",
                    "index2": "CCCCCCCC",
                    "sample_number": "S1",
                    "chemistry": "Nextera",
                    "disable_contamination_check": False,
                    "research": True,
                    "comments": "thisiscommentthree",
                    "orig_sample_name": "Sample3--Proj3",
                    "tags": "N501-N703",
                    "is_T_primer": False,
                },
                {
                    "sample": "Sample4",
                    "project": "Proj4",
                    "filename": "Sample4--Proj4_S2",
                    "index1": "GGGGGGGG",
                    "index2": "TTTTTTTT",
                    "sample_number": "S2",
                    "chemistry": "BreakingBad",
                    "disable_contamination_check": True,
                    "research": False,
                    "comments": "thisiscommentfour",
                    "orig_sample_name": "Sample4--Proj4",
                    "tags": "N501-N704",
                    "is_T_primer": False
                }
            ]
        )


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
CFE_SomeId_11-Jul-2014_N501-N701_Sample1_Proj1;CFE_SomeId_11-Jul-2014_N501-N701_Sample2_Proj2,Sample1_Proj1;Sample2_Proj2,11-Jul-2014_testing,N/A,ACGTACGT,TGCATGCA,11-Jul-2014_testing,Research:Sample1_Proj1:TRUE;Sample2_Proj2:FALSE Comments:Sample1_Proj1:thisiscommentone;Sample2_Proj2:thisiscommenttwo Disablecontamcheck:Sample1_Proj1:FALSE;Sample2_Proj2:TRUE,
CFE_SomeId_11-Jul-2014_N501-N702_Sample3_Proj3;CFE_SomeId_11-Jul-2014_N501-N702_Sample4_Proj4,Sample3_Proj3;Sample4_Proj4,11-Jul-2014_testing,N/A,AAAAGGGG,CCCCTTTT,11-Jul-2014_testing,Research:Sample3_Proj3:FALSE;Sample4_Proj4:FALSE Comments:Sample3_Proj3:thisiscommentthree;Sample4_Proj4:thisiscommentfour Chemistry:Sample3_Proj3:BreakingBad;Sample4_Proj4:MrWizard Disablecontamcheck:Sample3_Proj3:TRUE;Sample4_Proj4:TRUE,
"""
    clean_filenames = ["Sample1-Proj1-Sample2-Proj2_S1", "Sample3-Proj3-Sample4-Proj4_S2"]

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
             "Assay", "Description", "Chemistry", "Data", "DataSplit", "sample_sheet_version",
             "Reads"})

    def test_preamble_correct(self):
        """
        Test that all of the header stuff, as well as the sample sheet version, was set correctly.
        """
        self.assertEquals(self.ss["IEMFileVersion"], "3")
        self.assertEquals(self.ss["Investigator Name"], "RL")
        self.assertEquals(self.ss["Project Name"], "11-Jul-2014_v1multisamplestest")
        self.assertEquals(self.ss["Experiment Name"], "11-Jul-2014_v1multisamplestest")
        self.assertEquals(self.ss["Date"], "07/11/2014")
        self.assertEquals(self.ss["Workflow"], "GenerateFASTQ")
        self.assertEquals(self.ss["Assay"], "Nextera")
        self.assertEquals(self.ss["Description"], "Nextera")
        self.assertEquals(self.ss["Chemistry"], "Amplicon")
        self.assertEquals(self.ss["sample_sheet_version"], 1)

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


class MultipleSamplesV2Test(unittest.TestCase):
    stub_sample_sheet = """
[Header]
IEMFileVersion,3
Investigator Name,RL
Project Name,11-Jul-2014_v2multisamplestest
Experiment Name,11-Jul-2014_v2multisamplestest
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
CFE_SomeId_11-Jul-2014_N501-N701_Sample1--Proj1---CFE_SomeId_11-Jul-2014_N501-N701_Sample2--Proj2,Sample1--Proj1---Sample2--Proj2,11-Jul-2014_testing,N/A,ACGTACGT,TGCATGCA,11-Jul-2014_testing,Research:Sample1--Proj1--TRUE---Sample2--Proj2--FALSE Comments:Sample1--Proj1--thisiscommentone---Sample2--Proj2--thisiscommenttwo Disablecontamcheck:Sample1--Proj1--FALSE---Sample2--Proj2--TRUE,
CFE_SomeId_11-Jul-2014_N501-N702_Sample3--Proj3---CFE_SomeId_11-Jul-2014_N501-N702_Sample4--Proj4,Sample3--Proj3---Sample4--Proj4,11-Jul-2014_testing,N/A,AAAAGGGG,CCCCTTTT,11-Jul-2014_testing,Research:Sample3--Proj3--FALSE---Sample4--Proj4--FALSE Comments:Sample3--Proj3--thisiscommentthree---Sample4--Proj4--thisiscommentfour Chemistry:Sample3--Proj3--BreakingBad---Sample4--Proj4--MrWizard Disablecontamcheck:Sample3--Proj3--TRUE---Sample4--Proj4--TRUE,
"""
    clean_filenames = ["Sample1--Proj1---Sample2--Proj2_S1", "Sample3--Proj3---Sample4--Proj4_S2"]

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
             "Assay", "Description", "Chemistry", "Data", "DataSplit", "sample_sheet_version",
             "Reads"})

    def test_preamble_correct(self):
        """
        Test that all of the header stuff, as well as the sample sheet version, was set correctly.
        """
        self.assertEquals(self.ss["IEMFileVersion"], "3")
        self.assertEquals(self.ss["Investigator Name"], "RL")
        self.assertEquals(self.ss["Project Name"], "11-Jul-2014_v2multisamplestest")
        self.assertEquals(self.ss["Experiment Name"], "11-Jul-2014_v2multisamplestest")
        self.assertEquals(self.ss["Date"], "07/11/2014")
        self.assertEquals(self.ss["Workflow"], "GenerateFASTQ")
        self.assertEquals(self.ss["Assay"], "Nextera")
        self.assertEquals(self.ss["Description"], "Nextera")
        self.assertEquals(self.ss["Chemistry"], "Amplicon")
        self.assertEquals(self.ss["sample_sheet_version"], 2)

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
                        "orig_sample_name": "Sample1--Proj1---Sample2--Proj2",
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
                        "orig_sample_name": "Sample3--Proj3---Sample4--Proj4",
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
                    "filename": "Sample1--Proj1---Sample2--Proj2_S1",
                    "index1": "ACGTACGT",
                    "index2": "TGCATGCA",
                    "sample_number": "S1",
                    "chemistry": "Nextera",
                    "disable_contamination_check": False,
                    "research": True,
                    "comments": "thisiscommentone",
                    "orig_sample_name": "Sample1--Proj1---Sample2--Proj2",
                    "tags": "N501-N701",
                    "is_T_primer": False,
                },
                {
                    "sample": "Sample2",
                    "project": "Proj2",
                    "filename": "Sample1--Proj1---Sample2--Proj2_S1",
                    "index1": "ACGTACGT",
                    "index2": "TGCATGCA",
                    "sample_number": "S1",
                    "chemistry": "Nextera",
                    "disable_contamination_check": True,
                    "research": False,
                    "comments": "thisiscommenttwo",
                    "orig_sample_name": "Sample1--Proj1---Sample2--Proj2",
                    "tags": "N501-N701",
                    "is_T_primer": False,
                },
                {
                    "sample": "Sample3",
                    "project": "Proj3",
                    "filename": "Sample3--Proj3---Sample4--Proj4_S2",
                    "index1": "AAAAGGGG",
                    "index2": "CCCCTTTT",
                    "sample_number": "S2",
                    "chemistry": "BreakingBad",
                    "disable_contamination_check": True,
                    "research": False,
                    "comments": "thisiscommentthree",
                    "orig_sample_name": "Sample3--Proj3---Sample4--Proj4",
                    "tags": "N501-N702",
                    "is_T_primer": False
                },
                {
                    "sample": "Sample4",
                    "project": "Proj4",
                    "filename": "Sample3--Proj3---Sample4--Proj4_S2",
                    "index1": "AAAAGGGG",
                    "index2": "CCCCTTTT",
                    "sample_number": "S2",
                    "chemistry": "MrWizard",
                    "disable_contamination_check": True,
                    "research": False,
                    "comments": "thisiscommentfour",
                    "orig_sample_name": "Sample3--Proj3---Sample4--Proj4",
                    "tags": "N501-N702",
                    "is_T_primer": False
                }
            ]
        )


class OtherTest(unittest.TestCase):
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

        self.assertRaisesRegexp(
            ValueError,
            "sample sheet data header does not include Sample_Name",
            lambda: sample_sheet_parser(StringIO.StringIO(stub_sample_sheet))
        )

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
CFE_SomeId_10-Jul-2014_N501-N701_Sample1_Proj1,Sample1_Proj1,10-Jul-2014_testing,N/A,ACGTACGT,TGCATGCA,10-Jul-2014_testing,Research:Sample1_Proj1:TRUE Comments:Sample1_Proj1:thisiscommentone Disablecontamcheck:Sample1_Proj1:FALSE,
CFE_SomeId_10-Jul-2014_N501-N702_Sample2_Proj2,Sample2_Proj2,10-Jul-2014_testing,N/A,AAAAGGGG,CCCCTTTT,10-Jul-2014_testing,Research:Sample2_Proj2:FALSE Comments:Sample2_Proj2:thisiscommenttwo Chemistry:Sample2_Proj2:BreakingBad Disablecontamcheck:Sample2_Proj2:TRUE,
"""
 
        ss = sample_sheet_parser(StringIO.StringIO(stub_sample_sheet))
        self.assertEquals(ss["Experiment Name"], "10-Jul-2014")
