import unittest
import StringIO

import update_qai
from sample_sheet_parser import sample_sheet_parser
import json

class UpdateQaiTest(unittest.TestCase):
    def setUp(self):
        self.sample_sheet = sample_sheet_parser(StringIO.StringIO("""\
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
CFE_SomeId_10-Jul-2014_N501-N702_Sample2_Proj2;CFE_SomeId_10-Jul-2014_N501-N702_Sample3_Proj3,Sample2_Proj2;Sample3_Proj3,10-Jul-2014,N/A,AAAAGGGG,CCCCTTTT,10-Jul-2014,Research:Sample2_Proj2:FALSE;Sample3_Proj3:FALSE Comments:Sample2_Proj2:thisiscommenttwo;Sample3_Proj3:thisiscommentthree Disablecontamcheck:Sample2_Proj2:TRUE;Sample3_Proj3:TRUE,
"""))
        self.collated_counts_file = StringIO.StringIO("""\
sample_name,type,count
Sample1-Proj1_S1,raw,3000
Sample1-Proj1_S1,prelim R1-seed,1500
Sample1-Proj1_S1,remap R1-seed,2000
Sample1-Proj1_S1,prelim R2-seed,150
Sample1-Proj1_S1,remap R2-seed,200
Sample1-Proj1_S1,unmapped,800
Sample2-Proj2-Sample3-Proj3_S2,raw,3000
Sample2-Proj2-Sample3-Proj3_S2,prelim R2-seed,1500
Sample2-Proj2-Sample3-Proj3_S2,remap R2-seed,2000
Sample2-Proj2-Sample3-Proj3_S2,prelim R3-seed,150
Sample2-Proj2-Sample3-Proj3_S2,remap R3-seed,200
Sample2-Proj2-Sample3-Proj3_S2,unmapped,800
""")
        self.sequencings = [{'tag': 'N501-N701',
                             'id': 10001,
                             'target_project': 'Proj1'},
                            {'tag': 'N501-N702',
                             'id': 10002,
                             'target_project': 'Proj2'},
                            {'tag': 'N501-N702',
                             'id': 10003,
                             'target_project': 'Proj3'}]
        projectRegionsJson = """\
[
    {
        "project_name": "Proj1",
        "id": 20001,
        "seed_region_names": ["R1-seed"],
        "coordinate_region_name": "R1"
    },
    {
        "project_name": "Proj2",
        "id": 20002,
        "seed_region_names": ["R2-seed"],
        "coordinate_region_name": "R2"
    },
    {
        "project_name": "Proj3",
        "id": 20003,
        "seed_region_names": ["R3-seed"],
        "coordinate_region_name": "R3"
    },
    {
        "project_name": "Proj-all",
        "id": 20011,
        "seed_region_names": ["R1-seed"],
        "coordinate_region_name": "R1"
    },
    {
        "project_name": "Proj-all",
        "id": 20012,
        "seed_region_names": ["R2-seed"],
        "coordinate_region_name": "R2"
    },
    {
        "project_name": "Proj-all",
        "id": 20013,
        "seed_region_names": ["R3-seed"],
        "coordinate_region_name": "R3"
    }
]
"""
        self.project_regions = json.loads(projectRegionsJson)
        regionsJson = """\
[
    {
        "name": "R1",
        "id": 30001
    },
    {
        "name": "R1-seed",
        "id": 30011
    },
    {
        "name": "R2",
        "id": 30002
    },
    {
        "name": "R2-seed",
        "id": 30012
    },
    {
        "name": "R3",
        "id": 30003
    },
    {
        "name": "R3-seed",
        "id": 30013
    }
]
"""
        self.regions = json.loads(regionsJson)
                
    def test_on_target(self):
        coverage_file = StringIO.StringIO("""\
sample,project,region,seed,q.cut,min.coverage,which.key.pos,off.score,on.score
Sample1-Proj1_S1,Proj1,R1,R1-seed,15,1900,12,-3,4
""")
        expected_decisions = [{'sequencing_id': 10001,
                               'project_region_id': 20001,
                               'seed_region_id': 30011,
                               'sample_name': 'Sample1-Proj1_S1',
                               'score': 4,
                               'min_coverage': 1900,
                               'min_coverage_pos': 12,
                               'raw_reads': 3000,
                               'mapped_reads': 2000}]
        
        decisions = update_qai.build_review_decisions(coverage_file,
                                                      self.collated_counts_file,
                                                      self.sample_sheet,
                                                      self.sequencings,
                                                      self.project_regions,
                                                      self.regions)

        self.assertListEqual(expected_decisions, decisions)
        
    def test_off_target(self):
        coverage_file = StringIO.StringIO("""\
sample,project,region,seed,q.cut,min.coverage,which.key.pos,off.score,on.score
Sample1-Proj1_S1,Proj2,R2,R2-seed,15,1900,12,-3,4
""")
        expected_decisions = [{'sequencing_id': 10001,
                               'project_region_id': 20002,
                               'seed_region_id': 30012,
                               'sample_name': 'Sample1-Proj1_S1',
                               'score': -3,
                               'min_coverage': 1900,
                               'min_coverage_pos': 12,
                               'raw_reads': 3000,
                               'mapped_reads': 200}]
         
        decisions = update_qai.build_review_decisions(coverage_file,
                                                      self.collated_counts_file,
                                                      self.sample_sheet,
                                                      self.sequencings,
                                                      self.project_regions,
                                                      self.regions)
 
        self.assertListEqual(expected_decisions, decisions)
 
    def test_multiple_targets(self):
        coverage_file = StringIO.StringIO("""\
sample,project,region,seed,q.cut,min.coverage,which.key.pos,off.score,on.score
Sample2-Proj2-Sample3-Proj3_S2,Proj2,R2,R2-seed,15,1900,12,-3,4
Sample2-Proj2-Sample3-Proj3_S2,Proj3,R3,R3-seed,15,190,20,-2,3
""")
        expected_decisions = [{'sequencing_id': 10002,
                               'project_region_id': 20002,
                               'seed_region_id': 30012,
                               'sample_name': 'Sample2-Proj2-Sample3-Proj3_S2',
                               'score': 4,
                               'min_coverage': 1900,
                               'min_coverage_pos': 12,
                               'raw_reads': 3000,
                               'mapped_reads': 2000},
                              {'sequencing_id': 10003,
                               'project_region_id': 20003,
                               'seed_region_id': 30013,
                               'sample_name': 'Sample2-Proj2-Sample3-Proj3_S2',
                               'score': 3,
                               'min_coverage': 190,
                               'min_coverage_pos': 20,
                               'raw_reads': 3000,
                               'mapped_reads': 200}]
         
        decisions = update_qai.build_review_decisions(coverage_file,
                                                      self.collated_counts_file,
                                                      self.sample_sheet,
                                                      self.sequencings,
                                                      self.project_regions,
                                                      self.regions)
 
        self.assertListEqual(expected_decisions, decisions)
 
    def test_multiple_projects_one_on_target(self):
        coverage_file = StringIO.StringIO("""\
sample,project,region,seed,q.cut,min.coverage,which.key.pos,off.score,on.score
Sample1-Proj1_S1,Proj-all,R1,R1-seed,15,1900,12,-3,4
Sample1-Proj1_S1,Proj1,R1,R1-seed,15,1900,12,-3,4
""")
        expected_decisions = [{'sequencing_id': 10001,
                               'project_region_id': 20001,
                               'seed_region_id': 30011,
                               'sample_name': 'Sample1-Proj1_S1',
                               'score': 4,
                               'min_coverage': 1900,
                               'min_coverage_pos': 12,
                               'raw_reads': 3000,
                               'mapped_reads': 2000}]
         
        decisions = update_qai.build_review_decisions(coverage_file,
                                                      self.collated_counts_file,
                                                      self.sample_sheet,
                                                      self.sequencings,
                                                      self.project_regions,
                                                      self.regions)
 
        self.assertListEqual(expected_decisions, decisions)
