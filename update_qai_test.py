import unittest
import StringIO

import update_qai
from sample_sheet_parser import sample_sheet_parser

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
CFE_SomeId_10-Jul-2014_N501-N702_Sample2_Proj2,Sample2_Proj2,10-Jul-2014_testing,N/A,AAAAGGGG,CCCCTTTT,10-Jul-2014_testing,Research:Sample2_Proj2:FALSE Comments:Sample2_Proj2:thisiscommenttwo Chemistry:Sample2_Proj2:BreakingBad Disablecontamcheck:Sample2_Proj2:TRUE,
"""))
        self.collated_counts_file = StringIO.StringIO("""\
sample_name,type,count
Sample1-Proj1_S1,raw,3000
Sample1-Proj1_S1,prelim R1-seed,1500
Sample1-Proj1_S1,remap R1-seed,2000
Sample1-Proj1_S1,prelim R2-seed,150
Sample1-Proj1_S1,remap R2-seed,200
Sample1-Proj1_S1,unmapped,800
""")
        self.sequencings = [{'tag': 'N501-N701',
                             'id': 10001,
                             'target_project_regions': {'R1': {'id': 20001,
                                                               'seed': 'R1-seed'}}},
                            {'tag': 'N501-N702',
                             'id': 10002,
                             'target_project_regions': {'R2': {'id': 20002,
                                                               'seed': 'R2-seed'}}}]
        self.region_defaults = [{'name': 'R1',
                                 'default_project_region': {'id': 30001,
                                                            'seed': 'R1-seed'}},
                                {'name': 'R2',
                                 'default_project_region': {'id': 30002,
                                                            'seed': 'R2-seed'}}]
        
    def test_on_target(self):
        coverage_file = StringIO.StringIO("""\
sample,region,q.cut,min.coverage,which.key.pos,off.score,on.score
Sample1-Proj1_S1,R1,15,2000,12,-3,4
""")
        expected_decisions = [{'sequencing_id': 10001,
                               'project_region_id': 20001,
                               'score': 4,
                               'min_coverage': 2000,
                               'min_coverage_pos': 12,
                               'raw_reads': 3000,
                               'mapped_reads': 2000}]
        
        decisions = update_qai.build_review_decisions(coverage_file,
                                                      self.collated_counts_file,
                                                      self.sample_sheet,
                                                      self.sequencings,
                                                      self.region_defaults)

        self.assertListEqual(expected_decisions, decisions)
        
    def test_off_target(self):
        coverage_file = StringIO.StringIO("""\
sample,region,q.cut,min.coverage,which.key.pos,off.score,on.score
Sample1-Proj1_S1,R2,15,2000,12,-3,4
""")
        expected_decisions = [{'sequencing_id': 10001,
                               'project_region_id': 30002,
                               'score': -3,
                               'min_coverage': 2000,
                               'min_coverage_pos': 12,
                               'raw_reads': 3000,
                               'mapped_reads': 200}]
        
        decisions = update_qai.build_review_decisions(coverage_file,
                                                      self.collated_counts_file,
                                                      self.sample_sheet,
                                                      self.sequencings,
                                                      self.region_defaults)

        self.assertListEqual(expected_decisions, decisions)

#TODO: include seed regions in run sequencing summary and region defaults
#TODO: test multiple projects with same tags