from pathlib import Path

from micall.monitor.find_groups import find_groups, SampleGroup

BASIC_HEADER = """\
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
"""


def test_two_samples(tmpdir):
    sample_sheet_path = Path(tmpdir) / 'SampleSheet.csv'
    sample_sheet_path.write_text(BASIC_HEADER + """\
Sample_ID,Sample_Name,Sample_Plate,Sample_Well,index,index2,Sample_Project,Description,GenomeFolder
CFE_SomeId_10-Jul-2014_N501-N701_Sample1_Proj1,Sample1_Proj1,,,ACGTACGT,TGCATGCA,,,
CFE_SomeId_10-Jul-2014_N501-N702_Sample2_Proj2,Sample2_Proj2,,,AAAAGGGG,CCCCTTTT,,,
""")
    expected_groups = [SampleGroup('Sample1',
                                   ('Sample1-Proj1_S1_L001_R1_001.fastq.gz',
                                    None)),
                       SampleGroup('Sample2',
                                   ('Sample2-Proj2_S2_L001_R1_001.fastq.gz',
                                    None))]

    groups = list(find_groups(['Sample1-Proj1_S1_L001_R1_001.fastq.gz',
                               'Sample2-Proj2_S2_L001_R1_001.fastq.gz'],
                              sample_sheet_path))

    assert expected_groups == groups


def test_combine_midi(tmpdir):
    sample_sheet_path = Path(tmpdir) / 'SampleSheet.csv'
    sample_sheet_path.write_text(BASIC_HEADER + """\
Sample_ID,Sample_Name,Sample_Plate,Sample_Well,index,index2,Sample_Project,Description,GenomeFolder
CFE_SomeId_10-Jul-2014_N501-N701_Sample1_HCV,Sample1_HCV,,,ACGTACGT,TGCATGCA,,,
CFE_SomeId_10-Jul-2014_N501-N702_Sample1MIDI_MidHCV,Sample1MIDI_MidHCV,,,AAAAGGGG,CCCCTTTT,,,
""")
    expected_groups = [SampleGroup('Sample1',
                                   ('Sample1-HCV_S1_L001_R1_001.fastq.gz',
                                    'Sample1MIDI-MidHCV_S2_L001_R1_001.fastq.gz'))]

    groups = list(find_groups(['Sample1-HCV_S1_L001_R1_001.fastq.gz',
                               'Sample1MIDI-MidHCV_S2_L001_R1_001.fastq.gz'],
                              sample_sheet_path))

    assert expected_groups == groups


def test_midi_unmatched(tmpdir):
    sample_sheet_path = Path(tmpdir) / 'SampleSheet.csv'
    sample_sheet_path.write_text(BASIC_HEADER + """\
Sample_ID,Sample_Name,Sample_Plate,Sample_Well,index,index2,Sample_Project,Description,GenomeFolder
CFE_SomeId_10-Jul-2014_N501-N701_Sample1MIDI_MidHCV,Sample1MIDI_MidHCV,,,ACGTACGT,TGCATGCA,,,
""")
    expected_groups = [SampleGroup('Sample1MIDI',
                                   ('Sample1MIDI-MidHCV_S1_L001_R1_001.fastq.gz',
                                    None))]

    groups = list(find_groups(['Sample1MIDI-MidHCV_S1_L001_R1_001.fastq.gz'],
                              sample_sheet_path))

    assert expected_groups == groups


def test_unmatched_midi_file_not_found(tmpdir):
    sample_sheet_path = Path(tmpdir) / 'SampleSheet.csv'
    sample_sheet_path.write_text(BASIC_HEADER + """\
Sample_ID,Sample_Name,Sample_Plate,Sample_Well,index,index2,Sample_Project,Description,GenomeFolder
CFE_SomeId_10-Jul-2014_N501-N701_Sample1MIDI_MidHCV,Sample1MIDI_MidHCV,,,ACGTACGT,TGCATGCA,,,
CFE_SomeId_10-Jul-2014_N501-N702_Sample2_Proj2,Sample2_Proj2,,,AAAAGGGG,CCCCTTTT,,,
""")
    expected_groups = [SampleGroup('Sample2',
                                   ('Sample2-Proj2_S2_L001_R1_001.fastq.gz',
                                    None))]

    groups = list(find_groups(['Sample2-Proj2_S2_L001_R1_001.fastq.gz'],
                              sample_sheet_path))

    assert expected_groups == groups
