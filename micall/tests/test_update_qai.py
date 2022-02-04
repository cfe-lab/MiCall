from io import StringIO
from pathlib import Path
from textwrap import dedent

from micall.monitor.update_qai import build_conseqs
from micall.utils.sample_sheet_parser import sample_sheet_parser


def test_build_conseqs():
    sample_sheet_path = Path(__file__).parent / 'microtest' / 'SampleSheet.csv'
    with sample_sheet_path.open() as sample_sheet_file:
        sample_sheet = sample_sheet_parser(sample_sheet_file)
    conseq_csv = StringIO(dedent("""\
        sample,region,q-cutoff,consensus-percent-cutoff,offset,sequence
        2140A-HIV_S17,HIV1-B-FR-K03455-seed,15,MAX,2252,CCTCAGGTC
        2140A-HIV_S17,HIV1-B-FR-K03455-seed,15,0.010,2252,CCTCWGGTC
        2190A-SARSCOV2_S23,SARS-CoV-2-seed,15,MAX,13441,TCAGCTGAT
        """))
    expected_conseqs = [{'conseq_cutoff': 'MAX',
                         'ok_for_release': True,
                         'qcutoff': 15.0,
                         'region': 'HIV1-B-FR-K03455-seed',
                         'samplename': '2140A_HIV',
                         'seq': 'CCTCAGGTC',
                         'snum': 'S17',
                         'testcode': None},
                        {'conseq_cutoff': '0.010',
                         'ok_for_release': True,
                         'qcutoff': 15.0,
                         'region': 'HIV1-B-FR-K03455-seed',
                         'samplename': '2140A_HIV',
                         'seq': 'CCTCWGGTC',
                         'snum': 'S17',
                         'testcode': None},
                        {'conseq_cutoff': 'MAX',
                         'ok_for_release': True,
                         'qcutoff': 15.0,
                         'region': 'SARS-CoV-2-seed',
                         'samplename': '2190A_SARSCOV2',
                         'seq': 'TCAGCTGAT',
                         'snum': 'S23',
                         'testcode': None}]

    conseqs = build_conseqs(conseq_csv, '23-Jan-2014', sample_sheet)

    assert conseqs == expected_conseqs


def test_build_conseqs_multi_project():
    sample_sheet_path = Path(__file__).parent / 'microtest' / 'SampleSheet.csv'
    with sample_sheet_path.open() as sample_sheet_file:
        sample_sheet = sample_sheet_parser(sample_sheet_file)
    conseq_csv = StringIO(dedent("""\
        sample,region,q-cutoff,consensus-percent-cutoff,offset,sequence
        2100A-HCV-1337B-V3LOOP-PWND-HIV_S12,HCV-1a,15,MAX,6603,AGGAATACG
        2100A-HCV-1337B-V3LOOP-PWND-HIV_S12,HIV1-B-FR-K03455-seed,15,MAX,2776,TTTCAGAGA
        """))
    expected_conseqs = [{'conseq_cutoff': 'MAX',
                         'ok_for_release': True,
                         'qcutoff': 15.0,
                         'region': 'HCV-1a',
                         'samplename': '2100A_HCV;1337B_V3LOOP;PWND_HIV',
                         'seq': 'AGGAATACG',
                         'snum': 'S12',
                         'testcode': None},
                        {'conseq_cutoff': 'MAX',
                         'ok_for_release': True,
                         'qcutoff': 15.0,
                         'region': 'HIV1-B-FR-K03455-seed',
                         'samplename': '2100A_HCV;1337B_V3LOOP;PWND_HIV',
                         'seq': 'TTTCAGAGA',
                         'snum': 'S12',
                         'testcode': None}]

    conseqs = build_conseqs(conseq_csv, '23-Jan-2014', sample_sheet)

    assert conseqs == expected_conseqs


def test_build_conseqs_unknown_project(caplog):
    sample_sheet_file = StringIO("""\
[Header]
Assay,Amplicon
[Data]
Sample_ID,Sample_Name,Sample_Plate,Sample_Well,index,index2,Sample_Project,Description,GenomeFolder
CFE_MS1_23-Jan-2014_N512-N701_2100A_HCV_nopid;\
CFE_MS1_23-Jan-2014_N512-N701_1337B_Unknown_nopid,2100A_HCV;1337B_Unknown,\
23-Jan-2014,N/A,TTTAAAAA,AAATTTTT,23-Jan-2014,Research:2100A_HCV:FALSE;1337B_Unknown:FALSE \
Comments:2100A_HCV:;1337B_Unknown: Disablecontamcheck:2100A_HCV:FALSE;1337B_Unknown:FALSE,
""")
    sample_sheet = sample_sheet_parser(sample_sheet_file)
    conseq_csv = StringIO(dedent("""\
        sample,region,q-cutoff,consensus-percent-cutoff,offset,sequence
        2100A-HCV-1337B-Unknown_S1,HCV-1a,15,MAX,6603,AGGAATACG
        2100A-HCV-1337B-Unknown_S1,HIV1-B-FR-K03455-seed,15,MAX,2776,TTTCAGAGA
        """))
    expected_conseqs = [{'conseq_cutoff': 'MAX',
                         'ok_for_release': True,
                         'qcutoff': 15.0,
                         'region': 'HCV-1a',
                         'samplename': '2100A_HCV;1337B_Unknown',
                         'seq': 'AGGAATACG',
                         'snum': 'S1',
                         'testcode': None},
                        {'conseq_cutoff': 'MAX',
                         'ok_for_release': False,
                         'qcutoff': 15.0,
                         'region': 'HIV1-B-FR-K03455-seed',
                         'samplename': '2100A_HCV;1337B_Unknown',
                         'seq': 'TTTCAGAGA',
                         'snum': 'S1',
                         'testcode': None}]

    conseqs = build_conseqs(conseq_csv, '23-Jan-2014', sample_sheet)

    assert conseqs == expected_conseqs
    assert caplog.messages == []


def test_build_conseqs_bogus_project(caplog):
    sample_sheet_file = StringIO("""\
[Header]
Assay,Amplicon
[Data]
Sample_ID,Sample_Name,Sample_Plate,Sample_Well,index,index2,Sample_Project,Description,GenomeFolder
CFE_MS1_23-Jan-2014_N512-N701_2100A_HCV_nopid;\
CFE_MS1_23-Jan-2014_N512-N701_1337B_BOGUS_nopid,2100A_HCV;1337B_BOGUS,\
23-Jan-2014,N/A,TTTAAAAA,AAATTTTT,23-Jan-2014,\
Research:2100A_HCV:FALSE;1337B_BOGUS:FALSE Comments:2100A_HCV:;1337B_BOGUS: \
Disablecontamcheck:2100A_HCV:FALSE;1337B_BOGUS:FALSE,
""")
    sample_sheet = sample_sheet_parser(sample_sheet_file)
    conseq_csv = StringIO(dedent("""\
        sample,region,q-cutoff,consensus-percent-cutoff,offset,sequence
        2100A-HCV-1337B-BOGUS_S1,HCV-1a,15,MAX,6603,AGGAATACG
        2100A-HCV-1337B-BOGUS_S1,HIV1-B-FR-K03455-seed,15,MAX,2776,TTTCAGAGA
        """))
    expected_conseqs = [{'conseq_cutoff': 'MAX',
                         'ok_for_release': True,
                         'qcutoff': 15.0,
                         'region': 'HCV-1a',
                         'samplename': '2100A_HCV;1337B_BOGUS',
                         'seq': 'AGGAATACG',
                         'snum': 'S1',
                         'testcode': None},
                        {'conseq_cutoff': 'MAX',
                         'ok_for_release': False,
                         'qcutoff': 15.0,
                         'region': 'HIV1-B-FR-K03455-seed',
                         'samplename': '2100A_HCV;1337B_BOGUS',
                         'seq': 'TTTCAGAGA',
                         'snum': 'S1',
                         'testcode': None}]

    conseqs = build_conseqs(conseq_csv, '23-Jan-2014', sample_sheet)

    assert conseqs == expected_conseqs
    assert caplog.messages == [
        'Failed to load project seeds for 2100A-HCV-1337B-BOGUS_S1 in 23-Jan-2014.']
