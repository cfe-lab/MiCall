from io import StringIO
from pathlib import Path
from textwrap import dedent
import pytest
from unittest.mock import patch

from micall.monitor.update_qai import build_conseqs, build_review_decisions, upload_review_to_qai
from micall.utils.sample_sheet_parser import _sample_sheet_parser
import micall.monitor.qai_helper


@pytest.fixture
def sample_sheet():
    sample_sheet_path = Path(__file__).parent / 'microtest' / 'SampleSheet.csv'
    with sample_sheet_path.open() as sample_sheet_file:
        sample_sheet = _sample_sheet_parser(sample_sheet_file)
    return sample_sheet


def test_build_conseqs(sample_sheet):
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


def test_build_conseqs_multi_project(sample_sheet):
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
    sample_sheet = _sample_sheet_parser(sample_sheet_file)
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
    sample_sheet = _sample_sheet_parser(sample_sheet_file)
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


# noinspection DuplicatedCode
def test_build_review_decisions(sample_sheet):
    coverage_file = StringIO(dedent("""\
    sample,project,region,seed,q.cut,min.coverage,which.key.pos,off.score,on.score
    2140A-HIV_S17,HIVB,HIV1B-gag,HIV1-B-FR-K03455-seed,15,150,1,-3,4
    """))
    collated_counts_file = StringIO(dedent("""\
    sample,type,count,filtered_count,seed_dist,other_dist,other_seed
    2140A-HIV_S17,raw,20000
    """))
    cascade_file = StringIO(dedent("""\
    sample,demultiplexed,v3loop,g2p,prelim_map,remap,aligned
    2140A-HIV_S17,10000,0,0,10000,10000,10000
    """))

    sequencings = [{'tag': 'N505-N702', 'target_project': 'HIVB', 'id': 'HIVid'}]

    project_regions = [{'id': 'HIV project region id',
                        'project_name': 'HIVB',
                        'seed_region_names': 'HIV1-seed',
                        'coordinate_region_name': 'HIV1B-gag'}]
    regions = [{'id': 'HIV seed region id', 'name': 'HIV1-B-FR-K03455-seed'}]

    expected_decisions = [{'sequencing_id': 'HIVid',
                           'project_region_id': 'HIV project region id',
                           'seed_region_id': 'HIV seed region id',
                           'sample_name': '2140A-HIV_S17',
                           'score': 4,
                           'min_coverage': 150,
                           'min_coverage_pos': 1,
                           'raw_reads': 20000,
                           'mapped_reads': None}]

    decisions = build_review_decisions(coverage_file,
                                       collated_counts_file,
                                       cascade_file,
                                       sample_sheet,
                                       sequencings,
                                       project_regions,
                                       regions)

    assert decisions == expected_decisions


# noinspection DuplicatedCode
def test_build_review_decisions_sequencing_different_project(sample_sheet):
    coverage_file = StringIO(dedent("""\
    sample,project,region,seed,q.cut,min.coverage,which.key.pos,off.score,on.score
    2140A-HIV_S17,HIVB,HIV1B-gag,HIV1-B-FR-K03455-seed,15,150,1,-3,4
    """))
    collated_counts_file = StringIO(dedent("""\
    sample,type,count,filtered_count,seed_dist,other_dist,other_seed
    2140A-HIV_S17,raw,20000
    """))
    cascade_file = StringIO(dedent("""\
    sample,demultiplexed,v3loop,g2p,prelim_map,remap,aligned
    2140A-HIV_S17,10000,0,0,10000,10000,10000
    """))

    sequencings = [{'tag': 'N505-N702', 'target_project': 'HIVGHA', 'id': 'HIVid'}]

    project_regions = [{'id': 'HIV project region id',
                        'project_name': 'HIVB',
                        'seed_region_names': 'HIV1-seed',
                        'coordinate_region_name': 'HIV1B-gag'}]
    regions = [{'id': 'HIV seed region id', 'name': 'HIV1-B-FR-K03455-seed'}]

    expected_decisions = [{'sequencing_id': 'HIVid',
                           'project_region_id': 'HIV project region id',
                           'seed_region_id': 'HIV seed region id',
                           'sample_name': '2140A-HIV_S17',
                           'score': -3,
                           'min_coverage': 150,
                           'min_coverage_pos': 1,
                           'raw_reads': 20000,
                           'mapped_reads': None}]

    decisions = build_review_decisions(coverage_file,
                                       collated_counts_file,
                                       cascade_file,
                                       sample_sheet,
                                       sequencings,
                                       project_regions,
                                       regions)

    assert decisions == expected_decisions


# noinspection DuplicatedCode
def test_build_review_decisions_no_reads(sample_sheet):
    coverage_file = StringIO(dedent("""\
    sample,project,region,seed,q.cut,min.coverage,which.key.pos,off.score,on.score
    """))
    collated_counts_file = StringIO(dedent("""\
    sample,type,count,filtered_count,seed_dist,other_dist,other_seed
    2140A-HIV_S17,raw,20000
    """))
    cascade_file = StringIO(dedent("""\
    sample,demultiplexed,v3loop,g2p,prelim_map,remap,aligned
    2140A-HIV_S17,10000,0,0,10000,10000,10000
    """))

    sequencings = [{'tag': 'N505-N702', 'target_project': 'HIVB', 'id': 'HIVid'}]

    project_regions = [{'id': 'HIV project region id',
                        'project_name': 'HIVB',
                        'seed_region_names': 'HIV1-seed',
                        'coordinate_region_name': 'HIV1B-gag'}]
    regions = [{'id': 'HIV seed region id', 'name': 'HIV1-B-FR-K03455-seed'}]

    expected_decisions = [{'sequencing_id': 'HIVid',
                           'sample_name': '2140A-HIV_S17',
                           'raw_reads': 20000,
                           'mapped_reads': 0}]

    decisions = build_review_decisions(coverage_file,
                                       collated_counts_file,
                                       cascade_file,
                                       sample_sheet,
                                       sequencings,
                                       project_regions,
                                       regions)

    assert decisions == expected_decisions


# noinspection DuplicatedCode
def test_build_review_decisions_more_counts(sample_sheet):
    coverage_file = StringIO(dedent("""\
    sample,project,region,seed,q.cut,min.coverage,which.key.pos,off.score,on.score
    2140A-HIV_S17,HIVB,HIV1B-gag,HIV1-B-FR-K03455-seed,15,100,10,-2,3
    """))
    collated_counts_file = StringIO(dedent("""\
    sample,type,count,filtered_count,seed_dist,other_dist,other_seed
    2140A-HIV_S17,raw,200
    2140A-HIV_S17,prelim HIV1-B-FR-K03455-seed,200,200
    2140A-HIV_S17,remap-1 HIV1-B-FR-K03455-seed,200
    2140A-HIV_S17,remap-final HIV1-B-FR-K03455-seed,200
    2140A-HIV_S17,unmapped,0
    """))
    cascade_file = StringIO(dedent("""\
    sample,demultiplexed,v3loop,g2p,prelim_map,remap,aligned
    2140A-HIV_S17,100,0,0,100,100,100
    """))

    sequencings = [{'tag': 'N505-N702', 'target_project': 'HIVB', 'id': 'HIVid'}]

    project_regions = [{'id': 'HIV project region id',
                        'project_name': 'HIVB',
                        'seed_region_names': 'HIV1-seed',
                        'coordinate_region_name': 'HIV1B-gag'}]
    regions = [{'id': 'HIV seed region id', 'name': 'HIV1-B-FR-K03455-seed'}]

    expected_decisions = [{'sequencing_id': 'HIVid',
                           'project_region_id': 'HIV project region id',
                           'seed_region_id': 'HIV seed region id',
                           'sample_name': '2140A-HIV_S17',
                           'score': 3,
                           'min_coverage': 100,
                           'min_coverage_pos': 10,
                           'raw_reads': 200,
                           'mapped_reads': 200}]

    decisions = build_review_decisions(coverage_file,
                                       collated_counts_file,
                                       cascade_file,
                                       sample_sheet,
                                       sequencings,
                                       project_regions,
                                       regions)

    assert decisions == expected_decisions


def return_json_side_effect(url):
    if url.split('=')[0] == "/lab_miseq_project_regions?pipeline":
        return [{'id': 'HIV project region id',
                        'project_name': 'HIVB',
                        'seed_region_names': 'HIV1-seed',
                        'coordinate_region_name': 'HIV1B-gag'}]
    elif url == "/lab_miseq_regions":
        return [{'id': 'HIV seed region id', 'name': 'HIV1-B-FR-K03455-seed'}]
    elif url == "/lab_miseq_pipelines?version=7.15":
        return [{'id': 'pipelineID'}]
    else:
        return None


# noinspection DuplicatedCode
def test_upload_review_to_qai(sample_sheet):
    coverage_file = StringIO(dedent("""\
    sample,project,region,seed,q.cut,min.coverage,which.key.pos,off.score,on.score
    2140A-HIV_S17,HIVB,HIV1B-gag,HIV1-B-FR-K03455-seed,15,100,10,-2,3
    """))
    collated_counts_file = StringIO(dedent("""\
    sample,type,count,filtered_count,seed_dist,other_dist,other_seed
    2140A-HIV_S17,raw,200
    2140A-HIV_S17,prelim HIV1-B-FR-K03455-seed,200,200
    2140A-HIV_S17,remap-1 HIV1-B-FR-K03455-seed,200
    2140A-HIV_S17,remap-final HIV1-B-FR-K03455-seed,200
    2140A-HIV_S17,unmapped,0
    """))
    cascade_file = StringIO(dedent("""\
    sample,demultiplexed,v3loop,g2p,prelim_map,remap,aligned
    2140A-HIV_S17,100,0,0,100,100,100
    """))
    run = {'id': 'test_run_id', 'sequencing_summary': [{'tag': 'N505-N702', 'target_project': 'HIVB', 'id': 'HIVid'}]}
    conseqs = []

    expected_decision = {'sequencing_id': 'HIVid',
                         'project_region_id': 'HIV project region id',
                         'seed_region_id': 'HIV seed region id',
                         'sample_name': '2140A-HIV_S17',
                         'score': 3,
                         'min_coverage': 100,
                         'min_coverage_pos': 10,
                         'raw_reads': 200,
                         'mapped_reads': 200}

    with patch('micall.monitor.qai_helper.Session') as mock_session:
        mock_session.get_json.side_effect = return_json_side_effect
        upload_review_to_qai(coverage_file,
                             collated_counts_file,
                             cascade_file,
                             run,
                             sample_sheet,
                             conseqs,
                             mock_session,
                             pipeline_version='7.15')

        mock_session.post_json.assert_called_once_with("/lab_miseq_reviews",
                                                       {'runid': 'test_run_id',
                                                        'pipeline_id': 'pipelineID',
                                                        'lab_miseq_review_decisions': [expected_decision],
                                                        'lab_miseq_conseqs': conseqs})


def test_upload_unknown_pipeline():
    run = {'id': 'runid', 'sequencing_summary': [{'tag': 'N505-N702', 'target_project': 'HIVB', 'id': 'HIVid'}]}
    coverage_file = StringIO()
    collated_counts_file = StringIO()
    cascade_file = StringIO()
    sample_sheet = {}
    conseqs = []
    session = micall.monitor.qai_helper.Session()
    pipeline_version = "unknown"

    with patch('micall.monitor.qai_helper.Session.get_json') as mock_get_json:
        mock_get_json.return_value = []
        with pytest.raises(RuntimeError) as exception:
            upload_review_to_qai(coverage_file,
                                 collated_counts_file,
                                 cascade_file,
                                 run,
                                 sample_sheet,
                                 conseqs,
                                 session,
                                 pipeline_version)
        assert str(exception.value) == f"Unknown pipeline: {pipeline_version}"
