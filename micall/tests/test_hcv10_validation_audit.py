"""Independent source-tracked synthetic HCV10 validation math tests."""
from collections import Counter
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

import pytest


SOURCE=(Path(__file__).resolve().parents[2]/
        'docs/validation/hcv10_synthetic_truth_20261109.py')
SPEC=spec_from_file_location('hcv10_synthetic_truth_20261109',SOURCE)
assert SPEC and SPEC.loader
AUDIT=module_from_spec(SPEC)
SPEC.loader.exec_module(AUDIT)


def fresh_stats():
    return {
        'source_mapped':Counter(),
        'source_mapped_correct':Counter(),
        'source_mapped_other':Counter(),
        'alternate_contig_count':Counter(),
        'physical_pairs_seen':set(),
    }


def test_known_source_assignments_not_confused_with_competing_hypotheses():
    stats=fresh_stats()
    AUDIT.consume('c10_p00000000',{10,4},stats)
    AUDIT.consume('c2_p00000001',{2,10},stats)
    AUDIT.consume('c3_p00000002',{10},stats)
    assert stats['source_mapped_correct'][10] == 1
    assert stats['source_mapped_correct'][2] == 1
    assert stats['source_mapped_correct'][3] == 0
    assert stats['alternate_contig_count'][10] == 3
    assert stats['source_mapped_other'][10] == 1
    assert stats['source_mapped_other'][2] == 1
    assert stats['source_mapped_other'][3] == 1


def test_original_pair_alternatives_are_deduplicated_before_truth_scoring(
        tmp_path):
    folder=tmp_path/'example'
    folder.mkdir()
    (folder/'remap.csv').write_text(
        'qname,flag,rname\n'
        'c1_p00000000,99,1-HCV-1a\n'
        'c1_p00000000,147,1-HCV-1a\n'
        'c1_p00000000,355,10-HCV-1a\n'
        'c1_p00000000,403,10-HCV-1a\n'
        'c10_p00000001,99,10-HCV-1a\n'
        'c10_p00000001,147,10-HCV-1a\n')
    (folder/'remap_counts.csv').write_text(
        'type,count\nraw,4\n'
        'remap-final 1-HCV-1a,2\n'
        'remap-final 10-HCV-1a,4\n')
    truth={'total_original_pairs':2}
    stats=AUDIT.map_metrics(folder,truth)
    assert stats['mapped_original_pairs']==2
    assert stats['correct_source_presence']['10']==1
    assert stats['absent_or_minor_contig_alternate_support']==2
    assert stats['fixed_k_saturated_groups']==0


def test_corrupted_seed_bases_must_be_checked_against_original_haplotype(
        tmp_path):
    folder=tmp_path
    (folder/'conseq.csv').write_text(
        'region,q-cutoff,consensus-percent-cutoff,offset,sequence\n'
        '1-HCV-1a,15,MAX,0,ACGTA\n')
    truth={'true_original_pair_count_by_clone':{
        f'clone{i}':0 for i in range(1,11)},
        'reference_corrupted_seed_SNV_sites':{
            f'clone{i}':{} for i in range(1,11)}}
    truth['reference_corrupted_seed_SNV_sites']['clone1']={
        '2':{'original':'C','corrupted_seed':'T'},
        '4':{'original':'T','corrupted_seed':'A'}}
    source=['ACGTA']*10
    accuracy=AUDIT.consensus_metrics(folder,source,truth)
    assert accuracy['1']['Q15_MAX_correct']==5
    assert accuracy['1']['Q15_MAX_wrong']==0
    assert accuracy['1']['seed_snv_corrections']==2
    assert accuracy['2']['Q15_MAX_correct']==0


@pytest.mark.parametrize('source_id',[
    'read1','c0_p00000000','c10_pINVALID','c1_p'])
def test_nonconforming_truth_identifier_refused(source_id):
    with pytest.raises(ValueError,match='source label missing'):
        AUDIT.consume(source_id,{10},fresh_stats())
