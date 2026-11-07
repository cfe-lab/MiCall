"""Truth-free regression tests for the real-data report's accounting logic."""
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

import pytest

PATH = (Path(__file__).resolve().parents[2] / 'docs/validation' /
        'real_hiv_hcv_20261106_analyze.py')
SPEC = spec_from_file_location('real_hiv_hcv_20261106_analyze', PATH)
assert SPEC and SPEC.loader
ANALYSIS = module_from_spec(SPEC)
SPEC.loader.exec_module(ANALYSIS)


def test_base_calls_compare_using_reference_offsets():
    baseline = {'ref': (10, 'A-CGT'), 'second': (0, 'NN')}
    candidate = {'ref': (11, 'CGT'), 'second': (0, 'AA')}
    result = ANALYSIS.calc_consensus_delta(baseline, candidate)
    assert result['callable_upstream'] == 4
    assert result['callable_nonstealing'] == 5
    assert result['gained_callable'] == 3
    assert result['lost_callable'] == 2
    assert result['discordant_bases'] == 2


def test_derived_q15_depth_obeys_merged_cigar_and_overlap(tmp_path):
    (tmp_path / 'aligned.csv').write_text(
        'refname,qcut,rank,count,offset,seq\n'
        'A,15,0,1,1,AC-GT\n'
        'A,15,1,1,1,ACCGT\n'
        'A,10,2,1,1,ACCGT\n'
        'B,15,0,2,0,TTTT\n')
    depth, accepted, invalid = ANALYSIS.depth_arrays(
        tmp_path, {'A': 'TTTTTTT', 'B': 'AAAA'})
    assert accepted == 4
    assert invalid == 0
    assert list(depth['A']) == [0, 2, 2, 1, 2, 2, 0]
    assert list(depth['B']) == [2, 2, 2, 2]
    assert ANALYSIS.calc_coverage_delta(
        {'A': depth['A']}, {'A': depth['A']})[
            'depth_positions_decreased'] == 0


@pytest.mark.parametrize('missing', (False, True))
def test_sidecar_validates_all_reference_positions(tmp_path, missing):
    (tmp_path / 'depth.csv').write_text(
        'refname,position,max_depth,avg_depth\n'
        'A,1,1,0.5\n' + ('' if missing else
                          'A,2,0,0\n'))
    from array import array

    result = ANALYSIS.verify_sidecar(
        tmp_path / 'depth.csv', {'A': array('q', [1, 0])})
    assert result['sidecar_violations'] == int(missing)


def test_true_sequenced_pairs_separate_from_multiple_mapping_hypotheses(tmp_path):
    (tmp_path / 'remap.csv').write_text(
        'qname,flag,rname,pos\n'
        'one,99,A,2\none,147,A,9\n'
        'one,355,B,2\none,403,B,9\n'
        'two,99,A,2\ntwo,147,A,9\n')
    (tmp_path / 'remap_counts.csv').write_text(
        'type,count\nraw,4\nremap-final A,4\nremap-final B,2\n'
        'remap-unique,4\nunmapped,0\n')
    stats = ANALYSIS.remap_stats(tmp_path)
    assert stats['mapped_pairs'] == 2
    assert stats['selected_placements'] == 3
    assert stats['cross_contig_pairs'] == 1
    assert stats['original_pairs'] == 2
