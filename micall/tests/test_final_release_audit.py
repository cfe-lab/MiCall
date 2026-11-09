"""Regression tests for the independent final-release scoring and provenance gate."""
from importlib.util import module_from_spec, spec_from_file_location
import json
from pathlib import Path
import sys

import pytest

VALIDATION = (Path(__file__).resolve().parents[2] / 'docs/validation')
sys.path.insert(0, str(VALIDATION))
SOURCE = VALIDATION / 'final_release_audit_20261108.py'
SPEC = spec_from_file_location('final_release_audit_20261108', SOURCE)
assert SPEC is not None and SPEC.loader is not None
AUDIT = module_from_spec(SPEC)
SPEC.loader.exec_module(AUDIT)


def build_case(folder):
    folder.mkdir(parents=True)
    (folder / 'status.json').write_text(json.dumps({
        'ok': True, 'timings': {'total_s': 3.5, 'remap_s': 1.25},
        'provenance': {
            'sha256_r1': 'r1', 'sha256_r2': 'r2',
            'sha256_contigs': 'panel'},
    }))
    (folder / 'conseq.csv').write_text(
        'region,q-cutoff,consensus-percent-cutoff,offset,sequence\n'
        'A,15,MAX,1,CGT\n'
        'B,15,MAX,0,TGAA\n'
        'A,10,MAX,0,AAAA\n')


def test_truth_scores_correct_wrong_and_absent_bases(tmp_path, monkeypatch):
    build_case(tmp_path / 'trial')
    monkeypatch.setattr(AUDIT, 'known_truth',
                        lambda _: {'A': 'ACGT', 'B': 'TGCA'})
    result = AUDIT.score_truth(tmp_path / 'trial', 'TEST')
    assert result['totals'] == {
        'correct': 6, 'wrong': 1, 'uncallable': 1}
    assert result['per_clone']['A'] == {
        'correct': 3, 'wrong': 0, 'uncallable': 1}
    assert result['timing_seconds'] == 3.5
    assert result['sha']['sha256_contigs'] == 'panel'


def test_comparison_rejects_different_contig_or_fastq_sources():
    first = {'sha': {'sha256_r1': 'R1', 'sha256_r2': 'R2',
                     'sha256_contigs': 'REF'}}
    same = {'sha': first['sha'].copy()}
    AUDIT.check_matching_inputs(first, same, 'valid')
    same['sha']['sha256_contigs'] = 'different-ref'
    with pytest.raises(ValueError, match='different inputs'):
        AUDIT.check_matching_inputs(first, same, 'invalid')


def test_deliberately_wrong_reference_seed_is_scored_independently(
        tmp_path, monkeypatch):
    build_case(tmp_path / 'trial')
    study = tmp_path / 'study-TEST'
    study.mkdir()
    (study / 'perturbations.json').write_text(json.dumps([
        {'refname': 'A', 'position': 2, 'original': 'C', 'mutated': 'G'},
        {'refname': 'B', 'position': 3, 'original': 'C', 'mutated': 'A'},
        {'refname': 'A', 'position': 1, 'original': 'A', 'mutated': 'T'},
    ]))
    monkeypatch.setattr(AUDIT, 'PRIVATE', tmp_path)
    outcome = AUDIT.count_perturbations(tmp_path / 'trial', 'TEST')
    assert outcome == {'reverted': 1, 'wrong_seed_retained': 1,
                       'other_wrong': 0, 'uncallable': 1}
