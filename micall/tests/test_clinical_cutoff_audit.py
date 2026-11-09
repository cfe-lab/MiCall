"""Strict verification of paired MiCall clinical mixture cutoff comparisons."""
import json
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

import pytest


SOURCE = (Path(__file__).resolve().parents[2] / 'docs/validation' /
          'clinical_cutoff_audit_20261109.py')
SPEC = spec_from_file_location('clinical_cutoff_audit_20261109', SOURCE)
assert SPEC and SPEC.loader
AUDIT = module_from_spec(SPEC)
SPEC.loader.exec_module(AUDIT)


def synthetic_case(root, alias, changed=False):
    outputs = root / 'runs' / alias
    for arm in ('upstream','nonstealing'):
        folder=outputs / arm
        folder.mkdir(parents=True)
        status={
            'ok':True,
            'provenance':{'sha256_r1':'R1',
                          'sha256_r2':'R2',
                          'sha256_contigs':'C'},
        }
        (folder/'status.json').write_text(json.dumps(status))
        before='ACCTA'
        after='AYCTA' if changed and arm=='nonstealing' else before
        (folder/'conseq.csv').write_text(
            'region,q-cutoff,consensus-percent-cutoff,offset,sequence\n'
            f'1-HCV-1a,15,0.200,7,{after}\n'
            '1-HCV-1a,15,MAX,7,ACTTA\n')
    return outputs


def test_reports_all_clinical_threshold_changes_by_coordinate_and_symbol(
        tmp_path, monkeypatch):
    monkeypatch.setattr(AUDIT, 'STUDY', tmp_path)
    synthetic_case(tmp_path,'TEST',changed=True)
    d=AUDIT.compare_sample('TEST','HCV')
    assert len(d)==1
    assert d[0]['mixture_cutoff']=='0.200'
    assert d[0]['upstream_symbol']=='C'
    assert d[0]['candidate_symbol']=='Y'
    assert d[0]['one_based_contig_position']==9


def test_mismatched_source_fingerprints_fail_closed(tmp_path, monkeypatch):
    monkeypatch.setattr(AUDIT, 'STUDY', tmp_path)
    original=synthetic_case(tmp_path,'TEST',changed=True)
    p=original/'nonstealing/status.json'
    status=json.loads(p.read_text())
    status['provenance']['sha256_r1']='DIFFERENT ORIGINAL FASTQ'
    p.write_text(json.dumps(status))
    with pytest.raises(ValueError,match='source or fixed references mismatch'):
        AUDIT.compare_sample('TEST','HCV')


def test_missing_20pct_reporting_row_cannot_silently_pass(
        tmp_path, monkeypatch):
    monkeypatch.setattr(AUDIT, 'STUDY', tmp_path)
    folder=synthetic_case(tmp_path,'TEST',changed=False)
    target=folder/'nonstealing/conseq.csv'
    target.write_text(
        'region,q-cutoff,consensus-percent-cutoff,offset,sequence\n'
        '1-HCV-1a,15,MAX,7,ACTTA\n')
    with pytest.raises(ValueError,match='reporting thresholds'):
        AUDIT.compare_sample('TEST','HCV')
