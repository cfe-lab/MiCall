"""Check real-read source dilution provenance and paired-end analysis."""
from collections import Counter
import gzip
from importlib.util import module_from_spec, spec_from_file_location
import json
from pathlib import Path
import sys

import pytest

VALIDATION = (Path(__file__).resolve().parents[2] / 'docs/validation')
sys.path.insert(0, str(VALIDATION))


def load(name):
    path = VALIDATION / f'{name}.py'
    spec = spec_from_file_location(name, path)
    assert spec and spec.loader
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


MIX = load('real_hcv_dilution_20261108')
ORIGIN = load('hcv_origin_eval')
AUDIT = load('real_hcv_dilution_eval_20261108')


def fake_source(parent):
    parent.mkdir()
    label = parent.name
    names = ('L001_R1_001', 'L001_R2_001')
    for side, token in enumerate(names):
        path = parent / f'{label}_{token}.fastq.gz'
        with gzip.open(path, 'wt') as stream:
            for origin, count in [('HCV1a', 5), ('HCV2b', 3)]:
                for index in range(count):
                    seq = 'ACGT' if side == 0 else 'TGCA'
                    stream.write(
                        f'@{origin}|fragment{index}\n'
                        f'{seq}\n+\n{"I" * len(seq)}\n')


def test_dilution_preserves_real_paired_sequences_quality_and_exact_source_fraction(
        tmp_path):
    sample = tmp_path / 'source'
    fake_source(sample)
    output = tmp_path / 'experiment'
    summary = MIX.create(sample, output, 2, 41,
                         expected_major=5, expected_minor=3)
    assert summary['source_minor_pairs'] == 2
    assert summary['total_pairs'] == 7
    assert summary['actual_minor_source_fraction'] == pytest.approx(2/7)
    read_pairs = list(MIX.pairs(output))
    assert len(read_pairs) == 7
    assert Counter(row[0][0].split('|')[0] for row in read_pairs) == {
        '@HCV1a': 5, '@HCV2b': 2}
    assert all(left[0] == right[0] for left, right in read_pairs)
    assert all(left[1:] == ('ACGT\n', '+\n', 'IIII\n')
               for left, right in read_pairs)
    assert all(right[1:] == ('TGCA\n', '+\n', 'IIII\n')
               for left, right in read_pairs)


def test_dilution_raises_for_wrong_reference_pool(tmp_path):
    sample = tmp_path / 'source'
    fake_source(sample)
    with pytest.raises(ValueError, match='Source pair counts'):
        MIX.create(sample, tmp_path / 'new', 2, 41)


def test_origin_evaluator_does_not_call_primary_singleton_a_paired_fragment(
        monkeypatch, tmp_path):
    monkeypatch.setattr(ORIGIN, 'OUT', tmp_path)
    folder = tmp_path / 'case' / 'nonstealing'
    folder.mkdir(parents=True)
    (folder / 'status.json').write_text(
        json.dumps({'ok': True, 'timings': {'total_s': 3.0}}))
    (folder / 'remap.csv').write_text(
        'qname,flag,rname\n'
        'HCV1a|singleton,73,1-HCV-1a\n'
        'HCV1a|proper,99,1-HCV-1a\n'
        'HCV1a|proper,147,1-HCV-1a\n'
        'HCV2b|minor,99,2-HCV-2b\n'
        'HCV2b|minor,147,2-HCV-2b\n'
        'HCV1a|multimap,99,1-HCV-1a\n'
        'HCV1a|multimap,147,1-HCV-1a\n'
        'HCV1a|multimap,355,3-HCV-1b\n'
        'HCV1a|multimap,403,3-HCV-1b\n')
    (folder / 'conseq.csv').write_text(
        'region,q-cutoff,consensus-percent-cutoff,offset,sequence\n'
        '1-HCV-1a,15,MAX,0,ACGT\n'
        '2-HCV-2b,15,MAX,0,ACGT\n')
    results = ORIGIN.measure('case', 'nonstealing')
    majority = results['counts']['HCV1a']
    minority = results['counts']['HCV2b']
    assert majority['mapped'] == 3
    assert majority['paired_mapped'] == 2
    assert majority['single_end_only'] == 1
    assert majority['paired_correct'] == 2
    assert majority['paired_decoy'] == 1
    assert minority['paired_mapped'] == 1
    assert minority['paired_correct'] == 1


def test_fractional_depth_report_checks_decoy_consensus_invariants(tmp_path):
    (tmp_path / 'conseq.csv').write_text(
        'region,q-cutoff,consensus-percent-cutoff,offset,sequence\n'
        '3-HCV-1b,15,MAX,9,ACGT\n')
    depth = tmp_path / 'multimap_depth.csv'
    depth.write_text(
        'refname,position,max_depth,avg_depth\n'
        '3-HCV-1b,10,110,110\n3-HCV-1b,11,105,99\n'
        '3-HCV-1b,12,106,100\n3-HCV-1b,13,110,104\n')
    result = AUDIT.fractional_decoy_evidence(tmp_path)
    assert result == {'fractional_depth_at_least_100': 3,
                      'fractional_depth_below_100': 1}
    depth.write_text('refname,position,max_depth,avg_depth\n3-HCV-1b,10,110,110\n')
    with pytest.raises(ValueError, match='Missing decoy depth'):
        AUDIT.fractional_decoy_evidence(tmp_path)


def test_real_17_reference_audit_counts_original_fastq_pairs(tmp_path):
    audit = load('real17_panel_audit_20261108')
    fastq = tmp_path / 'test.fastq'
    fastq.write_text('@read1\nACGT\n+\nIIII\n'
                     '@read2\nTTTT\n+\nIIII\n')
    assert audit.count_fastq_pairs(fastq) == 2
    fastq.write_text(fastq.read_text() + '@truncated\n')
    with pytest.raises(ValueError, match='multiple of four'):
        audit.count_fastq_pairs(fastq)
