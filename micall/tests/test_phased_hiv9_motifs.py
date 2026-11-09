"""Tests for alignment-independent linked private-SNP evidence."""
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

import pytest

PATH = (Path(__file__).resolve().parents[2] / 'docs/validation' /
        'phased_hiv9_motifs.py')
SPEC = spec_from_file_location('phased_hiv9_motifs', PATH)
assert SPEC and SPEC.loader
MOD = module_from_spec(SPEC)
SPEC.loader.exec_module(MOD)


def write_fastq(path, reads, mate):
    path.write_text(''.join(
        f'@{name}/{mate}\n{bases}\n+\n{"I" * len(bases)}\n'
        for name, bases in reads))


def test_linked_snp_support_is_counted_once_per_original_pair(tmp_path):
    majority = 'A' * 50
    minor = majority[:14] + 'C' + majority[15:25] + 'G' + majority[26:]
    truth = {'major': majority, 'minor': minor}
    assert [pos for pos, _, _ in MOD.find_private_motifs(truth)['minor']] == [
        15, 26]
    first = [('majority', majority), ('minority', minor),
             ('single-marker', minor[:25] + 'A' * 25)]
    second = [(name, 'T' * len(bases)) for name, bases in first]
    r1, r2 = tmp_path / 'r1.fastq', tmp_path / 'r2.fastq'
    write_fastq(r1, first, '1')
    write_fastq(r2, second, '2')
    actual = MOD.evaluate(truth, r1, r2)
    assert actual['original_pairs'] == 3
    minor_counts = actual['per_clone']['minor']
    assert minor_counts['at_least_one'] == 2
    assert minor_counts['at_least_two'] == 1
    assert minor_counts['all_marker_motifs'] == 1


def test_pair_streams_must_be_synchronized(tmp_path):
    r1, r2 = tmp_path / 'r1.fastq', tmp_path / 'r2.fastq'
    write_fastq(r1, [('pair1', 'A' * 50)], '1')
    write_fastq(r2, [('pair2', 'T' * 50)], '2')
    with pytest.raises(ValueError, match='not synchronized'):
        MOD.evaluate({'clone': 'A' * 50}, r1, r2)


def test_reverse_complement_lookup_and_invalid_motif_width():
    assert MOD.reverse_complement('ACGTT') == 'AACGT'
    with pytest.raises(ValueError):
        MOD.find_private_motifs({'A': 'ACGT'}, width=20)
