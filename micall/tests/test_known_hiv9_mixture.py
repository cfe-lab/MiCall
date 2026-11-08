"""Known-proportion clone mixture: independent source truth and scoring tests."""
from collections import Counter
import csv
import gzip
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import pytest


VALIDATION = Path(__file__).resolve().parents[2] / 'docs/validation'


def import_script(name):
    filename = VALIDATION / (name + '.py')
    spec = spec_from_file_location(name, filename)
    assert spec and spec.loader
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


GENERATOR = import_script('known_hiv9_mixture_20261108')
AUDIT = import_script('known_hiv9_mixture_eval_20261108')


def make_truth(tmp_path):
    truth = tmp_path / 'truth.fasta'
    records = []
    for idx in range(9):
        seq = list('ACGT' * 125)
        seq[30 + idx] = next(b for b in 'ACGT' if b != seq[30 + idx])
        records.append(SeqRecord(Seq(''.join(seq)),
                                 id=f'{idx+1}-HIV1-B-FR-K03455-seed'))
    SeqIO.write(records, truth, 'fasta')
    return truth


@pytest.mark.parametrize('fraction,expected', [
    ('0', 0), ('0p2', 2), ('1', 10)])
def test_generator_has_exact_known_truth_pair_count_and_synchronized_mates(
        tmp_path, fraction, expected):
    truth = make_truth(tmp_path)
    dest = tmp_path / f'mixture-{fraction}'
    details = GENERATOR.generate(truth, dest, fraction, 1000, 400)
    assert details['input_pairs'] == 1000
    assert details['clone_pairs']['clone9'] == expected
    first = dest / (dest.name+'_L001_R1_001.fastq.gz')
    second = dest / (dest.name+'_L001_R2_001.fastq.gz')
    qnames1 = []
    qnames2 = []
    with gzip.open(first, 'rt') as a, gzip.open(second, 'rt') as b:
        for _ in range(1000):
            rows1 = [a.readline() for _ in range(4)]
            rows2 = [b.readline() for _ in range(4)]
            qnames1.append(rows1[0].strip())
            qnames2.append(rows2[0].strip())
            assert len(rows1[1].strip()) == 250
            assert len(rows2[1].strip()) == 250
            assert rows1[3] == rows2[3] == 'I' * 250 + '\n'
        assert not a.readline()
        assert not b.readline()
    assert qnames1 == qnames2
    assert len(set(qnames1)) == 1000
    assert Counter(name.split('_')[0] for name in qnames1)[
        '@clone9'] == expected


def test_pair_support_distinguishes_true_original_pairs_from_alternatives(
        tmp_path):
    folder = tmp_path / 'out'
    folder.mkdir()
    lines = [
        ['clone9_pair00001', '99', '9-HIV1-B-FR-K03455-seed'],
        ['clone9_pair00001', '147', '9-HIV1-B-FR-K03455-seed'],
        ['clone9_pair00001', '355', '8-HIV1-B-FR-K03455-seed'],
        ['clone9_pair00001', '403', '8-HIV1-B-FR-K03455-seed'],
        ['clone9_pair00002', '73', '9-HIV1-B-FR-K03455-seed'],
        ['clone8_pair00003', '99', '9-HIV1-B-FR-K03455-seed'],
        ['clone8_pair00003', '147', '9-HIV1-B-FR-K03455-seed'],
    ]
    with (folder / 'remap.csv').open('w') as f:
        w = csv.writer(f)
        w.writerow(['qname', 'flag', 'rname'])
        w.writerows(lines)
    result = AUDIT.pair_support(folder, {'clone9': 2, 'clone8': 1})
    assert result['clone9']['full_paired'] == 1
    assert result['clone9']['paired_to_correct'] == 1
    assert result['clone9']['multi_reference'] == 1
    assert result['clone9']['unpaired_only'] == 1
    assert result['clone9']['not_fully_mapped'] == 1
    assert result['clone8']['paired_exclusively_wrong'] == 1


def test_consensus_scoring_counts_false_callable_on_absent_reference(
        tmp_path):
    folder = tmp_path / 'out'
    folder.mkdir()
    (folder / 'conseq.csv').write_text(
        'region,q-cutoff,consensus-percent-cutoff,offset,sequence\n'
        '1-HIV1-B-FR-K03455-seed,15,MAX,0,ACGT\n'
        '9-HIV1-B-FR-K03455-seed,15,MAX,0,ACAA\n')
    ref = {'1-HIV1-B-FR-K03455-seed': 'ACGT',
           '9-HIV1-B-FR-K03455-seed': 'ACGT'}
    result = AUDIT.score_calls(ref, folder)
    assert result['clone1']['correct'] == 4
    assert result['clone1']['wrong'] == 0
    assert result['clone9']['correct'] == 2
    assert result['clone9']['wrong'] == 2


def test_invalid_source_qname_is_rejected():
    with pytest.raises(ValueError):
        AUDIT.pair_origin('some_unverified_read')
