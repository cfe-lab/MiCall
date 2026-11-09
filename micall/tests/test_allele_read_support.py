"""Orthogonal original-read allele evidence: errors, phase, and fragment counting."""
import gzip
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import random

import pytest


PATH = (Path(__file__).resolve().parents[2] / 'docs/validation' /
        'allele_read_support.py')
SPEC = spec_from_file_location('allele_read_support', PATH)
assert SPEC and SPEC.loader
MODULE = module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def truth():
    rng = random.Random(184)
    sequence = ''.join(rng.choice('ACGT') for _ in range(201))
    # Two linked positions six nucleotides apart.
    sequence = sequence[:94] + 'T' + sequence[95:100] + 'A' + sequence[101:]
    assert sequence[94] == 'T' and sequence[100] == 'A'
    return sequence


def fastqs(tmp_path, reads):
    first, second = tmp_path/'R1.fastq.gz', tmp_path/'R2.fastq.gz'
    with gzip.open(first, 'wt') as l, gzip.open(second, 'wt') as r:
        for index, (a, b, q1, q2) in enumerate(reads):
            name=f'@read{index}'
            l.write(f'{name}\n{a}\n+\n{q1}\n')
            r.write(f'{name}\n{b}\n+\n{q2}\n')
    return first, second


def test_original_pair_votes_are_unique_ambiguous_and_quality_filtered(tmp_path):
    sequence=truth()
    word=MODULE.allele_words(sequence,101,('A','C'),31)
    a,c=word['A'],word['C']
    reads=[
        (a,MODULE.reverse_complement(a),'I'*31,'I'*31),
        (c,MODULE.reverse_complement(c),'I'*31,'I'*31),
        (a,c,'I'*31,'I'*31),
        (a,a,'!'*31,'!'*31),
    ]
    f1,f2=fastqs(tmp_path,reads)
    scores=MODULE.original_pair_allele_support(
        sequence,101,('A','C'),f1,f2,length=31,quality=30)
    assert scores['original_fragments']==4
    assert scores['original_pair_A']==1
    assert scores['original_pair_C']==1
    assert scores['original_pairs_supporting_multiple_alleles']==1
    assert scores['original_pairs_no_exact_context']==1


def test_neighboring_snp_is_required_to_recover_phased_minor_allele(tmp_path):
    sequence=truth()
    shifted=sequence[:94]+'G'+sequence[95:]
    assert shifted[94]=='G' and shifted[100]=='A'
    minor=MODULE.allele_words(shifted,101,('C',),31)['C']
    r1,r2=fastqs(tmp_path,[(
        minor,MODULE.reverse_complement(minor),'I'*31,'I'*31)])
    unlinked=MODULE.original_pair_allele_support(
        sequence,101,('C',),r1,r2,length=31,quality=30)
    linked=MODULE.original_pair_allele_support(
        shifted,101,('C',),r1,r2,length=31,quality=30)
    assert not unlinked.get('original_pair_C',0)
    assert linked['original_pair_C']==1


def test_reference_mismatch_or_invalid_read_pair_rejected(tmp_path):
    sequence=truth()
    with pytest.raises(ValueError,match='odd length'):
        MODULE.allele_words(sequence,101,('A','C'),30)
    with pytest.raises(ValueError,match='flanking'):
        MODULE.allele_words(sequence,1,('A','C'),31)
    first,second=fastqs(tmp_path,[('ACGT','ACGT','IIII','IIII')])
    with gzip.open(second,'wt') as f:
        f.write('@out-of-sync\nACGT\n+\nIIII\n')
    with pytest.raises(ValueError,match='synchronization'):
        list(MODULE.paired_fastq(first,second))
