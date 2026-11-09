"""Independent original-read genotype-discriminating sequence evidence."""
import gzip
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

import pytest


SOURCE = (Path(__file__).resolve().parents[2] / 'docs/validation' /
          'original_read_unique_motifs.py')
SPEC = spec_from_file_location('original_read_unique_motifs', SOURCE)
assert SPEC and SPEC.loader
MOTIF = module_from_spec(SPEC)
SPEC.loader.exec_module(MOTIF)


def write_fastq(path, records, mate):
    with gzip.open(path, 'wt') as f:
        for name, sequence, quality in records:
            f.write(f'@{name}/{mate}\n{sequence}\n+\n{quality}\n')


def test_distinct_reference_motifs_are_not_just_repeated_identical_kmers():
    word = 'ACGT' * 7 + 'ACG'
    other = 'TGCA' * 7 + 'TGC'
    assert MOTIF.linked_distinct_motifs({(word, 1), (word, 150)}) is False
    assert MOTIF.linked_distinct_motifs({(word, 1), (other, 20)}) is False
    assert MOTIF.linked_distinct_motifs({(word, 1), (other, 150)}) is True


def test_reference_exclusive_motifs_require_both_sequence_and_quality(
        tmp_path):
    refs = {'major': 'A' * 60 + 'C' * 90 + 'G' * 60,
            'decoy': 'A' * 60 + 'T' * 90 + 'G' * 60}
    _, counts = MOTIF.unique_markers(refs)
    assert counts['major'] > 0 and counts['decoy'] > 0
    x = refs['major'][70:155]
    y = MOTIF.reverse_complement(x)
    fq1, fq2 = tmp_path/'R1.fq.gz',tmp_path/'R2.fq.gz'
    write_fastq(fq1, [('good', x, 'I' * len(x)),
                      ('bad_quality', x, '!' * len(x)),
                      ('reverse', y, 'I' * len(y))], 1)
    write_fastq(fq2, [('good', x, 'I' * len(x)),
                      ('bad_quality', x, '!' * len(x)),
                      ('reverse', y, 'I' * len(y))], 2)
    output=MOTIF.evaluate(refs, fq1, fq2)
    assert output['input_pairs'] == 3
    assert output['per_reference']['major']['at_least_one'] == 2
    assert not output['per_reference']['decoy'].get('at_least_one', 0)
    assert output['minimum_per_base_phred'] == 30


def test_missing_or_out_of_sync_fastq_mates_fail_closed(tmp_path):
    first,second=tmp_path/'R1.gz',tmp_path/'R2.gz'
    write_fastq(first, [('first','ACGT','IIII')],1)
    write_fastq(second,[('second','ACGT','IIII')],2)
    with pytest.raises(ValueError,match='mates do not match'):
        list(MOTIF.original_pairs(first,second))


def test_repeated_candidate_word_is_panel_ambiguous():
    shared = {'A':'A' * 100+'C'*70,'B':'A'*100+'T'*70}
    evidence, counts = MOTIF.unique_markers(shared)
    assert MOTIF.canonical('A'*31) not in evidence
    assert counts['A']>0 and counts['B']>0



def test_broad_seed_specificity_cannot_infer_genotype_from_three_contigs():
    from importlib.util import spec_from_file_location, module_from_spec

    source = (Path(__file__).resolve().parents[2] / 'docs/validation' /
              'original_read_unique_motifs_report.py')
    import sys
    sys.path.insert(0, str(source.parent))
    spec = spec_from_file_location('original_read_unique_motifs_report',
                                   source)
    assert spec and spec.loader
    report = module_from_spec(spec)
    spec.loader.exec_module(report)
    contigs = {'HCV-1b': 'C' * 70, 'HCV-2b': 'A' * 70}
    other_subtypes = {
        'HCV-1b': 'C' * 90,
        'HCV-2b': 'A' * 90,
        'HCV-3a': 'C' * 90,
    }
    check = report.validate_broader_hcv_subtype_panel(
        contigs, other_subtypes)
    assert check['decoy_panel_exclusive_motifs'] == 1
    assert check['other_full_reference_subtypes_checked'] == 2
    assert check['motifs_shared_with_another_full_reference_subtype'] == 1
    assert check['motifs_not_seen_in_other_subtype_seeds'] == 0
