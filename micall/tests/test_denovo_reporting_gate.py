"""De novo clinical reporting must fail closed on incomplete hit searches."""
from io import StringIO

import pytest

from micall.core.denovo_multimap import (
    require_complete_denovo_mapping_for_clinical_report)
from micall.drivers.sample import Sample


def check(text):
    return require_complete_denovo_mapping_for_clinical_report(
        StringIO(text))


def test_complete_mapping_counts_allow_post_assembly_reporting():
    check(
        'type,count\nraw,400\n'
        'remap-final 1-HCV-1a,250\n'
        'multimap-capped-pairs,0\n'
        'multimap-oversized-pairs,0\n'
        'unmapped,150\n')


@pytest.mark.parametrize('flag', [
    'multimap-capped-pairs', 'multimap-oversized-pairs'])
def test_incomplete_mapping_refuses_clinical_consensus(flag):
    with pytest.raises(RuntimeError, match='Clinical de novo consensus withheld'):
        check(f'type,count\nraw,400\n{flag},1\n')


def test_missing_or_corrupt_mapping_qc_fails_closed():
    with pytest.raises(ValueError,match='type/count headers'):
        check('unrecognized,data\nrow,2\n')
    with pytest.raises(ValueError,match='raw input count'):
        check('type,count\nmultimap-capped-pairs,0\n')
    with pytest.raises(ValueError,match='Negative'):
        check('type,count\nraw,200\nmultimap-capped-pairs,-1\n')


def test_clinical_sample_stops_before_generating_any_consensus(tmp_path):
    (tmp_path/'remap_counts.csv').write_text(
        'type,count\nraw,400\nmultimap-capped-pairs,200\n')
    sample=Sample(
        fastq1='example_L001_R1_001.fastq',
        scratch_path=str(tmp_path))
    with pytest.raises(RuntimeError,match='capped-pairs=200'):
        sample.process_post_assembly(
            prefix='',use_denovo=True,excluded_projects=())
    assert not (tmp_path/'aligned.csv').exists()
    assert not (tmp_path/'conseq.csv').exists()
    assert not (tmp_path/'cascade.csv').exists()


def test_mapped_mode_is_not_subject_to_denovo_guard(tmp_path):
    sample=Sample(
        fastq1='example_L001_R1_001.fastq',
        scratch_path=str(tmp_path))
    # A missing mapping-count file does not itself block mapped mode:
    # downstream mapped-mode processing will fail on a different missing
    # normal input, rather than the de novo QC guard.
    with pytest.raises(FileNotFoundError) as error:
        sample.process_post_assembly(
            prefix='',use_denovo=False,excluded_projects=())
    assert 'remap.csv' in str(error.value)


def test_unstitched_mapping_guard_runs_before_referencefull_stitching(
        tmp_path, monkeypatch):
    from micall.drivers import sample as sample_module

    (tmp_path/'merged_contigs.csv').write_text('contig\n')
    sample=Sample(
        fastq1='demo_L001_R1_001.fastq',
        scratch_path=str(tmp_path))
    reached=[]

    def fake_assembly(read1, read2, result, merged):
        result.write_text('>authentic-computational-assembly\nACGTACGT\n')

    def fake_csv_to_fasta(source, output):
        output.touch()

    def fake_fasta_to_csv(source, output, blast_csv):
        output.write(
            'ref,match,group_ref,contig\n'
            'HCV-1a,0.5,HCV-1a,ACGTACGT\n')
        blast_csv.write('query,reference\n')

    def fake_referenceless():
        (tmp_path/'stitched_contigs.fasta').write_text(
            '>authentic-computational-assembly\nACGTACGT\n')

    def fake_mapper(
            read1, read2, contigs, remap, counts, references,
            unmapped1, unmapped2, path, **kwargs):
        counts.write(
            'type,count\nraw,4\nmultimap-capped-pairs,2\n')
        remap.write('qname,flag,rname\n')

    def must_not_stitch(*args, **kwargs):
        reached.append('unsafe_referencefull_stitcher_called')

    monkeypatch.setattr(sample_module, 'denovo', fake_assembly)
    monkeypatch.setattr(sample_module, 'csv_to_fasta', fake_csv_to_fasta)
    monkeypatch.setattr(sample_module, 'fasta_to_csv', fake_fasta_to_csv)
    monkeypatch.setattr(sample_module, 'map_to_contigs', fake_mapper)
    monkeypatch.setattr(
        sample_module, 'referencefull_contig_stitcher',must_not_stitch)
    monkeypatch.setattr(
        sample, 'run_referenceless_stitcher',fake_referenceless)

    with pytest.raises(RuntimeError,match='multimap-capped-pairs=2'):
        sample.run_denovo(excluded_seeds=())
    assert not reached
    assert not (tmp_path/'contigs.csv').exists()
