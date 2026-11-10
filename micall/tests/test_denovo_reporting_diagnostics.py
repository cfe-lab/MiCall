"""Alignment-cap diagnostics remain visible without blocking de novo reports."""
import pytest

from micall.drivers.sample import Sample


@pytest.mark.parametrize('diagnostic', [
    'multimap-capped-pairs', 'multimap-oversized-pairs'])
@pytest.mark.parametrize('prefix', ['', 'unstitched_'])
def test_post_assembly_reports_continue_despite_mapping_diagnostics(
        tmp_path, monkeypatch, diagnostic, prefix):
    from micall.drivers import sample as sample_module

    (tmp_path / (prefix + 'remap_counts.csv')).write_text(
        f'type,count\nraw,400\n{diagnostic},2\n')
    (tmp_path / (prefix + 'remap.csv')).write_text(
        'qname,flag,rname,pos,mapq,cigar,rnext,pnext,tlen,seq,qual\n')
    (tmp_path / (prefix + 'remap_conseq.csv')).write_text(
        'region,sequence\n')

    sample = Sample(fastq1='example_L001_R1_001.fastq',
                    scratch_path=str(tmp_path))
    called = []

    class ReachedSam2Aln(Exception):
        pass

    def fake_sam2aln(*args, **kwargs):
        called.append('sam2aln')
        raise ReachedSam2Aln

    monkeypatch.setattr(sample_module, 'sam2aln', fake_sam2aln)
    with pytest.raises(ReachedSam2Aln):
        sample.process_post_assembly(prefix=prefix, use_denovo=True,
                                     excluded_projects=())
    assert called == ['sam2aln']


@pytest.mark.parametrize('diagnostic', [
    'multimap-capped-pairs', 'multimap-oversized-pairs'])
def test_unstitched_mapping_diagnostics_allow_referencefull_stitching(
        tmp_path, monkeypatch, diagnostic):
    from micall.drivers import sample as sample_module

    (tmp_path / 'merged_contigs.csv').write_text('contig\n')
    sample = Sample(fastq1='demo_L001_R1_001.fastq',
                    scratch_path=str(tmp_path))
    reached = []

    def fake_assembly(read1, read2, result, merged):
        result.write_text('>assembled\nACGTACGT\n')

    def fake_csv_to_fasta(source, output):
        output.touch()

    def fake_fasta_to_csv(source, output, blast_csv):
        output.write(
            'ref,match,group_ref,contig\n'
            'HCV-1a,0.5,HCV-1a,ACGTACGT\n')
        blast_csv.write('query,reference\n')

    def fake_referenceless():
        (tmp_path / 'stitched_contigs.fasta').write_text(
            '>assembled\nACGTACGT\n')

    def fake_mapper(read1, read2, contigs, remap, counts, references,
                    unmapped1, unmapped2, path, **kwargs):
        reached.append('remap')
        counts.write(f'type,count\nraw,4\n{diagnostic},2\n')
        remap.write('qname,flag,rname\n')

    def fake_referencefull(stitched, output, plot, original_counts):
        reached.append('referencefull')
        output.write(
            'ref,match,group_ref,contig\n'
            'HCV-1a,0.5,HCV-1a,ACGTACGT\n')

    monkeypatch.setattr(sample_module, 'denovo', fake_assembly)
    monkeypatch.setattr(sample_module, 'csv_to_fasta', fake_csv_to_fasta)
    monkeypatch.setattr(sample_module, 'fasta_to_csv', fake_fasta_to_csv)
    monkeypatch.setattr(sample_module, 'map_to_contigs', fake_mapper)
    monkeypatch.setattr(sample_module, 'referencefull_contig_stitcher',
                        fake_referencefull)
    monkeypatch.setattr(sample, 'run_referenceless_stitcher',
                        fake_referenceless)

    sample.run_denovo(excluded_seeds=())
    assert reached == ['remap', 'referencefull', 'remap']
    assert (tmp_path / 'contigs.csv').exists()
    assert (tmp_path / 'remap_counts.csv').read_text().count(diagnostic) == 1
