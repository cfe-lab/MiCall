"""Real Bowtie2 de novo multi-placement integration smoke test."""
import csv
from io import StringIO
import random
import shutil

import pytest

from micall.core.remap import map_to_contigs
from micall.core.sam2aln import sam2aln
from micall.utils.translation import reverse_and_complement


@pytest.mark.skipif(not (shutil.which('bowtie2-align-s') and
                         shutil.which('bowtie2-build-s')),
                    reason='Requires Bowtie2 executables')
def test_denovo_retains_equal_bowtie2_pair_placements(monkeypatch, tmp_path):
    """End-to-end FASTQ -> Bowtie2 -> MiCall remap -> sam2aln depths."""
    from micall.core import remap

    rng = random.Random(12345)
    dna = ''.join(rng.choices('ACGT', k=330))
    # The contigs differ outside the covered molecule: both placements
    # genuinely have the same best Bowtie2 alignment score.
    references = {'A': dna, 'B': dna[:320] + 'T' * 10}
    forward = dna[25:125]
    reverse = reverse_and_complement(dna[170:270])
    r1, r2 = tmp_path / 'R1.fastq', tmp_path / 'R2.fastq'
    for path, seq in ((r1, forward), (r2, reverse)):
        path.write_text('@pair1\n{}\n+\n{}\n'.format(seq, 'I' * len(seq)))

    monkeypatch.setattr(remap, 'read_contigs', lambda *_: references)

    def perform(contigs):
        monkeypatch.setattr(remap, 'read_contigs', lambda *_: contigs)
        mapping, counts = StringIO(), StringIO()
        map_to_contigs(str(r1), str(r2), StringIO(''), mapping, counts,
                       StringIO(), StringIO(), StringIO(), str(tmp_path))
        aligned, depth = StringIO(), StringIO()
        remap_csv = mapping.getvalue()
        sam2aln(StringIO(remap_csv), aligned, depth_csv=depth,
                depth_reference_csv=StringIO(
                    'region,sequence\n' + ''.join(
                        f'{ref},{seq}\n' for ref, seq in contigs.items())))
        return (list(csv.DictReader(StringIO(remap_csv))),
                list(csv.DictReader(StringIO(depth.getvalue()))),
                aligned.getvalue(),
                list(csv.DictReader(StringIO(counts.getvalue()))))

    old_rows, old_depth, _, old_counts = perform({'A': dna})
    new_rows, new_depth, aligned, new_counts = perform(references)
    assert not any(row['type'] == 'remap-unique' for row in old_counts)
    assert next(row['count'] for row in new_counts
                if row['type'] == 'remap-unique') == '2'
    assert len(old_rows) == 2
    assert len(new_rows) == 4
    assert [row['rname'] for row in new_rows] == ['A', 'A', 'B', 'B']
    assert {row['qname'] for row in new_rows} == {'pair1'}
    assert all(row['max_depth'] == '1' for row in new_depth
               if row['refname'] == 'A' and int(row['position']) == 26)
    assert {row['avg_depth'] for row in new_depth
            if row['position'] == '26'} == {'0.5'}
    assert next(row for row in old_depth
                if row['position'] == '26')['avg_depth'] == '1'
    assert sum(1 for row in new_depth if row['max_depth'] == '0') > 0
    assert {row['refname'] for row in csv.DictReader(
        StringIO(aligned))} == {'A', 'B'}


@pytest.mark.skipif(not (shutil.which('bowtie2-align-s') and
                         shutil.which('bowtie2-build-s')),
                    reason='Requires Bowtie2 executables')
def test_denovo_does_not_admit_a_strictly_worse_snv_alignment(
        monkeypatch, tmp_path):
    from micall.core import remap

    rng = random.Random(34455)
    dna = ''.join(rng.choices('ACGT', k=330))
    alternative = list(dna)
    alternative[55] = next(base for base in 'ACGT' if base != dna[55])
    references = {'A': dna, 'B': ''.join(alternative)}
    r1seq = dna[25:125]
    r2seq = reverse_and_complement(dna[170:270])
    for name, seq in (('R1', r1seq), ('R2', r2seq)):
        (tmp_path / f'{name}.fastq').write_text(
            '@pair1\n{}\n+\n{}\n'.format(seq, 'I' * len(seq)))
    monkeypatch.setattr(remap, 'read_contigs', lambda *_: references)

    remapped, depths = StringIO(), StringIO()
    map_to_contigs(str(tmp_path / 'R1.fastq'),
                   str(tmp_path / 'R2.fastq'),
                   StringIO(''), remapped, StringIO(), StringIO(),
                   StringIO(), StringIO(), str(tmp_path))
    rows = list(csv.DictReader(StringIO(remapped.getvalue())))
    # An SNV in the first mate strictly distinguishes A from B.
    assert len(rows) == 2
    assert {row['rname'] for row in rows} == {'A'}
    sam2aln(StringIO(remapped.getvalue()), StringIO(), depth_csv=depths)
    assert {row['refname'] for row in
            csv.DictReader(StringIO(depths.getvalue()))} == {'A'}


@pytest.mark.skipif(not (shutil.which('bowtie2-align-s') and
                         shutil.which('bowtie2-build-s')),
                    reason='Requires Bowtie2 executables')
def test_denovo_warns_on_more_equivalent_contigs_than_the_reporting_cap(
        monkeypatch, tmp_path, caplog):
    from micall.core import remap
    from micall.core.denovo_multimap import MAX_DENOVO_ALIGNMENTS

    rng = random.Random(1456)
    dna = ''.join(rng.choices('ACGT', k=330))
    references = {f'contig{i}': dna
                  for i in range(MAX_DENOVO_ALIGNMENTS + 1)}
    first = dna[25:125]
    second = reverse_and_complement(dna[170:270])
    for name, read in (('R1', first), ('R2', second)):
        (tmp_path / f'{name}.fastq').write_text(
            '@one\n{}\n+\n{}\n'.format(read, 'I' * len(read)))
    monkeypatch.setattr(remap, 'read_contigs', lambda *_: references)
    output = StringIO()
    map_to_contigs(str(tmp_path / 'R1.fastq'),
                   str(tmp_path / 'R2.fastq'),
                   StringIO(''), output, StringIO(), StringIO(),
                   StringIO(), StringIO(), str(tmp_path))

    rows = list(csv.DictReader(StringIO(output.getvalue())))
    assert len(rows) == MAX_DENOVO_ALIGNMENTS * 2
    assert len({row['rname'] for row in rows}) == MAX_DENOVO_ALIGNMENTS
    assert 'reporting limit' in caplog.text


@pytest.mark.parametrize('use_denovo,prefix', [
    (True, ''), (True, 'unstitched_'), (False, ''),
])
def test_sample_process_post_assembly_writes_only_denovo_depth_sidecars(
        monkeypatch, tmp_path, use_denovo, prefix):
    from micall.drivers import sample as sample_module
    from micall.drivers.sample import Sample
    from micall.tests.test_denovo_multimap import (
        placement, normalize, remap_csv)

    # Run MiCall's actual process_post_assembly entry point and real
    # sam2aln, while stubbing unrelated graphical and scoring pipelines.
    for name in ('aln2counts', 'coverage_plot', 'plot_genome_coverage',
                 'concordance_plot'):
        monkeypatch.setattr(sample_module, name, lambda *args, **kwargs: None)
    monkeypatch.setattr(sample_module.CascadeReport, 'generate',
                        lambda self: None)

    sample = Sample(scratch_path=str(tmp_path))
    text, _ = normalize(
        placement('pair1', 'A'),
        placement('pair1', 'B', secondary=True))
    (tmp_path / f'{prefix}remap.csv').write_text(remap_csv(text).getvalue())
    (tmp_path / f'{prefix}remap_conseq.csv').write_text(
        'region,sequence\nA,ACGTA\nB,ACGTA\n')
    for filename in ['g2p_aligned.csv', 'g2p_summary.csv',
                     f'{prefix}remap_counts.csv']:
        (tmp_path / filename).write_text('')
    if use_denovo:
        (tmp_path / f'{prefix}contigs.csv').write_text('')

    sample.process_post_assembly(use_denovo, (), prefix)
    aligned = (tmp_path / f'{prefix}aligned.csv').read_text()
    assert len(aligned.splitlines()) == 3
    expected_depth = tmp_path / f'{prefix}multimap_depth.csv'
    assert expected_depth.exists() is use_denovo
    if use_denovo:
        rows = list(csv.DictReader(StringIO(expected_depth.read_text())))
        assert len(rows) == 10
        assert {r['avg_depth'] for r in rows} == {'0.5'}
