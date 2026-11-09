"""Real Bowtie2 de novo multi-placement integration smoke test."""
import csv
from io import StringIO
import random
import shutil

import pytest

from micall.core.remap import map_to_contigs
from micall.core.sam2aln import sam2aln
from micall.utils.translation import reverse_and_complement


def test_final_unmapped_fastqs_include_rejected_and_singleton_mates(
        monkeypatch, tmp_path):
    """Unmapped FASTQs represent the *final* mapping, not initial Bowtie2."""
    from micall.core import remap
    from micall.tests.test_denovo_multimap import placement

    input_names = ('accepted', 'rejected', 'singleton')
    originals = {}
    for mate in (1, 2):
        lines = []
        for name in input_names:
            record = (f'@{name}/{mate} description\n'
                      f'ACGTA\n+\n{mate * 5}\n')
            lines.append(record)
            originals[name, mate] = record
        (tmp_path / f'R{mate}.fastq').write_text(''.join(lines))

    accepted = placement('accepted', 'A')
    rejected = placement('rejected', 'A')
    rejected[1][7] = '999'  # Inconsistent mate coordinates: reject pair.
    singleton = placement('singleton', 'A', flag1=73, flag2=133)
    singleton[1][2] = '*'
    singleton[1][5] = '*'
    pairs = (accepted, rejected, singleton)

    def fake_mapping(*args, **kwargs):
        samfile = args[4]
        new_counts = args[12]
        new_counts.update({'A': 5})
        with open(samfile, 'w') as output:
            output.write('@SQ\tSN:A\tLN:50\n')
            for pair in pairs:
                output.writelines('\t'.join(row) + '\n' for row in pair)
        return 0

    monkeypatch.setattr(remap, 'read_contigs', lambda *_: {'A': 'ACGTA' * 10})
    monkeypatch.setattr(remap, 'find_bowtie2',
                        lambda: (object(), object()))
    monkeypatch.setattr(remap, 'map_to_reference', fake_mapping)
    remapped, counts, unmapped1, unmapped2 = (
        StringIO(), StringIO(), StringIO(), StringIO())
    remap.map_to_contigs(
        str(tmp_path / 'R1.fastq'), str(tmp_path / 'R2.fastq'),
        StringIO(''), remapped, counts, StringIO(),
        unmapped1, unmapped2, str(tmp_path))

    assert unmapped1.getvalue() == originals['rejected', 1]
    assert unmapped2.getvalue() == (
        originals['rejected', 2] + originals['singleton', 2])
    summary = {row['type']: int(row['count'])
               for row in csv.DictReader(StringIO(counts.getvalue()))}
    assert summary['unmapped-after-selection'] == 3
    assert summary['remap-unique-fragments'] == 2
    assert {r['qname'] for r in
            csv.DictReader(StringIO(remapped.getvalue()))} == {
                'accepted', 'singleton'}


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
@pytest.mark.parametrize('extra_contigs,expect_cap', [
    (1, False), (9, False), (25, True),
])
def test_denovo_reports_all_small_panel_contigs_and_warns_on_large_cap(
        monkeypatch, tmp_path, caplog, extra_contigs, expect_cap):
    from micall.core import remap
    from micall.core.denovo_multimap import MAX_DENOVO_ALIGNMENTS

    rng = random.Random(1456)
    dna = ''.join(rng.choices('ACGT', k=330))
    references = {f'contig{i}': dna
                  for i in range(MAX_DENOVO_ALIGNMENTS + extra_contigs)}
    first = dna[25:125]
    second = reverse_and_complement(dna[170:270])
    for name, read in (('R1', first), ('R2', second)):
        (tmp_path / f'{name}.fastq').write_text(
            '@one\n{}\n+\n{}\n'.format(read, 'I' * len(read)))
    monkeypatch.setattr(remap, 'read_contigs', lambda *_: references)
    output, count_csv = StringIO(), StringIO()
    map_to_contigs(str(tmp_path / 'R1.fastq'),
                   str(tmp_path / 'R2.fastq'),
                   StringIO(''), output, count_csv, StringIO(),
                   StringIO(), StringIO(), str(tmp_path))

    count_rows = {r['type']: int(r['count']) for r in
                  csv.DictReader(StringIO(count_csv.getvalue()))}
    assert count_rows.get('multimap-capped-pairs', 0) == int(expect_cap)
    rows = list(csv.DictReader(StringIO(output.getvalue())))
    expected = MAX_DENOVO_ALIGNMENTS if expect_cap else len(references)
    assert len(rows) == expected * 2
    assert len({row['rname'] for row in rows}) == expected
    assert ('reporting limit' in caplog.text) is expect_cap
    # The limitation must be machine-readable as well as logged.
    # For a small all-placement panel no artificial saturation is recorded.



@pytest.mark.skipif(not (shutil.which('bowtie2-align-s') and
                         shutil.which('bowtie2-build-s')),
                    reason='Requires Bowtie2 executables')
def test_multithreaded_all_alignments_keep_qname_groups_together(
        monkeypatch, tmp_path):
    from collections import defaultdict
    from micall.core import remap

    rng = random.Random(867)
    dna = ''.join(rng.choices('ACGT', k=350))
    monkeypatch.setattr(remap, 'read_contigs',
                        lambda *_: {f'R{i}': dna for i in range(9)})
    first = dna[25:145]
    second = reverse_and_complement(dna[175:295])
    for name, read in (('R1', first), ('R2', second)):
        fastq = tmp_path / f'{name}.fastq'
        fastq.write_text(''.join(
            f'@pair{i}\n{read}\n+\n{"I" * len(read)}\n'
            for i in range(50)))
    output = StringIO()
    map_to_contigs(str(tmp_path / 'R1.fastq'),
                   str(tmp_path / 'R2.fastq'),
                   StringIO(''), output, StringIO(), StringIO(),
                   StringIO(), StringIO(), str(tmp_path))
    counts = defaultdict(lambda: defaultdict(int))
    for row in csv.DictReader(StringIO(output.getvalue())):
        counts[row['qname']][row['rname']] += 1
    assert len(counts) == 50
    assert all(len(regions) == 9 for regions in counts.values())
    assert all(all(mates == 2 for mates in regions.values())
               for regions in counts.values())


@pytest.mark.skipif(not (shutil.which('bowtie2-align-s') and
                         shutil.which('bowtie2-build-s')),
                    reason='Requires Bowtie2 executables')
def test_rejected_mapped_pair_counts_as_unmapped_original_ends(
        monkeypatch, tmp_path):
    from collections import Counter
    from micall.core import remap

    rng = random.Random(901)
    dna = ''.join(rng.choices('ACGT', k=350))
    monkeypatch.setattr(remap, 'read_contigs',
                        lambda *_: {'A': dna})
    for name, read in (('R1', dna[25:125]),
                       ('R2', reverse_and_complement(dna[170:270]))):
        (tmp_path / f'{name}.fastq').write_text(
            f'@one\n{read}\n+\n{"I" * len(read)}\n')

    def reject_all(source, destination, limit=None, diagnostics=None):
        for line in source:
            if line.startswith('@'):
                destination.write(line)
        return Counter()

    monkeypatch.setattr(remap, 'select_best_pairs', reject_all)
    counts = StringIO()
    remapped = StringIO()
    map_to_contigs(str(tmp_path / 'R1.fastq'),
                   str(tmp_path / 'R2.fastq'), StringIO(''),
                   remapped, counts, StringIO(), StringIO(),
                   StringIO(), str(tmp_path))
    summary = {row['type']: int(row['count'])
               for row in csv.DictReader(StringIO(counts.getvalue()))}
    assert summary['raw'] == 2
    assert summary['unmapped'] == 0  # legacy Bowtie2 classification
    assert summary['unmapped-after-selection'] == 2
    assert len(list(csv.DictReader(StringIO(remapped.getvalue())))) == 0


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
    coverage_destinations = []
    monkeypatch.setattr(sample_module, 'coverage_plot',
                        lambda *args, **kwargs: coverage_destinations.append(
                            kwargs['coverage_maps_path']))
    for name in ('aln2counts', 'plot_genome_coverage', 'concordance_plot'):
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
    for filename in ['g2p_aligned.csv', 'g2p_summary.csv']:
        (tmp_path / filename).write_text('')
    # Clinical de novo reporting must have valid remap QC metadata even
    # when other reporting integrations are stubbed for this focused test.
    (tmp_path / f'{prefix}remap_counts.csv').write_text(
        'type,count\nraw,2\nremap-final A,2\n')
    if use_denovo:
        (tmp_path / f'{prefix}contigs.csv').write_text('')

    sample.process_post_assembly(use_denovo, (), prefix)
    aligned = (tmp_path / f'{prefix}aligned.csv').read_text()
    assert len(aligned.splitlines()) == 3
    expected_depth = tmp_path / f'{prefix}multimap_depth.csv'
    assert expected_depth.exists() is use_denovo
    expected_originals = tmp_path / f'{prefix}unique_aligned_fragments.csv'
    assert expected_originals.exists() is use_denovo
    assert coverage_destinations == [str(tmp_path / f'{prefix}coverage_maps')]
    if use_denovo:
        rows = list(csv.DictReader(StringIO(expected_depth.read_text())))
        assert len(rows) == 10
        assert {r['avg_depth'] for r in rows} == {'0.5'}
        fragment_rows = list(csv.DictReader(
            StringIO(expected_originals.read_text())))
        assert fragment_rows == [{'unique_aligned_fragments': '1'}]
