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
                aligned.getvalue())

    old_rows, old_depth, _ = perform({'A': dna})
    new_rows, new_depth, aligned = perform(references)
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
