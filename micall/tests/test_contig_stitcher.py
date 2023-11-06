import pytest
from micall.core.contig_stitcher import stitch_contigs, GenotypedContig


def test_3():
    ref_seq = 'A' * 100 + 'C' * 100

    contigs = [
        GenotypedContig(name='a',
                        seq='A' * 50 + 'C' * 20,
                        ref_name='testref',
                        ref_seq=ref_seq,
                        matched_fraction=0.5,
                        ),
        GenotypedContig(name='b',
                        seq='A' * 20 + 'C' * 50,
                        ref_name='testref',
                        ref_seq=ref_seq,
                        matched_fraction=0.5,
                        ),
        ]

    result = list(stitch_contigs(contigs))
    assert 100 == sum(len(x.seq) for x in result)
    assert result[0].contig.name == 'a+overlap(a,b)+b'
