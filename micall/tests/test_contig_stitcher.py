import pytest
from micall.core.contig_stitcher import stitch_contigs, GenotypedContig


def test_1():
    contigs = [
        GenotypedContig(name='a',
                        seq='ACTGACTG' * 100,
                        ref_name='testref',
                        ref_seq='ACTGACTG' * 100,
                        matched_fraction=1.0,
                        ),
        ]

    result = list(stitch_contigs(contigs))
    assert sorted(map(lambda x: x.seq, contigs)) \
        == sorted(map(lambda x: x.seq, result))


def test_2():
    ref_seq = 'A' * 100

    contigs = [
        GenotypedContig(name='a',
                        seq=ref_seq,
                        ref_name='testref',
                        ref_seq=ref_seq,
                        matched_fraction=0.5,
                        ),
        GenotypedContig(name='b',
                        seq='C' * 100,
                        ref_name='testref',
                        ref_seq=ref_seq,
                        matched_fraction=0.5,
                        ),
        ]

    result = list(stitch_contigs(contigs))
    assert sorted(map(lambda x: x.seq, contigs)) \
        == sorted(map(lambda x: x.seq, result))


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
