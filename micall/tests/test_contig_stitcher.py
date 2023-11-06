import pytest
from micall.core.contig_stitcher import stitch_contigs, GenotypedContig
from micall.tests.utils import MockAligner


@pytest.fixture(autouse=True)
def mock_mappy_aligner(monkeypatch):
    monkeypatch.setattr('micall.core.contig_stitcher.Aligner', MockAligner)


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
    assert len(result) == 1

    result = result[0]

    assert 100 == len(result.seq)
    assert result.seq == 'A' * 50 + 'C' * 50
    assert result.contig.name == 'a+overlap(a,b)+b'


def test_4():
    ref_seq = 'A' * 100 + 'C' * 100 + 'T' * 100

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
        GenotypedContig(name='c',
                        seq='C' * 20 + 'T' * 50,
                        ref_name='testref',
                        ref_seq=ref_seq,
                        matched_fraction=0.5,
                        ),
        ]

    result = list(stitch_contigs(contigs))
    assert len(result) == 2

    assert 100 == len(result[0].seq)
    assert result[0].seq == 'A' * 50 + 'C' * 50
    assert result[0].contig.name == 'a+overlap(a,b)+b'

    assert result[1].contig == contigs[2]
