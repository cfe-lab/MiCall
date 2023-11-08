
import pytest
from micall.core.contig_stitcher import stitch_contigs, GenotypedContig
from micall.tests.utils import MockAligner


@pytest.fixture(autouse=True)
def mock_mappy_aligner(monkeypatch):
    monkeypatch.setattr('micall.core.contig_stitcher.Aligner', MockAligner)


def test_identical_stitching_of_one_contig():
    # Scenario: When stitching one contig, it remains the same.

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


def test_separate_stitching_of_non_overlapping_contigs():
    # Scenario: When stitching multiple non-overlapping contigs, the order doesn't matter.

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

    contigs = [
        GenotypedContig(name='b',
                        seq='C' * 100,
                        ref_name='testref',
                        ref_seq=ref_seq,
                        matched_fraction=0.5,
                        ),
        GenotypedContig(name='a',
                        seq=ref_seq,
                        ref_name='testref',
                        ref_seq=ref_seq,
                        matched_fraction=0.5,
                        ),
        ]

    result = list(stitch_contigs(contigs))
    assert sorted(map(lambda x: x.seq, contigs)) \
        == sorted(map(lambda x: x.seq, result))


def test_correct_stitching_of_two_partially_overlapping_contigs():
    # Scenario: Two partially overlapping contigs are stitched correctly into a single sequence.

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
    assert result.query.name == 'a+overlap(a,b)+b'


def test_correct_processing_of_two_overlapping_and_one_separate_contig():
    # Scenario: Two overlapping contigs are stitched together, the non-overlapping is kept separate.

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
    assert result[0].query.name == 'a+overlap(a,b)+b'

    assert result[1].query == contigs[2]


def test_stitching_of_all_overlapping_contigs_into_one_sequence():
    # Scenario: All contigs have some overlapping parts, resulting in one continuous sequence after stitching.

    ref_seq = 'A' * 100 + 'C' * 100 + 'T' * 100

    contigs = [
        GenotypedContig(name='a',
                        seq='A' * 50 + 'C' * 20,
                        ref_name='testref',
                        ref_seq=ref_seq,
                        matched_fraction=0.5,
                        ),
        GenotypedContig(name='b',
                        seq='A' * 20 + 'C' * 100 + 'T' * 20,
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
    assert len(result) == 1

    result = result[0]

    assert 200 == len(result.seq)
    assert result.seq == 'A' * 50 + 'C' * 100 + 'T' * 50
    assert result.query.name == 'a+overlap(a,b)+b+overlap(a+overlap(a,b)+b,c)+c'


def test_stitching_with_empty_contigs():
    # Scenario: The function is able to handle and ignore empty contigs.

    ref_seq = 'A' * 100

    contigs = [
        GenotypedContig(name='a',
                        seq=ref_seq,
                        ref_name='testref',
                        ref_seq=ref_seq,
                        matched_fraction=0.5,
                        ),
        GenotypedContig(name='b',
                        seq='',
                        ref_name='testref',
                        ref_seq=ref_seq,
                        matched_fraction=0.5,
                        ),
        ]

    result = list(stitch_contigs(contigs))
    assert sorted(map(lambda x: x.seq, contigs)) \
        == sorted(map(lambda x: x.seq, result))


def test_stitching_of_identical_contigs():
    # Scenario: The function correctly handles and avoids duplication when identical contigs are stitched together.

    ref_seq = 'A' * 100
    contigs = [
        GenotypedContig(name=name,
                        seq='ACTGACTG' * 100,
                        ref_name='testref',
                        ref_seq='ACTGACTG' * 100,
                        matched_fraction=1.0,
                        )
        for name in ["a", "b", "c"]]

    result = list(stitch_contigs(contigs))
    assert len(result) == 1
    assert result[0].query == contigs[2]


def test_stitching_of_zero_contigs():
    # Scenario: The function does not crash if no contigs given.

    contigs = []
    result = list(stitch_contigs(contigs))
    assert result == contigs


def test_correct_stitching_of_two_partially_overlapping_different_organism_contigs():
    # Scenario: Two partially overlapping contigs, but which come from different organism,
    # are not stitched into a single sequence.

    ref_seq = 'A' * 100 + 'C' * 100

    contigs = [
        GenotypedContig(name='a',
                        seq='A' * 50 + 'C' * 20,
                        ref_name='testref-1',
                        ref_seq=ref_seq,
                        matched_fraction=0.5,
                        ),
        GenotypedContig(name='b',
                        seq='A' * 20 + 'C' * 50,
                        ref_name='testref-2',
                        ref_seq=ref_seq,
                        matched_fraction=0.5,
                        ),
        ]

    result = list(stitch_contigs(contigs))
    assert len(result) == 2

    assert sorted(map(lambda x: x.seq, contigs)) \
        == sorted(map(lambda x: x.seq, result))


def test_correct_processing_complex_nogaps():
    # Scenario: There are two reference organisms.
    # Each with 4 contigs.
    # For each, three overlapping contigs are stitched together, the non-overlapping is kept separate.
    # This seems like the most general scenario if no gaps or complete goverage is involved.

    ref_seq = 'A' * 100 + 'C' * 100 + 'T' * 100 + 'G' * 100

    contigs = [[
        GenotypedContig(name='a',
                        seq='A' * 50 + 'C' * 20,
                        ref_name=ref_name,
                        ref_seq=ref_seq,
                        matched_fraction=0.5,
                        ),
        GenotypedContig(name='b',
                        seq='A' * 20 + 'C' * 50,
                        ref_name=ref_name,
                        ref_seq=ref_seq,
                        matched_fraction=0.5,
                        ),
        GenotypedContig(name='c',
                        seq='C' * 70 + 'T' * 20,
                        ref_name=ref_name,
                        ref_seq=ref_seq,
                        matched_fraction=0.5,
                        ),
        GenotypedContig(name='d',
                        seq='T' * 20 + 'G' * 50,
                        ref_name=ref_name,
                        ref_seq=ref_seq,
                        matched_fraction=0.5,
                        ),
        ] for ref_name in ['testref-1', 'testref-2']]

    contigs = sum(contigs, start=[])

    result = list(stitch_contigs(contigs))
    assert len(result) == 4

    assert 170 == len(result[0].seq)
    assert result[0].seq == 'A' * 50 + 'C' * 100 + 'T' * 20
    assert result[0].query.name == 'a+overlap(a,b)+b+overlap(a+overlap(a,b)+b,c)+c'
    assert result[0].query.ref_name == 'testref-1'

    assert 170 == len(result[1].seq)
    assert result[1].seq == 'A' * 50 + 'C' * 100 + 'T' * 20
    assert result[1].query.name == 'a+overlap(a,b)+b+overlap(a+overlap(a,b)+b,c)+c'
    assert result[1].query.ref_name == 'testref-2'

    assert result[2].query == contigs[3]
    assert result[3].query == contigs[7]
