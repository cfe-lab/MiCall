
import pytest
from micall.core.contig_stitcher import split_contigs_with_gaps, stitch_contigs, GenotypedContig, merge_intervals, find_covered_contig
from micall.tests.utils import MockAligner


@pytest.fixture()
def exact_aligner(monkeypatch):
    monkeypatch.setattr('micall.core.contig_stitcher.Aligner', MockAligner)


def test_identical_stitching_of_one_contig(exact_aligner):
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
    assert set(map(lambda x: x.seq, contigs)) \
        == set(map(lambda x: x.seq, result))


def test_separate_stitching_of_non_overlapping_contigs(exact_aligner):
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

    # No claims about the output order, so wrap into set()
    assert set(map(lambda x: x.seq, contigs)) \
        == set(map(lambda x: x.seq, result))

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

    # No claims about the output order, so wrap into set()
    assert set(map(lambda x: x.seq, contigs)) \
        == set(map(lambda x: x.seq, result))


def test_correct_stitching_of_two_partially_overlapping_contigs(exact_aligner):
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
    assert result.query.name == 'left(a)+overlap(a,b)+right(b)'


def test_correct_processing_of_two_overlapping_and_one_separate_contig(exact_aligner):
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
    assert result[0].query.name == 'left(a)+overlap(a,b)+right(b)'

    assert result[1].query == contigs[2]


def test_stitching_of_all_overlapping_contigs_into_one_sequence(exact_aligner):
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
    assert result.query.name == 'left(a)+overlap(a,b)+left(right(b))+overlap(left(a)+overlap(a,b)+right(b),c)+right(c)'


def test_stitching_with_empty_contigs(exact_aligner):
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
    assert set(map(lambda x: x.seq, contigs)) \
        == set(map(lambda x: x.seq, result))


def test_stitching_of_identical_contigs(exact_aligner):
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


def test_stitching_of_zero_contigs(exact_aligner):
    # Scenario: The function does not crash if no contigs given.

    contigs = []
    result = list(stitch_contigs(contigs))
    assert result == contigs


def test_correct_stitching_of_two_partially_overlapping_different_organism_contigs(exact_aligner):
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

    assert set(map(lambda x: x.seq, contigs)) \
        == set(map(lambda x: x.seq, result))


def test_correct_processing_complex_nogaps(exact_aligner):
    # Scenario: There are two reference organisms.
    # Each with 4 contigs.
    # For each, three overlapping contigs are stitched together, the non-overlapping is kept separate.

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
    assert result[0].query.name == 'left(a)+overlap(a,b)+left(right(b))+overlap(left(a)+overlap(a,b)+right(b),c)+right(c)'
    assert result[0].query.ref_name == 'testref-1'

    assert 170 == len(result[1].seq)
    assert result[1].seq == 'A' * 50 + 'C' * 100 + 'T' * 20
    assert result[1].query.name == 'left(a)+overlap(a,b)+left(right(b))+overlap(left(a)+overlap(a,b)+right(b),c)+right(c)'
    assert result[1].query.ref_name == 'testref-2'

    assert result[2].query == contigs[3]
    assert result[3].query == contigs[7]


def test_stitching_when_one_contig_completely_covered_by_another(exact_aligner):
    # Scenario: If one contig is completely covered by another contig,
    # the completely covered contig must be dropped.

    ref_seq = 'A' * 100 + 'C' * 100

    contigs = [
        GenotypedContig(name='a',
                        seq='A' * 20 + 'C' * 20,
                        ref_name='testref',
                        ref_seq=ref_seq,
                        matched_fraction=0.5,
                        ),
        GenotypedContig(name='b',
                        seq='A' * 50 + 'C' * 50,
                        ref_name='testref',
                        ref_seq=ref_seq,
                        matched_fraction=0.5,
                        ),
        ]

    result = list(stitch_contigs(contigs))
    assert len(result) == 1

    # Test to ensure that the final result contains the contig 'b' and
    # does not contain the completely covered contig 'a'.
    assert result[0].query.name == 'b'
    assert result[0].query == contigs[1]


def test_stitching_contig_with_big_noncovered_gap(exact_aligner):
    # Scenario: One contig has a big gap, which is however not covered by anything else.

    ref_seq = 'A' * 100 + 'C' * 100 + 'T' * 100

    contigs = [
        GenotypedContig(name='a',
                        seq= 'A' * 50 + 'T' * 50, # mind the C gap
                        ref_name='testref',
                        ref_seq=ref_seq,
                        matched_fraction=0.5,
                        ),
        ]

    result = list(stitch_contigs(contigs))

    assert set(map(lambda x: x.seq, contigs)) \
        == set(map(lambda x: x.seq, result))


def test_stitching_contig_with_big_noncovered_gap_2(exact_aligner):
    # Scenario: One contig has a big gap, which is however not covered by anything else.

    ref_seq = 'A' * 100 + 'C' * 100 + 'T' * 100 + 'G' * 100

    contigs = [
        GenotypedContig(name='a',
                        seq='A' * 50 + 'T' * 50, # mind the C gap
                        ref_name='testref',
                        ref_seq=ref_seq,
                        matched_fraction=0.5,
                        ),
        GenotypedContig(name='B',
                        seq='G' * 50,
                        ref_name='testref',
                        ref_seq=ref_seq,
                        matched_fraction=0.5,
                        ),
        ]

    result = list(stitch_contigs(contigs))

    assert set(map(lambda x: x.seq, contigs)) \
        == set(map(lambda x: x.seq, result))


def test_stitching_contig_with_big_covered_gap(exact_aligner):
    # Scenario: If one contig has a big gap covered by another contig.

    ref_seq = 'G' * 100 + 'A' * 100 + 'C' * 100 + 'T' * 100 + 'G' * 100

    contigs = [
        GenotypedContig(name='a',
                        seq='G' * 50 + 'A' * 50 + 'T' * 100, # mind the gap
                        ref_name='testref',
                        ref_seq=ref_seq,
                        matched_fraction=0.5,
                        ),
        GenotypedContig(name='b',
                        seq='A' * 100 + 'C' * 100 + 'T' * 100 + 'G' * 50,
                        ref_name='testref',
                        ref_seq=ref_seq,
                        matched_fraction=0.5,
                        ),
        ]

    contigs = [x.align_to_reference() for x in contigs]
    assert len(list(contigs[0].gaps())) == 1
    assert len(list(contigs[1].gaps())) == 0

    result = list(split_contigs_with_gaps(contigs))
    assert len(result) == 3
    assert all(list(contig.gaps()) == [] for contig in result)


def test_stitching_contig_with_small_covered_gap(exact_aligner):
    # Scenario: If one contig has a small gap covered by another contig.

    ref_seq = 'G' * 100 + 'A' * 9 + 'C' * 100 + 'T' * 100

    contigs = [
        GenotypedContig(name='a',
                        seq='G' * 100 + 'A' * 0 + 'C' * 100, # mind the gap
                        ref_name='testref',
                        ref_seq=ref_seq,
                        matched_fraction=0.5,
                        ),
        GenotypedContig(name='b',
                        seq='A' * 9 + 'C' * 100 + 'T' * 100,
                        ref_name='testref',
                        ref_seq=ref_seq,
                        matched_fraction=0.5,
                        ),
        ]

    contigs = [x.align_to_reference() for x in contigs]
    assert len(list(contigs[0].gaps())) == 1
    assert len(list(contigs[1].gaps())) == 0

    result = list(split_contigs_with_gaps(contigs))

    assert set(map(lambda x: x.seq, contigs)) \
        == set(map(lambda x: x.seq, result))



#  _   _       _ _     _            _
# | | | |_ __ (_) |_  | |_ ___  ___| |_ ___
# | | | | '_ \| | __| | __/ _ \/ __| __/ __|
# | |_| | | | | | |_  | ||  __/\__ \ |_\__ \
#  \___/|_| |_|_|\__|  \__\___||___/\__|___/
#

@pytest.mark.parametrize("intervals, expected", [
    ([], []),
    ([(1, 3)], [(1, 3)]),

    # Non-overlapping intervals
    ([(1, 3), (5, 6)], [(1, 3), (5, 6)]),

    # Directly overlapping intervals
    ([(1, 3), (2, 5)], [(1, 5)]),

    # Adjacent intervals that exactly touch each other
    ([(1, 2), (3, 4)], [(1, 4)]),

    # Nested intervals
    ([(1, 10), (2, 5)], [(1, 10)]),

    # Multiple merged intervals
    ([(1, 3), (2, 4), (6, 8), (10, 11), (11, 12)],
     [(1, 4), (6, 8), (10, 12)]),

    # Intervals out of initial order
    ([(4, 6), (1, 2)],
     [(1, 2), (4, 6)]),

    # Overlapping intervals with out of order inputs
    ([(1, 4), (3, 5), (2, 3), (7, 10), (9, 12)],
     [(1, 5), (7, 12)]),

    # Large set of intervals with various overlaps
    ([(1, 4), (2, 6), (5, 8), (7, 8), (10, 15), (11, 12), (13, 14), (17, 18)],
     [(1, 8), (10, 15), (17, 18)]),

    # Intervals where end is less than start should return as is or be handled explicitly depending on implementation
    ([(5, 3), (1, 2)],
     [(1, 2), (5, 3)]),

    # Intervals that are exactly one after the other in sequence / Intervals that are completely disjoint
    ([(1, 2), (4, 5), (7, 8)],
     [(1, 2), (4, 5), (7, 8)]),

    # Overlapping intervals that merge into one large interval
    ([(2, 6), (4, 10), (5, 15), (14, 20)],
     [(2, 20)]),

    # Same interval repeated multiple times
    ([(1, 5), (1, 5), (1, 5)],
     [(1, 5)]),

    # Single point intervals
    ([(1, 1), (5, 5), (3, 3)],
     [(1, 1), (3, 3), (5, 5)]),

    ([(1, 1), (5, 5), (3, 3), (1, 1), (1, 1)],
     [(1, 1), (3, 3), (5, 5)]),

    ([(1, 1), (2, 3)],
     [(1, 3)]),

    # Intervals that start with negative numbers
    ([(-5, 0), (-2, 3), (1, 7), (9, 12)],
     [(-5, 7), (9, 12)]),
])
def test_merge_intervals(intervals, expected):
    assert merge_intervals(intervals) == expected


class MockAlignedContig:
    def __init__(self, ref_name, r_st, r_ei, name="contig"):
        self.ref_name = ref_name
        self.alignment = MockAlignment(r_st, r_ei)
        self.name = name


class MockAlignment:
    def __init__(self, r_st, r_ei):
        self.r_st = r_st
        self.r_ei = r_ei


# Simple function to create mock AlignedContig objects for testing, including ref_name.
def create_mock_aligned_contig(ref_name, r_st, r_ei, name="contig"):
    return MockAlignedContig(ref_name, r_st, r_ei, name)


@pytest.mark.parametrize("contigs, expected_covered_name", [
    # No contigs are completely covered.
    ([('ref1', 0, 100), ('ref1', 101, 200)], None),
    ([('ref1', 0, 50), ('ref1', 51, 100)], None),

    # A single contig is completely covered by one other contig.
    ([('ref1', 0, 100), ('ref1', 0, 200)], 'contig1'),
    ([('ref1', 50, 150), ('ref1', 0, 200)], 'contig1'),

    # A single contig completely covers another, but with different reference names.
    ([('ref1', 0, 50), ('ref2', 0, 100)], None),

    # Single coverage with exact match.
    ([('ref1', 0, 100), ('ref1', 0, 100)], 'contig1'),

    # A single contig is completely covered at the beginning by one and at the end by another contig.
    ([('ref1', 0, 50), ('ref1', 50, 100), ('ref1', 25, 75)], 'contig3'),

    # Contigs overlap but none are completely covered.
    ([('ref1', 0, 50), ('ref1', 40, 90), ('ref1', 80, 120)], None),

    # Multiple contigs with some covered completely by a single other contig.
    ([('ref1', 0, 200), ('ref1', 10, 30), ('ref1', 170, 190)], 'contig2'),

    # Multiple contigs with complex overlaps and one completely covered.
    ([('ref1', 30, 60), ('ref1', 0, 50), ('ref1', 20, 70), ('ref1', 60, 90)], 'contig1'),

    # Edge case where a contig starts where another ends.
    ([('ref1', 0, 50), ('ref1', 50, 100)], None),

    # Contigs are completely covered in a nested fashion.
    ([('ref1', 0, 200), ('ref1', 50, 150), ('ref1', 100, 125)], 'contig2'),

    # Contigs are adjacent and cover each other completely.
    ([('ref1', 0, 100), ('ref1', 101, 200), ('ref1', 0, 200)], 'contig1'),

    # Single large contig covers several smaller non-adjacent contigs.
    ([('ref1', 0, 500), ('ref1', 50, 100), ('ref1', 200, 250), ('ref1', 300, 350)], 'contig2'),

    # Single large contig covers several smaller adjacent contigs.
    ([('ref1', 50, 100), ('ref1', 70, 300), ('ref1', 101, 199), ('ref1', 200, 350)], 'contig2'),

    # Single small contig is covered by several larger contigs.
    ([('ref1', 0, 250), ('ref1', 200, 300), ('ref1', 600, 800), ('ref1', 250, 700)], 'contig2'),

    # Complex case with multiple contigs and complete coverage by combinations.
    ([('ref1', 0, 100), ('ref1', 30, 130), ('ref1', 60, 160), ('ref1', 90, 190), ('ref1', 120, 220)], 'contig2'),

    # Contigs with same start but different end, where one is covered.
    ([('ref1', 0, 100), ('ref1', 0, 50)], 'contig2'),

    # Contigs with same end but different start, where one is covered.
    ([('ref1', 50, 100), ('ref1', 0, 100)], 'contig1'),

    # Contig covered by two overlapping contigs that don't individually cover the whole range.
    ([('ref1', 0, 75), ('ref1', 25, 100), ('ref1', 0, 100)], 'contig1'),

    # Two contigs are covered completely by one large contig.
    ([('ref1', 0, 300), ('ref1', 50, 100), ('ref1', 200, 250)], 'contig2'),

    # No contigs at all.
    ([], None),
])
def test_find_covered(contigs, expected_covered_name):
    mock_contigs = [create_mock_aligned_contig(ref_name, r_st, r_ei, f'contig{i+1}')
                    for i, (ref_name, r_st, r_ei) in enumerate(contigs)]
    covered = find_covered_contig(mock_contigs)
    if expected_covered_name is None:
        assert covered is None
    else:
        assert covered is not None
        assert covered.name == expected_covered_name

