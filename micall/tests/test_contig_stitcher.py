
import pytest
import random
from micall.core.contig_stitcher import split_contigs_with_gaps, stitch_contigs, GenotypedContig, merge_intervals, find_covered_contig, stitch_consensus, calculate_concordance, align_all_to_reference
from micall.tests.utils import MockAligner, fixed_random_seed


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
                        match_fraction=1.0,
                        ),
        ]

    results = list(stitch_contigs(contigs))
    assert set(map(lambda x: x.seq, contigs)) \
        == set(map(lambda x: x.seq, results))


def test_separate_stitching_of_non_overlapping_contigs(exact_aligner):
    # Scenario: When stitching multiple non-overlapping contigs, the order doesn't matter.

    ref_seq = 'A' * 100

    contigs = [
        GenotypedContig(name='a',
                        seq=ref_seq,
                        ref_name='testref',
                        ref_seq=ref_seq,
                        match_fraction=0.5,
                        ),
        GenotypedContig(name='b',
                        seq='C' * 100,
                        ref_name='testref',
                        ref_seq=ref_seq,
                        match_fraction=0.5,
                        ),
        ]

    results = list(stitch_contigs(contigs))

    # No claims about the output order, so wrap into set()
    assert set(map(lambda x: x.seq, contigs)) \
        == set(map(lambda x: x.seq, results))

    contigs = [
        GenotypedContig(name='b',
                        seq='C' * 100,
                        ref_name='testref',
                        ref_seq=ref_seq,
                        match_fraction=0.5,
                        ),
        GenotypedContig(name='a',
                        seq=ref_seq,
                        ref_name='testref',
                        ref_seq=ref_seq,
                        match_fraction=0.5,
                        ),
        ]

    results = list(stitch_contigs(contigs))

    # No claims about the output order, so wrap into set()
    assert set(map(lambda x: x.seq, contigs)) \
        == set(map(lambda x: x.seq, results))


def test_correct_stitching_of_two_partially_overlapping_contigs(exact_aligner):
    # Scenario: Two partially overlapping contigs are stitched correctly into a single sequence.

    ref_seq = 'A' * 100 + 'C' * 100

    contigs = [
        GenotypedContig(name='a',
                        seq='A' * 50 + 'C' * 20,
                        ref_name='testref',
                        ref_seq=ref_seq,
                        match_fraction=0.5,
                        ),
        GenotypedContig(name='b',
                        seq='A' * 20 + 'C' * 50,
                        ref_name='testref',
                        ref_seq=ref_seq,
                        match_fraction=0.5,
                        ),
        ]

    results = list(stitch_contigs(contigs))
    assert len(results) == 1

    result = results[0]

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
                        match_fraction=0.5,
                        ),
        GenotypedContig(name='b',
                        seq='A' * 20 + 'C' * 50,
                        ref_name='testref',
                        ref_seq=ref_seq,
                        match_fraction=0.5,
                        ),
        GenotypedContig(name='c',
                        seq='C' * 20 + 'T' * 50,
                        ref_name='testref',
                        ref_seq=ref_seq,
                        match_fraction=0.5,
                        ),
        ]

    results = list(stitch_contigs(contigs))
    assert len(results) == 2

    assert 100 == len(results[0].seq)
    assert results[0].seq == 'A' * 50 + 'C' * 50
    assert results[0].query.name == 'left(a)+overlap(a,b)+right(b)'

    assert results[1].query == contigs[2]


def test_stitching_of_all_overlapping_contigs_into_one_sequence(exact_aligner):
    # Scenario: All contigs have some overlapping parts, resulting in one continuous sequence after stitching.

    ref_seq = 'A' * 100 + 'C' * 100 + 'T' * 100

    contigs = [
        GenotypedContig(name='a',
                        seq='A' * 50 + 'C' * 20,
                        ref_name='testref',
                        ref_seq=ref_seq,
                        match_fraction=0.5,
                        ),
        GenotypedContig(name='b',
                        seq='A' * 20 + 'C' * 100 + 'T' * 20,
                        ref_name='testref',
                        ref_seq=ref_seq,
                        match_fraction=0.5,
                        ),
        GenotypedContig(name='c',
                        seq='C' * 20 + 'T' * 50,
                        ref_name='testref',
                        ref_seq=ref_seq,
                        match_fraction=0.5,
                        ),
        ]

    results = list(stitch_contigs(contigs))
    assert len(results) == 1

    result = results[0]

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
                        match_fraction=0.5,
                        ),
        GenotypedContig(name='b',
                        seq='',
                        ref_name='testref',
                        ref_seq=ref_seq,
                        match_fraction=0.5,
                        ),
        ]

    results = list(stitch_contigs(contigs))
    assert set(map(lambda x: x.seq, contigs)) \
        == set(map(lambda x: x.seq, results))


def test_stitching_of_identical_contigs(exact_aligner):
    # Scenario: The function correctly handles and avoids duplication when identical contigs are stitched together.

    ref_seq = 'A' * 100
    contigs = [
        GenotypedContig(name=name,
                        seq='ACTGACTG' * 100,
                        ref_name='testref',
                        ref_seq='ACTGACTG' * 100,
                        match_fraction=1.0,
                        )
        for name in ["a", "b", "c"]]

    results = list(stitch_contigs(contigs))
    assert len(results) == 1
    assert results[0].query == contigs[2]


def test_stitching_of_zero_contigs(exact_aligner):
    # Scenario: The function does not crash if no contigs given.

    contigs = []
    results = list(stitch_contigs(contigs))
    assert results == contigs


def test_correct_stitching_of_two_partially_overlapping_different_organism_contigs(exact_aligner):
    # Scenario: Two partially overlapping contigs, but which come from different organism,
    # are not stitched into a single sequence.

    ref_seq = 'A' * 100 + 'C' * 100

    contigs = [
        GenotypedContig(name='a',
                        seq='A' * 50 + 'C' * 20,
                        ref_name='testref-1',
                        ref_seq=ref_seq,
                        match_fraction=0.5,
                        ),
        GenotypedContig(name='b',
                        seq='A' * 20 + 'C' * 50,
                        ref_name='testref-2',
                        ref_seq=ref_seq,
                        match_fraction=0.5,
                        ),
        ]

    results = list(stitch_contigs(contigs))
    assert len(results) == 2

    assert set(map(lambda x: x.seq, contigs)) \
        == set(map(lambda x: x.seq, results))


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
                        match_fraction=0.5,
                        ),
        GenotypedContig(name='b',
                        seq='A' * 20 + 'C' * 50,
                        ref_name=ref_name,
                        ref_seq=ref_seq,
                        match_fraction=0.5,
                        ),
        GenotypedContig(name='c',
                        seq='C' * 70 + 'T' * 20,
                        ref_name=ref_name,
                        ref_seq=ref_seq,
                        match_fraction=0.5,
                        ),
        GenotypedContig(name='d',
                        seq='T' * 20 + 'G' * 50,
                        ref_name=ref_name,
                        ref_seq=ref_seq,
                        match_fraction=0.5,
                        ),
        ] for ref_name in ['testref-1', 'testref-2']]

    contigs = sum(contigs, start=[])

    results = list(stitch_contigs(contigs))
    assert len(results) == 4

    assert 170 == len(results[0].seq)
    assert results[0].seq == 'A' * 50 + 'C' * 100 + 'T' * 20
    assert results[0].query.name == 'left(a)+overlap(a,b)+left(right(b))+overlap(left(a)+overlap(a,b)+right(b),c)+right(c)'
    assert results[0].query.ref_name == 'testref-1'

    assert 170 == len(results[1].seq)
    assert results[1].seq == 'A' * 50 + 'C' * 100 + 'T' * 20
    assert results[1].query.name == 'left(a)+overlap(a,b)+left(right(b))+overlap(left(a)+overlap(a,b)+right(b),c)+right(c)'
    assert results[1].query.ref_name == 'testref-2'

    assert results[2].query == contigs[3]
    assert results[3].query == contigs[7]


def test_stitching_when_one_contig_completely_covered_by_another(exact_aligner):
    # Scenario: If one contig is completely covered by another contig,
    # the completely covered contig must be dropped.

    ref_seq = 'A' * 100 + 'C' * 100

    contigs = [
        GenotypedContig(name='a',
                        seq='A' * 20 + 'C' * 20,
                        ref_name='testref',
                        ref_seq=ref_seq,
                        match_fraction=0.5,
                        ),
        GenotypedContig(name='b',
                        seq='A' * 50 + 'C' * 50,
                        ref_name='testref',
                        ref_seq=ref_seq,
                        match_fraction=0.5,
                        ),
        ]

    results = list(stitch_contigs(contigs))
    assert len(results) == 1

    # Test to ensure that the final result contains the contig 'b' and
    # does not contain the completely covered contig 'a'.
    assert results[0].query.name == 'b'
    assert results[0].query == contigs[1]


def test_stitching_contig_with_big_noncovered_gap(exact_aligner):
    # Scenario: One contig has a big gap, which is however not covered by anything else.

    ref_seq = 'A' * 100 + 'C' * 100 + 'T' * 100

    contigs = [
        GenotypedContig(name='a',
                        seq= 'A' * 50 + 'T' * 50, # mind the C gap
                        ref_name='testref',
                        ref_seq=ref_seq,
                        match_fraction=0.5,
                        ),
        ]

    results = list(stitch_contigs(contigs))

    assert set(map(lambda x: x.seq, contigs)) \
        == set(map(lambda x: x.seq, results))


def test_stitching_contig_with_big_noncovered_gap_2(exact_aligner):
    # Scenario: One contig has a big gap, which is however not covered by anything else.

    ref_seq = 'A' * 100 + 'C' * 100 + 'T' * 100 + 'G' * 100

    contigs = [
        GenotypedContig(name='a',
                        seq='A' * 50 + 'T' * 50, # mind the C gap
                        ref_name='testref',
                        ref_seq=ref_seq,
                        match_fraction=0.5,
                        ),
        GenotypedContig(name='B',
                        seq='G' * 50,
                        ref_name='testref',
                        ref_seq=ref_seq,
                        match_fraction=0.5,
                        ),
        ]

    results = list(stitch_contigs(contigs))

    assert set(map(lambda x: x.seq, contigs)) \
        == set(map(lambda x: x.seq, results))


def test_stitching_contig_with_big_covered_gap(exact_aligner):
    # Scenario: If one contig has a big gap covered by another contig.

    ref_seq = 'G' * 100 + 'A' * 100 + 'C' * 100 + 'T' * 100 + 'G' * 100

    contigs = [
        GenotypedContig(name='a',
                        seq='G' * 50 + 'A' * 50 + 'T' * 100, # mind the gap
                        ref_name='testref',
                        ref_seq=ref_seq,
                        match_fraction=0.5,
                        ),
        GenotypedContig(name='b',
                        seq='A' * 100 + 'C' * 100 + 'T' * 100 + 'G' * 50,
                        ref_name='testref',
                        ref_seq=ref_seq,
                        match_fraction=0.5,
                        ),
        ]

    contigs = align_all_to_reference(contigs)
    assert len(list(contigs[0].alignment.gaps())) == 1
    assert len(list(contigs[1].alignment.gaps())) == 0

    results = list(split_contigs_with_gaps(contigs))
    assert len(results) == 3
    assert all(list(contig.alignment.gaps()) == [] for contig in results)


def test_stitching_contig_with_small_covered_gap(exact_aligner):
    # Scenario: If one contig has a small gap covered by another contig.

    ref_seq = 'G' * 100 + 'A' * 9 + 'C' * 100 + 'T' * 100

    contigs = [
        GenotypedContig(name='a',
                        seq='G' * 100 + 'A' * 0 + 'C' * 100, # mind the gap
                        ref_name='testref',
                        ref_seq=ref_seq,
                        match_fraction=0.5,
                        ),
        GenotypedContig(name='b',
                        seq='A' * 9 + 'C' * 100 + 'T' * 100,
                        ref_name='testref',
                        ref_seq=ref_seq,
                        match_fraction=0.5,
                        ),
        ]

    contigs = align_all_to_reference(contigs)
    assert len(list(contigs[0].alignment.gaps())) == 1
    assert len(list(contigs[1].alignment.gaps())) == 0

    results = list(split_contigs_with_gaps(contigs))

    assert all(x.seq == x.lstrip_query().rstrip_query().seq for x in results)

    assert set(map(lambda x: x.seq, contigs)) \
        == set(map(lambda x: x.seq, results))


def test_stitching_partial_align(exact_aligner):
    # Scenario: A single contig has a sequence that partially aligns to the reference sequence.

    contigs = [
        GenotypedContig(name='a',
                        seq='T' * 10 + 'C' * 20 + 'A' * 10,
                        ref_name='testref',
                        ref_seq='A' * 20 + 'C' * 20 + 'T' * 20,
                        match_fraction=0.3,
                        ),
        ]

    results = list(stitch_contigs(contigs))
    assert len(results) == len(contigs)
    for result in results:
        assert any(result.seq in contig.seq for contig in contigs)

    assert all(x.seq != x.lstrip_query().rstrip_query().seq for x in results)

    assert set(map(lambda x: x.seq, contigs)) \
        != set(map(lambda x: x.lstrip_query().rstrip_query().seq, results))


def test_partial_align_consensus(exact_aligner):
    # Scenario: A single contig partially aligns to the reference sequence, and a consensus sequence is being stitched.

    contigs = [
        GenotypedContig(name='a',
                        seq='T' * 10 + 'C' * 20 + 'A' * 10,
                        ref_name='testref',
                        ref_seq='A' * 20 + 'C' * 20 + 'T' * 20,
                        match_fraction=0.3,
                        ),
        ]

    results = list(stitch_consensus(contigs))
    assert len(results) == len(contigs)
    assert set(map(lambda x: x.seq, contigs)) \
        == set(map(lambda x: x.seq, results))


def test_stitching_partial_align_multiple_sequences(exact_aligner):
    # Scenario: Multiple contigs have sequences that partially align to the same reference sequence.

    ref_seq='A' * 20 + 'C' * 20 + 'T' * 20

    contigs = [
        GenotypedContig(name='a',
                        seq='T' * 10 + 'C' * 20 + 'A' * 10,
                        ref_name='testref',
                        ref_seq=ref_seq,
                        match_fraction=0.3,
                        ),
        GenotypedContig(name='b',
                        seq='C' * 20 + 'A' * 10 + 'G' * 10,
                        ref_name='testref',
                        ref_seq=ref_seq,
                        match_fraction=0.3,
                        ),
        ]

    results = list(stitch_contigs(contigs))
    assert len(results) == 1
    for result in results:
        assert any(result.seq in contig.seq for contig in contigs)

    assert set(map(lambda x: x.seq, contigs)) \
        != set(map(lambda x: x.lstrip_query().rstrip_query().seq, results))


def test_partial_align_consensus_multiple_sequences(exact_aligner):
    # Scenario: Multiple contigs partially align to the same reference sequence, and a consensus sequence is being stitched from them.

    ref_seq='A' * 20 + 'C' * 20 + 'T' * 20

    contigs = [
        GenotypedContig(name='a',
                        seq='A' * 20,
                        ref_name='testref',
                        ref_seq=ref_seq,
                        match_fraction=0.3,
                        ),
        GenotypedContig(name='b',
                        seq='T' * 20,
                        ref_name='testref',
                        ref_seq=ref_seq,
                        match_fraction=0.3,
                        ),
        ]

    results = list(stitch_consensus(contigs))
    assert len(results) == 1
    assert results[0].seq == contigs[0].seq + contigs[1].seq
    assert results[0].name == 'a+b'


def test_partial_align_consensus_multiple_overlaping_sequences(exact_aligner):
    # Scenario: Multiple contigs partially align to the same reference sequence, and a consensus sequence is being stitched from them.

    ref_seq='A' * 20 + 'C' * 20 + 'T' * 20

    contigs = [
        GenotypedContig(name='a',
                        seq='T' * 10 + 'A' * 5 + 'C' * 20 + 'A' * 10,
                        ref_name='testref',
                        ref_seq=ref_seq,
                        match_fraction=0.3,
                        ),
        GenotypedContig(name='b',
                        seq='C' * 20 + 'T' * 5 + 'A' * 10 + 'G' * 10,
                        ref_name='testref',
                        ref_seq=ref_seq,
                        match_fraction=0.3,
                        ),
        ]

    results = list(stitch_consensus(contigs))
    assert len(results) == 1
    assert results[0].seq == 'T' * 10 + 'A' * 5 + 'C' * 20 + 'T' * 5 + 'A' * 10 + 'G' * 10
    assert results[0].seq == contigs[0].seq[:-10] + contigs[1].seq[20:]
    assert results[0].name == 'left(a)+overlap(a,b)+right(b)'


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


def test_concordance_same_length_inputs():
    with pytest.raises(ValueError):
        calculate_concordance('abc', 'ab')

def test_concordance_completely_different_strings():
    result = calculate_concordance('a'*30, 'b'*30)
    assert all(n == 0 for n in result)

def generate_random_string_pair(length):
    left = ''.join(random.choice('ACGT') for _ in range(length))
    right = ''.join(random.choice('ACGT') for _ in range(length))
    return left, right

def generate_test_cases(num_cases):
    with fixed_random_seed(42):
        length = random.randint(1, 80)
        return [generate_random_string_pair(length) for _ in range(num_cases)]

concordance_cases = generate_test_cases(num_cases=100)


@pytest.mark.parametrize('left, right', concordance_cases)
def test_concordance_output_is_list_of_floats(left, right):
    result = calculate_concordance(left, right)
    assert isinstance(result, list), "Result should be a list"
    assert all(isinstance(n, float) for n in result), "All items in result should be float"


@pytest.mark.parametrize('left, right', concordance_cases)
def test_concordance_output_range(left, right):
    result = calculate_concordance(left, right)
    assert all(0 <= n <= 1 for n in result), "All values in result should be between 0 and 1"


@pytest.mark.parametrize('left, right', concordance_cases)
def test_concordance_higher_if_more_matches_added(left, right):
    # Insert exact matches in the middle
    matching_sequence = 'A' * 30
    insert_position = len(left) // 2
    new_left = left[:insert_position] + matching_sequence + left[insert_position + len(matching_sequence):]
    new_right = right[:insert_position] + matching_sequence + right[insert_position + len(matching_sequence):]

    old_conc = calculate_concordance(left, right)
    new_conc = calculate_concordance(new_left, new_right)
    old_average = sum(old_conc) / len(old_conc)
    new_average = sum(new_conc) / len(new_conc)
    assert old_average <= new_average


@pytest.mark.parametrize('left, right', concordance_cases)
def test_concordance_higher_in_matching_areas(left, right):
    # Insert exact matches in the middle
    matching_sequence = 'A' * 30
    insert_position = len(left) // 2
    new_left = left[:insert_position] + matching_sequence + left[insert_position + len(matching_sequence):]
    new_right = right[:insert_position] + matching_sequence + right[insert_position + len(matching_sequence):]

    concordance_scores = calculate_concordance(new_left, new_right)

    # Check concordance in the matching area
    matching_area_concordance = concordance_scores[insert_position:insert_position + len(matching_sequence)]

    # Calculate average concordance inside and outside the matching area
    average_inside = sum(matching_area_concordance) / len(matching_sequence)
    average_outside = (sum(concordance_scores) - sum(matching_area_concordance)) / (len(concordance_scores) - len(matching_sequence))

    # Assert that the concordance is indeed higher in the matching area
    assert average_inside > average_outside, "Concordance in matching areas should be higher than in non-matching areas"
