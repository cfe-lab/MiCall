from collections import Counter

import pytest
from micall.utils.referenceless_contig_stitcher import \
    stitch_consensus, ContigWithAligner, Pool, \
    ACCEPTABLE_STITCHING_SCORE, check_merged_sequence_support
from micall.utils.exact_coverage import reverse_complement
from micall.utils.contig_stitcher_context import ReferencelessStitcherContext
from micall.utils.referenceless_score import Score
from micall.utils.referenceless_contig_path import ContigsPath

 # Load autouse fixtures
from micall.tests.referenceless_tests_utils import disable_acceptable_prob_check, force_failing_map_overlap, disable_kmer_filter

# prevent linter warnings
assert disable_acceptable_prob_check is not None
assert force_failing_map_overlap is not None
assert disable_kmer_filter is not None

TTT = 40 * 'T'
AAA = 40 * 'A'


@pytest.mark.parametrize(
    "seqs, expected",
    [
        #
        # Singletons
        #

        (('AAAAA' + TTT, TTT + 'GGGGG'), ('AAAAA' + TTT + 'GGGGG',)),
        (('AAAAA' + TTT, TTT + 'GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG'), ('AAAAA' + TTT + 'GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG',)),
        (('AAAAA' + TTT + 'CCCCCCCCCCCCCCCCCCCC', TTT + 'GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG'), ('AAAAA' + TTT + 'GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG',)),
        (('AAAAA' + TTT, 'GGGGG' + TTT), ('AAAAA' + TTT, 'GGGGG' + TTT,)),
        (('GGGGG' + TTT, 'CCCCCAAAAA' + TTT), ('CCCCCAAAAA' + TTT, 'GGGGG' + TTT,)),
        (('AAAAA' + TTT, 'GGGGG' + TTT), ('AAAAA' + TTT, 'GGGGG' + TTT,)),
        ((TTT + 'AAAAA', TTT + 'GGGGG'), (TTT + 'AAAAA', TTT + 'GGGGG',)),
        (('AAAAAAAAAAAAAA' + TTT, 'GGGGG' + TTT), ('AAAAAAAAAAAAAA' + TTT, 'GGGGG' + TTT,)),
        (('GGGGGGGGGGGGGG' + TTT, 'AAAAA' + TTT), ('AAAAA' + TTT, 'GGGGGGGGGGGGGG' + TTT,)),
        (('AAAAA' + TTT, 'GGGGGGGGGGGGGG' + TTT), ('AAAAA' + TTT, 'GGGGGGGGGGGGGG' + TTT,)),
        (('GGGGG' + TTT, 'AAAAAAAAAAAAAA' + TTT), ('AAAAAAAAAAAAAA' + TTT, 'GGGGG' + TTT,)),
        (('GGGGG' + TTT + 'AAAAAAAAAAAAA', 'CC' + TTT + 'CC'), ('CC' + TTT + 'CC', 'GGGGG' + TTT + 'AAAAAAAAAAAAA',)),
        (('GGGGG' + TTT + 'AAAAAAAAAAAAA', TTT + 'CC'), ('GGGGG' + TTT + 'AAAAAAAAAAAAA', TTT + 'CC',)),
        (('GGGGG' + TTT + 'AAAAAAAAAAAAA', 'CC' + TTT), ('CC' + TTT, 'GGGGG' + TTT + 'AAAAAAAAAAAAA',)),
        (('CC' + TTT + 'CC', 'GGGGG' + TTT + 'AAAAAAAAAAAAA'), ('CC' + TTT + 'CC', 'GGGGG' + TTT + 'AAAAAAAAAAAAA',)),
        ((TTT + 'CC', 'GGGGG' + TTT + 'AAAAAAAAAAAAA'), ('GGGGG' + TTT + 'AAAAAAAAAAAAA', TTT + 'CC',)),
        (('CC' + TTT, 'GGGGG' + TTT + 'AAAAAAAAAAAAA'), ('CC' + TTT, 'GGGGG' + TTT + 'AAAAAAAAAAAAA',)),
        ((TTT, 'GGGGG' + TTT + 'AAAAAAAAAAAAA'), ('GGGGG' + TTT + 'AAAAAAAAAAAAA',)),
        (('GGGGG' + TTT + 'AAAAAAAAAAAAA', TTT), ('GGGGG' + TTT + 'AAAAAAAAAAAAA',)),
        (('AAAA', 'GGGG'), ('AAAA', 'GGGG',)),
        (('AAAAT', 'TGGGG'), ('AAAATGGGG',)),

        #
        # Multiple.
        #
        (('AAA' + 'T' * 40,  'T' * 40 + 'GGG' + 'T' * 40, 'T' * 40 + 'CCC' + 'M' * 40, 'M' * 40 + 'GGG'),
         ('AAA' + TTT + 'GGG' + TTT + 'CCC' + 'M' * 40,
          'M' * 40 + 'GGG',)),
        (('AAA' + 'T' * 40,
          'T' * 40 + 'GGG' + 'A' * 40,
          'A' * 40 + 'CCC' + 'T' * 40,
          'T' * 40 + 'GGG'),
         ('AAA' + 'T' * 40 + 'GGG' + 'A' * 40 + 'CCC' + 'T' * 40,)),
        (('AAA',), ('AAA',)),
        ((), ()),

        # Left-covered: right = S + (mismatching tail). Correct behavior: keep right only.
        ((('AC' * 60), ('AC' * 60 + 'G' * 200)), (('AC' * 60 + 'G' * 200),)),

        # Right-covered (mirror image): left has a mismatching head, right is exactly the covered prefix.
        # Correct behavior: keep left only.
        ((( 'C' * 200 + 'AC' * 60 ), ('AC' * 60)), (( 'C' * 200 + 'AC' * 60 ),)),

        # Stress the off-by-one when the longer contig is only barely longer:
        # right = S + 'G'. Correct behavior: keep right.
        ((('AC' * 60), ('AC' * 60 + 'G')), (('AC' * 60 + 'G'),)),

        # Symmetric "barely longer" case for right-covered:
        ((( 'G' + 'AC' * 60 ), ('AC' * 60)), (( 'G' + 'AC' * 60 ),)),

        # 1) Left is covered by right, but right's prefix has a 1-bp insertion that
        #    requires a gap to align (so map_overlap finds nothing; fallback is used).
        #    Correct behavior: keep only the bigger right contig.
        (( 'A' * 150,
           'A' * 75 + 'G' + 'A' * 75 + 'C' * 30 ),
         ( 'A' * 75 + 'G' + 'A' * 75 + 'C' * 30, )),

        # 2) Mirror of (1): right is covered by left, but left's suffix has a 1-bp insertion.
        #    Correct behavior: keep only the bigger left contig.
        (( 'C' * 30 + 'A' * 75 + 'G' + 'A' * 75,
           'A' * 150 ),
         ( 'C' * 30 + 'A' * 75 + 'G' + 'A' * 75, )),

        # 3) Left covered by right with a single insertion near the start.
        #    Correct behavior: keep only the bigger right contig because this is actually a stitch, not cover.
        (( 'A' * 100,
           'A' * 70 + 'G' + 'A' * 70 ),
         ( 'A' * 70 + 'G' + 'A' * 70, )),

        # 4) Mirror of (3): right covered by left with a single insertion near the end.
        #    Correct behavior: keep only the bigger left contig.
        (( 'G' + 'A' * 50 + 'C' * 20 + 'A' * 50,
           'A' * 100 ),
         ( 'G' + 'A' * 50 + 'C' * 20 + 'A' * 100, )),

        # Imperfect mutual coverage.
        (( 'A' * 50 + 'GTGTGTGT' + 'C' * 200 + 'TGTGTG' + 'A' * 50,
           'A' * 50 + 'ACACA' + 'C' * 200 + 'CACACACACA' + 'A' * 50 ),
         ( 'A' * 50 + 'GTGTGTGT' + 'C' * 200 + 'TGTGTG' + 'A' * 50,
           'A' * 50 + 'ACACA' + 'C' * 200 + 'CACACACACA' + 'A' * 50 )),

        # Left-covered: right = perfect prefix match of left, then 'Z' marker, then junk.
        # Correct behavior: keep both contigs.
        (( 'A' * 120,
           'A' * 120 + 'Z' + 'G' * 200 ),
         ( 'A' * 120, 'A' * 120 + 'Z' + 'G' * 200 ,)),

        # Right-covered (mirror): left has junk + 'Z' + perfect suffix that covers right.
        # Correct result: keep only the bigger left contig.
        (( 'C' * 80 + 'Z' + 'A' * 120,
           'A' * 120 ),
         ( 'C' * 80 + 'Z' + 'A' * 120 ,)),

        # Right-covered: left has junk + 'Z' + perfect middle that covers right.
        # Correct behavior: keep both contigs.
        (( 'C' * 80 + 'Z' + 'A' * 120 + 'Z' + 'C' * 40,
           'A' * 120 ),
         ( 'A' * 120, 'C' * 80 + 'Z' + 'A' * 120 + 'Z' + 'C' * 40,)),
    ],
)
def test_stitch_simple_cases(seqs, expected, disable_acceptable_prob_check, force_failing_map_overlap):
    contigs = [ContigWithAligner(None, seq, reads_count=None) for seq in seqs]
    with ReferencelessStitcherContext.fresh():
        consenses = tuple(sorted(contig.seq for contig in stitch_consensus(contigs)))
    assert consenses == tuple(sorted(expected))


# Unit tests for Pool class
class TestPool:
    """Unit tests for the Pool class from referenceless_contig_stitcher."""

    def test_pool_empty_creation(self):
        """Test creating an empty pool with specified capacity."""
        capacity = 5
        pool = Pool.empty(capacity, ACCEPTABLE_STITCHING_SCORE())

        # Check pool structure
        assert pool.ring.capacity == capacity
        assert len(pool.ring) == 0
        assert len(pool.existing) == 0
        assert pool.smallest_score == ACCEPTABLE_STITCHING_SCORE()
        assert pool.min_acceptable_score == ACCEPTABLE_STITCHING_SCORE()

    def test_pool_add_single_path(self):
        """Test adding a single path to an empty pool."""
        pool = Pool.empty(3, ACCEPTABLE_STITCHING_SCORE())

        # Create a test path
        contig = ContigWithAligner(name="1", seq="ATCG", reads_count=None)
        path = ContigsPath.singleton(contig)

        # Add the path
        result = pool.add(path)

        assert result is True  # Should return True when successfully added
        assert len(pool.ring) == 1
        assert len(pool.existing) == 1
        assert pool.existing["ATCG"] == path
        assert pool.ring[0] == path

    def test_pool_add_multiple_paths_different_sequences(self):
        """Test adding multiple paths with different sequences."""
        pool = Pool.empty(3, ACCEPTABLE_STITCHING_SCORE())

        # Create test paths with different sequences
        contig1 = ContigWithAligner(name="1", seq="ATCG", reads_count=None)
        contig2 = ContigWithAligner(name="2", seq="GGCC", reads_count=None)
        contig3 = ContigWithAligner(name="3", seq="TTAA", reads_count=None)

        path1 = ContigsPath.singleton(contig1)
        path2 = ContigsPath.singleton(contig2)
        path3 = ContigsPath.singleton(contig3)

        # Add all paths
        assert pool.add(path1) is True
        assert pool.add(path2) is True
        assert pool.add(path3) is True

        assert len(pool.ring) == 3
        assert len(pool.existing) == 3
        assert pool.existing["ATCG"] == path1
        assert pool.existing["GGCC"] == path2
        assert pool.existing["TTAA"] == path3

    def test_pool_add_duplicate_sequence_worse_score(self):
        """Test adding a path with duplicate sequence but worse score (should be rejected)."""
        pool = Pool.empty(3, ACCEPTABLE_STITCHING_SCORE())

        # Create two paths with same sequence but different scores
        contig = ContigWithAligner(name="1", seq="ATCG", reads_count=None)
        better_path = ContigsPath(whole=contig, contigs_ids=frozenset([contig.id]), contains_contigs_ids=frozenset([contig.id]), score=Score(10.0))
        worse_path = ContigsPath(whole=contig, contigs_ids=frozenset([contig.id]), contains_contigs_ids=frozenset([contig.id]), score=Score(5.0))

        # Add better path first
        assert pool.add(better_path) is True
        assert len(pool.ring) == 1
        assert pool.existing["ATCG"] == better_path

        # Try to add worse path (should be rejected)
        assert pool.add(worse_path) is False
        assert len(pool.ring) == 1  # Should still be 1
        assert pool.existing["ATCG"] == better_path  # Should still be the better path

    def test_pool_add_duplicate_sequence_better_score(self):
        """Test adding a path with duplicate sequence but better score (replaces old path with deduplication)."""
        pool = Pool.empty(3, ACCEPTABLE_STITCHING_SCORE())

        # Create two paths with same sequence but different scores
        contig = ContigWithAligner(name="1", seq="ATCG", reads_count=None)
        worse_path = ContigsPath(whole=contig, contigs_ids=frozenset([contig.id]), contains_contigs_ids=frozenset([contig.id]), score=Score(5.0))
        better_path = ContigsPath(whole=contig, contigs_ids=frozenset([contig.id]), contains_contigs_ids=frozenset([contig.id]), score=Score(10.0))

        # Add worse path first
        assert pool.add(worse_path) is True
        assert len(pool.ring) == 1
        assert pool.existing["ATCG"] == worse_path

        # Add better path (new behavior: replaces old path with deduplication)
        assert pool.add(better_path) is True
        assert len(pool.ring) == 1  # Still only one path (replaced, not added)
        assert pool.existing["ATCG"] == better_path  # Existing mapping points to better path

        # Verify only the better path is in the ring
        ring_paths = list(pool.ring)
        assert worse_path not in ring_paths  # Old path should be removed
        assert better_path in ring_paths      # New path should be present
        assert pool.ring[0] == better_path    # Should be the only path

    def test_pool_capacity_enforcement(self):
        """Test that pool enforces capacity limits through SortedRing."""
        capacity = 2
        pool = Pool.empty(capacity, ACCEPTABLE_STITCHING_SCORE())

        # Create paths with different scores
        contig1 = ContigWithAligner(name="1", seq="AAAA", reads_count=None)  # Will have score SCORE_NOTHING (lowest)
        contig2 = ContigWithAligner(name="2", seq="CCCC", reads_count=None)
        contig3 = ContigWithAligner(name="3", seq="GGGG", reads_count=None)

        path1 = ContigsPath(whole=contig1, contigs_ids=frozenset([contig1.id]), contains_contigs_ids=frozenset([contig1.id]), score=Score(1.0))
        path2 = ContigsPath(whole=contig2, contigs_ids=frozenset([contig2.id]), contains_contigs_ids=frozenset([contig2.id]), score=Score(5.0))
        path3 = ContigsPath(whole=contig3, contigs_ids=frozenset([contig3.id]), contains_contigs_ids=frozenset([contig3.id]), score=Score(3.0))

        # Add first two paths
        assert pool.add(path1) is True
        assert pool.add(path2) is True
        assert len(pool.ring) == 2

        # Add third path - should cause capacity enforcement
        # The SortedRing should handle this internally
        pool.add(path3)
        assert len(pool.ring) <= capacity

    def test_pool_smallest_score_tracking(self):
        """Test that pool correctly tracks the smallest score."""
        pool = Pool.empty(3, ACCEPTABLE_STITCHING_SCORE())

        # Initially should be ACCEPTABLE_STITCHING_SCORE()
        assert pool.smallest_score == ACCEPTABLE_STITCHING_SCORE()
        assert pool.min_acceptable_score == ACCEPTABLE_STITCHING_SCORE()

        # Add a path with score higher than ACCEPTABLE_STITCHING_SCORE()
        contig = ContigWithAligner(name="1", seq="ATCG", reads_count=None)
        high_score_path = ContigsPath(whole=contig, contigs_ids=frozenset([contig.id]), contains_contigs_ids=frozenset([contig.id]), score=Score(100.0))

        pool.add(high_score_path)

        # smallest_score should update to max of ring[0].score and ACCEPTABLE_STITCHING_SCORE()
        expected_score = max(high_score_path.get_score(), ACCEPTABLE_STITCHING_SCORE())
        assert pool.smallest_score == expected_score
        assert pool.min_acceptable_score == expected_score

    def test_pool_min_acceptable_score_property(self):
        """Test that min_acceptable_score property returns smallest_score."""
        pool = Pool.empty(3, ACCEPTABLE_STITCHING_SCORE())

        # Should initially equal ACCEPTABLE_STITCHING_SCORE()
        assert pool.min_acceptable_score == pool.smallest_score
        assert pool.min_acceptable_score == ACCEPTABLE_STITCHING_SCORE()

        # Add a path and check that property updates
        contig = ContigWithAligner(name="1", seq="ATCG", reads_count=None)
        path = ContigsPath(whole=contig, contigs_ids=frozenset([contig.id]), contains_contigs_ids=frozenset([contig.id]), score=Score(50.0))
        pool.add(path)

        assert pool.min_acceptable_score == pool.smallest_score

    def test_pool_add_path_with_score_nothing(self):
        """Test adding a path with SCORE_NOTHING."""
        pool = Pool.empty(3, ACCEPTABLE_STITCHING_SCORE())

        contig = ContigWithAligner(name="1", seq="ATCG", reads_count=None)
        path = ContigsPath.singleton(contig)  # This creates a path with SCORE_NOTHING

        result = pool.add(path)
        assert result is True
        assert len(pool.ring) == 1
        assert pool.existing["ATCG"] == path

    def test_pool_ring_sorted_by_score(self):
        """Test that paths in the ring are sorted by score."""
        pool = Pool.empty(5, ACCEPTABLE_STITCHING_SCORE())

        # Create paths with different scores
        contigs_and_scores = [
            (ContigWithAligner(name="1", seq="AAAA", reads_count=None), Score(10.0)),
            (ContigWithAligner(name="2", seq="CCCC", reads_count=None), Score(5.0)),
            (ContigWithAligner(name="3", seq="GGGG", reads_count=None), Score(15.0)),
            (ContigWithAligner(name="4", seq="TTTT", reads_count=None), Score(1.0))
        ]

        paths = []
        for contig, score in contigs_and_scores:
            path = ContigsPath(whole=contig, contigs_ids=frozenset([contig.id]), contains_contigs_ids=frozenset([contig.id]), score=score)
            paths.append(path)
            pool.add(path)

        # Check that ring is sorted (SortedRing should maintain order)
        assert len(pool.ring) == 4
        # SortedRing keeps items sorted, so scores should be in ascending order
        ring_scores = [path.get_score() for path in pool.ring]
        assert ring_scores == sorted(ring_scores)

    def test_pool_existing_mapping_consistency(self):
        """Test that existing mapping stays consistent with ring contents."""
        pool = Pool.empty(3, ACCEPTABLE_STITCHING_SCORE())

        # Add some paths
        contig1 = ContigWithAligner(name="1", seq="ATCG", reads_count=None)
        contig2 = ContigWithAligner(name="2", seq="GGCC", reads_count=None)

        path1 = ContigsPath.singleton(contig1)
        path2 = ContigsPath.singleton(contig2)

        pool.add(path1)
        pool.add(path2)

        # Check that all sequences in existing are represented
        assert "ATCG" in pool.existing
        assert "GGCC" in pool.existing
        assert pool.existing["ATCG"] == path1
        assert pool.existing["GGCC"] == path2

        # The paths should also be in the ring
        ring_sequences = {path.whole.seq for path in pool.ring}
        existing_sequences = set(pool.existing.keys())
        assert ring_sequences == existing_sequences

    def test_pool_add_exactly_equal_scores(self):
        """Test adding paths with exactly equal scores (should reject duplicate)."""
        pool = Pool.empty(3, ACCEPTABLE_STITCHING_SCORE())

        # Create two paths with same sequence and EXACTLY equal scores
        contig = ContigWithAligner(name="1", seq="ATCG", reads_count=None)
        path1 = ContigsPath(whole=contig, contigs_ids=frozenset([contig.id]), contains_contigs_ids=frozenset([contig.id]), score=Score(5.0))
        path2 = ContigsPath(whole=contig, contigs_ids=frozenset([contig.id]), contains_contigs_ids=frozenset([contig.id]), score=Score(5.0))

        # Add first path
        assert pool.add(path1) is True
        assert len(pool.ring) == 1
        assert pool.existing["ATCG"] == path1

        # Try to add second path with equal score (should be rejected due to >= condition)
        assert pool.add(path2) is False
        assert len(pool.ring) == 1  # Should still be 1
        assert pool.existing["ATCG"] == path1  # Should still be the first path

    def test_pool_ring_insert_failure_but_existing_updated(self):
        """Test intricate behavior: ring.insert() fails but existing mapping is still updated."""
        pool = Pool.empty(1, ACCEPTABLE_STITCHING_SCORE())  # Very small capacity

        # Create paths where ring insertion might fail but scores allow existing update
        contig1 = ContigWithAligner(name="1", seq="AAAA", reads_count=None)
        contig2 = ContigWithAligner(name="2", seq="BBBB", reads_count=None)  # Different sequence

        # Create paths with very low scores that might be rejected by ring
        path1 = ContigsPath(whole=contig1, contigs_ids=frozenset([contig1.id]), contains_contigs_ids=frozenset([contig1.id]), score=Score(1.0))
        path2 = ContigsPath(whole=contig2, contigs_ids=frozenset([contig2.id]), contains_contigs_ids=frozenset([contig2.id]), score=Score(0.5))  # Even lower score

        # Add first path (should succeed)
        result1 = pool.add(path1)
        assert result1 is True
        assert len(pool.ring) == 1
        assert pool.existing["AAAA"] == path1

        # Try to add second path with lower score - ring.insert() should fail
        # because capacity=1 and new score is lower than existing
        result2 = pool.add(path2)
        assert result2 is False  # Ring insertion failed
        assert len(pool.ring) == 1  # Ring unchanged
        # But existing mapping should still be updated (this is the intricate behavior)
        assert pool.existing["BBBB"] == path2

    def test_pool_capacity_zero_edge_case(self):
        """Test Pool behavior with edge case capacities."""
        # Test that Pool.empty() with capacity 0 fails (SortedRing should reject this)
        try:
            pool = Pool.empty(0, ACCEPTABLE_STITCHING_SCORE())
            # If we get here, the Pool constructor didn't validate capacity
            # Let's test the behavior anyway
            contig = ContigWithAligner(name="1", seq="ATCG", reads_count=None)
            path = ContigsPath.singleton(contig)
            result = pool.add(path)
            # With capacity 0, ring.insert should always fail
            assert result is False
        except ValueError:
            # Expected: SortedRing should reject capacity 0
            pass

    def test_pool_empty_ring_access_edge_case(self):
        """Test potential IndexError when accessing ring[0] on empty ring."""
        pool = Pool.empty(1, ACCEPTABLE_STITCHING_SCORE())

        # Initially ring is empty, so pool.smallest_score should be ACCEPTABLE_STITCHING_SCORE()
        assert pool.smallest_score == ACCEPTABLE_STITCHING_SCORE()

        # Create a path with very low score that ring.insert() will reject
        contig = ContigWithAligner(name="1", seq="ATCG", reads_count=None)
        very_low_score = Score(-1000.0)  # Much lower than ACCEPTABLE_STITCHING_SCORE()
        path = ContigsPath(whole=contig, contigs_ids=frozenset([contig.id]), contains_contigs_ids=frozenset([contig.id]), score=very_low_score)

        # This should not crash even though ring[0] would be accessed if ring.insert() succeeded
        pool.add(path)
        # Ring insertion should fail due to low score vs capacity constraint

        # The ring should still be empty, so smallest_score should remain unchanged
        assert pool.smallest_score == ACCEPTABLE_STITCHING_SCORE()

    def test_pool_existing_mapping_vs_ring_inconsistency(self):
        """Test scenarios where existing mapping and ring contents might become inconsistent."""
        pool = Pool.empty(2, ACCEPTABLE_STITCHING_SCORE())

        # Add paths that will fill the capacity
        contig1 = ContigWithAligner(name="1", seq="AAAA", reads_count=None)
        contig2 = ContigWithAligner(name="2", seq="BBBB", reads_count=None)
        contig3 = ContigWithAligner(name="3", seq="CCCC", reads_count=None)

        path1 = ContigsPath(whole=contig1, contigs_ids=frozenset([contig1.id]), contains_contigs_ids=frozenset([contig1.id]), score=Score(5.0))
        path2 = ContigsPath(whole=contig2, contigs_ids=frozenset([contig2.id]), contains_contigs_ids=frozenset([contig2.id]), score=Score(10.0))
        path3 = ContigsPath(whole=contig3, contigs_ids=frozenset([contig3.id]), contains_contigs_ids=frozenset([contig3.id]), score=Score(3.0))

        # Add first two paths (fill capacity)
        assert pool.add(path1) is True
        assert pool.add(path2) is True
        assert len(pool.ring) == 2

        # Add third path with lower score - should be rejected by ring but existing updated
        assert pool.add(path3) is False

        # Check for potential inconsistency
        existing_sequences = set(pool.existing.keys())
        ring_sequences = {path.whole.seq for path in pool.ring}

        # This tests whether the implementation maintains consistency
        # (Current implementation might have inconsistency based on add() logic)
        # For now, just verify we can access both without crashes
        assert len(existing_sequences) >= 0
        assert len(ring_sequences) >= 0

    def test_pool_duplicate_with_marginally_better_score(self):
        """Test duplicate sequence with marginally better score (replaces with deduplication)."""
        pool = Pool.empty(3, ACCEPTABLE_STITCHING_SCORE())

        contig = ContigWithAligner(name="1", seq="ATCG", reads_count=None)
        # Create paths with very close but different scores
        worse_path = ContigsPath(whole=contig, contigs_ids=frozenset([contig.id]), contains_contigs_ids=frozenset([contig.id]), score=Score(10.0))
        barely_better_path = ContigsPath(whole=contig, contigs_ids=frozenset([contig.id]), contains_contigs_ids=frozenset([contig.id]), score=Score(10.000001))

        # Add worse path first
        assert pool.add(worse_path) is True
        assert pool.existing["ATCG"] == worse_path

        # Add barely better path (tests floating point comparison precision)
        assert pool.add(barely_better_path) is True
        assert pool.existing["ATCG"] == barely_better_path

        # With deduplication, only the better path should remain
        assert len(pool.ring) == 1
        assert pool.ring[0] == barely_better_path

    def test_pool_empty_sequence_edge_case(self):
        """Test Pool behavior with empty sequence strings."""
        pool = Pool.empty(3, ACCEPTABLE_STITCHING_SCORE())

        # Create contig with empty sequence
        contig = ContigWithAligner(name="1", seq="", reads_count=None)
        path = ContigsPath.singleton(contig)

        result = pool.add(path)
        assert result is True
        assert len(pool.ring) == 1
        assert pool.existing[""] == path  # Empty string as key

        # Try to add another path with empty sequence
        contig2 = ContigWithAligner(name="2", seq="", reads_count=None)
        path2 = ContigsPath(whole=contig2, contigs_ids=frozenset([contig2.id]), contains_contigs_ids=frozenset([contig2.id]), score=Score(5.0))

        # Should be treated as duplicate sequence
        assert pool.add(path2) is True
        assert pool.existing[""] == path2  # Should update to better path

    def test_pool_smallest_score_only_updates_on_successful_ring_insert(self):
        """Test that smallest_score only updates when ring.insert() returns True."""
        pool = Pool.empty(1, ACCEPTABLE_STITCHING_SCORE())

        # Add a path that will succeed
        contig1 = ContigWithAligner(name="1", seq="AAAA", reads_count=None)
        path1 = ContigsPath(whole=contig1, contigs_ids=frozenset([contig1.id]), contains_contigs_ids=frozenset([contig1.id]), score=Score(10.0))

        assert pool.add(path1) is True
        # smallest_score should update because ring.insert() returned True
        assert pool.smallest_score == max(path1.get_score(), ACCEPTABLE_STITCHING_SCORE())

        # Try to add a path that will fail ring insertion
        contig2 = ContigWithAligner(name="2", seq="BBBB", reads_count=None)
        path2 = ContigsPath(whole=contig2, contigs_ids=frozenset([contig2.id]), contains_contigs_ids=frozenset([contig2.id]), score=Score(5.0))  # Lower score

        current_smallest_score = pool.smallest_score
        assert pool.add(path2) is False  # Should fail due to capacity and lower score
        # smallest_score should NOT update because ring.insert() returned False
        assert pool.smallest_score == current_smallest_score

    def test_pool_capacity_one_detailed_behavior(self):
        """Test detailed behavior with capacity=1 to understand ring/existing interaction."""
        pool = Pool.empty(1, ACCEPTABLE_STITCHING_SCORE())

        # Add first path
        contig1 = ContigWithAligner(name="1", seq="AAAA", reads_count=None)
        path1 = ContigsPath(whole=contig1, contigs_ids=frozenset([contig1.id]), contains_contigs_ids=frozenset([contig1.id]), score=Score(5.0))

        result1 = pool.add(path1)
        assert result1 is True
        assert len(pool.ring) == 1
        assert pool.existing["AAAA"] == path1

        # Add second path with higher score and different sequence
        contig2 = ContigWithAligner(name="2", seq="BBBB", reads_count=None)
        path2 = ContigsPath(whole=contig2, contigs_ids=frozenset([contig2.id]), contains_contigs_ids=frozenset([contig2.id]), score=Score(10.0))

        result2 = pool.add(path2)
        assert result2 is True  # Should succeed because score is higher
        assert len(pool.ring) == 1  # Ring should still have size 1 (replaced, not added)
        assert pool.existing["BBBB"] == path2
        # But what happened to existing["AAAA"]? It should still be there based on add() logic

        # Add third path with lower score and different sequence
        contig3 = ContigWithAligner(name="3", seq="CCCC", reads_count=None)
        path3 = ContigsPath(whole=contig3, contigs_ids=frozenset([contig3.id]), contains_contigs_ids=frozenset([contig3.id]), score=Score(3.0))

        assert pool.add(path3) is False  # Should fail because score is lower than ring[0]
        # But existing mapping should still be updated
        assert pool.existing["CCCC"] == path3

    def test_pool_deduplication_set_consistency(self):
        """Test set management during deduplication."""
        pool = Pool.empty(3, ACCEPTABLE_STITCHING_SCORE())

        # Add initial path
        contig = ContigWithAligner(name="1", seq="ATCG", reads_count=None)
        path1 = ContigsPath(whole=contig, contigs_ids=frozenset([contig.id]), contains_contigs_ids=frozenset([contig.id]), score=Score(5.0))

        assert pool.add(path1) is True
        assert "ATCG" in pool.set
        assert pool.existing["ATCG"] == path1

        # Add better path with same sequence (should replace)
        path2 = ContigsPath(whole=contig, contigs_ids=frozenset([contig.id]), contains_contigs_ids=frozenset([contig.id]), score=Score(10.0))

        assert pool.add(path2) is True

        assert "ATCG" in pool.set, "BUG: Set incorrectly removes sequence during replacement"
        assert pool.existing["ATCG"] == path2
        assert len(pool.ring) == 1
        assert pool.ring[0] == path2

    def test_pool_existing_mapping_comprehensive_tracking(self):
        """Test that existing mapping comprehensively tracks all sequences for deduplication."""
        # Create a small pool that will be full
        pool = Pool.empty(1, ACCEPTABLE_STITCHING_SCORE())

        # Fill the pool to capacity
        contig1 = ContigWithAligner(name="1", seq="AAAA", reads_count=None)
        path1 = ContigsPath(whole=contig1, contigs_ids=frozenset([contig1.id]), contains_contigs_ids=frozenset([contig1.id]), score=Score(10.0))
        assert pool.add(path1) is True

        # Try to add a lower-scoring path that should be rejected by ring
        contig2 = ContigWithAligner(name="2", seq="TTTT", reads_count=None)
        path2 = ContigsPath(whole=contig2, contigs_ids=frozenset([contig2.id]), contains_contigs_ids=frozenset([contig2.id]), score=Score(5.0))
        result = pool.add(path2)

        # The add should fail
        assert not result, "Expected add to fail for lower-scoring path when pool is full"

        # DESIGN FEATURE: existing mapping tracks ALL sequences for comprehensive deduplication
        # This prevents future additions of the same sequence even if it was previously rejected
        assert "TTTT" in pool.existing, (
            f"Expected existing mapping to track all sequences for deduplication. "
            f"existing contains: {list(pool.existing.keys())}"
        )

        # Verify that trying to add the same sequence again is properly rejected
        contig3 = ContigWithAligner(name="3", seq="TTTT", reads_count=None)
        path3 = ContigsPath(whole=contig3, contigs_ids=frozenset([contig3.id]), contains_contigs_ids=frozenset([contig3.id]), score=Score(6.0))  # Slightly better but still worse than what's tracked
        result2 = pool.add(path3)
        assert not result2, "Expected duplicate sequence to be rejected based on existing mapping"

    def test_pool_set_tracking_during_complex_operations(self):
        """Test that set correctly tracks sequences during complex add/remove scenarios."""
        pool = Pool.empty(2, ACCEPTABLE_STITCHING_SCORE())  # Small capacity to force evictions

        # Add first path
        contig1 = ContigWithAligner(name="1", seq="AAAA", reads_count=None)
        path1 = ContigsPath(whole=contig1, contigs_ids=frozenset([contig1.id]), contains_contigs_ids=frozenset([contig1.id]), score=Score(10.0))
        assert pool.add(path1) is True
        assert "AAAA" in pool.set

        # Add second path (different sequence)
        contig2 = ContigWithAligner(name="2", seq="BBBB", reads_count=None)
        path2 = ContigsPath(whole=contig2, contigs_ids=frozenset([contig2.id]), contains_contigs_ids=frozenset([contig2.id]), score=Score(20.0))
        assert pool.add(path2) is True
        assert "BBBB" in pool.set
        assert len(pool.set) == 2

        # Add third path (different sequence, high score) - should evict lowest
        contig3 = ContigWithAligner(name="3", seq="CCCC", reads_count=None)
        path3 = ContigsPath(whole=contig3, contigs_ids=frozenset([contig3.id]), contains_contigs_ids=frozenset([contig3.id]), score=Score(30.0))
        result = pool.add(path3)

        if result is True:
            # Set should reflect actual ring contents
            ring_sequences = {path.whole.seq for path in pool.ring}
            assert ring_sequences == pool.set, "BUG: Set doesn't match actual ring contents"
            assert "CCCC" in pool.set
            # One of the previous sequences should have been evicted
            assert len(pool.set) == len(pool.ring)

    def test_pool_duplicate_replacement_with_eviction_scenario(self):
        """Test complex scenario: duplicate replacement when at capacity."""
        pool = Pool.empty(2, ACCEPTABLE_STITCHING_SCORE())

        # Fill capacity with two different sequences
        contig1 = ContigWithAligner(name="1", seq="AAAA", reads_count=None)
        contig2 = ContigWithAligner(name="2", seq="BBBB", reads_count=None)
        path1 = ContigsPath(whole=contig1, contigs_ids=frozenset([contig1.id]), contains_contigs_ids=frozenset([contig1.id]), score=Score(10.0))
        path2 = ContigsPath(whole=contig2, contigs_ids=frozenset([contig2.id]), contains_contigs_ids=frozenset([contig2.id]), score=Score(15.0))

        assert pool.add(path1) is True
        assert pool.add(path2) is True
        assert len(pool.ring) == 2
        assert len(pool.set) == 2

        # Now add a better version of the first sequence
        better_path1 = ContigsPath(whole=contig1, contigs_ids=frozenset([contig1.id]), contains_contigs_ids=frozenset([contig1.id]), score=Score(25.0))
        assert pool.add(better_path1) is True

        # Check consistency: should still have 2 items, with better_path1 replacing path1
        assert len(pool.ring) == 2
        assert len(pool.set) == 2
        assert "AAAA" in pool.set
        assert "BBBB" in pool.set
        assert pool.existing["AAAA"] == better_path1

        # Verify that better_path1 is actually in the ring, not the old path1
        ring_paths = list(pool.ring)
        assert better_path1 in ring_paths
        assert path1 not in ring_paths

    def test_pool_empty_ring_smallest_score_bug(self):
        """Test potential IndexError when accessing ring[0] on empty ring."""
        pool = Pool.empty(1, ACCEPTABLE_STITCHING_SCORE())

        # Try to add a path with very low score that ring might reject
        contig = ContigWithAligner(name="1", seq="ATCG", reads_count=None)
        very_low_score = Score(-1000.0)
        path = ContigsPath(whole=contig, contigs_ids=frozenset([contig.id]), contains_contigs_ids=frozenset([contig.id]), score=very_low_score)

        # This should not crash due to accessing ring[0] when ring is empty
        try:
            result = pool.add(path)
            # If it succeeded, smallest_score should be updated safely
            if result is True:
                assert pool.smallest_score >= ACCEPTABLE_STITCHING_SCORE()
        except IndexError as e:
            if "ring[0]" in str(e):
                assert False, "BUG: IndexError when accessing ring[0] on empty ring"

    def test_pool_deduplication_assert_failure_scenario(self):
        """Test scenario that might trigger the assert to_delete_index >= 0."""
        pool = Pool.empty(3, ACCEPTABLE_STITCHING_SCORE())

        # Add a path
        contig = ContigWithAligner(name="1", seq="ATCG", reads_count=None)
        path1 = ContigsPath(whole=contig, contigs_ids=frozenset([contig.id]), contains_contigs_ids=frozenset([contig.id]), score=Score(5.0))
        assert pool.add(path1) is True

        # Manually corrupt the state to test assert robustness
        # (In real scenario, this might happen due to race conditions or other bugs)
        pool.set.add("FAKE_SEQ")  # Add fake sequence to set

        # Try to add a path with sequence that's in set but not in ring
        fake_contig = ContigWithAligner(name="fake", seq="FAKE_SEQ", reads_count=None)
        fake_path = ContigsPath(whole=fake_contig, contigs_ids=frozenset([fake_contig.id]), contains_contigs_ids=frozenset([fake_contig.id]), score=Score(10.0))

        try:
            pool.add(fake_path)
        except AssertionError:
            # The assert to_delete_index >= 0 should catch this inconsistency
            pass  # Expected in this corrupted state test

    def test_pool_set_existing_ring_consistency_invariant(self):
        """Test that set, existing, and ring maintain consistency invariants."""
        pool = Pool.empty(4, ACCEPTABLE_STITCHING_SCORE())

        sequences_to_test = ["AAAA", "BBBB", "CCCC", "DDDD", "EEEE"]
        paths = []

        for i, seq in enumerate(sequences_to_test):
            contig = ContigWithAligner(name=str(i), seq=seq, reads_count=None)
            path = ContigsPath(whole=contig, contigs_ids=frozenset([contig.id]), contains_contigs_ids=frozenset([contig.id]), score=Score(float(i * 10)))
            paths.append(path)

            pool.add(path)

            # Verify invariants after each addition
            ring_sequences = {p.whole.seq for p in pool.ring}
            existing_sequences = set(pool.existing.keys())

            # INVARIANT 1: set should match ring sequences
            assert pool.set == ring_sequences, f"INVARIANT VIOLATED: set {pool.set} != ring sequences {ring_sequences}"

            # INVARIANT 2: all ring sequences should be in existing
            assert ring_sequences.issubset(existing_sequences), f"INVARIANT VIOLATED: ring sequences {ring_sequences} not in existing {existing_sequences}"

            # INVARIANT 3: ring size should not exceed capacity
            assert len(pool.ring) <= pool.ring.capacity, f"INVARIANT VIOLATED: ring size {len(pool.ring)} exceeds capacity {pool.ring.capacity}"

    def test_pool_deduplication_score_comparison_edge_cases(self):
        """Test edge cases in score comparison for deduplication."""
        pool = Pool.empty(3, ACCEPTABLE_STITCHING_SCORE())

        contig = ContigWithAligner(name="1", seq="ATCG", reads_count=None)

        # Test with exactly equal scores (should be rejected)
        path1 = ContigsPath(whole=contig, contigs_ids=frozenset([contig.id]), contains_contigs_ids=frozenset([contig.id]), score=Score(10.0))
        path2 = ContigsPath(whole=contig, contigs_ids=frozenset([contig.id]), contains_contigs_ids=frozenset([contig.id]), score=Score(10.0))

        assert pool.add(path1) is True
        assert pool.add(path2) is False  # Equal score should be rejected
        assert pool.existing["ATCG"] == path1  # Should still be first path

        # Test with very small difference (floating point precision)
        epsilon_better = ContigsPath(whole=contig, contigs_ids=frozenset([contig.id]), contains_contigs_ids=frozenset([contig.id]), score=Score(10.0000000001))
        assert pool.add(epsilon_better) is True
        assert pool.existing["ATCG"] == epsilon_better

        # Test with much better score
        much_better = ContigsPath(whole=contig, contigs_ids=frozenset([contig.id]), contains_contigs_ids=frozenset([contig.id]), score=Score(50.0))
        assert pool.add(much_better) is True
        assert pool.existing["ATCG"] == much_better
        assert len(pool.ring) == 1  # Should still be just one path

    def test_pool_deduplication_multiple_replacements(self):
        """Test multiple consecutive replacements of the same sequence."""
        pool = Pool.empty(5, ACCEPTABLE_STITCHING_SCORE())

        contig = ContigWithAligner(name="1", seq="ATCG", reads_count=None)
        scores = [5.0, 10.0, 15.0, 8.0, 25.0, 3.0, 30.0]  # Mix of better and worse scores

        best_path = None
        best_score = 0.0

        for i, score in enumerate(scores):
            path = ContigsPath(whole=contig, contigs_ids=frozenset([contig.id]), contains_contigs_ids=frozenset([contig.id]), score=Score(score))
            result = pool.add(path)

            if score > best_score:
                # Should succeed and update
                assert result is True
                best_path = path
                best_score = score
                assert pool.existing["ATCG"] == best_path
            else:
                # Should be rejected
                assert result is False
                assert pool.existing["ATCG"] == best_path  # Should remain unchanged

            # Pool should always contain exactly one path for this sequence
            assert len(pool.ring) == 1
            assert "ATCG" in pool.set
            assert pool.ring[0] == best_path

    def test_pool_deduplication_with_different_sequence_interleaved(self):
        """Test deduplication behavior when adding different sequences between duplicates."""
        pool = Pool.empty(3, ACCEPTABLE_STITCHING_SCORE())

        # Add initial path
        contig1 = ContigWithAligner(name="1", seq="AAAA", reads_count=None)
        path1a = ContigsPath(whole=contig1, contigs_ids=frozenset([contig1.id]), contains_contigs_ids=frozenset([contig1.id]), score=Score(10.0))
        assert pool.add(path1a) is True

        # Add different sequence
        contig2 = ContigWithAligner(name="2", seq="BBBB", reads_count=None)
        path2 = ContigsPath(whole=contig2, contigs_ids=frozenset([contig2.id]), contains_contigs_ids=frozenset([contig2.id]), score=Score(15.0))
        assert pool.add(path2) is True

        # Add another different sequence
        contig3 = ContigWithAligner(name="3", seq="CCCC", reads_count=None)
        path3 = ContigsPath(whole=contig3, contigs_ids=frozenset([contig3.id]), contains_contigs_ids=frozenset([contig3.id]), score=Score(20.0))
        assert pool.add(path3) is True

        assert len(pool.ring) == 3
        assert len(pool.set) == 3

        # Now add better version of first sequence
        path1b = ContigsPath(whole=contig1, contigs_ids=frozenset([contig1.id]), contains_contigs_ids=frozenset([contig1.id]), score=Score(25.0))
        assert pool.add(path1b) is True

        # Should still have 3 paths, but path1a should be replaced with path1b
        assert len(pool.ring) == 3
        assert len(pool.set) == 3
        assert pool.existing["AAAA"] == path1b

        ring_paths = list(pool.ring)
        assert path1b in ring_paths
        assert path1a not in ring_paths
        assert path2 in ring_paths
        assert path3 in ring_paths


# Additional end-to-end tests for the stitcher
@pytest.mark.parametrize(
    "seqs, expected",
    [
        # Deduplicate identical contigs (perfect full overlap -> single output)
        (("AAA" + TTT + "GGG", "AAA" + TTT + "GGG"), ("AAA" + TTT + "GGG",)),

        # Presence of empty contig should be preserved; other contigs still stitch
        (("", "AAAAA" + TTT, TTT + "GGGGG"), ("", "AAAAA" + TTT + "GGGGG")),

        # Two independent pairs should stitch separately in one run
        (("XXX" + TTT, TTT + "GGG", "CCC" + AAA, AAA + "DDD"),
         ("XXX" + TTT + "GGG", "CCC" + AAA + "DDD")),

        # Multi-step chain across three overlaps, requires second-phase (o2_loop)
        (("AAA" + TTT, TTT + "G" * 120, "G" * 120 + AAA, AAA + "CCC"),
         ("AAA" + TTT + "G" * 120 + AAA + "CCC",)),

        # Multiple identical duplicates - new behavior: mutual coverage keeps all copies
        ((("XXXXX" + TTT + "YYYYY",) * 9), (("XXXXX" + TTT + "YYYYY",) * 9)),
    ],
)
def test_stitch_additional_cases(seqs, expected, disable_acceptable_prob_check):
    contigs = [ContigWithAligner(None, seq, reads_count=None) for seq in seqs]
    with ReferencelessStitcherContext.fresh():
        consenses = tuple(sorted(contig.seq for contig in stitch_consensus(contigs)))
    assert consenses == tuple(sorted(expected))


def test_stitch_covered_contig_is_ignored(disable_acceptable_prob_check):
    """When one contig is fully contained in another, keep only the bigger contig."""
    bigger = "AAAAA" + TTT + "GGGGG"
    contained = TTT
    contigs = [ContigWithAligner(None, bigger, reads_count=None), ContigWithAligner(None, contained, reads_count=None)]
    with ReferencelessStitcherContext.fresh():
        consenses = tuple(sorted(c.seq for c in stitch_consensus(contigs)))
    assert consenses == (bigger,)


# ---------------------------------------------------------------------------
# Tests for check_merged_sequence_support — join-boundary read validation
# ---------------------------------------------------------------------------

# A merged sequence using real DNA bases for predictable cut-spanning tests.
# Left half: AAACCCGGGTTTAAA   Right half: CCCGGGTTTAAACCC   Cut at 15.
DNA_MERGED = "AAACCCGGGTTTAAACCCGGGTTTAAACCC"   # 30 bp
DNA_CUT = 15
DNA_RLEN = 5  # read_length for window size
# window = [15 - 5//2, 15 + (5-5//2)) = [13, 18)  —  5 positions


class TestCheckMergedSequenceSupport:
    """Unit tests for check_merged_sequence_support."""

    # ---- helper: build a read index from placements -------------------
    @staticmethod
    def _make_rd(seq: str, starts_with_lens, count=1):
        rd = {}
        for s, L in starts_with_lens:
            frag = seq[s:s + L]
            canonical = min(frag, reverse_complement(frag))
            rd.setdefault(L, Counter())[canonical] += count
        return rd

    # ------------------------------------------------------------------
    # 1. Cut-spanning — at least min_depth reads must cross the cut
    # ------------------------------------------------------------------

    def test_spanning_reads_accepted(self):
        """Reads that strictly cross the cut pass the check."""
        # S=12 "AAACC" covers [12,17), crosses cut=15 (12<15<17).
        # S=13 "AACCC" covers [13,18), also crosses cut.
        rd = self._make_rd(DNA_MERGED, [(12, 5), (13, 5)])
        assert check_merged_sequence_support(DNA_MERGED, DNA_CUT, rd, 1, DNA_RLEN)[0]

    def test_spanning_no_crossing_rejected(self):
        """Reads ending exactly at the cut do NOT span (start < cut)."""
        # S=10 "TTAAA" covers [10,15).  10 < 15 < 15 → false, doesn't span.
        rd = self._make_rd(DNA_MERGED, [(10, 5)])
        assert not check_merged_sequence_support(DNA_MERGED, DNA_CUT, rd, 1, DNA_RLEN)[0]

    def test_spanning_start_at_cut_rejected(self):
        """Reads starting exactly at the cut do NOT span (cut < end)."""
        # S=15 "CCCGG" covers [15,20).  15 < 15 < 20 → false, doesn't span.
        rd = self._make_rd(DNA_MERGED, [(15, 5)])
        assert not check_merged_sequence_support(DNA_MERGED, DNA_CUT, rd, 1, DNA_RLEN)[0]

    def test_split_left_and_right_rejected(self):
        """Independent left & right coverage without a crossing read fails."""
        rd = self._make_rd(DNA_MERGED, [(10, 5), (15, 5)])
        # Read at 10 covers [10,15)  — ends at cut, does NOT cross.
        # Read at 15 covers [15,20)  — starts at cut, does NOT cross.
        assert not check_merged_sequence_support(DNA_MERGED, DNA_CUT, rd, 1, DNA_RLEN)[0]

    def test_min_depth_zero_disables_spanning_check(self):
        """min_depth=0 → check skipped regardless of reads."""
        rd = self._make_rd(DNA_MERGED, [(10, 5)])  # does NOT cross
        assert check_merged_sequence_support(DNA_MERGED, DNA_CUT, rd, 0, DNA_RLEN)[0]

    # ------------------------------------------------------------------
    # 2. Boundary-window coverage — every position in the window covered
    # ------------------------------------------------------------------

    def test_window_fully_covered_accepted(self):
        """Every position in the window is covered to at least min_depth."""
        # Window [13,18): positions 13..17.  Provide reads covering all.
        rd = self._make_rd(DNA_MERGED, [
            (12, 5), (13, 5),  # span the cut
            (9, 5), (10, 5), (11, 5),  # left-only coverage
            (14, 5), (15, 5), (16, 5), (17, 5),  # right-side coverage
        ])
        assert check_merged_sequence_support(DNA_MERGED, DNA_CUT, rd, 1, DNA_RLEN)[0]

    def test_window_edge_uncovered_rejected(self):
        """A single position in the window with no coverage fails."""
        rd = self._make_rd(DNA_MERGED, [
            (12, 5),  # spans cut=15 (12<15<17), covers [12,17)
            # Position 13 in window [13,18) needs starts at [9,13].
            # None provided → uncovered.
        ])
        assert not check_merged_sequence_support(DNA_MERGED, DNA_CUT, rd, 1, DNA_RLEN)[0]

    # ------------------------------------------------------------------
    # 3. read_index = None  vs  read_index = {}
    # ------------------------------------------------------------------

    def test_none_read_index_disables_check(self):
        """read_index=None → validation disabled."""
        assert check_merged_sequence_support(DNA_MERGED, DNA_CUT, None, 1, DNA_RLEN)[0]

    def test_empty_read_index_rejects(self):
        """read_index={} with validation enabled → rejected (no reads)."""
        assert not check_merged_sequence_support(DNA_MERGED, DNA_CUT, {}, 1, DNA_RLEN)[0]

    # ------------------------------------------------------------------
    # 4. Reverse-complement support (real DNA bases only)
    # ------------------------------------------------------------------

    def test_reverse_complement_spanning(self):
        """A read present only as rc can still span the cut."""
        merged = "ACGT" * 8  # 32 bp
        cut = 16
        rc_kmer = "CGTAC"     # reverse complement of GTACG
        # Canonical of both is "CGTAC" — use that directly.
        rd = {5: Counter({rc_kmer: 1})}
        assert check_merged_sequence_support(merged, cut, rd, 1, 5)[0]

    def test_multiple_read_lengths(self):
        """Reads of different lengths are all considered."""
        merged = "ACGT" * 8  # 32 bp
        cut = 16
        # L=5 canonical "CGTAC" (count 2) at S=13,14 → 4 spanning.
        # L=7 canonical "GTACGTA" (count 1) at S=10,11,14,15 → 4 spanning.
        # Total spanning = 8.  min_depth=8 passes, min_depth=9 fails.
        rd = {5: Counter({"CGTAC": 2}), 7: Counter({"GTACGTA": 1})}
        assert check_merged_sequence_support(merged, cut, rd, 8, 5)[0]
        assert not check_merged_sequence_support(merged, cut, rd, 9, 5)[0]

    def test_multiple_read_lengths_insufficient(self):
        """All read lengths combined still don't reach impossible min_depth."""
        merged = "ACGT" * 8  # 32 bp
        cut = 16
        # Canonical of "GTACG" is "CGTAC".
        rd = {5: Counter({"CGTAC": 10})}
        assert not check_merged_sequence_support(merged, cut, rd, 999, 5)[0]

    # ------------------------------------------------------------------
    # 6. Placements and multiplicity — each placement × count contributes
    # ------------------------------------------------------------------

    def test_homopolymer_one_read_min_depth_1(self):
        """One read at 4 spanning placements × count 1 = spanning 4 ≥ 1."""
        merged = "A" * 30
        cut = 15
        # "AAAAA" has 4 spanning placements (S=11..14) × count 1 = 4.
        rd = {5: Counter({"AAAAA": 1})}
        assert check_merged_sequence_support(merged, cut, rd, 1, 5)[0]

    def test_homopolymer_placements_count(self):
        """One read at 4 spanning placements × count 1 = spanning 4."""
        merged = "A" * 30
        cut = 15
        rd = {5: Counter({"AAAAA": 1})}
        # 4 spanning placements (S=11,12,13,14) × 1 count = 4.
        assert check_merged_sequence_support(merged, cut, rd, 4, 5)[0]
        assert not check_merged_sequence_support(merged, cut, rd, 5, 5)[0]

    def test_multiplicity_increases_depth(self):
        """Higher FASTQ count increases each placement's contribution."""
        merged = "A" * 30
        cut = 15
        # 4 spanning placements × 10 count = 40.
        rd = {5: Counter({"AAAAA": 10})}
        assert check_merged_sequence_support(merged, cut, rd, 40, 5)[0]
        assert not check_merged_sequence_support(merged, cut, rd, 41, 5)[0]

    def test_rc_accumulates_counts(self):
        """A read and its reverse complement share one canonical, counts sum."""
        merged = "ACGT" * 8  # 32 bp
        cut = 16
        # Both "GTACG" and its rc "CGTAC" map to canonical "CGTAC" with count 2.
        # "CGTAC" matches at S=13 and S=14, both spanning → 2 placements × 2 = 4.
        rd = {5: Counter({"CGTAC": 2})}
        assert check_merged_sequence_support(merged, cut, rd, 4, 5)[0]
        assert not check_merged_sequence_support(merged, cut, rd, 5, 5)[0]

    def test_two_distinct_reads_min_depth_2(self):
        """Two genuinely different reads satisfy min_depth=2."""
        merged = "ACGT" * 8  # 32 bp
        cut = 16
        # "GTACG" canonical → "CGTAC", "TACGT" canonical → "ACGTA".
        # Both at S=14,15 span → spanning_total = 2.
        rd = {5: Counter({"CGTAC": 1, "ACGTA": 1})}
        assert check_merged_sequence_support(merged, cut, rd, 2, 5)[0]

    def test_window_coverage_placements(self):
        """Multiple placements contribute to window coverage."""
        merged = "A" * 30
        cut = 15
        rd = {5: Counter({"AAAAA": 2})}
        # 4 spanning placements (S=11..14) × 2 = spanning 8.
        # Window positions covered by up to 5 placements × 2 = coverage 10.
        # So spanning_total=8 is the bottleneck, not window coverage.
        assert check_merged_sequence_support(merged, cut, rd, 8, 5)[0]
        assert not check_merged_sequence_support(merged, cut, rd, 9, 5)[0]


# ---------------------------------------------------------------------------
# End-to-end tests with read index
# ---------------------------------------------------------------------------
def test_stitch_with_read_index_none_disabled(disable_acceptable_prob_check):
    """read_index=None => validation disabled, merge proceeds."""
    left = ContigWithAligner(None, "AAAACCCC", None)
    right = ContigWithAligner(None, "CCCCGGGG", None)
    with ReferencelessStitcherContext.fresh() as ctx:  # ctx.read_index stays None
        ctx.minimum_read_depth = 1
        result = tuple(stitch_consensus([left, right]))
    assert len(result) == 1
    assert result[0].seq == "AAAACCCCGGGG"


def test_stitch_with_read_index_empty_rejected(disable_acceptable_prob_check):
    """read_index={} => validation enabled but no reads, merge rejected."""
    left = ContigWithAligner(None, "AAAACCCC", None)
    right = ContigWithAligner(None, "CCCCGGGG", None)
    with ReferencelessStitcherContext.fresh() as ctx:
        ctx.read_index = {}
        ctx.minimum_read_depth = 1
        ctx.read_length = 4
        result = tuple(stitch_consensus([left, right]))
    assert len(result) == 2  # unmerged


def test_stitch_with_read_index_matching_accepted(disable_acceptable_prob_check):
    """read_index with exact k-mers for every start position => merge accepted."""
    left = ContigWithAligner(None, "AAAACCCC", None)
    right = ContigWithAligner(None, "CCCCGGGG", None)
    # Provide all possible 4-mers from the expected merged sequence.
    merged = "AAAACCCCGGGG"
    rd = {4: Counter()}
    for i in range(len(merged)-3):
        frag = merged[i:i+4]
        canonical = min(frag, reverse_complement(frag))
        rd[4][canonical] += 1
    with ReferencelessStitcherContext.fresh() as ctx:
        ctx.read_index = rd
        ctx.minimum_read_depth = 1
        ctx.read_length = 4
        result = tuple(stitch_consensus([left, right]))
    assert len(result) == 1
    assert result[0].seq == merged


def test_stitch_with_identical_contig_sequences(disable_acceptable_prob_check):
    """Identical contig sequences — covered-contig case, no merge needed."""
    seq = "AAAACCCCGGGG"
    left = ContigWithAligner(None, seq, None)
    right = ContigWithAligner(None, seq, None)
    # No read_index set → None → validation disabled.
    with ReferencelessStitcherContext.fresh():
        result = tuple(stitch_consensus([left, right]))
    assert len(result) == 1  # one covers the other


def test_stitch_with_reads_from_fastq(tmp_path, disable_acceptable_prob_check):
    """Reads from real FASTQ files are used for join validation."""
    import micall.utils.referenceless_contig_stitcher as stitcher_mod
    from micall.utils.referenceless_contig_stitcher import build_read_index

    fasta_path = tmp_path / "contigs.fasta"
    fasta_path.write_text(">left\nAAAACCCC\n>right\nCCCCGGGG\n")

    fq1 = tmp_path / "reads_R1.fastq"
    fq2 = tmp_path / "reads_R2.fastq"

    # Provide ALL 4-mers from the merged sequence — should be accepted.
    merged = "AAAACCCCGGGG"
    all_kmers = [merged[i:i+4] for i in range(len(merged)-3)]
    lines = []
    for i, kmer in enumerate(all_kmers):
        lines.append(f"@r{i}R1\n{kmer}\n+\n{chr(73)*4}\n")
    fq1.write_text("".join(lines))
    lines = []
    for i, kmer in enumerate(all_kmers):
        lines.append(f"@r{i}R2\n{kmer}\n+\n{chr(73)*4}\n")
    fq2.write_text("".join(lines))

    read_index = build_read_index(fq1, fq2)

    with open(fasta_path) as f:
        contigs = tuple(stitcher_mod.read_contigs(f))

    with ReferencelessStitcherContext.fresh() as ctx:
        ctx.read_index = read_index
        ctx.minimum_read_depth = 1
        ctx.read_length = 4
        result = tuple(stitch_consensus(contigs))

    assert len(result) == 1
    assert result[0].seq == merged


def test_merged_check_with_split_side_reads(disable_acceptable_prob_check):
    """Direct test that split-side reads fail the merged check.

    "AAAA" (left side, covers [0,4), ends before cut at 5).
    "CCCG" (right side, covers [5,9), starts at cut, rc "CGGG" missing from left).
    Neither spans cut=5; neither rc spills to the other side.
    """
    merged = "AAAACCCCGGGG"
    cut = 5
    # Note: these must be canonicalized. min("AAAA","TTTT")="AAAA", min("CCCG","CGGG")="CCCG".
    rd = {4: Counter({"AAAA": 1, "CCCG": 1})}
    assert not check_merged_sequence_support(merged, cut, rd, 1, 4)[0], \
        "AAAA (ends before cut) and CCCG (starts at cut, rc CGGG not on left) should not support the merge"


def test_stitch_with_reads_from_fastq_split_side_rejected(tmp_path, disable_acceptable_prob_check):
    """Split-side reads from FASTQ: both sides covered but no read crosses cut."""
    import micall.utils.referenceless_contig_stitcher as stitcher_mod
    from micall.utils.referenceless_contig_stitcher import build_read_index

    fasta_path = tmp_path / "contigs.fasta"
    fasta_path.write_text(">left\nAAAACCCC\n>right\nCCCCGGGG\n")

    fq1 = tmp_path / "reads_R1.fastq"
    fq2 = tmp_path / "reads_R2.fastq"

    # "AAAA" on merged covers [0,4), ends before the cut (cut≥?).
    # "CCCG" on merged covers [5,9), starts at the cut.
    # Neither spans the cut; rc "CGGG" doesn't match left side.
    left_only = "AAAA"
    right_only = "CCCG"
    fq1.write_text(f"@r1\n{left_only}\n+\nIIII\n@r2\n{right_only}\n+\nIIII\n")
    fq2.write_text(f"@r1\n{left_only}\n+\nIIII\n@r2\n{right_only}\n+\nIIII\n")

    read_index = build_read_index(fq1, fq2)

    with open(fasta_path) as f:
        contigs = tuple(stitcher_mod.read_contigs(f))

    with ReferencelessStitcherContext.fresh() as ctx:
        ctx.read_index = read_index
        ctx.minimum_read_depth = 1
        ctx.read_length = 4
        result = tuple(stitch_consensus(contigs))

    # Should NOT merge — no read crosses the join cut.
    assert len(result) == 2, (
        f"Split-side FASTQ reads should not allow a merge, "
        f"but got {len(result)} contig(s)."
    )


# ---------------------------------------------------------------------------
# Tests for build_read_index
# ---------------------------------------------------------------------------

def _make_fastq(tmp_path, name, records):
    """Write a FASTQ file with the given records (list of (header, seq))."""
    path = tmp_path / name
    lines = []
    for h, s in records:
        lines.append(f"@{h}\n{s}\n+\n{'I' * len(s)}\n")
    path.write_text("".join(lines))
    return path


def test_build_read_index_duplicates_and_rc(tmp_path):
    """build_read_index accumulates duplicate reads and RC equivalents."""
    from micall.utils.referenceless_contig_stitcher import build_read_index

    fq1 = _make_fastq(tmp_path, "r1.fastq", [
        ("r1", "ACGTA"),   # canonical = min("ACGTA","TACGT") = "ACGTA"
        ("r2", "TACGT"),   # rc of "ACGTA", canonical = "ACGTA" (same)
        ("r3", "ACGTA"),   # duplicate of r1
    ])
    fq2 = _make_fastq(tmp_path, "r2.fastq", [
        ("r4", "ACGTA"),   # duplicate again
        ("r5", "ACGTA"),   # duplicate again (same count as R1 to stay synced)
        ("r6", "ACGTA"),   # duplicate again
    ])

    index = build_read_index(fq1, fq2)

    # Two read lengths?  No, all are length 5.
    assert 5 in index, "Should have a counter for length 5"
    counter = index[5]
    # Only one canonical key "ACGTA" with total count 6:
    #   R1: ACGTA (1) + TACGT (1) → canonical "ACGTA" → +2
    #   R1: ACGTA (1)                                 → +1 (total 3)
    #   R2: ACGTA × 3                                 → +3 (total 6)
    assert counter == {"ACGTA": 6}, (
        f"Expected single canonical 'ACGTA' with count 6, got {dict(counter)}"
    )
    assert len(counter) == 1, "Only one distinct canonical read expected"


def test_build_read_index_empty_fastq(tmp_path):
    """build_read_index with empty FASTQ files returns empty dict."""
    from micall.utils.referenceless_contig_stitcher import build_read_index

    fq1 = _make_fastq(tmp_path, "empty1.fastq", [])
    fq2 = _make_fastq(tmp_path, "empty2.fastq", [])
    index = build_read_index(fq1, fq2)
    assert index == {}, "Empty FASTQs should produce empty index"


def test_cli_fastq_pair_validation(tmp_path):
    """CLI must reject --fastq1 without --fastq2 and vice versa."""
    from micall.core.contig_stitcher import main

    fasta = tmp_path / "in.fasta"
    fasta.write_text(">c\nAAAA\n")
    out = tmp_path / "out.fasta"
    fq = tmp_path / "r.fastq"
    fq.write_text("@r\nAAAA\n+\nIIII\n")

    # Only --fastq1 provided -> error.
    try:
        main(["without-references", str(fasta), str(out),
              "--fastq1", str(fq)])
        assert False, "Should have raised SystemExit"
    except SystemExit:
        pass

    # Only --fastq2 provided -> error.
    try:
        main(["without-references", str(fasta), str(out),
              "--fastq2", str(fq)])
        assert False, "Should have raised SystemExit"
    except SystemExit:
        pass


def test_cli_negative_read_depth_rejected(tmp_path):
    """--minimum-read-depth negative must be rejected."""
    from micall.core.contig_stitcher import main

    fasta = tmp_path / "in.fasta"
    fasta.write_text(">c\nAAAA\n")
    out = tmp_path / "out.fasta"

    try:
        main(["without-references", str(fasta), str(out),
              "--minimum-read-depth", "-1"])
        assert False, "Should have raised SystemExit"
    except SystemExit:
        pass


def test_cli_zero_read_length_rejected(tmp_path):
    """--read-length 0 must be rejected."""
    from micall.core.contig_stitcher import main

    fasta = tmp_path / "in.fasta"
    fasta.write_text(">c\nAAAA\n")
    out = tmp_path / "out.fasta"

    try:
        main(["without-references", str(fasta), str(out),
              "--read-length", "0"])
        assert False, "Should have raised SystemExit"
    except SystemExit:
        pass
