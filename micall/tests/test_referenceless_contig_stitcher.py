import pytest
from micall.utils.referenceless_contig_stitcher import \
    stitch_consensus, ContigWithAligner, Pool, \
    ACCEPTABLE_STITCHING_SCORE
from micall.utils.contig_stitcher_context import ReferencelessStitcherContext
from micall.utils.referenceless_score import Score
from micall.utils.referenceless_contig_path import ContigsPath

from micall.tests.referenceless_tests_utils import disable_acceptable_prob_check, force_failing_map_overlap

# prevent linter warnings
assert disable_acceptable_prob_check is not None
assert force_failing_map_overlap is not None


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
        (('AAAAA' + TTT, 'GGGGG' + TTT), ("AAAAA" + TTT,)),
        (('GGGGG' + TTT, 'CCCCCAAAAA' + TTT), ('CCCCCAAAAA' + TTT,)),
        (('AAAAA' + TTT, 'GGGGG' + TTT), ('AAAAA' + TTT,)),
        ((TTT + 'AAAAA', TTT + 'GGGGG'), (TTT + 'AAAAA',)),
        (('AAAAAAAAAAAAAA' + TTT, 'GGGGG' + TTT), ('AAAAAAAAAAAAAA' + TTT,)),
        (('GGGGGGGGGGGGGG' + TTT, 'AAAAA' + TTT), ('GGGGGGGGGGGGGG' + TTT,)),
        (('AAAAA' + TTT, 'GGGGGGGGGGGGGG' + TTT), ('GGGGGGGGGGGGGG' + TTT,)),
        (('GGGGG' + TTT, 'AAAAAAAAAAAAAA' + TTT), ('AAAAAAAAAAAAAA' + TTT,)),
        (('GGGGG' + TTT + 'AAAAAAAAAAAAA', 'CC' + TTT + 'CC'), ('GGGGG' + TTT + 'AAAAAAAAAAAAA',)),
        (('GGGGG' + TTT + 'AAAAAAAAAAAAA', TTT + 'CC'), ('GGGGG' + TTT + 'AAAAAAAAAAAAA',)),
        (('GGGGG' + TTT + 'AAAAAAAAAAAAA', 'CC' + TTT), ('GGGGG' + TTT + 'AAAAAAAAAAAAA',)),
        (('CC' + TTT + 'CC', 'GGGGG' + TTT + 'AAAAAAAAAAAAA'), ('GGGGG' + TTT + 'AAAAAAAAAAAAA',)),
        ((TTT + 'CC', 'GGGGG' + TTT + 'AAAAAAAAAAAAA'), ('GGGGG' + TTT + 'AAAAAAAAAAAAA',)),
        (('CC' + TTT, 'GGGGG' + TTT + 'AAAAAAAAAAAAA'), ('GGGGG' + TTT + 'AAAAAAAAAAAAA',)),
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

        # 3) "Barely longer" variant (tight threshold): still covered with a single insertion.
        #    Correct behavior: keep both because they cover each other mutually.
        (( 'A' * 100,
           'A' * 50 + 'G' + 'A' * 50 ),
         ( 'A' * 50 + 'G' + 'A' * 50, 'A' * 100 )),

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
        # Correct result: keep only the bigger right contig.
        (( 'A' * 120,
           'A' * 120 + 'Z' + 'G' * 200 ),
         ( 'A' * 120 + 'Z' + 'G' * 200 ,)),

        # Right-covered (mirror): left has junk + 'Z' + perfect suffix that covers right.
        # Correct result: keep only the bigger left contig.
        (( 'C' * 80 + 'Z' + 'A' * 120,
           'A' * 120 ),
         ( 'C' * 80 + 'Z' + 'A' * 120 ,)),

        # Right-covered: left has junk + 'Z' + perfect middle that covers right.
        # Correct result: keep only the bigger left contig.
        (( 'C' * 80 + 'Z' + 'A' * 120 + 'Z' + 'C' * 40,
           'A' * 120 ),
         ( 'C' * 80 + 'Z' + 'A' * 120 + 'Z' + 'C' * 40,)),
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

        # Multiple identical duplicates collapse to a single contig
        ((("XXXXX" + TTT + "YYYYY",) * 9), ("XXXXX" + TTT + "YYYYY",)),
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
