"""
Pool for managing candidate contig paths during referenceless stitching.

This module encapsulates the logic for tracking a bounded set of the best
candidate contig paths while the referenceless contig stitcher explores
possible extensions. It provides:

- A stable, typed container (`Pool`) around a SortedRing of `ContigsPath`
- Deduplication by full sequence (keeps only the highest scoring path per seq)
- A fast check of the current minimum acceptable score for pruning

Design notes
------------
- The pool has a fixed capacity enforced by `SortedRing`, which stores
  the best-scoring paths and can evict worse ones as needed.
- The pool tracks the smallest acceptable score, which is initialized by
  the caller (stitcher) to its global acceptance threshold. As higher
  quality paths fill the pool, this threshold can only increase.
"""

from dataclasses import dataclass
from typing import MutableMapping, Set

from micall.utils.sorted_ring import SortedRing
from micall.utils.referenceless_contig_path import ContigsPath
from micall.utils.referenceless_score import Score


@dataclass
class Pool:
    """A bounded pool of best candidate `ContigsPath`.

    Attributes:
        ring: Sorted container of candidate paths honoring a fixed capacity.
        set: Tracks sequences currently present to support O(1) membership checks.
        existing: Maps sequence -> best `ContigsPath` with that sequence.
        smallest_score: The current minimum acceptable score for new paths.
    """

    ring: SortedRing[ContigsPath]
    set: Set[str]
    existing: MutableMapping[str, ContigsPath]
    smallest_score: Score

    @staticmethod
    def empty(capacity: int, min_acceptable_score: Score) -> "Pool":
        """Create an empty `Pool` with the given capacity and baseline threshold.

        Args:
            capacity: Maximum number of candidates to retain.
            min_acceptable_score: Initial minimum acceptable score. The pool
                will not consider combinations below this, and this threshold
                can only increase as better paths fill the pool.
        """
        ring: SortedRing[ContigsPath] = SortedRing(capacity=capacity)
        return Pool(ring, set(), {}, min_acceptable_score)

    @property
    def min_acceptable_score(self) -> Score:
        """Return the current minimum acceptable score for this pool."""
        return self.smallest_score

    def add(self, path: ContigsPath) -> bool:
        """Insert or update a path candidate if it improves the pool.

        Keeps only the best-scoring candidate for a given full sequence.
        Maintains the minimum acceptable score according to the current
        worst element retained by the `ring` and the initial baseline.

        Returns:
            True if the pool changed (inserted or improved an entry), False otherwise.
        """
        key = path.whole.seq
        alternative = self.existing.get(key)
        if alternative is not None and alternative.get_score() >= path.get_score():
            return False

        # Record/replace the best candidate for this sequence.
        self.existing[key] = path

        deleted_seq = None
        old_size = len(self.ring)

        if key in self.set:
            to_delete_index = -1
            for i, candidate in enumerate(self.ring):
                if candidate.whole.seq == key:
                    to_delete_index = i
                    break

            assert to_delete_index >= 0
            deleted_seq = key
            del self.ring[to_delete_index]
        else:
            if len(self.ring) > 0:
                deleted_seq = self.ring[0].whole.seq

        if self.ring.insert(path):
            new_size = len(self.ring)

            if (
                new_size == old_size
                and deleted_seq is not None
                and deleted_seq in self.set
            ):
                self.set.remove(deleted_seq)

            self.set.add(key)
            # Keep the threshold at least as high as the worst element retained.
            # Never let it drop below the initial baseline.
            self.smallest_score = max(self.ring[0].get_score(), self.smallest_score)
            return True

        return False
