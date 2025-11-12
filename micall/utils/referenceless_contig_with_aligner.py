"""
High-level utilities for referenceless contig overlap discovery and anchoring
with lightweight local alignment.

What this module provides
-------------------------
- A contig wrapper (`ContigWithAligner`) exposing multiple aligner views over the
  same sequence:
  - generic mapping anywhere (mappy)
  - anchored-to-left and anchored-to-right mappings (forward/reversed)
- A convolution-based scorer (`find_maximum_overlap`) that estimates the most
  likely relative placement (shift) between two contigs without a reference.
- A thin adapter (`map_overlap`) that prepares and queries the appropriate
  aligner depending on the stitching relation being tested.

Why this exists (in the stitcher)
---------------------------------
The referenceless stitcher needs two complementary capabilities:
1) Coarse placement: roughly where should two contigs sit relative to each other?
   We answer this with a fast convolution across softened one-hot encodings of
   nucleotides. The result is a shift maximizing expected match count, later
   turned into an overlap score.
2) Local, end-aware anchoring: once a candidate window is known, we need to
   decide how much of each side to trust and align. For this, we use mappy in a
   controlled way to act like "align to left end" or "align to right end". The
   stitcher uses the returned anchor points to compute conservative cutoffs of
   the overlap to align and score (see referenceless_contig_stitcher.py:
   cutoffs_left_* / cutoffs_right_* and friends).
"""

import numpy as np

from abc import ABC, abstractmethod
from typing import Iterator, Tuple, Mapping, Literal, NoReturn
from dataclasses import dataclass
from mappy import Aligner as OriginalMappyAligner
from functools import cached_property

from micall.utils.referenceless_score import Score
from micall.utils.contig_stitcher_contigs import Contig
from micall.utils.overlap_stitcher import \
    exp_dropoff_array, find_max_overlap_length
from micall.utils.find_maximum_overlap import \
    get_overlap_results, choose_convolution_method


# Relation requested by the stitcher when probing overlap anchors.
# - "left": anchor at the left end of this contig (we seek earliest start)
# - "right": anchor at the right end (we seek latest end)
# - "cover": allow mapping anywhere (used when one contig may be fully covered)
OverlapRelation = Literal["left", "right", "cover"]


class LocalAligner(ABC):
    @abstractmethod
    def map(self, query: str) -> Iterator[Tuple[int, int]]: ...


class MappyAligner(LocalAligner):
    """Unconstrained mapper backed by mappy.

    What: Allows the query to align anywhere on the contig.
    Why: Used when probing general coverage or when precise end anchoring is not
    required (relation == "cover"). We keep only primary placements as signals
    to the stitcher.
    """

    def __init__(self, seq: str) -> None:
        self.aligner = OriginalMappyAligner(seq=seq)

    def map(self, query: str) -> Iterator[Tuple[int, int]]:
        for x in self.aligner.map(query):
            if x.is_primary:
                yield (x.r_st, x.r_en)


# Padding used to force the mappy DP to respect a chosen end by giving it a
# non-biological, homogeneous context that cannot be mistaken for the real
# contig and prevents pathological sliding.
PADDING = 40


class ForwardAligner(LocalAligner):
    """Right-end anchorer.

    What: Encourages mappings whose right boundary is informative and clamps the
    left boundary with a synthetic homogeneous prefix.
    Why: The stitcher needs the latest end position when the right side is the
    trusted boundary (relation == "right"). We add either A* or C* padding on
    both the reference and query to stabilize the right-edge placement.

    Yielded coordinates: (-inf, end) in spirit; implemented as
    (-very_big_number, end - PADDING) so the consumer can take max(end).
    """

    def __init__(self, seq: str) -> None:
        # Pick a padding base distinct from the natural first base to avoid
        # accidental continuation through padding.
        if seq.startswith('A'):
            self.seq = 'C' * PADDING + seq
        else:
            self.seq = 'A' * PADDING + seq
        self.aligner = OriginalMappyAligner(seq=self.seq)

    def map(self, query: str) -> Iterator[Tuple[int, int]]:
        # Match the query padding to the chosen synthetic flank.
        if self.seq.startswith('A'):
            query = 'A' * PADDING + query
        else:
            query = 'C' * PADDING + query

        for x in self.aligner.map(query):
            if x.is_primary:
                end = x.r_en - PADDING
                if end > 0:
                    # Signal we only trust the right boundary.
                    yield (-99999999999, end)


class ReversedAligner(LocalAligner):
    """Left-end anchorer.

    What: Encourages mappings whose left boundary is informative and clamps the
    right boundary with a synthetic homogeneous suffix.
    Why: The stitcher needs the earliest start position when the left side is
    the trusted boundary (relation == "left"). We add A* or C* padding after
    the sequence and mirror that on the query to stabilize left-edge placement.

    Yielded coordinates: (start, +inf) in spirit; implemented as
    (start, +very_big_number) so the consumer can take min(start).
    """

    def __init__(self, seq: str) -> None:
        if seq.endswith('A'):
            self.seq = seq + 'C' * PADDING
        else:
            self.seq = seq + 'A' * PADDING
        self.aligner = OriginalMappyAligner(seq=self.seq)

    def map(self, query: str) -> Iterator[Tuple[int, int]]:
        if self.seq.endswith('A'):
            query = query + 'A' * PADDING
        else:
            query = query + 'C' * PADDING

        for x in self.aligner.map(query):
            if x.is_primary:
                start = x.r_st
                if start < len(self.seq) - PADDING:
                    # Signal we only trust the left boundary.
                    yield (start, -99999999999999999)


@dataclass(frozen=True)
class ContigWithAligner(Contig):
    """A contig bundled with cached alignment-friendly views.

    What: A lightweight wrapper around `Contig` that exposes cached aligners and
    vectorized representations of the sequence for overlap scoring.
    Why: The stitcher frequently needs to query the same contig in different
    ways; caching avoids repeated construction cost and ensures consistent
    behavior.
    """

    @cached_property
    def mappy_aligner(self) -> LocalAligner:
        return MappyAligner(seq=self.seq)

    @cached_property
    def forward_aligner(self) -> LocalAligner:
        return ForwardAligner(seq=self.seq)

    @cached_property
    def reversed_aligner(self) -> LocalAligner:
        return ReversedAligner(seq=self.seq)

    @staticmethod
    def make(contig: Contig) -> 'ContigWithAligner':
        return ContigWithAligner(name=contig.name, seq=contig.seq, reads_count=contig.reads_count)

    @staticmethod
    def empty() -> 'ContigWithAligner':
        return ContigWithAligner.make(Contig.empty())

    @cached_property
    def nucleotide_seq(self) -> np.ndarray:
        ret = np.frombuffer(self.seq.encode('utf-8'), dtype='S1')
        return ret

    @cached_property
    def alphabet(self) -> Tuple[str, ...]:
        return tuple(sorted(set(self.seq)))

    @cached_property
    def alignment_seqs(self) -> Mapping[str, np.ndarray]:
        """Softened one-hot arrays per nucleotide used in convolution scoring.

        What: For each symbol in the contig's alphabet, build an indicator
        vector with 1's at positions of that symbol and then apply exponential
        drop-off smoothing.
        Why: Smoothing makes the convolution-based overlap estimate more
        tolerant to small local disagreements and sequencing noise, improving
        the stability of the coarse shift detection.
        """
        def to_array(letter: str) -> np.ndarray:
            value = letter.encode('utf-8')
            ret = np.zeros(len(self.nucleotide_seq))
            ret[self.nucleotide_seq == value] = 1
            exp_dropoff_array(ret, factor=8)
            return ret

        return {x: to_array(x) for x in self.alphabet}


def find_maximum_overlap(
    left: ContigWithAligner, right: ContigWithAligner
) -> Tuple[int, float]:
    """Estimate the best relative placement of two contigs.

    What: Computes a cross-signal convolution between symbol-specific softened
    one-hot vectors from the left and the flipped right contig to produce a
    profile of expected matches across all possible shifts.

    Why: Provides the stitcher with a fast, reference-free estimate of where the
    overlap is most likely and a corresponding strength value. The shift is
    later converted into an initial overlap window and used to decide if a
    merge attempt is even worth pursuing.

    Returns:
      (shift, score_like):
        - shift: integer offset placing `right` relative to `left` (negative
          means right starts inside/left of `left`).
        - score_like: a float proportional to overlap quality (monotonic with
          the scoring model used by the stitcher). A value <= 0 means "no
          convincing overlap" and callers treat it as no placement.
    """

    total = np.zeros(len(left.seq) + len(right.seq) - 1)
    method = choose_convolution_method(len(left.seq), len(right.seq))
    keys = sorted(set(list(left.alphabet) + list(right.alphabet)))
    for key in keys:
        if key not in left.alignment_seqs or \
            key not in right.alignment_seqs:
            continue

        x = left.alignment_seqs[key]
        y = np.flip(right.alignment_seqs[key])
        total += method(x, y, mode='full')

    return get_overlap_results(total, len(left.seq), len(right.seq))


def map_overlap(self: ContigWithAligner,
                minimum_score: Score,
                relation: OverlapRelation,
                overlap: str,
                ) -> Iterator[Tuple[int, int]]:
    """Map an overlap probe against this contig under a stitching relation.

    What: Returns candidate (start, end) placements of the provided overlap
    string on this contig, using an aligner chosen by the relation:
      - "left": left-end anchored; we care about minimal start
      - "right": right-end anchored; we care about maximal end
      - "cover": unconstrained mapping; we care about both

    Why: The stitcher uses these placements to compute conservative cutoffs of
    the overlap region to align (see find_overlap_cutoffs in the stitcher).
    Before aligning, we optionally present a trimmed view of the contig to the
    aligner based on a theoretical upper bound on the maximum usable overlap
    length derived from `minimum_score`. This saves work and avoids infeasible
    searches when the required score cannot be attained with the full contig.

    Yields:
      Tuples of (start, end) in coordinates of the original contig. For anchored
      relations one of the bounds is a sentinel extreme as documented in the
      aligner classes; consumers pick min(start) or max(end) accordingly.
    """

    if relation == "left":
        aligner = self.reversed_aligner
    elif relation == "right":
        aligner = self.forward_aligner
    else:
        aligner = self.mappy_aligner
    shift = 0

    if relation != "cover":
        # Compute an upper bound on how much of this contig can contribute to a
        # valid overlap that reaches `minimum_score` given the number of matches
        # cannot exceed the probe length.
        optimistic_number_of_matches = len(overlap)
        max_length = find_max_overlap_length(M=optimistic_number_of_matches,
                                             X=minimum_score,
                                             L_high=len(self.seq),
                                             )

        assert max_length > 0
        assert max_length >= len(overlap)
        assert max_length <= len(self.seq)

        # If the bound is tighter than the full contig length, trim the view and
        # build an anchoring aligner over that view. Remember `shift` so we can
        # translate returned coordinates back to the original contig.
        if max_length < len(self.seq):
            if relation == "left":
                seq = self.seq[-max_length:]
                shift = len(self.seq) - max_length
                aligner = ReversedAligner(seq=seq)
            elif relation == "right":
                seq = self.seq[:max_length]
                shift = 0
                aligner = ForwardAligner(seq=seq)
            else:
                _x: NoReturn = relation

    for start, end in aligner.map(overlap):
        yield (start + shift, end + shift)
