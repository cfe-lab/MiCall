"""
Represents a maximal-overlap placement between two contigs.

Attributes:
    shift: Relative shift to place `left` vs `right` at their best overlap.
        The value is negative when `right` starts before the end of `left`.
        A value of 0 indicates no overlap (we never construct Overlap with 0).
    size: Length (in bases) of the overlapping window induced by `shift`.
"""

from typing import NamedTuple

Overlap = NamedTuple("Overlap", [("shift", int), ("size", int)])
