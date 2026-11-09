"""Memory-efficient per-position coverage for de novo alternative mappings.

The counters stay sparse for isolated positions and switch to packed arrays
once coverage is sufficiently dense. This avoids storing two Python dicts
per covered base across hundreds of contigs.
"""
from array import array
from collections import Counter
from csv import DictReader, DictWriter
import os


class _DepthPositions:
    def __init__(self):
        self.sparse = {}
        self.maximum = None
        self.average = None
        self.highest = 0

    def add(self, position, fraction):
        maximum = self.maximum
        if maximum is not None and position <= len(maximum) * 2 + 1024:
            assert self.average is not None
            if position > len(maximum):
                old_length = len(maximum)
                gap = position - old_length
                maximum.extend(array('Q', [0]) * gap)
                self.average.extend(array('d', [0.0]) * gap)
                # A distant position may have been stored sparsely before
                # subsequent reads filled the intervening positions.
                # Incorporate those values as the packed span grows.
                for old_pos in range(old_length + 1, position + 1):
                    previous = self.sparse.pop(old_pos, None)
                    if previous is not None:
                        maximum[old_pos - 1] = previous[0]
                        self.average[old_pos - 1] = previous[1]
            maximum[position - 1] += 1
            self.average[position - 1] += fraction
            return

        count, mean = self.sparse.get(position, (0, 0.0))
        self.sparse[position] = (count + 1, mean + fraction)
        if maximum is not None:
            return

        self.highest = max(self.highest, position)
        if len(self.sparse) < 256 or 4 * len(self.sparse) < self.highest:
            return
        self.maximum = array('Q', [0]) * self.highest
        self.average = array('d', [0.0]) * self.highest
        for pos, (value, avg) in self.sparse.items():
            self.maximum[pos - 1] = value
            self.average[pos - 1] = avg
        self.sparse.clear()

    def value(self, position):
        if self.maximum is not None and position <= len(self.maximum):
            assert self.average is not None
            return self.maximum[position - 1], self.average[position - 1]
        return self.sparse.get(position, (0, 0.0))

    def positions(self):
        if self.maximum is None:
            return sorted(self.sparse)
        covered = (i + 1 for i, value in enumerate(self.maximum) if value)
        return sorted([*covered, *self.sparse])


class AlternativeDepths:
    """Per-base max and fractional Q15 depth from accepted read-pair merges."""

    def __init__(self):
        self.references = {}

    def add_fragment(self, placements):
        """Add one physical fragment across its alternative Q15 placements.

        Each reference position gets at most one maximum-depth contribution
        per fragment; fractional depth is the fraction of alternatives
        covering that position. Partially overlapping placements therefore
        retain their unique coverage without doubling their shared span.
        """
        placements = list(placements)
        if not placements:
            return
        hits: Counter[tuple[str, int]] = Counter()
        for reference, merged_sequence in placements:
            for position, base in enumerate(merged_sequence, start=1):
                if base in 'ACGT':
                    hits[reference, position] += 1
        for (reference, position), count in hits.items():
            self.references.setdefault(reference, _DepthPositions()).add(
                position, count / len(placements))

    def write(self, output, reference_csv=None):
        writer = DictWriter(output, ['refname', 'position', 'max_depth',
                                     'avg_depth'], lineterminator=os.linesep)
        writer.writeheader()
        lengths = {}
        if reference_csv is not None:
            lengths = {row['region']: len(row['sequence'])
                       for row in DictReader(reference_csv)}
        for reference in sorted(self.references.keys() | lengths.keys()):
            bins = self.references.get(reference, _DepthPositions())
            positions = (range(1, lengths[reference] + 1)
                         if reference in lengths else bins.positions())
            for position in positions:
                maximum, average = bins.value(position)
                writer.writerow({'refname': reference, 'position': position,
                                 'max_depth': maximum,
                                 'avg_depth': format(average, '.12g')})
