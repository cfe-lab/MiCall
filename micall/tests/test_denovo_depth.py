"""Regression tests for compact de novo maximum/fractional coverage."""
from collections import defaultdict
from csv import DictReader
from io import StringIO
import random
import tracemalloc

from micall.core.denovo_depth import AlternativeDepths, _DepthPositions


def test_sparse_outlier_is_reconciled_when_dense_span_grows():
    positions = _DepthPositions()
    for pos in range(1, 301):
        positions.add(pos, 0.25)
    assert positions.maximum is not None
    positions.add(8000, 0.25)
    assert positions.sparse[8000] == (1, 0.25)
    for pos in (600, 1200, 2000, 4000, 8000):
        positions.add(pos, 0.25)
    assert positions.value(8000) == (2, 0.5)
    assert 8000 not in positions.sparse
    assert positions.positions() == sorted({
        *range(1, 301), 600, 1200, 2000, 4000, 8000})


def test_isolated_large_offset_does_not_allocate_dense_prefix():
    positions = _DepthPositions()
    tracemalloc.start()
    try:
        positions.add(1_000_000, 0.5)
        _, peak = tracemalloc.get_traced_memory()
    finally:
        tracemalloc.stop()
    assert positions.maximum is None
    assert positions.value(1_000_000) == (1, 0.5)
    assert peak < 100_000


def test_compact_counters_match_reference_for_random_sparse_dense_positions():
    rng = random.Random(617)
    data = _DepthPositions()
    reference = defaultdict(lambda: (0, 0.0))
    for _ in range(5000):
        pos = rng.choice([rng.randint(1, 1500), rng.randint(3000, 9000)])
        weight = rng.choice([0.125, 0.25, 0.5, 1.0])
        count, avg = reference[pos]
        reference[pos] = count + 1, avg + weight
        data.add(pos, weight)
    assert data.positions() == sorted(reference)
    for pos in range(1, 10_001):
        assert data.value(pos) == reference[pos]


def test_alternative_depths_preserve_fractional_and_zero_reference_positions():
    depths = AlternativeDepths()
    depths.add('A', '---ACGT', 0.5)
    depths.add('A', '---A-GT', 0.5)
    depths.add('B', '-' * 15 + 'AC', 1.0)
    output = StringIO()
    depths.write(output, StringIO('region,sequence\nA,ACGTACGTAC\n'
                                  'B,ACGTACGTACGTACGTAC\nC,ACGT\n'))
    rows = list(DictReader(StringIO(output.getvalue())))
    assert len(rows) == 10 + 18 + 4
    values = {(r['refname'], int(r['position'])):
              (int(r['max_depth']), float(r['avg_depth'])) for r in rows}
    assert values['A', 4] == (2, 1.0)
    assert values['A', 5] == (1, 0.5)
    assert values['A', 6] == (2, 1.0)
    assert values['B', 16] == (1, 1.0)
    assert values['C', 1] == (0, 0.0)
    assert values['A', 10] == (0, 0.0)
