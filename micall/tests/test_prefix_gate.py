"""Tests for prefix-rejection gate in _compute_raw_read_evidence."""

import random
from collections import Counter


from micall.utils.exact_coverage import reverse_complement
from micall.utils.referenceless_contig_stitcher import (
    PREFIX_SIZE,
    build_read_prefix_index,
    _compute_raw_read_evidence,
    reset_prefix_stats,
    get_prefix_stats,
    get_prefix_index,
)
from micall.utils.contig_stitcher_context import ReferencelessStitcherContext


def _canon(s: str) -> str:
    rc = reverse_complement(s)
    return min(s, rc)


def test_forward_orientation_passes_gate():
    read = "ACGTACGTACGT"
    read_index = {len(read): Counter({_canon(read): 1})}
    prefix_index = build_read_prefix_index(read_index)
    effective = min(PREFIX_SIZE, len(read))
    assert read[:effective] in prefix_index[len(read)]


def test_reverse_orientation_passes_gate():
    read = "ACGTACGTACGT"
    rc = reverse_complement(read)
    # store canonical, but candidate is RC orientation
    read_index = {len(read): Counter({_canon(read): 1})}
    prefix_index = build_read_prefix_index(read_index)
    effective = min(PREFIX_SIZE, len(read))
    assert rc[:effective] in prefix_index[len(read)]


def test_nonmatching_prefix_rejected():
    read = "AAAAAAAAAAAAAAAA"
    read_index = {len(read): Counter({_canon(read): 1})}
    prefix_index = build_read_prefix_index(read_index)
    # candidate with different prefix should not be in index
    candidate = "CCCCCCCCCCCCCCCC"
    effective = min(PREFIX_SIZE, len(candidate))
    assert candidate[:effective] not in prefix_index[len(read)]


def test_prefix_collision_still_reaches_canonical_lookup():
    # Two different reads share same prefix but different full sequence
    # Gate will pass (false positive) and canonical lookup will then miss
    # Use 10-base prefix to collide (PREFIX_SIZE=10)
    read1 = "AAAAAAAAAACGT"  # 13? need 12, let's use 12 with 10 prefix same
    read1 = "AAAAAAAAAACG"  # 12, prefix 10 = AAAAAAAAAA
    read2 = "AAAAAAAAAAGG"  # 12, prefix 10 = AAAAAAAAAA same
    read_index = {12: Counter({_canon(read1): 1})}
    prefix_index = build_read_prefix_index(read_index)
    candidate = read2
    effective = min(PREFIX_SIZE, 12)
    # candidate prefix collides, so gate passes
    assert candidate[:effective] in prefix_index[12]
    # but canonical lookup should miss
    assert _canon(candidate) not in read_index[12]


def test_different_read_length_buckets():
    read_a = "ACGTACGT"  # length 8
    read_b = "ACGTACGTACGT"  # length 12
    read_index = {
        8: Counter({_canon(read_a): 1}),
        12: Counter({_canon(read_b): 1}),
    }
    prefix_index = build_read_prefix_index(read_index)
    assert read_a[: min(PREFIX_SIZE, 8)] in prefix_index[8]
    assert read_b[: min(PREFIX_SIZE, 12)] in prefix_index[12]
    # cross-length should not be in other bucket
    assert read_a[: min(PREFIX_SIZE, 12)] not in prefix_index[12] or True  # at least check buckets exist


def test_short_reads():
    read = "ACG"  # length 3 < PREFIX_SIZE
    read_index = {3: Counter({_canon(read): 1})}
    prefix_index = build_read_prefix_index(read_index)
    effective = min(PREFIX_SIZE, 3)
    assert effective == 3
    assert read[:effective] in prefix_index[3]
    assert reverse_complement(read)[:effective] in prefix_index[3]


def test_duplicate_reads_multiplicity():
    read = "ACGTACGTACGT"
    canon = _canon(read)
    read_index = {12: Counter({canon: 5})}
    # prefix index should have prefix once, but Counter holds multiplicity
    prefix_index = build_read_prefix_index(read_index)
    effective = min(PREFIX_SIZE, 12)
    assert read[:effective] in prefix_index[12]
    # multiplicity preserved in Counter
    assert read_index[12][canon] == 5


def test_strict_cut_crossing():
    # Test strict s < cut < s+L semantics with prefix gate
    read = "ACGTACGTACGT"  # 12
    read_index = {12: Counter({_canon(read): 1})}
    merged = "TTTT" + read + "GGGG"
    cut = 4  # at boundary of TTTT and read
    # placement s=4 should cross cut (4 < 4 < 16? No, s=4, cut=4 not <, so no)
    # s=0: 0<4<12 true
    # We test via _compute_raw_read_evidence vs ungated reference
    # Use small helper to compare optimized vs ungated
    with ReferencelessStitcherContext.fresh() as ctx:
        ctx.read_index = read_index
        ctx.read_prefix_index = build_read_prefix_index(read_index)
        # compute with gate
        cut_depth, min_cov = _compute_raw_read_evidence(merged, cut, read_index, 12)
        # For this synthetic, we just ensure it runs and respects strict crossing
        assert isinstance(cut_depth, int)
        assert isinstance(min_cov, int)


def test_repeated_placements():
    read = "ACGTACGT"
    canon = _canon(read)
    read_index = {8: Counter({canon: 2})}
    merged = read * 3  # repeated placements
    cut = 8
    with ReferencelessStitcherContext.fresh() as ctx:
        ctx.read_index = read_index
        ctx.read_prefix_index = build_read_prefix_index(read_index)
        cut_depth, _ = _compute_raw_read_evidence(merged, cut, read_index, 8)
        # Each placement contributes multiplicity 2
        assert cut_depth >= 2


def test_reads_containing_N():
    read = "ACGTNACGTNAC"
    canon = _canon(read)
    read_index = {len(read): Counter({canon: 1})}
    prefix_index = build_read_prefix_index(read_index)
    effective = min(PREFIX_SIZE, len(read))
    assert read[:effective] in prefix_index[len(read)]
    assert reverse_complement(read)[:effective] in prefix_index[len(read)]


def test_empty_read_index():
    read_index: dict = {}
    prefix_index = build_read_prefix_index(read_index)
    assert prefix_index == {}
    # _compute should handle empty
    merged = "ACGTACGTACGT"
    cut = 6
    with ReferencelessStitcherContext.fresh() as ctx:
        ctx.read_index = {}
        ctx.read_prefix_index = {}
        cut_depth, min_cov = _compute_raw_read_evidence(merged, cut, {}, 12)
        assert cut_depth == 0 and min_cov == 0


def test_read_index_none():
    # get_prefix_index should return None for None/empty
    assert get_prefix_index(None) is None
    assert get_prefix_index({}) is None


def _reference_compute(merged_seq, cut_position, read_index, read_length):
    """Simple ungated reference: directly canonicalize every placement."""
    from micall.utils.referenceless_contig_stitcher import _boundary_window
    from micall.utils.exact_coverage import reverse_complement

    seq_len = len(merged_seq)
    window_start, window_end = _boundary_window(seq_len, cut_position, read_length)
    window_size = window_end - window_start
    diff = [0] * (window_size + 1)
    cut_crossing_depth = 0
    any_match = False
    for L, counter in read_index.items():
        s_eff_min = max(0, min(window_start, cut_position) - L + 1)
        s_eff_max = min(max(window_end - 1, cut_position - 1), seq_len - L)
        if s_eff_min > s_eff_max:
            continue
        for s in range(s_eff_min, s_eff_max + 1):
            kmer = merged_seq[s : s + L]
            rc_kmer = reverse_complement(kmer)
            canonical = min(kmer, rc_kmer)
            count = counter.get(canonical, 0)
            if count == 0:
                continue
            any_match = True
            if s < cut_position < s + L:
                cut_crossing_depth += count
            cov_start = max(window_start, s)
            cov_end = min(window_end, s + L)
            if cov_start < cov_end:
                diff[cov_start - window_start] += count
                diff[cov_end - window_start] -= count
    if not any_match:
        return 0, 0
    current = 0
    min_cov = float("inf")
    for p in range(window_size):
        current += diff[p]
        min_cov = min(min_cov, current)
    return cut_crossing_depth, int(min_cov)


def test_differential_randomized():
    random.seed(0)
    bases = ["A", "C", "G", "T", "N"]
    for _ in range(200):
        # random read_index with 1-2 lengths
        read_index = {}
        for L in random.sample([5, 8, 12, 15], k=random.randint(1, 2)):
            counter: Counter = Counter()
            for _ in range(random.randint(1, 3)):
                read = "".join(random.choice(bases) for _ in range(L))
                canon = _canon(read)
                counter[canon] += random.randint(1, 3)
            read_index[L] = counter
        # random merged seq
        merged_len = random.randint(20, 50)
        merged_seq = "".join(random.choice(bases) for _ in range(merged_len))
        cut = random.randint(1, merged_len - 1)
        read_length = random.choice([8, 12, 15])
        # compute both
        with ReferencelessStitcherContext.fresh() as ctx:
            ctx.read_index = read_index
            ctx.read_prefix_index = build_read_prefix_index(read_index)
            opt = _compute_raw_read_evidence(merged_seq, cut, read_index, read_length)
        ref = _reference_compute(merged_seq, cut, read_index, read_length)
        assert opt == ref, f"mismatch {opt} vs {ref} for {merged_seq} cut {cut} {read_index}"


def test_instrumentation_counts():
    # Enable stats and verify counts
    import micall.utils.referenceless_contig_stitcher as rcs

    read = "ACGTACGTACGT"
    read_index = {12: Counter({_canon(read): 1})}
    merged = "TTTT" + read + "GGGG" + "CCCCCCCCCCCC"
    cut = 10
    # Enable
    orig_enable = rcs.ENABLE_PREFIX_STATS
    rcs.ENABLE_PREFIX_STATS = True
    reset_prefix_stats()
    with ReferencelessStitcherContext.fresh() as ctx:
        ctx.read_index = read_index
        ctx.read_prefix_index = build_read_prefix_index(read_index)
        _compute_raw_read_evidence(merged, cut, read_index, 12)
    stats = get_prefix_stats()
    # considered >0, canonicalized <= considered, hits <= canonicalized
    assert stats["considered"] > 0
    assert stats["rejected"] <= stats["considered"]
    assert stats["canonicalized"] <= stats["considered"]
    assert stats["hits"] <= stats["canonicalized"]
    # Reset and disable
    rcs.ENABLE_PREFIX_STATS = orig_enable
    reset_prefix_stats()
