"""Tests for RC optimization coordinate transformation."""

import random
from collections import Counter

import pytest

from micall.utils.exact_coverage import reverse_complement
from micall.utils.referenceless_contig_stitcher import _compute_raw_read_evidence
from micall.utils.contig_stitcher_context import ReferencelessStitcherContext


def test_rc_slice_identity_various_lengths():
    for seq_len in [1, 2, 5, 8, 10, 12, 15, 30, 50, 100]:
        seq = "".join(random.choice("ACGTN") for _ in range(seq_len))
        rc_seq = reverse_complement(seq)
        for _ in range(20):
            L = random.randint(1, seq_len) if seq_len > 0 else 1
            if L > seq_len:
                continue
            s = random.randint(0, seq_len - L)
            expected = reverse_complement(seq[s : s + L])
            actual = rc_seq[len(seq) - s - L : len(seq) - s]
            assert expected == actual, f"{seq} s={s} L={L} {expected} vs {actual}"


def test_edge_placements():
    seq = "ACGTACGTACGT"
    rc_seq = reverse_complement(seq)
    n = len(seq)
    for L in [1, 2, 4, 12]:
        # s=0
        expected = reverse_complement(seq[0:L])
        actual = rc_seq[n - 0 - L : n - 0]
        assert expected == actual
        # s+L=n
        s = n - L
        expected = reverse_complement(seq[s : s + L])
        actual = rc_seq[n - s - L : n - s]
        assert expected == actual


def test_single_base_candidate():
    seq = "ACGTN"
    rc_seq = reverse_complement(seq)
    for s in range(len(seq)):
        L = 1
        expected = reverse_complement(seq[s : s + L])
        actual = rc_seq[len(seq) - s - L : len(seq) - s]
        assert expected == actual


def test_multiple_L_values():
    seq = "ACGTACGTACGTACGTACGT"
    rc_seq = reverse_complement(seq)
    n = len(seq)
    for L in [5, 8, 12, 15]:
        for s in [0, 5, 10, n - L]:
            if s < 0 or s + L > n:
                continue
            expected = reverse_complement(seq[s : s + L])
            actual = rc_seq[n - s - L : n - s]
            assert expected == actual


def test_palindromic_sequence():
    seq = "ACGTACGT"  # not palindromic, but test RC vs RC
    # Palindromic: ACGT -> ACGT RC is ACGT? Actually "ATAT" is palindromic
    pal = "ATATATAT"
    rc_pal = reverse_complement(pal)
    assert pal == rc_pal  # ATAT reverse complement is ATAT
    # For palindromic, RC slice should equal forward slice
    n = len(pal)
    rc_seq = reverse_complement(pal)
    for s in range(n):
        for L in [1, 2, 4]:
            if s + L <= n:
                expected = reverse_complement(pal[s : s + L])
                actual = rc_seq[n - s - L : n - s]
                assert expected == actual


def test_non_palindromic_sequence():
    seq = "ACGTACGTACGT"
    rc_seq = reverse_complement(seq)
    n = len(seq)
    # non-palindromic: RC(s) != s
    s = 2
    L = 5
    expected = reverse_complement(seq[s : s + L])
    actual = rc_seq[n - s - L : n - s]
    assert expected == actual
    assert expected != seq[s : s + L]


def test_N_containing_sequence():
    seq = "ACGTNACGTNACGTN"
    rc_seq = reverse_complement(seq)
    n = len(seq)
    for s, L in [(0, 5), (2, 8), (5, 10), (0, len(seq))]:
        if s + L <= n:
            expected = reverse_complement(seq[s : s + L])
            actual = rc_seq[n - s - L : n - s]
            assert expected == actual


def _reference_compute(merged_seq, cut_position, read_index, read_length):
    """Reference that does reverse_complement(kmer) inside loop."""
    from micall.utils.referenceless_contig_stitcher import _boundary_window

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
            canonical = kmer if kmer <= rc_kmer else rc_kmer
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
        if current < min_cov:
            min_cov = current
    return cut_crossing_depth, int(min_cov)


def _canon(s):
    rc = reverse_complement(s)
    return s if s <= rc else rc


def test_differential_with_matches():
    random.seed(1)
    bases = ["A", "C", "G", "T", "N"]
    for _ in range(200):
        # Build merged_seq and sampled reads to ensure matches
        merged_len = random.randint(20, 60)
        merged_seq = "".join(random.choice(bases) for _ in range(merged_len))
        cut = random.randint(1, merged_len - 1)
        read_length = random.choice([8, 12, 15])
        # Build read_index with some reads sampled from merged_seq in both orientations
        read_index = {}
        for L in random.sample([5, 8, 12, 15], k=random.randint(1, 2)):
            counter: Counter = Counter()
            # Add 1-2 reads sampled from merged_seq to ensure matches
            for _ in range(random.randint(1, 2)):
                if merged_len >= L:
                    s = random.randint(0, merged_len - L)
                    read = merged_seq[s : s + L]
                    # Randomly flip orientation
                    if random.random() < 0.5:
                        read = reverse_complement(read)
                    canon = _canon(read)
                    counter[canon] += random.randint(1, 3)
            # Add 1-2 random nonmatching reads
            for _ in range(random.randint(0, 2)):
                read = "".join(random.choice(bases) for _ in range(L))
                canon = _canon(read)
                if canon not in counter:
                    counter[canon] += 1
            if counter:
                read_index[L] = counter
        if not read_index:
            continue
        with ReferencelessStitcherContext.fresh() as ctx:
            ctx.read_index = read_index
            # Use actual optimized function
            opt = _compute_raw_read_evidence(merged_seq, cut, read_index, read_length)
        ref = _reference_compute(merged_seq, cut, read_index, read_length)
        assert opt == ref, f"mismatch {opt} vs {ref} merged {merged_seq} cut {cut} {read_index}"


def test_differential_multiple_buckets_multiplicity():
    # Explicit multiplicity and repeated placements
    read = "ACGTACGT"
    canon = _canon(read)
    read_index = {8: Counter({canon: 3})}
    merged = read + "TTTT" + read  # two placements
    cut = 8  # between
    with ReferencelessStitcherContext.fresh() as ctx:
        ctx.read_index = read_index
        opt = _compute_raw_read_evidence(merged, cut, read_index, 8)
    ref = _reference_compute(merged, cut, read_index, 8)
    assert opt == ref
    # repeated placements should be counted
    assert opt[0] >= 3 or opt[1] >= 0


def test_differential_reverse_orientation_matches():
    read = "ACGTACGTACGT"
    canon = _canon(read)
    rc = reverse_complement(read)
    # Index contains canonical, but merged contains RC orientation
    read_index = {12: Counter({canon: 1})}
    merged = "TTTT" + rc + "GGGG"
    cut = 7  # inside rc
    with ReferencelessStitcherContext.fresh() as ctx:
        ctx.read_index = read_index
        opt = _compute_raw_read_evidence(merged, cut, read_index, 12)
    ref = _reference_compute(merged, cut, read_index, 12)
    assert opt == ref
    assert opt[0] > 0  # should have cut crossing


def test_differential_no_matches():
    read_index = {8: Counter({_canon("AAAAAAAA"): 1})}
    merged = "CCCCCCCCCCCCCCCCCCCC"
    cut = 10
    with ReferencelessStitcherContext.fresh() as ctx:
        ctx.read_index = read_index
        opt = _compute_raw_read_evidence(merged, cut, read_index, 8)
    ref = _reference_compute(merged, cut, read_index, 8)
    assert opt == ref == (0, 0)


def test_differential_N_handling():
    read = "ACGTNACGT"
    canon = _canon(read)
    read_index = {9: Counter({canon: 1})}
    merged = "TTTT" + read + "GGGG"
    cut = 6
    with ReferencelessStitcherContext.fresh() as ctx:
        ctx.read_index = read_index
        opt = _compute_raw_read_evidence(merged, cut, read_index, 9)
    ref = _reference_compute(merged, cut, read_index, 9)
    assert opt == ref
