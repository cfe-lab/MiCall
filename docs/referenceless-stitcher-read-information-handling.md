# Specification: Read-Supported Overlap Validation for the Referenceless Contig Stitcher

## 1. Summary

Add read-coverage validation to the referenceless contig stitcher. Before
accepting a merge between two contigs, verify that every base in the overlap
region (plus a one-read-length buffer on each side) is covered by at least
`minimum_read_depth` reads from the original FASTQ files.

---

## 2. Motivation

The current stitcher validates overlaps purely on sequence similarity
(`calculate_overlap_score`). This can accept false merges from:

- Repetitive / low-complexity regions that score high by chance.
- Chimeric contigs produced by the de-novo assembler.
- Cross-contamination where unrelated sequences happen to share similarity.

Read support provides independent experimental evidence. If no read covers a
position in the overlap vicinity, the merge is rejected regardless of sequence
score.

---

## 3. Design

### 3.1 Configuration parameters

| Parameter | Type | Default | Scope | Description |
|---|---|---|---|---|
| `--fastq1` | `Path` | — | CLI | Forward reads FASTQ (plain or .gz) |
| `--fastq2` | `Path` | — | CLI | Reverse reads FASTQ (plain or .gz) |
| `--minimum-read-depth` | `int` | `1` | CLI + Context | Min reads per position. `0` = disable |
| `--read-length` | `int` | `150` | CLI + Context | Read length for boundary margin |

Coverage validation is **disabled** (backward compatible) when neither FASTQ is
provided, or when `--minimum-read-depth 0`.

### 3.2 CLI changes (`micall/core/contig_stitcher.py`)

Add four arguments to the `without-references` subparser (lines 37-40):

```python
# After existing contigs/stitched_contigs arguments, add:
parser.add_argument('--fastq1', type=Path, default=None,
                    help='Forward reads FASTQ for coverage validation.')
parser.add_argument('--fastq2', type=Path, default=None,
                    help='Reverse reads FASTQ for coverage validation.')
parser.add_argument('--minimum-read-depth', type=int, default=1,
                    help='Minimum read coverage to accept an overlap. '
                         '0 disables coverage validation.')
parser.add_argument('--read-length', type=int, default=150,
                    help='Read length used for the boundary margin. '
                         'Coverage is checked overlap ± read_length.')
```

When FASTQs are provided, coverage is computed before calling the stitcher.
The `referenceless_contig_stitcher` entry point gains optional keyword
arguments that populate the context.

### 3.3 Context changes (`micall/utils/contig_stitcher_context.py`)

Add coverage data and parameters to `ReferencelessStitcherContext` (after
the existing caches at line 87):

```python
# Coverage data for read-support validation.
# Maps contig sequence (str) → per-position coverage (np.ndarray[int32]).
# Coverage is computed once for all input contigs before stitching.
coverage_data: Dict[str, np.ndarray] = {}

# Copy of CLI parameters; set by referenceless_contig_stitcher entry point.
minimum_read_depth: int = 1
read_length: int = 150
```

These fields are populated once during `referenceless_contig_stitcher` (or
`referenceless_contig_stitcher_with_ctx`) and are read-only during the
algorithm.

### 3.4 Coverage precomputation

Before any stitching, compute exact coverage once for every input contig:

```
for each (contig_name, contig_seq) in input_contigs:
    coverage_array = np.zeros(len(contig_seq), dtype=np.int32)

for each read in FASTQ(s):
    for each exact match of read in any contig (forward or rc):
        for pos in [match_start, match_end):
            coverage_array[pos] += 1
```

Implementation notes:

- Reuse the k-mer hashing logic from `micall/utils/exact_coverage.py` (the
  `build_kmer_index_for_size` / `find_exact_matches` functions).
- Count both R1 and R2 reads independently. A read pair yields two
  independent counts.
- Use `coverage0` semantics (full read extent, NOT trimmed by
  `overlap_size`). The stitcher's own boundary logic handles edge effects.
- Coverage is computed for **original input contigs only**. Newly stitched
  contigs do not receive coverage arrays — the check is skipped for them
  (see §3.6).

Alternatively, call `calculate_exact_coverage` from `exact_coverage.py`
with `overlap_size=0` and use the returned `coverage0` dictionary.

### 3.5 Coverage validation function

New function (added to `referenceless_contig_stitcher.py`):

```python
def check_read_support(
    left: ContigWithAligner,
    right: ContigWithAligner,
    left_cutoff: int,
    right_cutoff: int,
    ctx: ReferencelessStitcherContext,
) -> bool:
    """Return True if the overlap region has sufficient read coverage.

    Checks every position in the overlap window plus a boundary margin
    of `read_length` on each side on both contigs.

    When coverage data is unavailable for a contig (it was produced by
    a previous merge), the check for that contig is skipped.
    """
    min_depth = ctx.minimum_read_depth
    if min_depth == 0:
        return True  # disabled

    read_len = ctx.read_length
    cov = ctx.coverage_data

    # --- Check left contig ---
    left_cov = cov.get(left.seq)
    if left_cov is not None:
        left_start = max(0, left_cutoff - read_len)
        left_end = min(len(left.seq), left_cutoff + right_cutoff + read_len)
        if np.any(left_cov[left_start:left_end] < min_depth):
            return False

    # --- Check right contig ---
    right_cov = cov.get(right.seq)
    if right_cov is not None:
        right_start = 0  # overlap begins at position 0 on the right contig
        right_end = min(len(right.seq), right_cutoff + read_len)
        if np.any(right_cov[right_start:right_end] < min_depth):
            return False

    return True
```

#### Boundary formulas

For the **left contig**, the overlap window occupies:
`[left_cutoff, left_cutoff + right_cutoff)`.

- Left boundary: one read-length before the overlap starts, clamped to 0.
- Right boundary: one read-length after the overlap ends, clamped to the
  contig length.

```
check_start = max(0, left_cutoff - read_length)
check_end   = min(len(left.seq), left_cutoff + right_cutoff + read_length)
```

For the **right contig**, the overlap window occupies `[0, right_cutoff)`.

- Left boundary: 0 (clamped; cannot go negative).
- Right boundary: one read-length after the overlap ends, clamped to the
  contig length.

```
check_start = 0
check_end   = min(len(right.seq), right_cutoff + read_length)
```

#### Geometric interpretation

Extending one read-length beyond the overlap on each side checks that at
least one full read-length of flanking sequence has coverage. This means:

- A read starting before the overlap and extending into it is counted.
- A read starting inside the overlap and extending past it is counted.
- A read that exactly spans the junction and slightly beyond is counted.

If a contig is too short to accommodate the full margin on one side, the
boundary is clamped to the contig's edge (no artificial extension).

### 3.6 Integration point in the merge algorithm

The check is called inside `try_combine_contigs` (line 756 of
`referenceless_contig_stitcher.py`), **after** the score threshold passes
and **before** the merge decision:

```python
    # Line 756 (existing):
    if result_score < minimum_base_score:
        return None

    # NEW: coverage check (insert here)
    if not check_read_support(
        left, right, left_cutoff, right_cutoff,
        ReferencelessStitcherContext.get(),
    ):
        return None

    # Lines 759-769 (existing): covered-contig handling
    if covered is not None:
        ...
```

The same check applies to both the "covered" and "non-covered" cases. In
the covered case (`covered is not None`), the check is performed on the
**original** `left` and `right` contigs (not `covered`/`bigger`), because
the coverage arrays index by original contig sequence.

### 3.7 When coverage data is unavailable

Coverage arrays are computed for original input contigs only. When a
previously stitched contig participates in a merge (e.g., during the
`o2_loop` greedy phase), its sequence won't be in `coverage_data`. In this
case the check for that specific contig is **skipped** (the other contig
is still checked if it has data).

This means:

- All original-to-original merges are validated.
- Original-to-stitched merges are partially validated (original side only).
- Stitched-to-stitched merges are skipped entirely.

This is a safe default: the most critical merges are the initial ones
between original contigs. As the stitcher builds longer paths, each
component has already been individually validated.

---

## 4. File-by-file changes

| File | Change |
|---|---|
| `micall/core/contig_stitcher.py` | Add `--fastq1`, `--fastq2`, `--minimum-read-depth`, `--read-length` args to `without-references` subparser. When FASTQs provided, compute coverage and pass to stitcher. |
| `micall/utils/contig_stitcher_context.py` | Add `coverage_data`, `minimum_read_depth`, `read_length` fields to `ReferencelessStitcherContext.__init__`. |
| `micall/utils/referenceless_contig_stitcher.py` | Add `check_read_support()` function. Call it inside `try_combine_contigs` after score check. Accept optional coverage parameters in entry points. |
| `micall/utils/exact_coverage.py` | No changes needed. Reuse existing `find_exact_matches` / `build_kmer_index_for_size` or the full `calculate_exact_coverage` with `overlap_size=0`. |
| `micall/drivers/sample.py` | Pass `self.trimmed1_fastq` / `self.trimmed2_fastq` / `minimum_read_depth` / `read_length` to `referenceless_contig_stitcher` call (line 440). |

---

## 5. Edge cases

### 5.1 No FASTQ provided (backward compat)

When neither `--fastq1` nor `--fastq2` is given, coverage data is empty and
`minimum_read_depth` defaults to 0 automatically (or is explicitly set to 0).
`check_read_support` returns `True` immediately. Old behaviour preserved.

### 5.2 `minimum_read_depth = 0` (explicit disable)

The function short-circuits at the top when `min_depth == 0`. No coverage
computation is performed even if FASTQs are provided.

### 5.3 Empty contig

A contig with `len(seq) == 0` has no coverage array. `cov.get("")` returns
`None`, so the check for that contig is skipped. This matches current
behaviour where empty contigs pass through.

### 5.4 Contig shorter than read_length

The boundary formula clamps to `len(contig.seq)`. If the contig is shorter
than `read_length`, the checked region is the entire contig.

### 5.5 Overlap at contig boundary

If the overlap sits flush against one end of a contig (e.g., `left_cutoff == 0`
or `left_cutoff + right_cutoff == len(left.seq)`), the boundary margin on
that side is clamped to 0 or `len(seq)` respectively. No out-of-range access.

### 5.6 No exact matches found

If no read from FASTQ maps exactly to any contig, all coverage arrays are
zero. With `minimum_read_depth ≥ 1`, every merge fails the check and the
stitcher produces no merges (only original contigs are output). This is
correct behaviour: if the reads don't match the contigs at all, any merge
would be speculative.

### 5.7 Stitched contigs from o2_loop

The greedy `o2_loop` phase creates merged contigs not present in the
original input. These have no coverage array. The check for these contigs
is skipped. If both sides of a candidate pair in `o2_loop` are stitched
contigs, the check is fully skipped.

---

## 6. Testing strategy

### 6.1 Unit tests for `check_read_support`

```
test_read_support_both_covered():
    coverage_data = {"ACGT"*10: np.array([5]*40)}
    ctx = ReferencelessStitcherContext()
    ctx.coverage_data = coverage_data
    ctx.minimum_read_depth = 2
    ctx.read_length = 10
    left = ContigWithAligner(None, "ACGT"*10, None)
    right = ContigWithAligner(None, "ACGT"*10, None)
    assert check_read_support(left, right, 20, 20, ctx)  # overlap in middle

test_read_support_insufficient_coverage():
    # coverage_array has zeros in the overlap region
    ...

test_read_support_min_depth_zero():
    # returns True regardless of coverage
    ...

test_read_support_no_coverage_data():
    # contig seq not in coverage_data -> skip
    ...

test_read_support_boundary_clamping():
    # overlap at very start of contig -> left boundary clamped to 0
    ...
```

### 6.2 Integration tests for the stitcher with coverage

```
test_stitch_rejects_uncovered_overlap():
    # Create two contigs that overlap in sequence
    # Provide FASTQ reads that cover everything EXCEPT the overlap region
    # Verify they are NOT stitched

test_stitch_accepts_covered_overlap():
    # Same contigs, reads covering everything including overlap
    # Verify they ARE stitched

test_stitch_backward_compat_no_fastq():
    # No FASTQs provided -> stitching proceeds as before
```

### 6.3 Existing test suites

The existing exact-match log tests in
`test_referenceless_contig_stitcher_exact.py` are not affected when FASTQs
are not provided (default behaviour unchanged).

---

## 7. Implementation order

1. Add `coverage_data`, `minimum_read_depth`, `read_length` to
   `ReferencelessStitcherContext.__init__`.
2. Write `check_read_support()` and its unit tests.
3. Integrate `check_read_support()` into `try_combine_contigs`.
4. Add CLI arguments to `contig_stitcher.py`.
5. Wire coverage precomputation into the `without-references` path.
6. Add integration tests.
7. Update `sample.py` to pass FASTQ paths and parameters.
