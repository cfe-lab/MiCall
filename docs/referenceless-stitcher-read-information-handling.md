# Specification: Read-Supported Join Validation for the Referenceless Contig Stitcher

## 1. Summary

Add read-validation to the referenceless contig stitcher. Before accepting a
merge between two contigs, check that the **merged sequence** around the **join
cut** is supported by reads from the original FASTQ files.

Two conditions must hold:

1. **Cut-spanning support** — at least ``minimum_read_depth`` reads have an
   exact-match interval that strictly crosses the cut position:

   ```
   start < cut_position < end
   ```

   A read that ends *at* the cut or starts *at* the cut does **not** support
   the join.  This prevents split-side coverage (reads independently covering
   left and right) from being accepted.

2. **Boundary-window coverage** — every base in a read-length-sized window
   around the cut has at least ``minimum_read_depth`` exact read coverage:

   ```
   left_radius   = read_length // 2
   right_radius  = read_length - left_radius
   window_start  = max(0, cut_position - left_radius)
   window_end    = min(len(merged), cut_position + right_radius)
   ```

   This ensures not just that *a* read crosses the cut, but that the
   neighbourhood around the join is well supported.

Both checks evaluate support in **coordinates of the proposed merged sequence**,
not in coordinates of the original source contigs.  Independent coverage on the
left and right contigs is insufficient — a read must actually span the junction
in the merged sequence.

---

## 2. Motivation

The current stitcher validates overlaps purely on sequence similarity
(`calculate_overlap_score`). This can accept false merges from:

- Repetitive / low-complexity regions that score high by chance.
- Chimeric contigs produced by the de-novo assembler.
- Cross-contamination where unrelated sequences happen to share similarity.

A naive per-contig coverage check (validating each original contig's overlap
region independently) is insufficient: two independently well-covered contig
ends do **not** prove that a single read spans the join. The check must
operate on the merged sequence itself and require actual cut-spanning reads.

---

## 3. Design

### 3.1 Configuration parameters

| Parameter | Type | Default | Scope | Description |
|---|---|---|---|---|
| `--fastq1` | `Path` | — | CLI | Forward reads FASTQ (plain or .gz) |
| `--fastq2` | `Path` | — | CLI | Reverse reads FASTQ (plain or .gz) |
| `--minimum-read-depth` | `int` | `1` | CLI + Context | Min reads per position. `0` = disable |
| `--read-length` | `int` | `150` | CLI + Context | Read length for window size |

Validation is **disabled** when `read_index=None` (no FASTQs provided) or
`--minimum-read-depth 0`. `--fastq1` and `--fastq2` must be provided together;
passing only one is a CLI error.  `--minimum-read-depth` must be ≥ 0 and
`--read-length` must be ≥ 1.

### 3.2 Read index (`build_read_index`)

Before stitching, reads from both FASTQ files are collected into a
**read index**: a `Dict[int, Counter[str]]` mapping read_length → Counter of
canonical read sequences with multiplicity.  Each read from R1 and R2 is
canonicalized via ``min(seq, reverse_complement(seq))`` and accumulated under
the canonical key.  Multiple FASTQ occurrences increase the count; a read and
its reverse complement accumulate under one key.

```python
def build_read_index(fastq1_path, fastq2_path) -> Dict[int, Counter[str]]:
    ...
```

The index is built **once** and reused for all merge candidates.  The Counter
values reflect the total number of observations of each canonical read in both
FASTQ files.

`read_index=None` (no FASTQs provided) means validation is disabled.
`read_index={}` (FASTQs provided but empty) means validation is enabled but
no reads were indexed — merges will be rejected.

### 3.3 Validation function (`check_merged_sequence_support`)

```python
def check_merged_sequence_support(
    merged_seq: str,
    cut_position: int,
    read_index: Optional[Dict[int, Counter[str]]],
    min_depth: int,
    read_length: int,
) -> Tuple[bool, int, int]:
```

The function proceeds in three stages:

1. **Early exit** — if `read_index is None` or `min_depth == 0`, return True
   (validation disabled).  If `not read_index`, return False (no reads to check).

2. **Coverage computation** — for each read length ``L`` in the index and each
   start position ``s`` whose span can touch the window or cross the cut, the
   k-mer ``merged[s:s+L]`` is extracted and canonicalized via
   ``min(kmer, reverse_complement(kmer))``.  If the canonical key exists in the
   index, its multiplicity (the Counter value) contributes to:

   * **cut-spanning**: the multiplicity is added to the spanning total if
     ``s < cut_position < s + L``.
   * **window coverage**: the multiplicity is added to a difference array for
     the window zone, producing per-position coverage.

   A single read sequence that matches at multiple start positions contributes
   its multiplicity at each placement.  Canonicalization means a read and its
   reverse complement accumulate under one key.

3. **Cut-spanning check** — the total weighted spanning count is compared to
   ``min_depth``.  If insufficient, the join is rejected.

4. **Window-coverage check** — the minimum weighted coverage across all
   positions in ``[window_start, window_end)`` is compared to ``min_depth``.

### 3.4 Integration in `try_combine_contigs`

After score threshold passes and the covered-contig case is handled:

```python
result_seq, overlap_size, join_boundary = merge_by_concordance(...)

ctx = ReferencelessStitcherContext.get()
passed, span_cnt, min_win_cov = check_merged_sequence_support(
    result_seq, join_boundary,
    ctx.read_index, ctx.minimum_read_depth, ctx.read_length,
)
if not passed:
    if ctx.is_debug2:
        log(events.ReadSupportRejected(
            left.unique_name, right.unique_name,
            join_boundary, ctx.minimum_read_depth,
            span_cnt, min_win_cov,
        ))
    return None
```

Covered-contig cases (one contig fully inside the overlap of another) do
not create a join boundary, so the check is skipped for them.

### 3.5 Canonicalization and multiplicity

Read identity is canonicalized via ``canonical = min(seq, reverse_complement(seq))``.
This means:

- A read and its reverse complement share one canonical key; their FASTQ counts
  are summed.
- Multiplicity is preserved: if a read appears 3 times in the FASTQ, its
  canonical key has count 3.
- Each exact-match placement of a canonical read on the merged sequence
  contributes the full count.  A single read sequence that matches at multiple
  start positions contributes at each placement.

In repetitive regions, one read may match many positions, each contributing
its multiplicity to spanning and coverage.

### 3.6 Limitations: exact matching and placements

The validation uses **exact matching only**: a read must match a substring
of the merged sequence perfectly, with zero mutations, insertions, or
deletions.  This means:

- Reads with sequencing errors in the overlap region are not counted.
- Reads from divergent quasispecies variants are not counted.
- In repetitive regions, one read sequence may be placed at many start
  positions, each contributing its full count.  This can inflate coverage
  estimates in low-complexity sequence.

However, the primary target is **spurious joins** from coincidental sequence
similarity (chimeric assemblies, contamination).  These will have **zero**
spanning reads, which the check correctly rejects.

---

## 4. File-by-file changes

| File | Change |
|---|---|
| `micall/core/contig_stitcher.py` | Add `--fastq1`, `--fastq2`, `--minimum-read-depth`, `--read-length` args.  Validate FASTQ pair and argument ranges.  Build read index or pass None. |
| `micall/utils/contig_stitcher_context.py` | `read_index` field (``Optional[Dict[int, Counter[str]]]``, default None).  `minimum_read_depth`, `read_length`. |
| `micall/utils/referenceless_contig_stitcher.py` | `build_read_index()` returning Counter of canonical reads.  `check_merged_sequence_support()` with cut-spanning + window-coverage, counting placements × multiplicity.  Restructured `try_combine_contigs`.  `merge_by_concordance` returns `join_boundary`.  `ReadSupportRejected` event for debug2. |
| `micall/drivers/sample.py` | Use `build_read_index`, pass read index to stitcher. |

---

## 5. Edge cases

### 5.1 No FASTQ provided (backward compat)

`read_index=None` → validation disabled → all merges proceed as before.

### 5.2 `minimum_read_depth = 0`

Check short-circuits.  No reads examined.

### 5.3 Empty FASTQ files

`read_index={}` → validation enabled but no reads → all merges rejected.
This is intentional: if reads were provided but none exist, there is no
evidence to support any merge.

### 5.4 Cut at position 0

Impossible in practice (a join requires both left and right contigs).
If it occurs, the spanning check fails because no read can satisfy
`start < 0 < end`.

### 5.5 Reverse-complement matching

For a k-mer extracted from the merged sequence, both the k-mer itself and
its reverse complement are looked up in the index.  Palindromic k-mers
are counted once to avoid double-counting.

### 5.6 Multiple read lengths

All read lengths in the index contribute to both the spanning and window
checks.  A 150-mer and a 75-mer (from trimmed reads) both count.

---

## 6. Test coverage

### Unit tests (`TestCheckMergedSequenceSupport`)

| Test | What it verifies |
|---|---|
| `test_spanning_reads_accepted` | Reads strictly crossing the cut pass |
| `test_spanning_no_crossing_rejected` | Read ending at cut does not span |
| `test_spanning_start_at_cut_rejected` | Read starting at cut does not span |
| `test_split_left_and_right_rejected` | Independent L/R coverage fails |
| `test_min_depth_zero_disables_spanning_check` | min_depth=0 skips check |
| `test_window_fully_covered_accepted` | All window positions covered |
| `test_window_edge_uncovered_rejected` | One uncovered window position fails |
| `test_none_read_index_disables_check` | read_index=None disables |
| `test_empty_read_index_rejects` | read_index={} rejects |
| `test_reverse_complement_spanning_with_real_bases` | RC reads support the cut |
| `test_multiple_read_lengths` | Multiple lengths contribute |
| `test_multiple_read_lengths_insufficient` | Combined still below min_depth |

### End-to-end tests

| Test | What it verifies |
|---|---|
| `test_stitch_with_read_index_none_disabled` | None → merge proceeds |
| `test_stitch_with_read_index_empty_rejected` | {} → merge rejected |
| `test_stitch_with_read_index_matching_accepted` | Matching reads → merge accepted |
| `test_stitch_with_identical_contig_sequences` | Covered-contig case works |
| `test_stitch_with_reads_from_fastq` | Real FASTQ files → merge accepted |

### CLI validation tests

| Test | What it verifies |
|---|---|
| `test_cli_fastq_pair_validation` | --fastq1 without --fastq2 → error |
| `test_cli_negative_read_depth_rejected` | --minimum-read-depth -1 → error |
| `test_cli_zero_read_length_rejected` | --read-length 0 → error |

---

## 7. Acceptance criteria

- A candidate join is **not** accepted unless exact matching placements
  crossing the join cut (`start < cut < end`), weighted by FASTQ multiplicity,
  reach `minimum_read_depth`.
- The read-length-sized window around the cut must have per-base exact
  placement coverage, weighted by FASTQ multiplicity, of at least
  `minimum_read_depth`.
- The check runs against the **proposed merged sequence**, not the original
  contig sequences.
- A single read sequence matching at multiple start positions contributes its
  multiplicity at each placement (placements × count = total support).
- A read and its reverse complement share one canonical key; their counts are summed.
- Split left/right coverage without a cut-spanning read is **rejected**.
- Empty read index (`{}`) does **not** silently disable validation when
  validation was intended (FASTQs provided).
- `read_index=None` disables validation (backward compatible).
