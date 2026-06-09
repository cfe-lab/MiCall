# Specification: Read-Supported Join Validation for the Referenceless Contig Stitcher

## 1. Summary

Add read-validation to the referenceless contig stitcher. Before accepting a
merge between two contigs, check that the **merged sequence** around the **join
boundary** is supported by reads from the original FASTQ files.  Specifically,
every position within ±`read_length` of the join boundary must be covered by
at least `minimum_read_depth` reads that **exactly match** substrings of the
merged sequence.

This validates the actual join, not just the original contigs independently.
A merge is rejected if no read spans the boundary, even when both source
contigs have local coverage.

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
operate on the merged sequence itself.

---

## 3. Design

### 3.1 Configuration parameters

| Parameter | Type | Default | Scope | Description |
|---|---|---|---|---|
| `--fastq1` | `Path` | — | CLI | Forward reads FASTQ (plain or .gz) |
| `--fastq2` | `Path` | — | CLI | Reverse reads FASTQ (plain or .gz) |
| `--minimum-read-depth` | `int` | `1` | CLI + Context | Min reads per position. `0` = disable |
| `--read-length` | `int` | `150` | CLI + Context | Read length for boundary window |

Coverage validation is **disabled** when neither FASTQ is provided, or when
`--minimum-read-depth 0`. `--fastq1` and `--fastq2` must be provided together;
passing only one is a CLI error.

### 3.2 CLI changes (`micall/core/contig_stitcher.py`)

```python
# without-references subparser:
parser.add_argument('--fastq1', type=Path, default=None,
                    help='Forward reads FASTQ for join validation.')
parser.add_argument('--fastq2', type=Path, default=None,
                    help='Reverse reads FASTQ for join validation.')
parser.add_argument('--minimum-read-depth', type=int, default=1,
                    help='Minimum reads per position near join boundary. '
                         '0 disables validation.')
parser.add_argument('--read-length', type=int, default=150,
                    help='Read length for boundary window. '
                         'Coverage is checked join_boundary +/- read_length.')
```

Validation:

```python
if (args.fastq1 is None) != (args.fastq2 is None):
    parser.error("--fastq1 and --fastq2 must be provided together.")
if args.minimum_read_depth < 0:
    parser.error("--minimum-read-depth must be non-negative.")
if args.read_length < 1:
    parser.error("--read-length must be positive.")
```

### 3.3 Read index (`build_read_index`)

Before stitching, reads from both FASTQ files are collected into a
**read index**: a `Dict[int, Counter[str]]` mapping read_length → Counter of
read sequences with their counts.

```python
def build_read_index(fastq1_path: Path, fastq2_path: Path) -> Dict[int, Counter[str]]:
    read_index: Dict[int, Counter[str]] = {}
    with open_fastq(fastq1) as fq1, open_fastq(fastq2) as fq2:
        while True:
            header1 = fq1.readline()
            if not header1:
                break
            seq1 = fq1.readline().strip()
            fq1.readline(); fq1.readline()
            fq2.readline()
            seq2 = fq2.readline().strip()
            fq2.readline(); fq2.readline()
            for seq in (seq1, seq2):
                if not seq:
                    continue
                read_index.setdefault(len(seq), Counter())[seq] += 1
    return read_index
```

This index is built **once** and reused for all merge candidates throughout
the stitching process.

### 3.4 Context changes (`micall/utils/contig_stitcher_context.py`)

```python
# Read index for join-boundary validation.
# Maps read_length -> Counter of read sequences with their counts.
self.read_index: Dict[int, Counter[str]] = {}

# Coverage validation parameters.
self.minimum_read_depth: int = 1
self.read_length: int = 150
```

### 3.5 Validation function (`check_merged_sequence_support`)

After the merged sequence `M` and join boundary position `B` are determined
(by `merge_by_concordance`), the check validates every position in
`[B - read_length, B + read_length)` on `M`:

```python
def check_merged_sequence_support(
    merged_seq: str,
    join_boundary: int,
    read_index: Dict[int, Counter[str]],
    min_depth: int,
    read_length: int,
) -> bool:
    if min_depth == 0 or not read_index:
        return True

    window_start = max(0, join_boundary - read_length)
    window_end = min(len(merged_seq), join_boundary + read_length)

    # Pick the read length closest to the configured value.
    lookup_len = min(read_index.keys(), key=lambda x: abs(x - read_length))
    counter = read_index[lookup_len]

    for p in range(window_start, window_end):
        # Reads covering position p start in [p-L+1, p].
        s_min = max(0, p - lookup_len + 1)
        s_max = min(p, len(merged_seq) - lookup_len)
        total = 0
        for s in range(s_min, s_max + 1):
            kmer = merged_seq[s:s + lookup_len]
            total += counter.get(kmer, 0)
        if total < min_depth:
            return False
    return True
```

#### Geometric interpretation

For each position `p` near the join boundary, we count every read that
*exactly matches* a substring of the merged sequence `M` and whose span
`[s, s+L)` covers `p`.  At least `min_depth` such reads must exist.

A read that crosses the join boundary at position `B` will cover positions
on both sides of `B`, providing evidence that the merge is correct.  Reads
that only match one original contig do **not** count towards the other
side's coverage — the check fails if any position in the window lacks
sufficient coverage.

### 3.6 Integration in the merge algorithm

The check is inserted into `try_combine_contigs` after the merged sequence
has been constructed (the `merge_by_concordance` call is moved before the
check):

```python
# After score threshold check and covered-contig handling...

result_seq, overlap_size, join_boundary = merge_by_concordance(
    aligned_1, aligned_2, left_remainder, right_remainder
)

# Validate the merged sequence around the join boundary against reads.
ctx = ReferencelessStitcherContext.get()
if not check_merged_sequence_support(
    result_seq, join_boundary,
    ctx.read_index, ctx.minimum_read_depth, ctx.read_length,
):
    return None

# Build ContigWithAligner and return.
```

Covered-contig cases (one contig fully inside the overlap of another) do
not create a join boundary, so the check is skipped for them.

### 3.7 Limitations: exact matching

The validation uses **exact matching only**: a read must match a substring
of the merged sequence perfectly, with zero mutations, insertions, or
deletions.  This means:

- Reads with sequencing errors in the overlap region are not counted.
- Reads from divergent quasispecies variants are not counted.
- The check is conservative — it undercounts true biological support.

However, the primary target is **spurious joins** from coincidental sequence
similarity (low-complexity repeats, chimeric assemblies, contamination).
These will have **zero** spanning reads, which the check correctly rejects.

---

## 4. File-by-file changes

| File | Change |
|---|---|
| `micall/core/contig_stitcher.py` | Add `--fastq1`, `--fastq2`, `--minimum-read-depth`, `--read-length` args to `without-references` subparser. Validate FASTQ pair. Build read index and pass to stitcher. |
| `micall/utils/contig_stitcher_context.py` | Replace `coverage_data` with `read_index`. Keep `minimum_read_depth`, `read_length`. Remove numpy import. |
| `micall/utils/referenceless_contig_stitcher.py` | Add `build_read_index()` (replaces `compute_coverage_from_fastqs`). Add `check_merged_sequence_support()` (replaces `check_read_support`). Restructure `try_combine_contigs` to merge before checking. Modify `merge_by_concordance` to return `join_boundary`. Remove numpy and unused exact_coverage imports. |
| `micall/drivers/sample.py` | Use `build_read_index` instead of `compute_coverage_from_fastqs`. Pass read index to stitcher. |

---

## 5. Edge cases

### 5.1 No FASTQ provided (backward compat)

When neither `--fastq1` nor `--fastq2` is given, `read_index` stays empty.
`check_merged_sequence_support` returns `True` immediately when
`not read_index`. Old behaviour preserved.

### 5.2 `minimum_read_depth = 0` (explicit disable)

The function short-circuits at `min_depth == 0`. No read index is needed.

### 5.3 Empty or very short contig

A contig shorter than `read_length` is fully covered by the window.
The check still works: start/end are clamped to sequence boundaries.

### 5.4 Boundary at contig edge

If `join_boundary - read_length < 0`, the window starts at 0.
If `join_boundary + read_length >= len(merged_seq)`, the window ends at
`len(merged_seq)`.

### 5.5 Duplicate contig sequences

The read index is keyed by read sequence (not contig sequence), so
duplicate or identical contig sequences do not cause coverage confusion.
Each contig is a separate `ContigWithAligner` object with its own identity;
the read index works independently of contig identity.

### 5.6 Read lengths vary

The function picks the read length from the index closest to the configured
`--read-length`. If reads have very different lengths (e.g., 50bp and 250bp),
the closest match is used. Reads of other lengths are not considered.

---

## 6. Testing strategy

### 6.1 Unit tests for `check_merged_sequence_support`

```
TestCheckMergedSequenceSupport:
  test_sufficient_reads_span_boundary      → pass
  test_no_reads_span_boundary_rejected     → fail (correctly)
  test_min_depth_zero_passes               → pass
  test_empty_read_index_passes             → pass
  test_boundary_at_sequence_start          → pass (degenerate)
  test_boundary_at_sequence_end            → pass (degenerate)
  test_partial_coverage_around_boundary    → fail (position near boundary uncovered)
```

### 6.2 End-to-end tests with read index

```
test_stitch_with_read_index:
  Bad case: reads cover original contigs in overlap,
  but no read spans boundary → merge rejected (2 contigs remain)

  Good case: reads also span boundary → merge accepted (1 contig)

test_stitch_with_reads_from_fastq:
  Same as above but reads come from real FASTQ files via build_read_index.
```

### 6.3 CLI validation tests

```
test_cli_fastq_pair_validation:
  --fastq1 without --fastq2 → SystemExit (error)
  --fastq2 without --fastq1 → SystemExit (error)

test_cli_negative_read_depth_rejected:
  --minimum-read-depth -1 → SystemExit (error)

test_cli_zero_read_length_rejected:
  --read-length 0 → SystemExit (error)
```

### 6.4 Existing test suites

The existing exact-match log tests in
`test_referenceless_contig_stitcher_exact.py` are unaffected when FASTQs
are not provided (default behaviour unchanged).
