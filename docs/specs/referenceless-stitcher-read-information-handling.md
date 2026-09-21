# Read-Supported Join Validation for the Referenceless Contig Stitcher

## 1. Purpose

Reduce false joins in referenceless contig stitching by requiring exact
read-level evidence at the join cut.  Sequence-similarity-based overlap
scoring can accept merges from coincidental similarity, low-complexity
repeats, or chimeric assemblies.  Read evidence from the original FASTQ
files provides an independent check: if no read matches a substring of
the merged sequence that strictly crosses the join cut, the merge is
rejected regardless of the sequence score.

## 2. Scope

This feature applies to the referenceless (de-novo) contig stitching
path in `micall/utils/referenceless_contig_stitcher.py`.  The stitcher
validates candidate joins after `merge_by_concordance()` computes the
merged sequence and join boundary.  Read evidence comes from the paired
FASTQ files that were used for the original assembly.

In the standard MiCall denovo pipeline (`micall/drivers/sample.py`),
validation is enabled by default using trimmed FASTQ files.  CLI users
can also enable it with `--fastq1`/`--fastq2`.

## 3. Definitions

**merged sequence** — the concatenation of two contigs joined at a
cut position determined by concordance analysis of their overlapping
alignment.

**join cut** — the position in the merged sequence where the left
contig ends and the right contig begins (half-open interval).  The base
immediately left of the cut is `merged[cut-1]`; the base immediately
right is `merged[cut]`.

**cut-spanning placement** — an exact match of a read at start
position `s` in the merged sequence such that
`s < cut < s + len(read)`.  A read ending at the cut or starting at the
cut does not span.

**boundary window** — a read-length-sized window centred on the cut:

```
left_radius   = read_length // 2
right_radius  = read_length - left_radius
window_start  = max(0, cut - left_radius)
window_end    = min(len(merged), cut + right_radius)
```

Every base in `[window_start, window_end)` must have coverage at least
`minimum_read_depth`.

**read index** — a `Dict[int, Counter[str]]` mapping read length to a
Counter of canonical read sequences and their FASTQ multiplicities.
Built once by `build_read_index()` before stitching.

**canonical read** — `min(seq, reverse_complement(seq))`.  A read and
its reverse complement share one canonical key.  FASTQ multiplicity is
preserved under that key.

**validation disabled** — when `read_index is None` (no FASTQs
provided) or `minimum_read_depth == 0`, all candidate merges are
accepted without read checking.  When `read_index` is a non-`None`
empty dict `{}`, validation is enabled but no reads were indexed, so
all merges are rejected.

## 4. Configuration and API

### CLI (`micall/core/contig_stitcher.py`, mode `without-references`)

| Argument | Type | Default | Description |
|---|---|---|---|
| `--fastq1` | `Path` | — | Forward reads FASTQ (plain or .gz) |
| `--fastq2` | `Path` | — | Reverse reads FASTQ (plain or .gz) |
| `--minimum-read-depth` | `int` | `1` | Minimum exact-placement depth at join cut and across validation window. `0` disables. |
| `--read-length` | `int` | `150` | Read length for the centred coverage window |

`--fastq1` and `--fastq2` must be provided together.  `--minimum-read-depth`
must be ≥0.  `--read-length` must be ≥1.

### Entry point

```python
def referenceless_contig_stitcher(
    input_fasta,
    output_fasta,
    *,
    read_index=None,
    minimum_read_depth=1,
    read_length=150,
)
```

`read_index=None` → validation disabled.  `read_index={}` → validation
enabled, no reads, merges rejected.

## 5. Required behavior

1. **Cut-spanning support.**  Sum the multiplicity of every canonical
   read whose exact-match placement strictly crosses the join cut
   (`start < cut < end`).  If the total is below `minimum_read_depth`,
   reject the merge.

2. **Boundary-window coverage.**  For every base in the boundary window,
   sum the multiplicity of every canonical read whose exact-match
   placement covers that base.  If any base is below
   `minimum_read_depth`, reject the merge.

3. **Placement × multiplicity counting.**  A canonical read sequence
   that matches at multiple start positions in the merged sequence
   contributes its multiplicity at each placement.  This is
   intentional: it keeps the check efficient (direct counter lookups
   from k-mers at candidate start positions) and preserves FASTQ
   multiplicity.  In repetitive sequence, a single observed read may
   have multiple valid placements; each contributes.

4. **Read index construction.**  `build_read_index()` reads both R1 and
   R2, canonicalizes each sequence, and increments the per-length
   counter.  R1/R2 record counts must match; mismatched files raise
   `ValueError`.  Reads with empty sequences are skipped.

5. **Validation disabled states.**
   - `read_index is None` → all merges accepted (backward compatible).
   - `minimum_read_depth == 0` → all merges accepted.
   - `read_index == {}` and `minimum_read_depth > 0` → all merges
     rejected (validation enabled but no reads found).

## 6. Algorithm

```
check_merged_sequence_support(merged_seq, cut, read_index, min_depth, read_length):

1.  Early exit if disabled (read_index is None or min_depth == 0).
2.  Early reject if enabled-but-empty (not read_index and min_depth > 0).
3.  Compute boundary window coordinates.
4.  For each read length L in read_index:
      For each start s that could affect the window or span the cut:
        kmer = merged_seq[s : s+L]
        canonical = min(kmer, reverse_complement(kmer))
        count = read_index[L].get(canonical, 0)
        if count > 0:
          if s < cut < s+L: add count to cut-crossing total.
          Add count to difference array at [cov_start, cov_end).
5.  If no match found at any start position, reject.
6.  If cut-crossing total < min_depth, reject.
7.  Integrate difference array to compute per-base coverage in window.
    If any base < min_depth, reject.
8.  Accept.
```

### Efficiency

The algorithm scans O(read_length) start positions per read length per
merge candidate.  Each start position does a constant-time k-mer
extraction, canonicalization, and counter lookup.  No read-by-read
alignment or per-read index scan is performed inside the merge loop.

The read index is built once per stitching run.  Memory cost is
proportional to the number of distinct canonical read sequences.
More selective on-demand indexing could be considered later if
profiling shows this to be a bottleneck.

## 7. Default pipeline behavior

In `micall/drivers/sample.py`, the normal denovo path builds a read
index from trimmed FASTQ files and passes it to the referenceless
stitcher with `minimum_read_depth=1`.  This means:

- Read-supported join validation is enabled by default in the standard
  pipeline.
- A merge is rejected unless at least one exact placement (weighted by
  FASTQ multiplicity) crosses the join cut.
- This reduces false joins from coincidental overlap similarity,
  chimeric assembly, or contamination.
- Users who want to disable validation can set `minimum_read_depth=0`
  (or skip the read-index construction in custom pipeline code).

## 8. Limitations

### Exact matching

The check uses exact matching only: a read must match a substring of
the merged sequence perfectly.  Reads with sequencing errors or
biological variation relative to the contig are not counted.  This
undercounts true support in diverse samples.

### Placement × multiplicity overcounting

In repetitive or low-complexity sequence, one canonical read may have
multiple equally valid exact-match placements in the merged sequence.
Each placement contributes its multiplicity, so total support can be
inflated.  This is an accepted tradeoff: the algorithm is efficient,
deterministic, and does not attempt to disambiguate multi-mapping
reads.

The check is still effective because completely unsupported joins
(chimeric, contamination-driven) should have **zero** exact
cut-spanning placements regardless of multiplicity.

### Not a biological validation

Read-supported join validation reduces false joins from unsupported
overlap-only merges.  It does not guarantee the merged sequence is
biologically correct.  It does not replace manual review of
assemblies.

## 9. Test expectations

| Test | Expected behavior |
|---|---|
| `read_index=None` and `min_depth>0` | Validation disabled, merge proceeds |
| `read_index={}` and `min_depth>0` | Validation enabled, no reads, merge rejected |
| `read_index` with reads and `min_depth=0` | Validation disabled, merge proceeds |
| Split-side reads covering left and right independently | Merge rejected (no read crosses cut) |
| At least one read crossing the cut, window covered | Merge accepted |
| One read at 4 spanning placements × count 1 = spanning 4 | `min_depth=4` passes, `5` fails |
| Multiplicity: same canonical at count 10, 4 placements = 40 | `min_depth=40` passes, `41` fails |
| RC-equivalent reads: counts accumulate under canonical key | `min_depth=sum` passes |
| Multiple read lengths all contribute | Total spanning includes all lengths |
| Build read index: duplicates and RC accumulate | Canonical key has correct total count |
| Build read index: empty FASTQs | Returns `{}` |
| `--fastq1` without `--fastq2` | CLI error |
| `--minimum-read-depth -1` | CLI error |
| Sample pipeline: read validation enabled by default | Test verifies non-None read index |

## 10. File structure

| File | Role |
|---|---|
| `micall/utils/exact_coverage.py` | `open_fastq`, `read_fastq_pairs_strict`, `reverse_complement` — shared FASTQ utilities |
| `micall/utils/contig_stitcher_context.py` | `ReferencelessStitcherContext.read_index` field |
| `micall/utils/referenceless_contig_stitcher.py` | `build_read_index`, `check_merged_sequence_support`, integration in `try_combine_contigs` |
| `micall/utils/referenceless_contig_stitcher_events.py` | `ReadSupportRejected` event for debug2 |
| `micall/core/contig_stitcher.py` | CLI arguments for without-references mode |
| `micall/drivers/sample.py` | Pipeline integration with validation enabled by default |
