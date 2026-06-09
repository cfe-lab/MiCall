# Referenceless Stitcher Improvement: Read-Supported Overlap Validation

## Problem

The referenceless contig stitcher validates candidate overlaps between contigs
purely through sequence alignment statistics. The scoring function
`calculate_overlap_score(L, M)` = `(4M - L) * L^-0.90` measures how unlikely
the observed match count `M` in an overlap of length `L` would be under a
correlated 4-letter random model.

This statistical approach can be fooled:

- **Repetitive / low-complexity regions** — homopolymers, dinucleotide repeats,
  and conserved motifs can produce high scores by chance even when the two
  contigs do not actually belong together.
- **Chimeric assemblies** — IVA can join unrelated fragments. The junction might
  have coincidental sequence similarity without any supporting reads.
- **Reference contamination** — contigs derived from different genomes may share
  enough sequence to pass the statistical filter.

No read data enters the picture today. The only evidence for a merge is "do the
two sequences look similar in the overlap region?".

## Existing Tool: `exact_coverage.py`

Located at `micall/utils/exact_coverage.py` (566 lines). It:

1. Takes **paired FASTQ files** (the original reads) + a **FASTA/CSV of contigs**
2. Builds **lazy k-mer indices** of contigs (keyed by read length; the full read
   is the k-mer)
3. For each read (forward and reverse-complement), does a hash-table lookup for
   exact matches in any contig
4. Maintains two numpy arrays per contig:
   - `coverage0[pos]`: counts every exact match (full read extent)
   - `coverage[pos]`: only the inner portion (trims `overlap_size` bases from
     each read end, default = 70)
5. Returns `(coverage0, coverage, contigs)` — per-position integer coverage counts

### Key characteristics

- **"Exact"** means zero mutations/indels — read must match contig perfectly
- Each read occurrence counted once; identical reads from different FASTQ
  entries are independent
- R1 and R2 counted separately
- Lazy index building: indices created per read-length only when first
  encountered

## What Would Change

### Data flow

Currently the stitcher takes `(input_fasta, output_fasta)` only. Coverage data
must reach the inner merge logic. The most natural integration point is the
pipeline driver `micall/drivers/sample.py:438`:

```python
# Current:
referenceless_contig_stitcher(combined_contigs_fasta, stitched_contigs_fasta)

# Proposed:
coverage_arrays = calculate_exact_coverage(
    trimmed1_fastq, trimmed2_fastq, combined_contigs_fasta
)
referenceless_contig_stitcher(
    combined_contigs_fasta, stitched_contigs_fasta, coverage=coverage_arrays
)
```

Coverage arrays would be stored in `ReferencelessStitcherContext` (which already
holds per-run caches) and queried during `try_combine_contigs`.

### Where to insert the check

Inside `try_combine_contigs` (lines 672-790 of
`referenceless_contig_stitcher.py`), **after** the alignment and score pass but
**before** accepting the merge:

```
 1. Upper bound check     → fail early if impossible
 2. Overlap placement     → convolution-based shift
 3. Overlap cutoffs       → mappy-anchored boundary refinement
 4. Alignment + scoring   → global alignment, concordance
 5. Score threshold check → compare against pool minimum
 6. COVERAGE CHECK        ← NEW
 7. Merge or keep covered → return merged contig or single bigger one
```

### What to check

For every position in the overlap region **plus a boundary margin** (one
read-length on each side), verify that at least `minimum_read_depth` reads
cover the position in the original FASTQ files.

### Configuration

New CLI parameters (for `contig_stitcher without-references` mode):

| Parameter | Type | Default | Description |
|---|---|---|---|
| `--fastq1` | Path | — | Forward reads FASTQ (gzip allowed) |
| `--fastq2` | Path | — | Reverse reads FASTQ (gzip allowed) |
| `--minimum-read-depth` | int | 1 | Min reads per position in overlap+margin. 0 = disable |
| `--read-length` | int | 150 | Read length used for boundary margin |

If neither `--fastq1` nor `--fastq2` is provided, coverage validation is
skipped entirely (backward compatible).

## Consequences

### Benefits

**Rejects spurious overlaps from repetitive regions.** The statistical model
can be fooled by sequence that scores high by chance. Read support catches
these because no actual spanning reads exist at the junction.

**Validates against chimeric assemblies.** A chimeric junction from IVA would
have zero read support even if flanking regions are well-covered.

**Empirical grounding.** The stitcher currently trusts an abstract z-score model.
Read support is direct experimental evidence.

**Composability.** Coverage check is an orthogonal safety net. This may allow
lowering `MIN_MATCHES` in the future, accepting weaker sequence matches when
read support is strong.

### Risks

**Computational cost.** Coverage precomputation requires a single pass over all
reads with k-mer lookups. For typical viral samples (100K–1M reads) this is
fast (seconds), but adds measurable latency.

**False rejections on low-coverage samples.** A threshold of even 1× coverage
may reject samples near detection limits, or regions where RNA secondary
structure causes RT drop-off.

**"Exact" matching is conservative.** Genuine quasispecies diversity within the
overlap region will not be counted. Reads with mutations (substitutions, indels)
relative to the contig are invisible to `exact_coverage.py`. The coverage check
therefore undercounts true biological support.

**Boundary sensitivity.** The margin around the overlap (one read-length) is a
heuristic. If reads are systematically shorter or longer than the configured
`--read-length`, the check may be too strict or too loose.

**Integration surface.** Coverage data must be threaded through the stitcher's
entry points and stored in `ReferencelessStitcherContext`. This is a moderate
invasive change to what is currently a pure-function, stateless algorithm.

## Implementation Strategy

### Phase 1 — Precomputed coverage in context

1. Add `coverage_data: Dict[str, np.ndarray]` and coverage parameters to
   `ReferencelessStitcherContext`
2. Compute coverage before stitching using `exact_coverage.py`'s k-mer logic
   (with `overlap_size=0` to count full read extents)
3. In `try_combine_contigs`, add a coverage check that validates the overlap
   window ± `read_length` on each input contig
4. When one or both contigs have no coverage data (e.g., previously stitched
   contigs), skip the check

### Phase 2 — CLI integration

1. Add `--fastq1`, `--fastq2`, `--minimum-read-depth`, `--read-length` to
   `contig_stitcher.py`'s `without-references` subparser
2. When FASTQs are provided, compute coverage and populate context before
   stitching
3. Wire through the pipeline driver `sample.py`

### Phase 3 — Testing

1. Unit tests for the coverage check function with synthetic coverage arrays
2. Parametrized tests on known-good and known-bad overlaps
3. Integration test running the full stitcher with reads

## Design Decisions

**Default `minimum_read_depth = 1`** — Not 2. The stricter threshold (2×) was
considered in the initial draft but rejected because it would reject too many
legitimate merges in low-coverage samples. A value of 1 means "at least one
read supports every position in the overlap vicinity", which catches the case
of zero-support chimeric joins while being permissive for real data.

**Boundary = one read-length** — The margin of exactly one read-length on each
side checks that at least one full read extends into the overlap from each
direction. This is the minimal geometric guarantee that a read could span
the junction. If `read_length > contig_length`, the boundary is clamped to the
contig boundaries.

**Coverage computed on original contigs only** — When a merge involves a
previously-stitched contig (from `o2_loop`), coverage data may not exist for
the synthetic sequence. In this case the check is skipped. The most critical
merges (original-to-original) are always validated.
