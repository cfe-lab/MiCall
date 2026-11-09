---
title: De novo multi-mapping
description: How MiCall selects contig placements and counts depth without stealing reads
---

## What problem does this solve?

An assembled sample can contain overlapping contigs, similar haplotypes, or
repeated sequence. A sequenced fragment may align equally well to more than
one contig or position. Choosing **one** of those alignments arbitrarily
transfers its support away from other equally plausible locations.

In the de novo pipeline, MiCall retains the best-scoring valid **paired
placements**, rather than forcing each physical read pair onto one contig.
It then distinguishes two useful depth questions:

- **Maximum depth**: how many different original fragments *could* support
  this base, considering all of their accepted mappings?
- **Fractional depth**: how much of those fragments' support is allocated to
  this base when their accepted mapping alternatives share the evidence?

Both measurements use Q15 merged reads, not raw Bowtie2 read-end counts.
They are separate from the usual nucleotide and amino-acid counts.

This behavior applies only to **de novo** remapping, both before and after
contig stitching. The reference-based iterative remap path is unchanged.

## From FASTQs to selected alignments

1. `micall/core/remap.py:map_to_contigs` maps trimmed paired FASTQs against
   the de novo contigs using Bowtie2 local paired-end alignment. Unlike the
   older single-placement remap, it requests up to **12** paired
   placements per fragment (`-k 12`), for every contig panel.
   Bowtie2 uses four threads with `--reorder`, keeping each fragment's
   records adjacent. Reaching the limit is a diagnostic, not a failure.
2. `micall/core/denovo_multimap.py:select_best_pairs` examines each
   consecutive query-name group. It requires correctly paired R1/R2 SAM
   records, consistent mate coordinates, valid alignment flags and CIGARs,
   and usable alignment scores. The pair score is the sum of both mates'
   `AS:i:` values.
3. Only placements with the **highest pair score** are retained. Identical
   placements are deduplicated. Distinct placements can overlap on the
   same contig: a shifted alignment can cover bases that another
   alignment does not. Original query names are preserved.
4. If no valid proper pair is selected, a valid primary singleton or
   discordant pair may be retained for the pre-existing
   `MixedReferenceSplitter` rescue path. Fully scoreless records are also
   supported for older SAM inputs; mixed scored/unscored candidates are
   not silently treated as equivalent.
5. The selected alignments are written to `remap.csv` in SAM-column CSV
   format, with two rows per retained paired placement. Multiple
   placements are different hypotheses for **one fragment**, not
   additional sequenced molecules.

The selector caps buffering at 8,192 SAM records per original fragment.
Excessively large groups are discarded rather than partially treated as
complete. Reaching Bowtie2's reporting limit can also mean some alignments
have not been reported, although these may all be worse than the selected
best placement. Both conditions are recorded in `remap_counts.csv` as
**diagnostics**, not reasons to abort clinical reporting. Like MiCall's
earlier single-placement mapping (`-k 1`), this workflow does not require
exhaustive Bowtie2 search. Bowtie2 is a heuristic aligner, so the absence
of a warning does not establish that every possible alignment was found.

## How depth is counted

`micall/core/sam2aln.py` merges both mates of every selected placement
and evaluates Q15 bases. `micall/core/denovo_depth.py` accumulates coverage
per **original read pair**, grouping alternatives by their original query
name. At reference position `p`, suppose a fragment has `N` successfully
merged placements and `k` of them contain a Q15 A/C/G/T at `p`:

- Its contribution to **`avg_depth`** is `k / N`.
- Its contribution to **`max_depth`** is **1 if `k > 0`, otherwise 0**.

Contributions are summed over distinct input fragments. For example:

| One fragment's accepted placements | max depth at A:101 | avg depth at A:101 | avg depth at B:101 |
| --- | ---: | ---: | ---: |
| One placement on A | 1 | 1 | 0 |
| One placement on A, one on B | 1 | 0.5 | 0.5 |
| Two placements on A both cover 101 | 1 | 1 | 0 |
| Two placements on A; only one covers 101 | 1 | 0.5 | 0 |

Thus **max depth never counts an original fragment twice at the same
reference position**, even when its alternative placements overlap.
`avg_depth` expresses ambiguity but need not be an integer.

For the ordinary `aligned.csv` path, bases already represented by another
placement of that fragment at the same contig position and quality cutoff
are masked with lowercase `n`, denoting **no extra coverage**. They must
not be represented as `-`, because a dash is actual deletion evidence.
Nonoverlapping parts of each alternative remain available for downstream
consensus reporting. The original, unmasked accepted merged placements
are used to calculate the separate depth table.

The depth sidecar `multimap_depth.csv` has:

| Column | Meaning |
| --- | --- |
| `refname` | Contig/reference name |
| `position` | One-based reference nucleotide position |
| `max_depth` | Number of distinct fragments with at least one merged Q15 base here |
| `avg_depth` | Sum of fractional Q15 contributions of those fragments |

When a reference sequence is supplied, uncovered positions are included
with zero depths. The depth sidecar is usually written in the sample's
scratch directory; it is not a substitute for `nuc.csv`.

## Fragment counts, unmapped reads, and safeguards

A fragment with multiple mappings must count **once** in the pipeline
cascade. `remap_counts.csv` therefore distinguishes per-contig read-end
counts (`remap <contig>`, `remap-final <contig>`) from
`remap-unique` (distinct mapped read ends) and
`remap-unique-fragments` (distinct mapped fragments, including mapped
singletons). These additional rows appear when the extra information is
needed. `sam2aln` also writes `unique_aligned_fragments.csv` (one integer
field), which counts original fragments with at least one accepted Q15
merged placement. The de novo `cascade.csv` uses these distinct fragment
counts rather than the number of alternative placements.

The `unmapped1.fastq` and `unmapped2.fastq` files contain original
read ends **not accepted in the final remap**, including those rejected by
placement selection. A disk-backed mapped-mate index reconciles the final
output with the input FASTQs. The longstanding `unmapped` row retains
Bowtie2's earlier accounting; `unmapped-after-selection` is emitted when
that differs from the final number of unmapped read ends.

Two QC rows describe limited alignment reporting or discarded evidence:

- `multimap-capped-pairs`: reporting limit reached for one or more
  original fragments; further placements could be missing.
- `multimap-oversized-pairs`: too many SAM records for a fragment;
  its mapping group was discarded.

These are **informational diagnostics**: neither value prevents
reference-guided stitching, the second remap, or clinical consensus
generation. A report can proceed despite a saturated Bowtie2 search,
as it could with the preceding single-placement implementation.

## Where to look in the code

- `micall/core/denovo_multimap.py` — SAM validation, pair scoring,
  alternative selection, alignment limits, and QC diagnostics.
- `micall/core/remap.py` — Bowtie2 invocation, selected `remap.csv`,
  mapped-fragment counting, and final unmapped FASTQs.
- `micall/core/sam2aln.py` — merging, per-fragment grouping, and overlap
  masking for the ordinary aligned read output.
- `micall/core/denovo_depth.py` — maximum and fractional Q15 depths.
- `micall/core/cascade_report.py` and `micall/drivers/sample.py` —
  fragment-based cascade reporting and de novo pipeline integration.

For exact output column definitions, see [Steps and input/output files](../steps.md).
