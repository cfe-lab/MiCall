---
title: De novo multi-mapping
description: Paired-read placement selection and depth accounting on assembled contigs
---

MiCall maps paired-end reads to the contigs produced by de novo assembly.
A read pair can have several equally well-scoring placements, either on
different contigs or at different positions of the same contig. The de novo
remapper keeps those placements together as possible locations of a **single
sequenced fragment**.

This document describes placement selection and the resulting coverage
accounting. For output file fields and their place in the wider workflow,
see [Steps and input/output files](../steps.md).

## Terminology

- **Fragment**: one original sequenced read pair (R1 and R2), identified by
  its query name. A fragment remains one observation regardless of its number
  of alignments.
- **Placement**: an alignment of that fragment's two mates to a particular
  reference location. Its two SAM records belong to the same placement.
- **Alternative placements**: distinct locations at which the same fragment
  has the highest accepted paired-alignment score.
- **Reference position**: a one-based nucleotide position on a contig used
  for a particular remapping pass.

Placement alternatives describe uncertainty about *location*, not additional
reads or independent evidence.

Maximum depth (`max_depth`) is a counterfactual upper bound, not independent
molecular observations and not observed unique-molecule depth: one sequenced
fragment contributes one unit of `max_depth` to *each* compatible placement,
so summing `max_depth` across alternative contigs (or positions) overcounts
molecules. Downstream reports must not misstate `max_depth` as observed
unique-molecule depth; use `avg_depth` (fractional `1/N` attribution per
fragment) when a molecule-conserving summary is needed, and seek
biological/human review before clinical use.

## Mapping and placement selection

The de novo pipeline remaps trimmed paired FASTQs twice: first against the
assembled/unstitched contigs, then against the contigs produced by stitching.
Each pass selects alignments and calculates its own outputs.

1. `micall/core/remap.py:map_to_contigs` runs Bowtie2 in local paired-end
   mode, with four threads, `--reorder`, and `-k 12`. Bowtie2 can report up
   to 12 paired placements per fragment. Records for each query name are
   processed together.
2. `micall/core/denovo_multimap.py:select_best_pairs` checks that each
   candidate comprises the correct R1 and R2 records of a concordant pair:
   consistent reference and mate positions, alignment flags, CIGAR strings,
   and read/quality lengths. Unusable candidate pairs are excluded.
3. The score of a candidate pair is the sum of the two SAM `AS:i:` alignment
   scores. Both mates must have parseable scores for the placement to be
   considered. The selector retains every **distinct valid placement tied for
   the highest pair score** among the reported candidates. Repeated records for
   the same placement are deduplicated; distinct placements can overlap.
4. When there is no valid concordant placement, the primary alignment may
   still provide a mapped singleton or a discordant pair for the
   `MixedReferenceSplitter` remapping path. These cases do not constitute
   additional paired-placement alternatives.
5. The retained SAM records are written to `remap.csv` (SAM fields in CSV
   form), using the original query name for every placement.

Selection operates on the placements Bowtie2 reports. Bowtie2 uses a
heuristic search and the `-k 12` limit bounds reported placements; the
selected set is not a claim that every possible alignment was enumerated.

## One fragment's contribution to depth

`micall/core/sam2aln.py` merges the mates of each selected paired
placement and evaluates bases at quality cutoff Q15.
`micall/core/denovo_depth.py:AlternativeDepths` groups the accepted Q15
merged placements by original fragment before accumulating per-position
depth.

For a given fragment, let:

- `N` be its number of successfully merged Q15 placements;
- `k` be the number of those placements containing a Q15 `A`, `C`, `G`,
  or `T` at a particular reference position.

That fragment contributes **1** to `max_depth` if `k > 0` (otherwise
**0**), and contributes **`k / N`** to `avg_depth`.

The values in `multimap_depth.csv` sum these contributions over original
fragments. Maximum depth therefore counts at most one contribution from
a fragment at any given position. Average depth distributes that fragment's
contribution across its accepted placement alternatives. Both depths are
computed separately for each contig position; average depth can be
fractional.

For example, consider one fragment and positions 101 on contigs A and B:

| Its accepted Q15 placements | A:101 max | A:101 avg | B:101 max | B:101 avg |
| --- | ---: | ---: | ---: | ---: |
| One placement covering A:101 | 1 | 1 | 0 | 0 |
| One covering A:101, one covering B:101 | 1 | 0.5 | 1 | 0.5 |
| Two both covering A:101 | 1 | 1 | 0 | 0 |
| Two on A, only one covering A:101 | 1 | 0.5 | 0 | 0 |

These are contributions from **one fragment**, not the total depth of a
sample. A fragment with no accepted Q15 placement contributes nothing.

### Coverage in `aligned.csv`

The nucleotide and amino-acid reporting path consumes `aligned.csv`.
For each original fragment, contig, and quality cutoff, `sam2aln`
counts each reference position's nucleotide or deletion observation at
most once. When one accepted placement has a Q-qualified nucleotide and
another has a deletion at the same position, the nucleotide takes
precedence, independently of placement order. Repeated observations
are masked with lowercase `n`, meaning *no additional coverage or
deletion evidence*. Leading and trailing alignment padding (`-`) is not
treated as a deletion. Positions uniquely supported by an alternative
placement remain available for reporting.

The separate `multimap_depth.csv` calculation uses the accepted merged
Q15 placements **before** this overlap masking. Thus its maximum and
fractional depths describe the full set of accepted placement hypotheses,
while the aligned-read output avoids counting the same fragment's
overlapping bases repeatedly.

The `clipping.csv` diagnostic also uses fragment-level accounting:
soft-clipped positions from all placements of one fragment are combined
by contig and reference position, giving at most one clipping observation
per fragment at each position.

## Output files and counts

| Output | Interpretation |
| --- | --- |
| `remap.csv` | Retained alignments, potentially multiple paired placements per original fragment, with original query names |
| `aligned.csv` | Merged read placements used for consensus and other downstream reporting; repeated coverage of a fragment is masked with lowercase `n` |
| `multimap_depth.csv` | Q15 per-position `refname,position,max_depth,avg_depth` |
| `unique_aligned_fragments.csv` | One row containing the number of original fragments with at least one accepted Q15 merged placement |
| `remap_counts.csv` | Mapping counts and diagnostics, with the distinctions below |
| `unmapped1.fastq`, `unmapped2.fastq` | Original read ends without a selected mapped placement, including rejected pairs and unmapped mates |
| `cascade.csv` | Sample-level counts of distinct fragments through remapping and alignment |

In `multimap_depth.csv`, `position` is one-based. When reference
sequences are available, the file also includes positions with zero
coverage. The file is a depth sidecar, not a replacement for `nuc.csv`.

The rows of `remap_counts.csv` have different counting units:

- `raw`: number of input **read ends**.
- `remap <contig>` and `remap-final <contig>`: mapped **read ends per
  contig**. Alternatives can appear against multiple contigs, so these
  counts are not distinct-fragment totals.
- `remap-fragments`: **distinct mapped fragments**, including those with
  only one mapped mate and those mapped solely to an excluded-seed contig.
  This row is emitted in every de novo remap, including when its count is
  zero. The cascade uses it directly.
- `unmapped`: read ends without a mapped placement **after** placement
  selection and mate rescue. The count matches the two final unmapped
  FASTQs; mapped ends on excluded-seed contigs are not included.
- `multimap-capped-pairs`: original fragments reaching Bowtie2's reporting
  limit. Written for every de novo remap, including a zero value, so the
  saturation status is machine-readable without inspecting logs.

MiCall tracks mapped mates by original query name when producing the
final unmapped FASTQs. Alignments on excluded-seed contigs are omitted
from `remap.csv` and downstream consensus reporting but retain their
mapped classification for fragment counts and unmapped FASTQs.
One fragment with multiple placements still counts once as mapped, and
once as aligned if at least one reported merged Q15 placement was accepted.

## Alignment limits and diagnostics

`remap_counts.csv` may contain two additional diagnostic rows:

- `multimap-capped-pairs`: number of original fragments for which Bowtie2
  reached the `-k 12` reporting limit. Other placements may exist; this
  does not establish that any omitted placement ties for the best score.
- `multimap-oversized-pairs`: number of query-name groups containing more
  than 8,192 SAM records. The entire group is excluded from selection.

Both remapping passes record this diagnostic in their own
`remap_counts.csv` (the unstitched pass uses the `unstitched_` prefix).
When the capped-pair count is nonzero, some placements may be absent from
the retained alignments. If omitted placements tie for best score, support
for their reference positions may be understated, including in downstream
coverage and variant reporting. The count alone cannot determine whether
that happened: a saturated search can also report a unique best alignment
with additional lower-scoring candidates. A zero count does not prove
exhaustive mapping either, because Bowtie2's search is heuristic.

These diagnostics are informational and do not suppress clinical reporting
or label specific variant calls as incomplete. For coverage interpretation,
inspect the relevant remapping pass's diagnostic count alongside the
reported alignments and depth; do not infer that omitted best-scoring
placements either exist or do not exist solely from saturation.

## Validation status and risks (follow-up note)

Docs-only clarification on top of `e0073cb5` (PR #1734, merged 2026-11-10;
6/6 check-runs success at the merge head: build, main, docker-test,
singularity-test, deploy, report-build-status). No code, mapped (non-denovo)
behavior, or reference/threshold change in this follow-up, so no runtime
impact beyond the base commit.

- Cheap focused re-runs in this follow-up clone: `test_denovo_depth.py`
  4 passed; `test_denovo_multimap.py` 80 passed with 2 failures at collection
  from missing `gotoh` (`test_callback_counts_unique_pairs...`,
  `test_every_denovo_search...`) — environment gap, not a code regression.
  `test_denovo_reporting_diagnostics.py` was not runnable here for the same
  missing-`gotoh` reason.
- Board-reported (not re-run here) merged-line evidence: full regression
  ~4415 passed / 5 skipped / 54 xfailed with Ruff and mypy clean; 198/198
  archived fixed-contig remaps passing the fail-closed guard; real HIV/HCV
  A/B comparisons showing recovered shared-depth positions but with
  allele-retention caveats (reported single-allele flip and mixture-symbol
  changes at clinical cutoffs; gains are nonmonotone across contigs).
- Risk: `max_depth` consensus is a counterfactual upper bound (see above);
  false-positive coverage/allele risk remains where alternatives are
  correlated. Clinical production remains NO-GO pending independent typed
  controls, calibrated minority-fraction work, and human/scientific review.
- Coordination: Ivanka `ottojung/ivanka#147` (Beta gate: validate molecular
  accounting in conditional HCV remapping) is still open with no comments at
  time of writing; related MiCall history (#141, #264, #290, #200, #719,
  #1428) is tracked in the Antonina issue, not re-litigated here.

## Implementation map

- `micall/core/remap.py` — Bowtie2 invocation, candidate selection
  integration, remap counts, mate rescue, and unmapped FASTQs.
- `micall/core/denovo_multimap.py` — validation, scoring, and selection
  of paired placements.
- `micall/core/sam2aln.py` — merged reads, grouping by original query
  name, and per-fragment overlap masking.
- `micall/core/denovo_depth.py` — Q15 maximum and fractional depth.
- `micall/core/cascade_report.py` — distinct mapped/aligned fragment
  accounting.
- `micall/drivers/sample.py` — both de novo remapping passes and
  downstream reporting.
