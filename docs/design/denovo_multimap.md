# De novo non-stealing remapping and alternative depth

MiCall's de novo pipeline aligns original trimmed read pairs against
assembled contigs before and after stitching. For panels of up to
16 contigs, and also for 17–32 contigs whose combined reference size is at most 65,536 bases, Bowtie2 uses `-a` to return **all alignments found**
rather than imposing an arbitrary eight-placement cutoff. When those resource ceilings are exceeded it retains bounded `-k 8` output as a resource-safety
fallback. The selector retains equally best-scoring complete concordant
paired placements. De novo mapping uses four Bowtie2 alignment threads
and `--reorder` to keep QNAMEs grouped for MiCall read-pair processing.
Ordinary mapped-mode analyses retain their earlier mapping parameters.

**Neither path certifies complete biological alignment discovery.** The
finite fallback can silently miss alternative contigs if its cap is
saturated; Bowtie2 itself uses a heuristic search even with `-a`.
These are explicit remaining production-release limitations.

## Why

With ordinary single-placement Bowtie2 mapping, a shared fragment is assigned
to just one of several similar contigs. Adding a new contig can therefore
reduce the original contig's read depth and mask an otherwise supported
consensus base.

Multi-mapping means that the read's alignment is compatible with several
positions. It **does not** mean that several independent molecules were
observed. This implementation deliberately exposes both interpretations
rather than silently passing them off as equivalent.

## Selection and pairing

The de novo remapping stage:

1. Uses `-a` for 1–16 input contigs, or 17–32 contigs within the 65,536-base reference budget; larger panels retain bounded
   `-k 8`. A warning identifies reads reaching the finite bound.
   The four-thread mapper uses `--reorder` to guarantee contiguous
   original read-pair QNAME groups. If a QNAME group becomes pathologically
   large (over 8,192 SAM records), it is discarded with a warning rather
   than partially counted or buffered indefinitely.
2. Requires proper, mapped, non-supplementary, non-duplicate alignments whose
   two mates have the same QNAME, reference and reciprocal mate positions.
3. Keeps placements tied for the best **combined R1+R2 alignment score**.
   Suboptimal hits are excluded, even if Bowtie2 reported them.
4. Selects best-scoring placements in a deterministic reference/position
   order. Removes identical placements and overlapping alternative placements
   **on the same reference**. This prevents one molecule being counted twice
   at one position of a single contig; it also conservatively discards the
   non-overlapping tail of a rejected overlapping placement. Distinct
   non-overlapping sites remain separate possible placements.
5. Retains the original QNAME for every placement and writes each R1/R2 pair
   *adjacently*, including SAM secondary flags. This matters because
   `sam2aln.matchmaker()` pairs successive rows by QNAME.
6. Preserves the existing mixed-reference and discordant-primary rescue
   paths when no valid proper pair is available. It also retains a properly
   represented **primary singleton** (one mapped mate, one unmapped mate)
   instead of dropping that original read pair. The singleton remains
   in the remap records and cascade accounting but does **not** become
   Q15 consensus support: MiCall's `sam2aln` requires two valid mate
   CIGARs to construct a merged read. Supplementary, secondary, malformed,
   or contradictory singleton records are rejected. It does not resurrect
   structurally contradictory pairs simply because they were primary.

Bowtie2 `-a` may emit unequal numbers of R1/R2 alternatives on a
*discordant* pair, so the QNAME group can legitimately have an odd
number of SAM records. The selector preserves the valid primary
discordant/singleton pair for legacy rescue and does not falsely report
that case as a malformed SAM group. Odd groups with unresolved or
extraneous alignment records still produce diagnostics.

Malformed alignments must not be recombined across placements. A warning
summarizes mapped read pairs with no usable or safe placement; genuinely
unmapped pairs are not reported as malformed. Only the
de novo `map_to_contigs` calls use this selection logic. Existing
single-placement mapped-mode behaviour remains unchanged.

## Depth definitions

In de novo mode, `multimap_depth.csv` and
`unstitched_multimap_depth.csv` contain one row per contig reference
position (1-based), including positions with zero coverage:

| Field | Definition |
| --- | --- |
| `refname` | Remapping contig name |
| `position` | 1-based contig reference coordinate |
| `max_depth` | Number of accepted original read pairs with a compatible placement covering this position, each given weight 1 on each compatible contig |
| `avg_depth` | Sum of fractional contributions `1/N`, where `N` is the number of accepted placements **for that original read pair** |

Both are calculated from the **same accepted Q15 merged sequence** that
MiCall writes to `aligned.csv`, without a second independent CIGAR merge.
Overlapping R1 and R2 bases count once per placement; low-quality `N`
bases and deleted reference positions do not contribute; inserted bases
have no reference coordinate and do not appear in this depth table.

For example, if one read pair maps to A and B, and another maps only to A:

- A: `max_depth = 2`, `avg_depth = 1.5`
- B: `max_depth = 1`, `avg_depth = 0.5`

**Maximum depth is the default used for the de novo consensus.** The existing
`sam2aln`/`aln2counts` pipeline receives an aligned fragment per accepted
placement. `avg_depth` is a separate audit output, not a substitute for
max-based consensus and not a change in the existing coverage threshold.

**Caution:** Maximum depth is an *upper bound on possible molecule
assignment* at a contig position, not a count of independent physical
molecules across contigs. The sum of max depths across contigs can exceed
the number of original read pairs. Downstream users must not use max
depth as a measure of independent evidence or haplotype abundance.
Fractional depth is likewise conditional on the limited set of reported,
accepted alignments and is not a calibrated posterior probability.
No phasing information is inferred from shared reads.

## Pipeline statistics and progress

The `remap-final <contig>` rows in `remap_counts.csv` intentionally count
alternative placements; their sum can exceed the number of original reads.
When that happens (including when excluded contigs are filtered out), a
`remap-unique` row records the number of original read pairs represented in
the final remap CSV, in the file's existing **read-count** units (two per
pair). `cascade.csv` uses this row to report the distinct physically mapped
pairs, not the sum of hypotheses. The `aligned` cascade column and per-contig
coverage remain placement-based for consensus and are not estimates of
independent molecules across contigs.

Unique mapped read-end accounting is calculated as SAM records are written,
using **constant additional memory with respect to the number of input
read pairs**. Within each QNAME group, the mapped R1 bit and mapped R2 bit
are each counted at most once; a fully mapped fragment counts as two read
ends while a properly represented singleton counts as one. Secondary
placement copies never add to the total. Distinct main and single-reference
rescue streams contain disjoint original read pairs. This assumes the MiCall
input invariant of unique QNAMEs for independent pairs and contiguous
Bowtie2 QNAME groups (the four-thread mode uses `--reorder`).

This is **not PCR/optical duplicate removal** and should not be interpreted
as a count of distinct original biological molecules.

The GUI remapping progress callback uses original read-pair counts rather
than the larger number of SAM alignment records from `-a` or `-k`. Warnings identify
reporting-limit saturation, mapped pairs rejected because no safe paired
placement was found, and odd-length QNAME groups with orphaned SAM records.

## Validation

- Pure-Python tests cover equal and unequal Bowtie2 scores, exact and
  overlapping alternatives, repeat locations, pairing and secondary flags,
  per-read denominators, quality filtering, and CIGAR insertions/soft clips.
- Bowtie2 integration tests run paired FASTQ -> de novo remapping ->
  `sam2aln`, compare one and two near-identical contigs, and check that a
  strictly worse SNV-containing alternative is **not** multi-counted.
- MiCall `Sample.process_post_assembly()` integration tests check both the
  stitched and unstitched de novo outputs, alongside the unchanged mapped mode.
- Existing MiCall remap and `sam2aln` tests remain regression controls.
- In review, validate authentic clinical research controls for false
  genotype calls, minority alleles, raw-read identity, contig masks, runtime,
  and non-denovo invariance; do not treat successful unit tests as biological
  validation. Assess exhaustive-search heuristics and finite-cap behavior
  on repetitive references before a production release.

The alignment stage emits one record pair per retained placement, so runtime
and intermediate I/O can grow approximately linearly with the number of
alternatives even when Bowtie2 itself is fast. The two depth CSVs include
zero-coverage positions and therefore also grow with total contig length.

Depth accumulation uses a hybrid sparse/packed representation: isolated
positions do not allocate large zero-filled prefixes, while densely covered
regions use contiguous 64-bit integer and floating-point arrays instead of
two Python dictionaries per base. Sparse outliers are reconciled if later
reads extend the packed region over them. In an illustrative Python 3.12
memory probe involving 500 contigs and 1,500 covered positions each, this
reduced peak depth-accumulator memory from approximately 104 MiB to 12 MiB.
The tradeoff is additional Python bookkeeping per base; both runtime and
peak memory still depend on alternative alignments and dataset structure.
These costs must be measured on authentic datasets before deployment.

Integration tests cover all nine identical contigs under the
small-panel `-a` path and the resource-capped fallback using
seventeen contigs. Multi-threaded read-pair order is verified with fifty
original pairs across nine equivalently matching references.

The revised all-alignment candidate gained 30 additional truth-correct
callable bases relative to the earlier eight-hit candidate on the
five-strain independently verified HIV control, with zero incorrect calls.
On 59,770 real paired reads, wall time decreased from about 448 to
200 seconds; upstream MiCall was about 174 seconds on that comparison.
The tradeoff remains dataset-dependent and must be validated on
large, repetitive panels without giving up read accounting integrity.

Related history: #141, #200, #264, #290, #719.

## Independent genotype and minority-haplotype controls

Real HIV reference controls from separately resequenced virus strains
and Sanger-confirmed clone mixtures are evaluated alongside these
algorithmic tests. See `docs/validation/independent_truth_20261107_scientific.md`
for actual scientific outcomes and release blockers.

A hard `-k 8` cutoff was experimentally shown to omit
equally best-scoring alternatives present in an independently verified
nine-haplotype HIV panel. The small-panel `-a` path resolves this
specific avoidable truncation, but does not provide a mathematical
guarantee of exhaustive mapping. Even when all reported alignments
are retained, per-base maximum depth remains conditional and cannot
establish genotype presence or true minority-haplotype abundance.

## Extended real-reference large-panel stress

An additional constructed reference-panel stress used 1,200 authentic HIV9
paired reads and 17 authentic reference contigs (nine Sanger-verified HIV
clone sequences, five independently sequenced HIV strains and three
archived HCV contigs). The greater-than-16 fallback reached the eight-hit
limit on **1,167** pairs. Relative to a deeper Bowtie2 all-hit search,
**188** pairs lost at least one reference alternative (**290** missing
alternatives in total). Bowtie2 all-hit search took 10.0 seconds versus
4.1 seconds for the bounded run on this 1,200-pair subset.

This is concrete evidence that the bounded large-panel fallback is **not
yet fit for clinical genotype inference**. Merely increasing a fixed hit
limit is not proven safe or fast at production scale. The reference panel
is deliberately constructed from existing real sequences: it is not an
actual 17-contig patient assembly. The reproducible study script is
docs/validation/real17_panel_audit_20261108.py, and the release verdict
is in docs/validation/final_release_audit_20261108.md.

## Resource-aware 17-contig follow-up and controlled minority truth

The new search policy permits Bowtie2 all-hit reporting on the real
61,951-base, 17-contig competition panel and eliminates the arbitrary
8-hit cutoff in that condition. In real HIV9 A/B testing this produced
1,521 -> 1,537 Sanger-correct callable Q15 MAX bases, zero wrong calls
in either arm, and 1,167 -> 0 capped pair groups. However, the result
is not monotone per clone: A_50 lost 67 callable bases, while H_0.4
and J_0.2 gained 75 and eight. Processing time rose from 10.5 to
16.2 seconds on the shared host. The reference panel is constructed
from authentic biological sequences, not one patient assembly.
See docs/validation/real17_postpolicy_20261108.md.

The resource ceilings are operational, not a mathematical guarantee
about Bowtie2 search completeness or production runtime. Some larger
panels remain bounded by eight hits and can omit genuine alternative
references. Clinical genotype use requires resolving that limitation.

An independently controlled but synthetic mixture of nine Sanger-verified
HIV clones was also tested with 20,000 error-free paired fragments.
With exactly 40 original clone9 pairs (0.2%), both the upstream and
candidate mappers recovered all 40 correct clone assignments but neither
produced any callable clone9 Q15 MAX bases. At 1% (200 original pairs),
upstream produced 450 callable bases and mapped 196 true source pairs;
the candidate produced 453 callable bases and mapped all 200. A 0%
clone9 control gave no callable clone9 bases in either algorithm.
These are favorable synthetic results, NOT independently measured
biological minority frequencies or full-length haplotype phase.
See docs/validation/known_hiv9_mixture_20261108.md.
