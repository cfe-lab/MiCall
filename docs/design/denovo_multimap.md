# De novo non-stealing remapping and alternative depth

MiCall's de novo pipeline aligns the original trimmed read pairs against
assembled contigs twice: before and after stitching. These stages now ask
Bowtie2 for up to eight placements per read pair (`-k 8`) and retain all
**equally best-scoring, complete, concordant paired** placements. Mapped
(non-denovo) analyses retain their previous Bowtie2 parameters.

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

1. Obtains at most eight candidate paired placements from Bowtie2. This cap
   bounds runtime; it does *not* prove all possible placements were found.
   Logs a warning when read pairs reach that reporting limit. Adding more
   competing contigs can still displace a previously reported placement
   if the limit is saturated; the non-stealing behavior is only guaranteed
   among the actual equally best reported placements.
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
   Q15 consensus support: MiCall's \`sam2aln\` requires two valid mate
   CIGARs to construct a merged read. Supplementary, secondary, malformed,
   or contradictory singleton records are rejected. It does not resurrect
   structurally contradictory pairs simply because they were primary.

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
Bowtie2 QNAME groups (one Bowtie2 thread preserves the ordering).

This is **not PCR/optical duplicate removal** and should not be interpreted
as a count of distinct original biological molecules.

The GUI remapping progress callback uses original read-pair counts rather
than the larger number of SAM alignment records from `-k`. Warnings identify
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
  validation. Assess the impact of the eight-placement reporting cap
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

Integration tests explicitly cover nine otherwise identical contigs, which
saturate the eight-placement cap and emit a warning. It is not legitimate
to treat all nine contigs as confirmed equally supported in that case.

Related history: #141, #200, #264, #290, #719.

## Independent genotype and minority-haplotype controls

Real HIV reference controls from separately resequenced virus strains
and Sanger-confirmed clone mixtures are evaluated alongside these
algorithmic tests. See \`docs/validation/independent_truth_20261107_scientific.md\`
for actual scientific outcomes and release blockers.

A critical limitation of \`-k 8\` is now experimentally observed on real
nine-plasmid HIV read pairs: it omits some equally best-scoring placements
that \`-a\` enumerates, including in a panel with nine known biological
clone references. Mapping cap saturation must never be treated as
confirmation that all genotype hypotheses were covered.
