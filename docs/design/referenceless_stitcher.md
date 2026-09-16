---
title: Referenceless Contig Stitching in MiCall
---

This document describes the **referenceless** contig stitcher
(`micall/utils/referenceless_contig_stitcher.py`,
invoked as `micall contig_stitcher without-references`).
It is one of two stitching algorithms in MiCall; see
[Contig Stitching in MiCall](stitcher.md) for the overview and
[Referencefull stitcher](referencefull_stitcher.md) for the other
algorithm.

Unless noted otherwise, "the stitcher" below means the referenceless
stitcher.

## 1. Design objective

The referenceless stitcher is a post-de-novo-assembly refinement step
that attempts to combine contigs only when the relationship can be
supported without reference-derived structural information.

Reference independence is intentional.

The stitcher must not decide that contigs belong together merely
because:

* they align near each other on a reference;
* they have the expected reference ordering;
* they have the expected reference orientation;
* joining them would make the result look more like a canonical
  genome;
* a subtype or reference label suggests that they should form one
  genome.

The implementation may carry metadata around (contig names, read
counts where available, per-run caches), but reference-derived
biological expectations must not be the evidence used to establish a
join.

In the current code this restriction is structural: the
referenceless path takes FASTA contigs
(`Contig` / `ContigWithAligner` in
`micall/utils/contig_stitcher_contigs.py` and
`micall/utils/referenceless_contig_with_aligner.py`), never a
reference sequence or reference coordinates. Ordering, overlap
windows, alignments, and scores are all computed from the contig
sequences and from short-read evidence. The standard denovo pipeline
(`micall/drivers/sample.py`) feeds this stitcher the combined
assembler FASTA and writes a stitched FASTA; the referencefull CSV
fields (`ref`, `group_ref`, `match`) do not exist on this path.

## 2. Why this constraint exists

A sample can genuinely contain structure that differs from the
canonical reference, including things such as:

* large deletions;
* inversions;
* rearrangements;
* duplications;
* recombinant or otherwise noncanonical structure.

Such structure may itself be biologically important.

A reference-guided assembly can sometimes improve contiguity by
imposing reference-derived order, but for analyses where structural
fidelity matters this can also be undesirable: the output may look
complete while silently normalizing away the unusual structure.

The referenceless stitcher therefore deliberately asks a narrower
question:

> What additional assembly structure is supported by the sample
> itself?

When evidence is insufficient, leaving contigs separate can be
preferable to inventing a relationship.

## 3. Error asymmetry

A false negative generally means:

```text
two truly related contigs remain separate
```

The original sequence evidence remains visible. A downstream user or
tool can still see both pieces.

A false positive means:

```text
two contigs that should remain separate are fused
```

That can destroy or obscure biological structure. The fused sequence
asserts a junction the sample did not support.

Therefore the algorithm is intentionally conservative. Several
safeguards below — the minimum-agreement score, the independent
shared-k-mer check, the perfect-match rule for contained contigs,
and read validation around proposed joins — all raise the bar for
accepting a merge rather than lowering it.

## 4. Inputs, outputs, and overall flow

**Input:** a FASTA file of de novo contigs. Each record becomes a
`ContigWithAligner` (sequence plus cached aligner views; `reads_count`
is currently `None` on this file path).

**Output:** a FASTA file of refined contigs. Some outputs combine
several input contigs; others pass through unchanged when no
supported join was found.

The top-level flow (`stitch_consensus` in
`micall/utils/referenceless_contig_stitcher.py`) has two phases:

1. **Overlap-path stitching** (`stitch_consensus_overlaps`):
   iteratively select the most probable compatible path through the
   remaining contigs, emit its merged sequence, remove its members
   from consideration, and repeat.
2. **Greedy pairwise cleanup** (`o2_loop` / `try_combine_1`): try
   every unordered pair once per round and merge the first acceptable
   pair found, repeating until no acceptable pair remains.

All per-run state lives in `ReferencelessStitcherContext`
(`micall/utils/contig_stitcher_context.py`): overlap, k-mer,
alignment, cutoff, and read-evidence caches plus read-validation
parameters. This keeps repeated pairwise checks cheap without
changing the acceptance rules.

## 5. Evidence the stitcher uses

A candidate join must survive every applicable check below. The
checks are conjunctive safeguards, not alternative theories of
relatedness:

* **Terminal overlap placement** — a coarse convolution estimate of
  where two contigs would sit relative to each other, reduced to a
  terminal overlap window (section 6).
* **Overlap alignment** — a global pairwise alignment of the two
  overlap windows (`align_queries` in
  `micall/utils/overlap_stitcher.py`).
* **Overlap scoring** — a statistical score of the alignment that
  must clear a minimum-agreement threshold derived from
  `MIN_MATCHES = 40` (section 6).
* **Shared k-mers** — an independent exact-match requirement
  (`KMER_SIZE = 30`) that rejects statistically plausible overlaps
  with no shared exact sequence (section 7).
* **Containment handling** — a separate perfect-match rule when one
  contig is fully covered by another (section 8).
* **Raw-read support** — an independent check that the proposed
  junction is crossed by sample reads (section 9 and
  [Read-Supported Join Validation](../specs/referenceless-stitcher-read-information-handling.md)).
* **Path competition** — individually plausible edges compete for
  membership in a bounded set of candidate paths; only winners
  survive (section 10).
* **Concordance-based construction** — the merged sequence itself is
  cut where local agreement is strongest (section 11).

No step consults a reference genome, reference coordinates, or
expected gene order.

## 6. Overlap discovery and scoring

### 6.1 Coarse placement

The stitcher first needs a hypothesis for *where* two contigs
overlap. `find_maximum_overlap` (in
`micall/utils/referenceless_contig_with_aligner.py`, built on
`micall/utils/find_maximum_overlap.py` and
`micall/utils/overlap_stitcher.py`) answers this with a fast
convolution:

* each contig is expanded into per-symbol indicator vectors;
* each vector is smoothed with an exponential drop-off
  (`exp_dropoff_array`, factor 8), so near-misses still contribute
  weakly and the estimate tolerates small local disagreements;
* cross-correlating the softened vectors across all shifts yields an
  expected-match profile;
* each shift is scored with the same statistical overlap model used
  later (`calculate_overlap_score`), and the best shift becomes the
  candidate placement.

A non-positive best value means "no convincing overlap": the pair is
abandoned (`shift == 0` in `get_overlap`). Otherwise the shift is
converted to a terminal overlap window (`Overlap(shift, size)` in
`micall/utils/referenceless_contig_stitcher_overlap.py`;
`compute_overlap_size`, `normalize_orientation`,
`initial_overlap_windows`).

The model here is deliberately coarse. It proposes a window worth
aligning; it does not itself accept a merge.

### 6.2 End-aware anchoring and cutoffs

Before aligning, the stitcher trims the problem to the part of each
contig it is willing to trust. `map_overlap` queries lightweight
`mappy`-backed views of a contig under a stitching relation:

* `"left"` — anchor at the left end (keep the earliest start);
* `"right"` — anchor at the right end (keep the latest end);
* `"cover"` — unconstrained mapping (used when one contig may be
  fully covered).

End anchoring is implemented with synthetic homogeneous padding
(`ForwardAligner` / `ReversedAligner`) so the underlying mapper
respects the chosen edge rather than sliding to an interior repeat.

The returned anchors become cutoffs
(`compute_overlap_cutoffs` / `find_overlap_cutoffs`, with the
`cutoffs_left_*` / `cutoffs_right_*` helpers) delimiting the overlap
region to align and score. A theoretical upper bound
(`find_max_overlap_length`) can additionally narrow the contig view
presented to the aligner when the required score cannot use the full
length. Cutoffs are cached per contig pair; because they are
monotonic in the acceptance threshold, a cutoff computed for a lower
threshold remains valid for a higher one.

### 6.3 Alignment and concordance

The trimmed overlap windows are globally aligned with Biopython's
`PairwiseAligner` in global mode with penalized end gaps
(`align_queries`). Matches, mismatches, and indels in that alignment
are the evidence for or against the join.

From the alignment the stitcher derives two things:

* a **concordance** profile (`calculate_concordance`): a sliding
  average of per-position agreement, accumulated forward and
  backward with a square-root weighting so sustained runs of matches
  score higher than isolated matches. The eventual merge point is
  chosen where concordance is strongest
  (`sort_concordance_indexes`), with ties broken toward the middle
  of the overlap so cuts stay far from disagreements
  (`merge_by_concordance`; see section 11).
* an **overlap score** (`calculate_overlap_score` with
  `score_alignment`): a z-like rarity score over a four-letter
  alphabet, generalized for correlated genomic sequence with an
  exponent (`alpha = -0.60`, i.e. standard deviation growing as
  `L^0.8` rather than `sqrt(L)`). Higher means more unexpected under
  the null model and therefore stronger evidence. The scored length
  includes a small bonus (`+1` for ordinary overlaps, `+2` for
  covering overlaps) expressing that the overlap is flanked by
  non-matching context.

### 6.4 Minimum agreement and fast rejection

A minimum amount of sequence agreement is required. The raw
threshold is set by `MIN_MATCHES = 40`:
`ACCEPTABLE_STITCHING_SCORE` is the transformed score of an overlap
just above that size, and every candidate edge must ultimately reach
at least the pool's minimum acceptable score (section 10).

To avoid wasted alignments, `try_combine_contigs` / 
`precheck_and_prepare_overlap` applies optimistic upper bounds first:
if even a perfect overlap of the available lengths
(`max_possible_overlap_score`) or of the discovered window
(`optimistic_overlap_score`) cannot reach the needed score, the pair
is rejected before alignment. The transformed score
(`calculate_referenceless_overlap_score`) monotonically amplifies the
raw score (a `999 + (999 * base)^2` shaping) and keeps genuine
scores far from the `SCORE_EPSILON = 1` sentinel used for
covered-contig bookkeeping, so scoring and containment signalling can
never be confused.

## 7. Shared-k-mer requirement

Statistical similarity alone is not always sufficient evidence of a
meaningful overlap: repeats, low-complexity sequence, and smoothed
convolution estimates can all produce plausible-looking scores for
unrelated contigs.

The stitcher therefore applies an independent shared-k-mer check
(`get_kmers`, `does_share_kmers`, `get_overlap`):

* every contig yields the set of its exact k-mers with
  `KMER_SIZE = 30`;
* if both contigs are at least 30 bases long and their k-mer sets
  are disjoint, the pair is rejected before any alignment, however
  good its statistical score would have been;
* k-mer sets are cached per sequence in the stitching context.

Requiring shared exact sequence provides additional specificity: a
genuine terminal overlap of sufficient length should normally share
at least one 30-mer, while coincidental similarity often shares
none.

Its scope is deliberately narrow:

* contigs shorter than 30 bases are exempt (they cannot contain a
  full k-mer to share);
* sharing a k-mer does **not** by itself prove an overlap — it only
  permits the statistical and read checks to proceed;
* failing to share a k-mer rejects the candidate merge but does not
  prove the contigs are biologically unrelated.

## 8. Covered-contig handling

When one contig is fully covered by another
(`calculate_covered`: one sequence length is at most the overlap
size), the stitcher does not perform a normal concordance merge.
Instead it applies a strict perfect-match rule in
`try_combine_contigs`:

* the covered sequence and the corresponding window of the larger
  contig are aligned;
* if every base of the overlap matches (`number_of_matches ==
  overlap.size`), the larger contig is kept and the smaller one is
  recorded as contained (returned with `SCORE_EPSILON` and
  `covered_input` marking which side was absorbed);
* any mismatch means no merge at all: the pair is rejected.

The conservative rationale is:

> An imperfect contained contig may represent error or redundancy,
> but it may also encode real variation. Without sufficient
> evidence, silently absorbing it would destroy that uncertainty.

This does not claim that all imperfect contained contigs are
biologically important. The point is that the algorithm intentionally
refuses to assume that they are not. Exact duplicates collapse
safely; near-duplicates are left alone for downstream analysis
rather than fused on statistical grounds.

Containment is tracked separately from path membership
(`ContigsPath.contigs_ids` versus `contains_contigs_ids`), so a
perfectly covered contig is remembered as explained without
contributing a second copy of its sequence to the merged result.

## 9. Read support

Even a join that passes overlap, k-mer, and containment checks still
proposes a new junction — a sequence that neither input contig
contained on its own. The raw sample reads provide independent
evidence about whether that junction is supported.

At a high level (`check_merged_sequence_support` and its cached
caller in `try_combine_contigs`):

* the candidate merged contig and its join boundary (`join_boundary`
  from `merge_by_concordance`) define a cut position;
* the stitcher requires exact placements of sample reads that
  strictly cross the cut, plus exact coverage of every base in a
  read-length-sized window centred on the cut;
* placements are canonicalized (`min(seq, reverse_complement(seq))`)
  so either strand counts, weighted by FASTQ multiplicity, with each
  valid placement contributing;
* if support is below `minimum_read_depth`, the merge is rejected
  (emitting `ReadSupportRejected` in debug2).

Full contracts — cut-spanning definitions, window geometry,
counting model including the accepted placement-times-multiplicity
overcounting tradeoff, disabled states (`read_index is None` or
`minimum_read_depth == 0` accepts; enabled-but-empty `{}` rejects),
CLI flags (`--fastq1` / `--fastq2`, `--minimum-read-depth`,
`--read-length`), and pipeline defaults (enabled with trimmed FASTQs
at depth 1 in `micall/drivers/sample.py`) — belong to the
implementation spec and are not repeated here. See:

* [Read-Supported Join Validation](../specs/referenceless-stitcher-read-information-handling.md)

## 10. Path and candidate competition

The stitcher does not accept every individually plausible overlap
independently. Candidate relationships compete, and compatible joins
form paths.

The mechanism (`ContigsPath` in
`micall/utils/referenceless_contig_path.py`, `Pool` in
`micall/utils/referenceless_contig_stitcher_pool.py`,
`calculate_all_paths` / `extend_by_1` / `calc_extension`):

* every remaining contig starts as a singleton seed path with
  `SCORE_NOTHING = 0`, seeds sorted longest-first;
* each cycle tries to extend every retained path with every
  remaining contig via `try_combine_contigs`, scoring extensions by
  summing edge scores (`combine_scores`);
* a bounded `Pool` (a `SortedRing` plus sequence deduplication)
  keeps only the best paths: same merged sequence keeps only its
  highest score, and the pool's minimum acceptable score only rises,
  pruning progressively weaker extensions;
* capacity per cycle is set by
  `intrapolate_number_of_alternatives` (`999 / max(1, n - 2)`,
  clamped to `[1, 999]`), bounding total work while still exploring
  alternatives when few contigs remain;
* the best surviving path (`find_most_probable_path`) is emitted,
  its members (including contained ones) are removed from the
  remaining set, and the loop repeats;
* if the best path is a singleton, the stitcher gives up on further
  path extension (`GiveUp`) and emits the rest unchanged;
* the later `o2_loop` performs a final greedy pairwise pass for
  leftovers.

The important consequence is:

> Evidence is evaluated locally at candidate edges, while a final
> multi-contig result can be produced through a chain of supported
> relationships.

A final component therefore asserts a chain of pairwise-supported
joins, not an all-pairs guarantee about every member. Two contigs at
opposite ends of an emitted component were never directly compared;
they are joined because each link in the chain cleared the
thresholds.

## 11. How merged sequence is constructed

When a non-covering pair is accepted,
`merge_by_concordance` builds the output from the global alignment
of the two overlap windows:

* the alignment's concordance profile selects the best split index;
* the left part of the left alignment and the right part of the
  right alignment (dashes removed) become the overlap contribution;
* outer remainders (`left_remainder`, `right_remainder`) are
  prepended and appended unchanged;
* the boundary between the left-derived and right-derived overlap
  chunks is recorded as `join_boundary` for read validation.

The merged contig sums input read counts only when both are
available; on the file path both are currently `None`, so the result
carries `None`. The merged sequence then participates in further
extension rounds as an ordinary contig.

## 12. Limitations

Reference independence does not mean that the algorithm can always
recover biological truth.

If two distinct molecules contain a long, highly similar region and
the sample-intrinsic evidence available to the algorithm does not
phase that region to distinguishing sequence — no conclusive overlap
placement, no shared k-mer that anchors the true junction, no
cut-spanning reads — the relationship may be fundamentally ambiguous
to the stitcher. It will leave the contigs separate, even if a
reference would have suggested an order.

That restriction is about the evidence the implementation is allowed
to use, not an absolute claim about all reference-independent or all
short-read methods. Additional sample-intrinsic evidence such as
longer reads or stronger linkage could, in principle, resolve such
ambiguities without using a reference. The defining restriction is:

> Do not use external reference-derived structural assumptions to
> resolve the ambiguity.

Other limits follow from the conservative design:

* exact read matching undercounts true support when reads carry
  errors or variation relative to the contigs; see the spec for the
  accepted tradeoffs;
* in repetitive sequence, one read may contribute at multiple
  placements, inflating support counts without creating support out
  of nothing;
* the minimum-agreement threshold (`MIN_MATCHES = 40`) and k-mer
  size (30) will miss true short overlaps — a deliberate price for
  specificity;
* greedy path selection and the bounded pool can in principle prefer
  a locally strong chain over a globally better one; capacity tuning
  bounds the search rather than guaranteeing optimality.

## 13. Relationship to referencefull

The [referencefull stitcher](referencefull_stitcher.md) is allowed
to use information the referenceless stitcher deliberately excludes:
a reference sequence, reference coordinates, and reference-derived
ordering and adjacency.

Referencefull may therefore intentionally resolve cases that
referenceless leaves unresolved — for example, placing two
non-overlapping contigs in reference order, or bridging a gap with
no sample-supported overlap. That is not automatically a failure of
either algorithm. One trades structural caution for completeness;
the other trades completeness for reference independence. See
[Contig Stitching in MiCall](stitcher.md) for guidance on which
question each output answers.
