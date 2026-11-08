# Real HIV/HCV validation of MiCall de novo non-stealing remapping

This folder contains the technical evidence accompanying the proposed MiCall de novo mapper change. **The observed coverage recovery is not a claim of improved clinical accuracy.**

## Production readiness and independent truth

**Current release decision: NO-GO.** See the
[latest independent release audit](final_release_audit_20261108.md)
for the tested clinical and engineering blockers. The
[7 November production release-gate assessment](production_gate_20261107.md)
is retained as a historical snapshot.

The additional validation against independently sequenced HIV clones
uses the [independent genotype and linked-haplotype report](independent_truth_20261107_scientific.md),
the [short release verdict](independent_truth_20261107_executive.md),
and an alignment-independent exact private-21mer phase audit
(phased_hiv9_motifs.py). These supplement, rather than replace,
the fixed-contig archived HIV/HCV A/B results below.

The revised mapper now uses all Bowtie2-reported alignments for panels
with at most 16 contigs, and for 17–32 contigs with at most 65,536
combined reference bases. All larger/longer panels retain the capped
search. Both have explicit limitations. Known clinical genotype/mixture
truth and full diagnostic pipeline sign-off are still outstanding.

## Read the results

The most recent independent-release evidence is in
[the final-candidate genotype and minority-haplotype audit](final_release_audit_20261108.md)
and its [aggregate machine-readable evidence](final_release_audit_20261108_metrics.json).
It covers separately verified HIV5/HIV9 reference truth, a deliberately
corrupted-reference challenge, alignment-independent paired-marker linkage,
and source-tracked HCV minority/absent-decoy controls.

A separate [known-source, Sanger-reference HIV9 minority challenge](known_hiv9_mixture_20261108.md)
tests zero, 0.2% and 1% clone9 abundance with exactly controlled
error-free paired fragments. Its
[aggregate truth measurements](known_hiv9_mixture_20261108_metrics.json)
are independent of the mapper's read assignment, but do not establish
clinical sensitivity. The
[post-policy 17-reference A/B report](real17_postpolicy_20261108.md)
documents recovered HIV9 callable bases and clone-level regressions
using genuine HIV reads.

Further original-read tests are documented in the
[0.2%, 0.4% and 0.99% source dilution study](real_hcv_dilution_20261108.md),
its [aggregate evidence](real_hcv_dilution_20261108_metrics.json), and the
[17-authentic-reference large-panel stress script](real17_panel_audit_20261108.py).
The new dilution studies use only observed read pairs, not generated viral
sequences; the 17-panel experiment combines already sequenced viral contigs
to test the alignment cap. In this engineered 17-panel challenge, 188 of
1,200 original read pairs missed at least one compatible reference under
the capped eight-placement search compared with Bowtie2 all-hit reporting.

The source-dilution and reference-stress scripts require private inputs
and emit original reads, names, reference sequences and SAM only inside an
ignored local directory. Do not commit those experimental intermediates.

**Current clinical production decision: NO-GO.** This is an intentionally
strict assessment. Correct genotype presence/absence, independently calibrated
minority fractions, full MiCall downstream reporting and large-panel resource
behavior have not all been validated.

- [Executive report](real_hiv_hcv_20261106_executive.md): concise interpretation, major changes, and review recommendation
- [Scientific report](real_hiv_hcv_20261106_scientific.md): study design, per-case results, coverage and allele accounting, performance caveats, limitations
- [Standalone dashboard](real_hiv_hcv_20261106_dashboard.html): self-contained HTML with embedded charts and per-case comparison
- [Machine-readable results](real_hiv_hcv_20261106_metrics.csv): pseudonymous per-case aggregate metrics
- [Figures](figures/): editable SVG charts for callable positions, mapping ambiguity and elapsed time

**Primary finding:** in 15 completed upstream/non-stealing paired comparisons of real HIV/HCV sequencing data, the candidate gained 3,088 Q15 MAX callable contig positions, lost three, and changed one callable nucleotide. Of the gained positions, 2,123 had fractional depth below MiCall's 100× threshold. The implementation is therefore **experimental and requires independent variant truth**.

## Experiment design

The private source corpus contains 198 archived clinical-research samples. We deliberately selected 16 case-conditions from 14 specimens with genuine paired-end FASTQ and archived assembly contigs: nine HCV-focused and seven HIV-focused conditions, including two HCV specimens tested with both archived IVA and multi-contig Ivanka assembly outputs.

Both arms used *exactly the same* input FASTQ pairs and fixed contig CSVs. The controlled chain is:

```
real paired FASTQ + fixed assembled contigs
    → MiCall map_to_contigs
    → MiCall sam2aln
    → MiCall aln2counts
    → Q15 MAX consensus and per-base depth comparison
```

Control: upstream MiCall commit `cba1898a`. Candidate: `2bc94b59366b0ff22a8a673bdb961c28ba75487a`. Both used Python 3.12 and Bowtie2 2.2.8. Changes upstream of this chain (assembler, primer trimming, G2P) are intentionally not part of this comparison. Independent patient genotype truth was **not** available.

An unmodified upstream arm failed on one sample due to a malformed SAM-to-CSV row; it is excluded from comparative numerical totals. All successful run intermediates passed a binary NUL-contamination audit.

## Reproduce one pair without exposing private inputs

You need authorized access to the original archived paired FASTQs and the corresponding fixed MiCall `unstitched_contigs_*.csv`, or to a specifically archived Ivanka `contigs.csv`. Raw inputs and the private alias mapping are **not distributed** with the repository.

Create an ignored local workspace and make both source trees available:

```sh
mkdir -p .validation-work/upstream-cba1898a \\
  .validation-work/candidate-2bc94b59
git archive cba1898a | tar -x -C .validation-work/upstream-cba1898a

# Reproduce using the exact candidate source commit, not an evolving branch.
git archive 2bc94b59366b0ff22a8a673bdb961c28ba75487a \
  | tar -x -C .validation-work/candidate-2bc94b59
```

Alternatively, use a detached Git worktree pinned to that commit. Install the project's Python dependencies and put the tested Bowtie2, Bowtie2-build and downstream aligner executables on `PATH`.

Use privately supplied values for `SOURCE_ROOT` and `SAMPLE_ID`; use aliases that do not disclose specimen identities.

```sh
SOURCE_ROOT=/path/to/private/archived/samples
SAMPLE_ID='private-sample-folder'
ALIAS='case-01'
COMMIT=2bc94b59366b0ff22a8a673bdb961c28ba75487a
RUNNER=docs/validation/real_hiv_hcv_20261106_runner.py

PYTHONPATH="$PWD/.validation-work/upstream-cba1898a" \
  .venv/bin/python "$RUNNER" --source-root "$SOURCE_ROOT" \
  --sample "$SAMPLE_ID" --alias "$ALIAS" --arm upstream \
  --output ".validation-work/ab/$ALIAS/upstream" \
  --candidate-commit "$COMMIT"

PYTHONPATH="$PWD/.validation-work/candidate-2bc94b59" \
  .venv/bin/python "$RUNNER" --source-root "$SOURCE_ROOT" \
  --sample "$SAMPLE_ID" --alias "$ALIAS" --arm nonstealing \
  --output ".validation-work/ab/$ALIAS/nonstealing" \
  --candidate-commit "$COMMIT"
```

For a preserved Ivanka contig panel rather than the original IVA contig CSV, pass `--contigs /path/to/the/archived/contigs.csv` to **both** invocations. Use the same archive and version-specific reference panel in both arms. All inputs are SHA-256 fingerprinted in ignored `provenance.json` files. The original compressed FASTQ is decompressed into an ignored cache using checksum-validated, atomic file staging.

The comparison for a full cohort requires an ignored `.validation-work/private-manifest.json` with entries of the form:

```json
[
  {
    "alias": "case-01",
    "sample": "PRIVATE_SAMPLE_FOLDER",
    "cohort": "HCV",
    "panel": "archived_IVA"
  }
]
```

Use `"panel": "archived_Ivanka"` and provide an additional `"contigs"` absolute path for such cases. Once the corresponding A/B outputs have finished:

```sh
.venv/bin/python docs/validation/real_hiv_hcv_20261106_analyze.py
.venv/bin/python docs/validation/real_hiv_hcv_20261106_details.py
.venv/bin/python docs/validation/real_hiv_hcv_20261106_report.py
```

The analysis and report generators create pseudonymous aggregate CSV, Markdown, SVG and HTML outputs. Do **not** commit `.validation-work/` or distribute the private analysis JSON files: they contain specimen correspondence and coordinate-level allele observations.

## Internal consistency and review

- Maximum depth is independently reconstructed from accepted Q15 merged-read records and compared with *every* candidate depth sidecar row and expected reference position.
- The reported unique mapped read-pair count is checked from original QNAME/R1/R2 records, independently of repeated alignment hypotheses.
- Gained, lost and discordant bases are compared on the same contig and reference-relative coordinate, accounting for consensus offsets.
- Independently repeated experiment runs, known-truth genotype controls, minority-haplotype spike-ins, infection-mixture controls, and complete MiCall end-to-end reanalysis are **not** part of this evidence and remain necessary before clinical release.
- Runtime measurements were made concurrently on a shared host with up to seven A/B workers; they are exploratory and **not** controlled performance benchmarks.
- The `-k 8` reporting cap can truncate equally compatible placements; maximum depth is an upper-bound hypothesis count, not independent physical molecule depth. Fractional depth is conditional on the returned placement set, not a calibrated posterior.

The reports should be reviewed as evidence of software behavior and its limitations, not as evidence that any downstream clinical diagnosis, genotype, or drug-resistance call is correct.
