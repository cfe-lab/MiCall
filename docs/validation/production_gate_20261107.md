# Production readiness assessment — MiCall de novo non-stealing mapper

**Decision: NO-GO for clinical production as of 7 November 2026.**

This decision is based on executed control experiments, not the presence or absence of unit tests alone. The proposed mapper is a promising research implementation, and it performs better against several independently verified HIV reference panels. It is not yet demonstrated safe for clinical genotype-presence and rare-haplotype reporting.

## Evidence and release gates

| Gate | Current evidence | Gate |
|---|---|---|
| Basic MiCall regression | 4,362 passing tests, 5 skipped, 54 expected failures; IVA/installation prerequisites excluded | Passed in tested configuration |
| Independently verified HIV5 base truth | Revised four-thread all-alignment candidate: 40,772 correct, zero incorrect across 48,381 truth bases; 200s vs upstream 174s | Favorable, limited scope |
| Deliberately wrong reference challenge | Earlier bounded candidate corrected 47/50 planted erroneous seeds, versus 46/50 upstream; zero wrong callable bases | Favorable, must repeat on final candidate |
| Independently verified HIV9 base truth | Earlier bounded candidate: 8,647 correct and 14 wrong vs upstream 7,751 correct and 15 wrong | Favorable, final adaptive candidate pending rerun |
| Complete competing-reference discovery | 1,200 genuine HIV9 pairs: eight-hit cap missed best-reference alternatives for 78 groups versus deeper search. The new all-placement path removes the explicit eight-hit restriction for panels up to 16 contigs, but Bowtie2 remains heuristic. Panels above 16 remain capped. | **Fail / not established** |
| Read-pair completeness and integrity | Primary singleton-mate read regression repaired; remapping counts now distinguish legacy initial-unmapped from final unmapped-after-selection; capped and oversized groups reported | Passed for tested cases; larger panels require safety validation |
| Minority-haplotype identity and quantitative abundance | Nine verified reference sequences and linked private 21-mers support clone-specific molecules, but observed raw-source motif frequencies contradict the nominal 0.2% and 0.4% clone labels. This experiment cannot verify sensitivity at those abundances. | **Fail / unresolved truth** |
| Absent-genotype specificity | In real HCV read-source challenge, a third genotype with zero source-library reads nevertheless accumulated 26 (zero-minority) or 32 (1% minority) Q15 MAX callable bases, with all compared mappers | **Fail / not established** |
| Runtime under production constraints | HIV5 now 200s versus 174s upstream; previous HIV9 3,243s and all-placement HIV9 full rerun still needs evaluation; 16-contig HCV stress currently underway | **Open** |
| Complete MiCall diagnostic output | Controlled studies cover real FASTQ → map_to_contigs → sam2aln → aln2counts, **not** a complete de novo assembly, trimming, G2P, clinical report, and reproducible drug resistance interpretation | **Open** |
| Known biological truth across intended HIV/HCV panel | Independent HIV5/HIV9 nucleotide references, but no independent final true HCV genotype or abundance adjudication for every candidate allele | **Open** |

## Exact scientific interpretation

**Maximum depth is a per-contig compatibility upper bound.** Retaining a physical read pair across several contigs does not create new molecules, cannot itself establish infection by multiple genotypes, and must not be used as an independent evidence count for a minority haplotype.

**Fractional average depth is not a posterior genotype probability.** It assigns one unit per original fragment across the **reported** candidate alignments. Truncation of candidates, alignment heuristic behavior, and failure to see an alternative still bias the denominator.

**Linked exact private markers demonstrate local sequence linkage, not calibrated 0.2% mixture sensitivity.** In 239,006 unaligned real read pairs, the putative 0.2%-labelled clone had at least one clone-private exact 21-mer in 96,433 pairs, and at least two in 35,122 pairs. This is incompatible with interpreting the provided nominal label as a verified 0.2% frequency in the sequenced library.

**End-to-end clinical sign-off is separate from technical integration.** Passing unit tests and showing more true callable bases on reference-matched controls is insufficient evidence that a clinical MiCall run preserves diagnostic specificity, drug resistance interpretation, haplotype phase, and safe computational limits.

## Conditions required to change the decision

1. Establish independently mixed/quantified HIV/HCV controls with true genotype presence/absence, clone frequencies and known phase; predefine sensitivity and false-positive thresholds.
2. Demonstrate the final adaptive mapper on those controls, including complex panels above 16 contigs, with an explicit, tested treatment of search-cap saturation and oversized read groups.
3. Run complete MiCall clinical de novo workflows (including assembler, stitching, G2P, trimming and report generation) on the same controls and confirm output equivalence/correctness.
4. Repeat performance and peak memory tests on large real cohorts under the actual per-sample production time and parallelism budget, not just shared-host pilot timings.
5. Obtain an independent scientific/code review and preserve provenance of the exact merged candidate commit.

Related evidence is in [the independent truth report](independent_truth_20261107_scientific.md) and [the earlier real HIV/HCV cohort report](real_hiv_hcv_20261106_scientific.md). The private raw controls and alignment artifacts remain in Git-ignored working directories.
