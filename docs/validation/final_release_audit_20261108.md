# Final-candidate independent genotype and minority-haplotype audit

**Clinical production decision: NO-GO; release requirements remain unmet.**

## Real HIV reference truth — identical inputs in all arms

| Condition | Algorithm | Correct Q15 MAX | Wrong Q15 MAX | Uncallable | Total runtime |
|:--|:--|--:|--:|--:|--:|
| HIV5-truth | upstream | 39,944 | 0 | 8,437 | 173.9s |
| HIV5-truth | bounded_candidate | 40,742 | 0 | 7,639 | 447.7s |
| HIV5-truth | final_candidate | 40,772 | 0 | 7,609 | 199.8s |
| HIV9-truth | upstream | 7,751 | 15 | 1,900 | 1808.8s |
| HIV9-truth | bounded_candidate | 8,647 | 14 | 1,005 | 3243.1s |
| HIV9-truth | final_candidate | 8,674 | 15 | 977 | 1776.8s |
| HIV5-perturbed | upstream | 39,929 | 0 | 8,452 | 176.1s |
| HIV5-perturbed | bounded_candidate | 40,829 | 0 | 7,552 | 464.2s |
| HIV5-perturbed | final_candidate | 40,816 | 0 | 7,565 | 162.8s |

### Reference-error correction challenge

| Implementation | Planted errors corrected | Wrong seed retained | Other wrong | Uncallable |
|:--|--:|--:|--:|--:|
| upstream | 46/50 | 0 | 0 | 4 |
| bounded_candidate | 47/50 | 0 | 0 | 3 |
| final_candidate | 47/50 | 0 | 0 | 3 |

The five-strain truth FASTA derives from individually sequenced laboratory isolates; the nine-haplotype truth FASTA derives from separately Sanger-confirmed plasmids. However, the known-truth references are used directly as *mapping seeds* for the truth panel. The deliberately corrupted-reference panel is a partial control against circular seed agreement.

Scientific sources: the five-virus control references are described by Di Giallonardo and colleagues (2014), https://doi.org/10.1093/nar/gku537; the nine-clone Sanger-confirmed reference mixture by Knyazev and colleagues (2021), https://doi.org/10.1093/nar/gkab576. The original published nominal clone percentages have not been confirmed as measured frequencies in the selected read set.

## HCV real-read source-origin and absent-genotype negative controls

| Condition | Mapper | Source HCV2b read pairs mapped correctly | Major-library pairs on absent HCV1b decoy | Decoy callable positions |
|:--|:--|--:|--:|--:|
| HCV-origin-zero | upstream | 0 | 225 | 26 |
| HCV-origin-zero | nonstealing | 0 | 230 | 26 |
| HCV-origin-zero | nonstealing-all-p4 | 0 | 230 | 26 |
| HCV-origin-1pct | upstream | 428 | 225 | 32 |
| HCV-origin-1pct | nonstealing | 428 | 230 | 32 |
| HCV-origin-1pct | nonstealing-all-p4 | 428 | 230 | 32 |

Fractional-depth specificity check: in the final candidate, all 26 Q15 MAX callable positions on the absent HCV1b decoy in the zero-minority control, and all 32 in the 0.99% mixture, reached fractional depth 100 as well as maximum depth 100. Fractional depth alone therefore does not establish genotype absence.

These controls mix real reads from identified source *libraries*. The identity of the original read library is known, but the biological genotype of each library was not independently adjudicated. Apparent calls on an intentionally absent contig are a warning about specificity, not proof of a clinically reportable false-positive result. These apparent calls already occurred with upstream MiCall.

## Real-read source dilutions below one percent

Additional same-input A/B controls resampled authentic archived HCV reads at verified source-library fractions 0%, 0.2%, 0.4% and 0.99%. Both independently seeded 0.2% resamples yielded 88/100 correctly assigned source pairs but no callable minority-contig bases. Two 0.4% resamples produced 1 and 22 callable bases despite 176/200 and 178/200 correct source assignments. The absent-source decoy retained 26 callable bases. These are not certified clinical genotype fractions. See the dilution report real_hcv_dilution_20261108.md.

## Large-panel bounded-alignment completeness: real-read stress

The final >16-contig fallback was challenged with 1,200 real paired HIV9 reads and a constructed panel of 17 authentic reference sequences (nine Sanger-confirmed HIV9 plasmids, five HIV5 isolates and three archived HCV contigs). 1,167 paired-read groups reached the -k 8 reporting limit. 188 read pairs omitted at least one reference found by Bowtie2 all-alignment mode, with 290 missing alternatives in total. This is a deeper heuristic search, not a completeness oracle. The reference panel is an engineered validation panel of real contigs, not an actual 17-contig patient assembly. Runtime was 4.1s capped versus 10.0s deeper. This prevents endorsing the large-panel fallback for clinical use without an adaptive, bounded correctness strategy.

## Alignment-independent phased-minority signals

Private 21-mers around each Sanger-confirmed, clone-private SNP were searched against the **unaligned** original paired FASTQ reads, including reverse-complement matches. Two private motifs occurring in one pair demonstrate local linkage, not a full-length phased haplotype or calibrated clone abundance.

| Clone | Private 21-mers | ≥1 per pair | ≥2 linked per pair | ≥3 linked per pair |
|:--|--:|--:|--:|--:|
| A_50 | 0 | 0 | 0 | 0 |
| B_25 | 2 | 274 | 0 | 0 |
| C_12.5 | 0 | 0 | 0 | 0 |
| D_6.3 | 2 | 342 | 0 | 0 |
| E_3.2 | 2 | 117 | 0 | 0 |
| F_1.6 | 6 | 189 | 1 | 0 |
| G_0.8 | 4 | 142 | 3 | 1 |
| H_0.4 | 3 | 39,057 | 16,180 | 4,461 |
| J_0.2 | 3 | 96,433 | 35,122 | 6,536 |

The clone named **J_0.2** has 96,433 pairs matching at least one private motif and 35,122 with two or more linked motifs in 239,006 sampled paired reads. It is therefore scientifically unjustified to claim demonstrated 0.2% sensitivity from this particular sequenced read set merely because the reference FASTA carries a 0.2% nominal label. Composition, amplification, sequence error and library provenance need separate adjudication.

## Production gate and correctness limitations

- **PASS in tested configuration:** 4,373 broader MiCall tests passed, with five skips, one deselection, and 54 expected failures. Installation/IVA tests depend on unavailable prerequisites.
- **PASS on reference subsets:** the revised candidate recovers more correct callable HIV bases with no overall increase in wrong calls against upstream on the fixed-known-truth panels; some clones and positions still have incorrect calls.
- **NOT ESTABLISHED:** known negative-control genotype-absence specificity, minority mixture abundance and full haplotype phase.
- **NOT ESTABLISHED:** complete detection of competing alignments on Bowtie2 heuristic search; panels over 16 still use a k=8 cap.
- **NOT ESTABLISHED:** end-to-end production MiCall clinical outputs (assembly, stitcher, trimming, G2P, HIV resistance reports), real 198-sample throughput and peak memory budgets.
- **NOT ESTABLISHED:** independently adjudicated full biological ground truth for all real HCV libraries, despite their verified source-library identities.
- **Fixed engineering diagnostic:** legitimate discordant Bowtie2 SAM groups with unequal mate-alternative counts no longer trigger false malformed-pair warnings; rescue behavior is unit-tested.

The clinical production release remains blocked. An increase in per-contig maximum depth is not independent molecule evidence. Never use the alternative-depth sidecar as a calibrated genotype frequency or clinical diagnosis by itself.

### Provenance

Upstream source: cba1898a. Bounded candidate controls were run on earlier revisions. All-placement, four-thread controls used the working tree implementing the same alignment mode as local commit 379e80e9; exact commit-pinned reruns and independent review remain mandatory before final release. Original reads, haplotype sequences and private run artifacts remain in git-ignored local storage.

Generated 2026-11-08 06:32 UTC
