# Independent genotype and minority-haplotype: release verdict

**NOT PRODUCTION READY**. The PR has credible genotype-correctness benefits but important unresolved safety, completeness and runtime issues.

The five-strain HIV laboratory mixture produced 798 additional correct callable bases and zero wrong calls versus independently sequenced references. After 50 deliberate wrong-seed mutations, it corrected 47/50 wrong bases versus 46/50 upstream, and gained 900 correct calls.

The 8-placement cap missed best-scoring reference alternatives for 78/1,200 real nine-clone read pairs compared with a deeper alignment search. Even a 16-placement window did not remove every difference.

At 0.99% real-source HCV minor-library admixture, both mappers retained 428/500 true-source minor reads on the minor contig; but an intentionally absent third-genotype decoy had callable bases in both mappers. MAX depth is not proof of independent molecular infection evidence.

A real-data read-accounting regression was fixed: preserve genuine single-mapped-mate primary SAM pairs and count one mapped read end, without manufacturing Q15 consensus support.

The revised small-panel all-placement candidate scored 40,772 correct HIV5 bases with zero wrong, in 200 seconds (vs 448 seconds previously). HIV9 full-run base accuracy also improved by 896 correct bases, but nominal rare-clone fractions conflict with direct unaligned-read motif counts. Full-pipeline clinical reporting and calibrated haplotype abundance remain unverified.

See independent_truth_20261107_scientific.md for complete data, the nine-haplotype status, methods, and gate matrix.
