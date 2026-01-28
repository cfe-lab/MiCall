---
title: Steps and their input / output files
subtitle: Where the data goes
---

Individual files are described after the list of steps.

* `micall_basespace.summarize_run`: extract phiX error rates from MiSeq data.
  * in - InterOp folder
  * in - read lengths from RunInfo.xml or BaseSpace run parameters
  * quality.csv - QC error rate data, grouped by tile
* `filter_quality`: Summarize phiX quality data, and decide which tiles and
    cycles to censor.
  * in - quality.csv
  * bad_cycles.csv - List of tiles and cycles rejected for poor quality
* `trim_fastqs`: Censor reads based on phiX quality data, and also trim adapter
    sequences.
  * in - original1.fastq
  * in - original2.fastq
  * in - bad_cycles.csv - List of tiles and cycles rejected for poor quality')
  * trimmed1.fastq - uncompressed FASTQ containing trimmed reads (read 1)
  * trimmed2.fastq - uncompressed FASTQ containing trimmed reads (read 2)
* `fastq_g2p`: use pairwise alignment to map reads to the V3LOOP region,
    then use the geno2pheno algorithm to
    translate genotype to phenotype and predict whether a sample is X4 or R5.
    HIV-1 particles use coreceptors to enter cells. Different particles can use
    CXCR4, CCR5, or both. This analysis looks at the sequences from a virus
    population and reports either "X4" (able to use CXCR4) or "R5" (only able
    to use CCR5). See the [design page][fastq_g2p_design] for more details.
  * in - fastq1
  * in - fastq2
  * g2p.csv - downloaded - calls each individual sequence X4, R5, or error.
  * g2p_summary.csv - downloaded - calls the entire sample X4 or R5.
  * g2p_aligned.csv - reads that mapped to V3LOOP, aligned to the HIV seed
  * not_v3_r1.fastq - reads that did not map to V3LOOP (read 1)
  * not_v3_r2.fastq - reads that did not map to V3LOOP (read 2)
  * merged_contigs.csv - contig sequences generated from amplicon reads that
    show up as very common insert lengths after merging read pairs
* `prelim_map`: map reads to all references. (Takes reads that did not map to V3LOOP.)
  * in - fastq1
  * in - fastq2
  * prelim.csv - SAM file format
* `denovo`: assemble contigs from individual reads. (Replaces `prelim_map` in
    the denovo version of MiCall.)
  * in - fastq1
  * in - fastq2
  * in - merged_contigs.csv
  * unstitched_contigs.csv - the assembled contigs, plus any merged contigs, including
    the best blast results
  * contigs.csv - stitched version of `unstitched_contigs`
  * blast.csv - multiple blast results for each contig
* `remap`: iteratively use consensus from previous mapping as reference to try
    and map more reads. See [remap design] for more details. (The denovo version
    of MiCall just maps once onto the assembled contigs.)
  * in - fastq1
  * in - fastq2
  * in - prelim.csv
  * remap.csv - SAM file format with commas instead of tabs.
  * remap_counts.csv - downloaded - how many reads mapped to each reference at
    each stage.
  * remap_conseq.csv - downloaded - consensus sequence that reads were mapped to
    on the final iteration
  * unstitched_conseq.csv - downloaded - consensus sequence that reads were
    mapped to the unstitched contigs.
  * unmapped1.fastq - FASTQ format (unstructured text) reads that didn't map to
    any of the final references.
  * unmapped2.fastq - FASTQ 
* `sam2aln`: extract the aligned portion of each read from a CSV SAM file. Also
    combines duplicates.
  * in - remap.csv
  * aligned.csv - reads aligned to consensus sequence
  * conseq_ins.csv - insertions relative to consensus sequence
  * failed_read.csv - downloaded - reads that fail to merge
  * clipping.csv - count of soft-clipped reads at each position
* `aln2counts`: take the aligned reads, and count how many of each nucleotide
    or amino acid appear at each position.
  * in - aligned.csv
  * in - g2p_aligned.csv
  * in - clipping.csv
  * in - conseq_ins.csv
  * nuc.csv - downloaded - nucleotide counts at each position
  * amino.csv - downloaded - amino counts at each position
  * concordance.csv - downloaded - concordance measures for all regions, 
    relative to the coordinate reference, cumulative
  * concordance_seed.csv - downloaded - concordance measures for all regions, 
    relative to the seed reference, for each contig
  * conseq.csv - downloaded - consensus sequence, minimum coverage (read depth) 100
  * conseq_all.csv - downloaded - consensus for all seeds (whole and region-wise), no minimum coverage
  * conseq_region.csv - downloaded - consensus for all regions
  * conseq_stitched.csv - downloaded - stitched whole genome consensus
  * insertions.csv - downloaded - insertions for all regions, at different mixture cutoffs
  * failed_align.csv - downloaded - any consensus that failed to align to its ref
  * nuc_variants.csv - downloaded - top nucleotide variants for HLA
  * nuc_detail.csv - nucleotide counts split up by contig
  * amino_detail.csv - amino counts split up by contig
  * genome_coverage.csv - coverage counts in full-genome coordinates
  * concordance_detailed.csv - coordinate concordance for each contig, averaged over a window size of 20,
    for each genome position
* `coverage_plots`: convert amino counts to a coverage graph for each gene
    region.
  * in - amino.csv
  * coverage_maps.tar - binary file &rarr; untar in results/coverage_maps
  * coverage_scores.csv - downloaded - a score for each region based on the
    coverage at key positions.
* `plot_genome_coverage`: convert [genome coverage] to a graph
  * in - genome_coverage.csv
  * genome_coverage.svg
  * genome_concordance.svg, if `use_concordance` is set
* `concordance_plot`: convert concordance counts to a concordance graph for each gene
    region.
  * in - concordance_detailed.csv
  * output in coverage_maps.tar as well
* `cascade_report`: summarize how many reads made it through each step of the
    pipeline.
* `resistance`: Use the Stanford HIVdb algorithm rules to call a sample's resistance
    levels to HIV drugs, or other rules for Hepatitis C drugs. See
    [resistance design]
  * in - amino.csv
  * in - amino_midi.csv - amino counts from the HCV-MIDI sample (ignored if
    it's the same file as amino.csv)
  * resistance.csv - resistance level to each HIV drug
  * mutations.csv - significant mutations found
  * resistance_fail.csv - failure reasons
* `genreport`: Generate the resistance report
  * in - resistance.csv
  * in - mutations.csv
  * report.pdf - the report, ready to print

[fastq_g2p_design]: http://cfe-lab.github.io/MiCall/design/fastq_g2p
[remap design]: http://cfe-lab.github.io/MiCall/design/remap
[resistance design]: http://cfe-lab.github.io/MiCall/design/resistance
[genome coverage]: genome_coverage.md

## File descriptions ##

* aligned.csv and g2p_aligned.csv
  * refname - seed reference the reads mapped to
  * qcut - minimum Phred quality score to include a nucleotide
  * rank - sequences are sorted by count, with the most common at rank 0
  * count - how many times exact copies of this sequence and soft clip counts
    were found
  * offset - offset of the read within the consensus, or number of dashes to
    add at the start
  * seq - the mapped sequence of the read, aligned to the consensus
* amino.csv
  * seed - seed reference the reads mapped to
  * region - coordinate reference for reporting against, usually a gene
  * q-cutoff - minimum Phred quality score to include a nucleotide
  * query.aa.pos - the 1-based index of the amino acid in the consensus sequence
    that came from this set of counts
  * refseq.aa.pos - the 1-based index of the amino acid in the coordinate reference
  * A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y - counts for the amino acids at this
    position
  * `*` - count of stop codons at this position
  * X - count of codons with at least one nucleotide with a Phred quality score
    below the cutoff
  * partial - count of codons with a partial deletion within a read (partial
    codons at either end of a read are not counted at all)
  * del - count of codons with a deletion relative to the coordinate reference
  * ins - count of codons with an insertion in or after them relative to the
    coordinate reference
  * clip - count of reads with soft clipping that would have mapped at this
    codon
  * g2p_overlap - count of reads at this position that overlapped with the
    V3LOOP G2P mapping and were ignored.
* amino_detail.csv - same fields as amino.csv
* blast.csv
  * contig_name - name that appears in amino_detail.csv and nuc_detail.csv
  * ref_name - reference name that BLAST found matches in
  * score - calculated by BLAST
  * match - fraction of the sequence that matched the reference, negative for
    reverse-complemented matches
  * pident - percentage of the matching part that was identical to the reference
  * start - start position of the matching part in the contig sequence
  * end - end position of the matching part in the contig sequence
  * ref_start - start position of the matching part in the reference sequence
  * ref_end - end position of the matching part in the reference sequence  
* cascade.csv - number of read pairs that flow through the pipeline steps
  * demultiplexed - count from the raw FASTQ
  * v3loop - aligned with V3LOOP
  * g2p - valid reads to count in G2P
  * prelim_map - mapped to other references on first pass
  * remap - mapped to other references after remapping
  * aligned - aligned with a reference and merged with mate
* conseq.csv
  * region - the name of the contig. Includes the name of the reference seed, plus an optional prefix, which is a number that makes the name unique.
  * q-cutoff - minimum quality score
  * consensus-percent-cutoff - to be included in a mixture, a variant must make
    up at least this fraction of the total valid counts
  * offset - using the seed reference's coordinate system, this is the 1-based
    index of the first character in the consensus sequence.
  * sequence - the consensus sequence, aligned to the codon reading frame.
    Adding the number of dashes in offset will align the first base in the
    consensus sequence with its corresponding base in the seed reference. The
    whole consensus sequence may not be aligned with the seed reference because
    of insertions and deletions. Mixtures above the cutoff are displayed as
    [IUPAC nucleotide codes]. If deletions are present above the cutoff, the
    nucleotide code is lower case. If coverage is below 100 or if the only thing
    present above the cutoff is deletions, the nucleotide code is "x".
* conseq_all.csv: same columns as conseq.csv, with two different offsets:
  * seed-offset - offset of the seed relative to the reference
  * region-offset - offset of the seed relative to the region
* conseq_region.csv: same columns as conseq.csv
* conseq_stitched.csv: our best guess at a whole genome consensus, stitched
  together from the individual regions' consensus sequences, at different mixture levels.
  Insertions are included and gaps are not filled. In areas where multiple 
  regions overlap, we prioritize translated regions (because the amino acid alignment is 
  more precise), and we prioritize information from the centre of the region over information 
  from the edge of the region (again because the alignment in the centre is more reliable).
  * seed - reference name
  * q-cutoff - minimum quality score
  * consensus-percent-cutoff - to be included in a mixture, a variant must make
    up at least this fraction of the total valid counts
  * offset - offset relative to the reference
  * sequence - stitched whole genome consensus
* conseq_ins.csv
  * qname - query name from the read
  * fwd_rev - F for forward reads, R for reverse
  * refname - seed reference the read mapped to
  * pos - 1-based position in the consensus sequence that this insertion follows
  * insert - the nucleotide sequence that was inserted
  * qual - the Phred quality scores for the inserted sequence
* unstitched_contigs.csv
  * ref - the reference name with the best BLAST result
  * match - the fraction of the contig that matched in BLAST, negative for
    reverse-complemented matches
  * group_ref - the reference name chosen to best match all of
    the contigs in a sample
  * contig - the nucleotide sequence of the assembled contig
* contigs.csv
  Same as `unstitched_contigs.csv`, but contigs are stitched by `micall/core/contig_stitcher.py`.
* stitcher_plot_svg
  An SVG visualization showing how contigs were stitched together, including coverage information
  and operations performed during stitching. Only generated for denovo assembly with the referencefull
  stitcher.
* coverage_scores.csv
  * project - the project this score is defined by
  * region - the region being displayed
  * seed - the seed mapped to
  * q.cut - quality score cutoff
  * min.coverage - minimum coverage at all key positions
  * which.key.pos - which key position had the minimum
  * off.score - score to use when off target
  * on.score - score to use when on target
* clipping.csv
  * refname - seed reference the reads mapped to
  * pos - one-based nucleotide position within the consensus sequence
  * count - number of read pairs that soft clipped instead of mapping to this
    position
* failed.csv
  * qname - query name from the read
  * cause - reason the reads failed to merge
* failed_align.csv -
  * seed - seed the reads aligned to
  * region - where the consensus was trying to align
  * qcut - the quality score cutoff
  * queryseq - the consensus sequence
  * refseq - the reference sequence for the region
* g2p.csv - calls individual sequences as either R5 or X4
  * rank - ranks each sequence within the sample, starting at 1 for the most
    common.
  * count - how many times this sequence appeared in the sample. (Note that the
    same amino acid sequence can be produced by different nucleotide sequences,
    and these are counted separately.)
  * g2p - the score from the geno2pheno algorithm
  * fpr - the false positive rate
  * call - blank if there was an error, R5 if fpr > 3.5, otherwise X4.
  * seq - the amino acid sequence that was used
  * aligned - how the amino acid sequence aligned with the reference sequence
  * error - reason for rejecting a sequence
    * `low quality` - more than half of the sequence had read quality with phred
        score < 33
    * `notdiv3` - the sequence did not end on a codon boundary
    * `zerolength` - none of the sequence mapped to the expected portion of the
        reference sequence
    * `cysteines` - the start or end of the amino acid sequence was not a cycsteine
    * `> 2 ambiguous` - more than one amino acid was ambiguous, or one position
        had more than two possible amino acids
    * `stop codons` - the sequence contained a stop codon
    * `length` - the amino acid sequence length was outside of the range 32 to 40
    * `failed to align` - the nucleotide sequence failed to align with the reference
    * `count < 3` - this is a summary line reporting how many reads were rejected
        for low counts
  * comment - a warning for sequences that had a problem but were good enough
    to score
    * `ambiguous` - one amino acid had two possible values
* g2p_summary.csv - summarizes the sample and calls it either R5 or X4
  * mapped - the number of reads that mapped to V3LOOP
  * valid - the number of those reads that were valid sequences and not errors
  * X4calls - the number of those reads that were called X4
  * X4pct - X4calls as a percentage of valid
  * final - the final decision: blank if valid is 0, X4 if X4pct >= 2, otherwise
    R5
* genome_coverage.csv
  * contig - the name of the contig that appears in amino_detail.csv and
    nuc_detail.csv, or the name of the seed reference used for mapping
  * coordinates - the name of the coordinate reference the query was aligned to
  * query_nuc_pos - the one-based nucleotide position within the contig or the remap
    consensus
  * refseq_nuc_pos - the one-based nucleotide position within the reference, used for display - please note that this 
  is *NOT* the reference position that corresponds to the query positions in the alignment.
  * dels - number of deletions reported at this position
  * coverage - number of reads that aligned to this position, including
    deletions
  * link - type of link between the contig position and the reference position: 'M' for a match, 'D' for a deletion, 
  'I' for an insertion, and 'U' for unknown (a section of the contig that didn't map to the reference)
* insertions.csv
  * seed - the name of the contig
  * mixture_cutoff - to be included in a mixture, a variant must make
    up at least this fraction of the total valid counts
  * region - region where the insertion was found
  * ref_region_pos - 1-based position of the nucleotide before the insertion, 
    relative to the reference region
  * ref_genome_pos - 1-based position of the nucleotide before the insertion,
    relative to the full reference genome
  * query_pos - 1-based position of the nucleotide before the insertion,
    relative to the contig
  * insertion - sequence of the insertion (as nucleotides), computed in the same
    way as for conseq.csv (e.g. mixture cutoffs, deletions and low coverage).
    If there was not an insertion in the reads at this position, this is counted
    as a deletion relative to the insertion.
* merged_contigs.csv
  * contig - the consensus sequence from all merged reads of that insert length
* mutations.csv
  * drug_class - the drug class code from the HIVdb rules, like NRTI
  * mutation - the wild-type amino, position, and resistant amino, like Q80K
  * prevalence - the fraction of coverage that contained this mutation
  * genotype - the HCV genotype
* nuc.csv
  * seed - seed reference or contig the reads mapped to
  * region - coordinate reference for reporting against, usually a gene
  * q-cutoff - minimum Phred quality score to include a nucleotide
  * query.nuc.pos - the 1-based index of the base in the consensus sequence that
    came from this set of counts, stretching across the entire seed reference or
    contig that the reads were mapped to
  * refseq.nuc.pos - the 1-based index of the base in the coordinate reference
    for a single gene region on non-coding region, named in the region column
  * genome.pos - the absolute 1-based index of a nucleotide (in other words,
    relative to the beginning of the entire reference sequence). This is the same
    as refseq.nuc.pos, but not specific to a single gene or
    non-coding region.
  * A,C,G,T - counts for the nucleotides at this position
  * N - count of reads with Phred quality score below the cutoff
  * del - count of reads with a deletion at this position
  * ins - count of reads with an insertion after this position, relative to the
    coordinate reference
  * clip - count of reads with soft clipping that would have mapped at this
    position
  * g2p_overlap - count of reads at this position that overlapped with the
    V3LOOP G2P mapping and were ignored.
* nuc_detail.csv - same fields as nuc.csv
* quality.csv and bad_cycles.csv
  * tile
  * cycle
  * errorrate - as a percentage
* remap_conseq.csv
  * region - the region mapped to
  * sequence - the consensus sequence used
* unstitched_conseq.csv
  * region - the region mapped to
  * sequence - the consensus sequence used
* unstitched_cascade.csv - number of read pairs that flow through the pipeline steps
  * demultiplexed - count from the raw FASTQ
  * v3loop - aligned with V3LOOP
  * g2p - valid reads to count in G2P
  * prelim_map - mapped to other references on first pass
  * remap - mapped to other references after remapping
  * aligned - aligned with a reference and merged with mate
* resistance.csv
  * region - the region code, like PR or RT
  * drug_class - the drug class code from the HIVdb rules, like NRTI
  * drug - the drug code from the HIVdb rules, like ABC
  * drug_name - the drug name from the HIVdb rules, like abacavir
  * level - resistance level as a number
  * level_name - resistance level's name, like "Susceptible"
  * score - numeric score calculated by the HIVdb rules
  * genotype - string describing the HCV genotype
* run_quality.csv
  * q30_fwd - portion of tiles and cycles with quality score of at least 30
  for forward reads
  * q30_rev - portion of tiles and cycles with quality score of at least 30
  for reverse reads
  * cluster_density - average cluster density for all tiles and cycles K/mm2
  * pass_rate - portion of clusters passing filters over all tiles and cycles
  * error_rate_fwd - average error rate over tiles and cycles in forward reads
  (phiX error count/tile/cycle)
  * error_rate_rev - average error rate over tiles and cycles in reverse reads
  (phiX error count/tile/cycle)
  * avg_quality - average Phred score over all clusters and cycles
  * avg_coverage - average coverage across the best region for each sample

[IUPAC nucleotide codes]: https://www.bioinformatics.org/sms/iupac.html
