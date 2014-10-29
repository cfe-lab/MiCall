# MiseqPipeline #
Pipeline for processing FASTQ data from an Illumina MiSeq.
Monitoring system regularly checks network attached storage for unprocessed runs, transfers FASTQ.gz files to the cluster and executes the pipeline.

## Steps and their input / output files ##

* prelim_map
  * in - fastq1
  * in - fastq2
  * prelim.csv - SAM file format
* remap
  * in - fastq1
  * in - fastq2
  * in - prelim.csv
  * remap.csv - SAM file format.
  * remap_counts.csv - collated as collated_counts.csv
  * remap_conseq.csv - consensus sequence
  * unmapped1.fastq - FASTQ format (unstructured text) &rarr; results/unmapped
  * unmapped2.fastq - FASTQ &rarr; results/unmapped
* sam2aln
  * in - remap.csv
  * aligned.csv - reads aligned to consensus sequence
  * conseq_ins.csv - collated - insertions relative to consensus sequence
  * failed_read.csv - collated - reads that fail to merge
* aln2counts
  * in - aligned.csv
  * nuc.csv - collated as nucleotide_frequencies.csv - nucleotide counts
  * amino.csv - collated as amino_frequencies.csv - amino counts
  * coord_ins.csv - collated - insertions relative to coordinate reference
  * conseq.csv - collated as collated_conseqs.csv - consensus sequence
  * failed_align.csv - collated
* aln2nuc
  * in - aligned.csv
  * nuc_variants.csv - collated - top nucleotide variants for HLA
* sam_g2p
  * in - remap.csv
  * g2p.csv - collated
* coverage_plots
  * in - amino.csv
  * in - nuc.csv
  * coverage_maps.tar - binary file &rarr; untar in results/coverage_maps
  * coverage_scores.csv - collated
