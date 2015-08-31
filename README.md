# MiCall #
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
  * failed_align.csv - collated - any consensus that failed to align to its ref
  * nuc_variants.csv - collated - top nucleotide variants for HLA
* sam_g2p
  * in - remap.csv
  * g2p.csv - collated
* coverage_plots
  * in - amino.csv
  * in - nuc.csv
  * coverage_maps.tar - binary file &rarr; untar in results/coverage_maps
  * coverage_scores.csv - collated

## File descriptions ##
* conseq.csv
  * consensus-percent-cutoff - to be included in a mixture, a variant must make
    up at least this fraction of the total valid counts
  * offset - using the seed reference's coordinate system, this is the 1-based
    index of the first character in sequence that is not a dash. For example, if
    the seed was `ACTAGTCC` and the consensus sequence is `-AGTC`, then the
    offset would be 4.
    ** change to ** the number of dashes that are not shown at the start of the
    sequence
  * sequence - the consensus sequence, aligned to the codon reading frame.
    ** once offset is changed, add the following **
    Adding the number of dashes in offset will align the first base in the
    consensus sequence with its corresponding base in the seed reference. The
    whole consensus sequence may not be aligned with the seed reference because
    of insertions and deletions
* nuc.csv
  * query.nuc.pos - the 1-based index of the base in the consensus sequence that
    came from this set of counts
