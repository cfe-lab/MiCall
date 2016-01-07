# MiCall #
## Processing FASTQ data from an Illumina MiSeq ##
Maps all the reads from a sample against a set of reference sequences, then
stitches all the reads into consensus sequences and coverage maps.

A monitoring system regularly checks the file system for unprocessed runs,
transfers FASTQ.gz files to the cluster and executes the pipeline.

## Dual Licensing ##
Copyright (C) 2016, University of British Columbia

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as published
by the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, visit [gnu.org][gnu]. The source code for
this program is available from [github.com][github].

The program is also available for a fee under a more permissive license. For
example, if you want to run a changed version of the program on a network server
without publishing the changed source code, [contact us][contact] about
purchasing a license.

[gnu]: http://www.gnu.org/licenses/
[github]: https://github.com/cfe-lab/MiCall
[contact]: http://www.google.com/recaptcha/mailhide/d?k=01yIEKHPqcNIG1ecF-Khum4g==&c=SnES_swjOvb0d22WN4Q30vx9gKyzHNDkTxbY6Y1os4w=

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
