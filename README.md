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

* `micall_basespace.summarize_run`: choose tiles and cycles to censor for high
    error rates.
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
    to use CCR5).
  * in - fastq1
  * in - fastq2
  * g2p.csv - downloaded - calls each individual sequence X4, R5, or error.
  * g2p_summary.csv - downloaded - calls the entire sample X4 or R5.
  * g2p_aligned.csv - reads that mapped to V3LOOP, aligned to the HIV seed
  * not_v3_r1.fastq - reads that did not map to V3LOOP (read 1)
  * not_v3_r2.fastq - reads that did not map to V3LOOP (read 2)
* `prelim_map`: map reads to all references. (Takes reads that did not map to V3LOOP.)
  * in - fastq1
  * in - fastq2
  * prelim.csv - SAM file format
* `remap`: iteratively use consensus from previous mapping as reference to try
    and map more reads.
  * in - fastq1
  * in - fastq2
  * in - prelim.csv
  * remap.csv - SAM file format with commas instead of tabs.
  * remap_counts.csv - downloaded - how many reads mapped to each reference at
    each stage.
  * remap_conseq.csv - downloaded - consensus sequence
  * unmapped1.fastq - FASTQ format (unstructured text) reads that didn't map to
    any of the final references.
  * unmapped2.fastq - FASTQ 
  * Choose which seeds to use as follows:
    1. Count all the reads that the prelim_map step mapped to all the seeds,
      skipping any shorter than 50 bases, because they are often primers.
    2. Look at how many reads mapped to each seed, and drop any that mapped
      fewer than ten.
    3. If more than one seed from the same seed group (HCV-1a, 1b, 2, 3, 4, etc.)
      pass the threshold, select the one with the most reads. In case of a tie,
      choose the first seed, alphabetically.
    4. Each iteration of remapping, look at the consensus sequence for each seed.
      If it has drifted closer to one of the other seeds than to its original seed,
      drop it. Drop one at a time, and recalculate distances, because two seeds can
      drift toward each other and both be candidates for dropping.
* `sam2aln`: extract the aligned portion of each read from a CSV SAM file. Also
    combines duplicates.
  * in - remap.csv
  * aligned.csv - reads aligned to consensus sequence
  * conseq_ins.csv - downloaded - insertions relative to consensus sequence
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
  * coord_ins.csv - downloaded - insertions relative to coordinate reference
  * conseq.csv - downloaded - consensus sequence
  * failed_align.csv - downloaded - any consensus that failed to align to its ref
  * nuc_variants.csv - downloaded - top nucleotide variants for HLA
* `coverage_plots`: convert amino counts to a coverage graph for each gene
    region.
  * in - amino.csv
  * coverage_maps.tar - binary file &rarr; untar in results/coverage_maps
  * coverage_scores.csv - downloaded - a score for each region based on the
    coverage at key positions.

## File descriptions ##
* quality.csv and bad_cycles.csv
  * tile
  * cycle
  * errorrate - as a percentage
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
* aligned.csv and g2p_aligned.csv
  * refname - seed reference the reads mapped to
  * qcut - minimum Phred quality score to include a nucleotide
  * rank - sequences are sorted by count, with the most common at rank 0
  * count - how many times exact copies of this sequence and soft clip counts
    were found
  * offset - offset of the read within the consensus, or number of dashes to
    add at the start
  * seq - the mapped sequence of the read, aligned to the consensus
* conseq_ins.csv
  * qname - query name from the read
  * fwd_rev - F for forward reads, R for reverse
  * refname - seed reference the read mapped to
  * pos - 1-based position in the consensus sequence that this insertion follows
  * insert - the nucleotide sequence that was inserted
  * qual - the Phred quality scores for the inserted sequence
* clipping.csv
  * refname - seed reference the reads mapped to
  * pos - one-based nucleotide position within the consensus sequence
  * count - number of read pairs that soft clipped instead of mapping to this
    position
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
  * seed - seed reference the reads mapped to
  * region - coordinate reference for reporting against, usually a gene
  * q-cutoff - minimum Phred quality score to include a nucleotide
  * query.nuc.pos - the 1-based index of the base in the consensus sequence that
    came from this set of counts
  * refseq.nuc.pos - the 1-based index of the base in the coordinate reference
  * A,C,G,T - counts for the nucleotides at this position
  * N - count of reads with Phred quality score below the cutoff
  * del - count of reads with a deletion at this position
  * ins - count of reads with an insertion after this position, relative to the
    coordinate reference
  * clip - count of reads with soft clipping that would have mapped at this
    position
  * g2p_overlap - count of reads at this position that overlapped with the
    V3LOOP G2P mapping and were ignored.
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
* coord_ins.csv - insertions in consensus sequence, relative to coordinate
    reference.
  * seed - seed reference the reads mapped to
  * region - coordinate reference for reporting against, usually a gene
  * qcut - minimum Phred quality score to include a nucleotide
  * left - the one-based position within the consensus sequence of the first
    nucleotide in the insertion
  * insert - the insertion sequence of amino acids
  * count - the number of times the insertion occurred
  * before - the one-based position within the coordinate reference that it
    was inserted before
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
    * low quality - more than half of the sequence had read quality with phred
        score < 33
    * notdiv3 - the sequence did not end on a codon boundary
    * zerolength - none of the sequence mapped to the expected portion of the
        reference sequence
    * cysteines - the start or end of the amino acid sequence was not a cycsteine
    * > 2 ambiguous - more than one amino acid was ambiguous, or one position
        more than two possible amino acids
    * stop codons - the sequence contained a stop codon
    * length - the amino acid sequence length was outside of the range 32 to 40
    * failed to align - the nucleotide sequence failed to align with the reference
    * count < 3 - this is a summary line reporting how many reads were rejected
        for low counts
  * comment - a warning for sequences that had a problem but were good enough
    to score
    * ambiguous - one amino acid had two possible values
* g2p_summary.csv - summarizes the sample and calls it either R5 or X4
  * mapped - the number of reads that mapped to V3LOOP
  * valid - the number of those reads that were valid sequences and not errors
  * X4calls - the number of those reads that were called X4
  * X4pct - X4calls as a percentage of valid
  * final - the final decision: blank if valid is 0, X4 if X4pct >= 2, otherwise
    R5

## The MiCall Monitor ##

MiCall Monitor (or just Monitor) handles the automated processing of new MiSeq data
through the MiCall pipeline.  It periodically scans the `RAW_DATA` folder for data,
and when new data appears it interfaces with Kive to start the processing.
This folder is populated outside of MiCall:

* Runs get uploaded by the MiSeq to `RAW_DATA`.
* The `watch_runs.rb` script in that folder watches for the files to finish
    copying, and then creates a file named `needsprocessing` in the folder.
* The [MiseqQCReport][] scripts upload the QC data to QAI, and then create a 
`qc_uploaded` file.

The Monitor looks for folders with both of these flag files, and ignores ones
without.

[MiseqQCReport]: https://github.com/cfe-lab/MiSeqQCReport/tree/master/modules

### Hourly scan for new MiSeq runs ###

Every hour, Monitor looks for new data to process with the following procedure.

* Scan for folders that have a `needsprocessing` flag.  Any such folders get added 
to a list (in memory) of all run folders.  This distinguishes other random stuff or 
stuff that isn't ready to go from actual run folders.  The newest of these folders 
(going by the date in the folder name, *not* the filesystem creation date) is 
earmarked as The Newest Run, regardless of any processing that has already been 
done by older versions of MiCall.  Folders are visited in reverse chronological 
order based on the date in the folder name.

* Iterate through the run folder list in memory.  For each run folder, look for a 
results subfolder corresponding to the current version of MiCall (located in
`Results/version_X.Y`, where `X.Y` is the current version number) with a 
`doneprocessing` flag in it.  If this is identified, skip this folder (regardless 
of whether or not it's The Newest Run).

* When a run folder is found that does *not* have such a results folder, one of 
two things happens.  If this folder is The Newest Run, Monitor gets/creates a MiCall 
pipeline run for each sample, all at once (see "Get/Create Run" below).  All other 
folders have their samples added to an in-memory list of "Samples That Need Processing".
   
    * Get/Create Run: Monitor looks for the existence of the required datasets on Kive 
    (by both MD5 and filename) and creates them if they don't exist.  Then, Monitor 
    looks to see if this data is already being processed through the current version of 
    MiCall.  If not, it starts the processing.  Either way, the Kive processing task
    is added to an in-memory list of "Samples In Progress".
    
### Periodic load monitoring and adding new processing jobs ###

Every 30 seconds, Monitor checks on all Samples In Progress.

* If a MiCall sample is finished, and that is the last one from a given MiSeq run to 
finish, then all results for that run are downloaded into a subfolder specific to the 
current version of MiCall (located in `Results/version_X.Y` where `X.Y` is the current 
version number) in that run's corresponding folder in `RAW_DATA`, and a `doneprocessing` 
flag is added to that subfolder.  Post-processing (e.g. uploading stuff to QAI) occurs.

* Monitor keeps a specified lower limit of samples to keep active at any given time.  
If a MiCall sample finishes processing and the number of active samples dips below 
that limit, Monitor looks at its list of Samples That Need Reprocessing and starts 
the next one, moving it from that list to Samples In Progress.

### Ways to manipulate the order in which the Monitor processes data ###
Sometimes, you want to do unusual things. Here are a few scenarios we've run into.

#### You need to force Monitor to reprocess a given run ####
First, stop Monitor.  Remove the results subfolder for the current version of 
MiCall.  Restart Monitor.  On the next hourly scan, Monitor will handle this 
folder as if it had never been processed through the current version of MiCall.  That 
is, if it's The Newest Run, its samples will be immediately started and added to 
Samples In Progress, and if it's not The Newest Run, its samples will be added to 
Samples That Need Processing.

#### You need Monitor to skip a folder and handle an older one first ####
Stop Monitor.  Add an `errorprocessing` flag to the run's folder in `RAW_DATA`.  
This will make Monitor's next hourly scan believe that it's failed and should be 
skipped, and Monitor will move on to the next one.  Note though that this has no 
effect on which folder is The Newest Run: even if you're skipping The Newest Run, 
the next one down does not inherit the mantle.

* DO NOT delete the `needsprocessing` flag from a folder to try and keep Monitor from 
handling it.  This will cause Conan's scripts to reupload data to that folder.  

Restart Monitor.  After all samples from this run have finished processing, 
remove the fake `errorprocessing` flag you set (no need to stop and restart Monitor).

#### You need to stop what's currently running and handle an older run ####
Stop Monitor.  In Kive, stop the processing tasks that you need to clear out, 
but don't remove them; the progress you've already made can be reused later when 
revisiting these tasks later.  Now, do the manipulations you need to do as in the
above case to make Monitor deal with your desired run first.  Restart Monitor.

As in the above case, when you are ready to process the run you previously stopped,
you can remove the fake `errorprocessing` flag you created for that run, and Monitor
will then restart those processing tasks on its next hourly scan.  Kive will be able 
reuse the progress already made when you stopped them.
