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
  * quality.csv
* `prelim_map`: map reads to all references.
  * in - fastq1
  * in - fastq2
  * prelim.csv - SAM file format
* `remap`: iteratively use consensus from previous mapping as reference to try
    and map more reads. Eliminate a consensus if it has changed into a
    different reference.
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
* `sam2aln`: extract the aligned portion of each read from a CSV SAM file. Also
    combines duplicates.
  * in - remap.csv
  * aligned.csv - reads aligned to consensus sequence
  * conseq_ins.csv - downloaded - insertions relative to consensus sequence
  * failed_read.csv - downloaded - reads that fail to merge
* `aln2counts`: take the aligned reads, and count how many of each nucleotide
    or amino acid appear at each position.
  * in - aligned.csv
  * nuc.csv - downloaded - nucleotide counts at each position
  * amino.csv - downloaded - amino counts at each position
  * coord_ins.csv - downloaded - insertions relative to coordinate reference
  * conseq.csv - downloaded - consensus sequence
  * failed_align.csv - downloaded - any consensus that failed to align to its ref
  * nuc_variants.csv - downloaded - top nucleotide variants for HLA
* `sam_g2p`: use a position-specific scoring matrix on the V3LOOP region to
    translate genotype to phenotype and predict whether a sample is X4 or R5.
    HIV-1 particles use coreceptors to enter cells. Different particles can use
    CXCR4, CCR5, or both. This analysis looks at the sequences from a virus
    population and reports either "X4" (able to use CXCR4) or "R5" (only able
    to use CCR5).
  * in - remap.csv
  * in - nuc.csv
  * g2p.csv - downloaded - calls each individual sequence X4, R5, or error.
  * g2p_summary.csv - downloaded - calls the entire sample X4 or R5. 
* `coverage_plots`: convert amino counts to a coverage graph for each gene
    region.
  * in - amino.csv
  * coverage_maps.tar - binary file &rarr; untar in results/coverage_maps
  * coverage_scores.csv - downloaded - a score for each region based on the
    coverage at key positions.

## File descriptions ##
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
* coord_ins.csv - insertions in consensus sequence, relative to coordinate
    reference.
  * seed - seed reference the reads mapped to
  * region - coordinate reference for reporting against, usually a gene
  * qcut - minimum Phred quality score to include a nucleotide
  * left - the one-based position within the consensus sequence of the first
    amino acid in the insertion
  * insert - the insertion sequence of amino acids
  * count - the number of times the insertion occurred
  * before - the one-based position within the coordinate reference that it
    was inserted before
  

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
