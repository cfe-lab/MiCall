These sample files exercise certain features of the pipeline with tiny FASTQ
files that are very quick to process. To build them, we took reference
sequences from the projects.json file, and either took snippets from them
or introduced deletions, insertions, and variants.

This folder also contains some useful scripts for manipulating these sequences.

To process these samples, copy the FASTQ files, sample sheet, and quality file
into the working folder, then run the sample_pipeline or run_processor launch
configurations.

The scenarios that each file tests are:

* 1234A-V3LOOP - nine identical reads, plus one with a single change in the
  fifth base
* 2000A-V3LOOP - missing coverage of the first codon
* 2010A-V3LOOP - five reads at the beginning of the reference and five reads at
  the end, with a one-codon gap between
* 2020A-GP41 - starts at second codon and inserts extra codon after first five
* 2030A-V3LOOP - garbage reads, all different
* 2040-HLA-B - exon2 of HLA-B, half have a changed base at position 7 to test
  mixtures in the consensus
* 2050-V3LOOP - low quality
* 2060A-V3LOOP - clean read tests g2p with three codon overlap between forward
  and reverse
* 2070A-PR - seven reads with two deletions and three with one deletion
* 2080A-V3LOOP - PR contamination below remap threshold
* 2090A-HCV - mixed infection of HCV 1A and 1B, in the NS5b region.
* 2100A-HCV-1337B-V3LOOP - two samples using same tags on different regions
* 2110A-V3LOOP - 5 clean codons, then 1 stop codon in codon 6, 2 low-quality
  codons in codon 7, 3 partial deletions in codon 8, 4 insertions in codon 9,
  3 soft clipped reads in codons 10 and 11, 1 deletion in codon 12, followed by
  6 clean codons with only 7 reads covering them
* 2120A-PR - 5 clean codons, then 1 stop codon in codon 6, 2 low-quality
  codons in codon 7, 3 partial deletions in codon 8, 4 insertions in codon 9,
  3 soft clipped reads in codons 10 and 11, 1 deletion in codon 12, followed by
  6 clean codons with only 7 reads covering them (soft clipped and insertions not working)
