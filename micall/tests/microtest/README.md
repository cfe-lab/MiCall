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
* 2010A-V3LOOP - ten reads at the beginning of the reference and ten reads at
  the end, with a two-codon gap between
* 2020A-GP41 - starts at second codon and inserts extra codon after first ten
* 2030A-V3LOOP - garbage reads, all different
* 2040A-HLA-B - exon2 of HLA-B, 20% have a changed base at nucleotide 7 to test
  mixtures in the consensus
* 2050-V3LOOP - low quality
* 2060A-V3LOOP - clean read tests g2p with four codon overlap between forward
  and reverse
* 2070A-PR - three reads with two deletions and twelve with one deletion
* 2080A-V3LOOP - PR contamination below remap threshold
* 2100A-HCV-1337B-V3LOOP-PWND-HIV - three samples using same tags on different
  regions: Nextera HIV-RT and HCV-NS5a, plus V3LOOP amplicon.
* 2110A-V3LOOP - 5 clean codons, then 1 read with a stop codon in codon 6, 2
  low-quality reads in codon 7, 3 reads with partial deletions in codon 8, 4
  reads with insertions after codon 9,
  4 reads with a deletion in codon 10, followed by 6 clean codons.
  There are also 10 clean reads that trigger the amplicon detector.
* 2120A-PR - 29 codons with no coverage, clean codons from 30-46, then 1 read
  with stop codons in codons 47-49, 2 low-quality reads in codons 50-52, 3
  reads with partial deletions in codons 53 and 54, 4 reads with insertions
  after codon 56, 3 soft clipped reads in codons 59-64, 1 read with a deletion
  in codons 65-67, followed by 18 clean codons with only 10 reads covering them.
  There are also 10 reads that trigger the amplicon detector.
* 2130A-HCV - full coverage of the whole genome portion of NS5b, with a couple
  of mutations. See the `make_sample.py` script for details.
* 2130AMIDI-MidHCV - full coverage of the MIDI portion of NS5b, with a couple
  of mutations. See the `make_sample.py` script for details.
* 2140A-HIV - full coverage of PR, with a mutation. See the `make_sample.py`
  script for details.
* 2160A-HCV - enough coverage of the whole genome portion of NS5b for IVA to
  assemble a contig, with a couple of mutations. See the `make_sample.py`
  script for details.
* 2160AMIDI-MidHCV - MIDI region for 2160A. See the `make_sample.py` script for
  details.
* 2170A-HCV - a mixed infection of HCV-1a and HCV-2a, each about 3000 bases
  long, with enough coverage to assemble two contigs.
* 2180A-HIV - Nextera reads from GP120 and V3LOOP, to show how GP120 is now
  reported based on remapping or assembly, and V3LOOP is reported based on
  pairwise alignment to a reference.
* 2190A-SARSCOV2 - two lengths of SARS-CoV-2 reads, to test that similar amplicons
  of slightly different lengths can be combined.
* 2200A-SARSCOV2 - SARS-CoV-2 amplicon that gets primers trimmed off the left,
  so the coverage starts on amino 27 of nsp1, instead of 20.
* 2210A-NFLHIVDNA - Nextera reads that cover the primer sections we check for
  in HIV proviral samples.
