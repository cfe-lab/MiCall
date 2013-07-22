#!/bin/bash

# Generate matrix of different quality rules

# Already done for Luke - A43J1: DATAPATH="/Users/emartin/Desktop/130524_M01841_0004_000000000-A43J1/Data/Intensities/BaseCalls/remap_sams/"
DATAPATH="/Users/emartin/Desktop/MiSeq/MiSeq_Runs/130621_M01841_0006_000000000-A3RWC/Data/Intensities/BaseCalls"
exit

# Unzip fastq.gz files
echo gunzip $DATAPATH/*.gz
exit

# Run pipeline to generate <sample>.<regionCode>.remap.sam files
python 1_pipeline3.py $DATAPATH/
exit

# Generate fasta files from remap.sam files: censor low quality bases with an 'N'
python 2_sam2fasta_with_censoring.py $DATAPATH 20
exit

# Generate seq files
python 3_compress_fasta_with_min_count.py 3 $DATAPATH*remap.sam.[0-9]*.fasta
exit

# Extract V3 out of env sequences (qCutoff determines the glob of files to load, but not processing behavior)
python 4_extractV3_with_G2P.py 20 $DATAPATH &
exit

# Determine the proportion X4 per sample (3, 3.5, 4, 5, 5.75, 7.5)
rm results
python 5_determine_proportion_X4.py 3.0 $DATAPATH >> results
exit
