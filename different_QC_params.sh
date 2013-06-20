#!/bin/bash

# Generate a matrix of different QC cutoffs

# 1) Generate fasta files with different base censoring rules
python pipeline3_resume_from_sam2fasta.py ../../130524_M01841_0004_000000000-A43J1/Data/Intensities/BaseCalls/remap_sams/ 0
python pipeline3_resume_from_sam2fasta.py ../../130524_M01841_0004_000000000-A43J1/Data/Intensities/BaseCalls/remap_sams/ 10
python pipeline3_resume_from_sam2fasta.py ../../130524_M01841_0004_000000000-A43J1/Data/Intensities/BaseCalls/remap_sams/ 15
python pipeline3_resume_from_sam2fasta.py ../../130524_M01841_0004_000000000-A43J1/Data/Intensities/BaseCalls/remap_sams/ 20

# 2) Focus only on HIV1B-env fasta files
rm /Users/emartin/Desktop/130524_M01841_0004_000000000-A43J1/Data/Intensities/BaseCalls/remap_sams/*HLA*.fasta
rm /Users/emartin/Desktop/130524_M01841_0004_000000000-A43J1/Data/Intensities/BaseCalls/remap_sams/*pol*.fasta

# 3) Generate different min count seq files from fasta (Ex: F00142.HIV1B-env.remap.sam.15.fasta.3.seq)
python compress_fasta.py 1 ../../130524_M01841_0004_000000000-A43J1/Data/Intensities/BaseCalls/remap_sams/*remap.sam.[0-9]*.fasta
python compress_fasta.py 3 ../../130524_M01841_0004_000000000-A43J1/Data/Intensities/BaseCalls/remap_sams/*remap.sam.[0-9]*.fasta
python compress_fasta.py 50 ../../130524_M01841_0004_000000000-A43J1/Data/Intensities/BaseCalls/remap_sams/*remap.sam.[0-9]*.fasta

# 4) Extract V3 out of env sequences (Ex: F00142.HIV1B-env.remap.sam.15.fasta.3.seq.v3)
# 5) Feed into G2P with FPR cutoffs:	3.0  3.5  4.0  5.0  5.75  7.5
