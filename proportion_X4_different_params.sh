#!/bin/bash

# Generate matrix of different quality rules

DATAPATH="/Users/emartin/Desktop/130524_M01841_0004_000000000-A43J1/Data/Intensities/BaseCalls/remap_sams/"

#python determine_proportion_X4.py 3.0 /Users/emartin/Desktop/130524_M01841_0004_000000000-A43J1/Data/Intensities/BaseCalls/remap_sams/

exit # ONLY ACTIVATE THIS IF YOU REALLY WANT THIS TO GO...!!

# 1) Starting with *.remap.sam files, generate .fasta files (with base censoring rules)
python pipeline3_resume_from_sam2fasta.py $DATAPATH 0
python pipeline3_resume_from_sam2fasta.py $DATAPATH 10
python pipeline3_resume_from_sam2fasta.py $DATAPATH 15
python pipeline3_resume_from_sam2fasta.py $DATAPATH 20

# 2) Focus only on HIV1B-env fasta files
rm $DATAPATH*HLA*.fasta
rm $DATAPATH*pol*.fasta

# 3) Generate different min count seq files from fasta (Ex: F00142.HIV1B-env.remap.sam.15.fasta.3.seq)
python compress_fasta.py 1 $DATAPATH*remap.sam.[0-9]*.fasta
python compress_fasta.py 3 $DATAPATH*remap.sam.[0-9]*.fasta
python compress_fasta.py 50 $DATAPATH*remap.sam.[0-9]*.fasta
exit

# 4) Extract V3 out of env sequences (Ex: Check *.HIV1B-env.remap.sam.10.fasta.*.seq files, generate .v3 files)
# FIXME: NEEDS TO BE FORMALLY PARALLELIZED, WITH AMPERSAND?
python extractV3_with_G2P.py 0 $DATAPATH &
python extractV3_with_G2P.py 10 $DATAPATH &
python extractV3_with_G2P.py 15 $DATAPATH &
python extractV3_with_G2P.py 20 $DATAPATH &

exit

# 5) Feed into G2P with FPR cutoffs:	3.0  3.5  4.0  5.0  5.75  7.5
rm results
python determine_proportion_X4.py 3.0 $DATAPATH >> results
python determine_proportion_X4.py 3.5 $DATAPATH >> results
python determine_proportion_X4.py 4.0 $DATAPATH >> results
python determine_proportion_X4.py 5.0 $DATAPATH >> results
python determine_proportion_X4.py 5.75 $DATAPATH >> results
python determine_proportion_X4.py 7.5 $DATAPATH >> results

#python determine_proportion_X4.py 3.0 /Users/emartin/Desktop/130524_M01841_0004_000000000-A43J1/Data/Intensities/BaseCalls/remap_sams/ >> results
#python determine_proportion_X4.py 3.5 /Users/emartin/Desktop/130524_M01841_0004_000000000-A43J1/Data/Intensities/BaseCalls/remap_sams/ >> results
#python determine_proportion_X4.py 4.0 /Users/emartin/Desktop/130524_M01841_0004_000000000-A43J1/Data/Intensities/BaseCalls/remap_sams/ >> results
#python determine_proportion_X4.py 5.0 /Users/emartin/Desktop/130524_M01841_0004_000000000-A43J1/Data/Intensities/BaseCalls/remap_sams/ >> results
#python determine_proportion_X4.py 5.75 /Users/emartin/Desktop/130524_M01841_0004_000000000-A43J1/Data/Intensities/BaseCalls/remap_sams/ >> results
#python determine_proportion_X4.py 7.5 /Users/emartin/Desktop/130524_M01841_0004_000000000-A43J1/Data/Intensities/BaseCalls/remap_sams/ >> results
