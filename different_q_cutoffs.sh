#!/bin/bash

# Use absolute paths
# If referencing folders, do not put a slash at the end of the path
# This shell script must be in the same directory as the pipeline

date
python pipeline3_resume_from_sam2fasta.py ../../130524_M01841_0004_000000000-A43J1/Data/Intensities/BaseCalls/remap_sams/ 0
python pipeline3_resume_from_sam2fasta.py ../../130524_M01841_0004_000000000-A43J1/Data/Intensities/BaseCalls/remap_sams/ 5
python pipeline3_resume_from_sam2fasta.py ../../130524_M01841_0004_000000000-A43J1/Data/Intensities/BaseCalls/remap_sams/ 10
python pipeline3_resume_from_sam2fasta.py ../../130524_M01841_0004_000000000-A43J1/Data/Intensities/BaseCalls/remap_sams/ 15
python pipeline3_resume_from_sam2fasta.py ../../130524_M01841_0004_000000000-A43J1/Data/Intensities/BaseCalls/remap_sams/ 20
python pipeline3_resume_from_sam2fasta.py ../../130524_M01841_0004_000000000-A43J1/Data/Intensities/BaseCalls/remap_sams/ 25
date
