#!/usr/local/bin/python

import csv,os,shutil,sys

if len(sys.argv) != 3:
	print "Syntax: {} [PATH_TO_AMINO_CSV] [PATH_TO_IMAGES_FOLDER]".format(sys.argv[0])
	sys.exit()

amino_freq_csv = sys.argv[1]		# Ex: /data/miseq/140207_M01841_0055_000000000-A64ED/amino_frequencies.csv
images_folder = sys.argv[2]			# Ex: /media/RAW_DATA/MiSeq/runs/140207_M01841_0055_000000000-A64ED/Results/version_5/coverage_maps

# Load amino_freqs.csv into memory
my_rows = []
with open(amino_freq_csv, "r") as f:
	amino_freq_csv = csv.DictReader(f)
	for row in amino_freq_csv:
		my_rows.append(row)

# Make images folder if it doesn't already exist
if not os.access(images_folder, os.F_OK):
	os.mkdir(images_folder)

# Get all coordinates for one group of data (A particular sample, region, q-cutoff combination)
for sample in set([x["sample"] for x in my_rows]):
	for region in set([x["region"] for x in my_rows if x["sample"] == sample]):
		for q in set([x["q-cutoff"] for x in my_rows if x["sample"] == sample and x["region"] == region]):

			dataset = [x for x in my_rows if x["sample"] == sample and x["region"] == region and x["q-cutoff"] == q]
			dataset.sort(key=lambda row: int(row['refseq.aa.pos']))
			csv_tmp_file = "{}_{}_{}.csv".format(sample,region,q)

			# Write temp CSV file for R for this (sample, region, q-cutoff)
			with open(csv_tmp_file,"wb") as f_out:
				f_out.write("refseq.aa.pos,coverage\n")
				for row in dataset:
					coverage = 0
					for aa in "ACDEFGHIKLMNPQRSTVWY*":
						coverage += int(row[aa])
					f_out.write("{},{}\n".format(row['refseq.aa.pos'],coverage))

			# Call the R script on the temp csv, then move the png into place and remove the temp csv
			png_file = csv_tmp_file.replace(".csv", ".png")
			os.system("/usr/bin/env Rscript coverage_plot.R {} {} {} {} {}".format(csv_tmp_file, sample, region, q, png_file))
			os.remove(csv_tmp_file)
			shutil.move(png_file, "{}/{}".format(images_folder.rstrip("/"),png_file))
