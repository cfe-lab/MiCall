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
	amino_freq_csv_reader = csv.DictReader(f)
	for row in amino_freq_csv_reader:
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

			# Write temp CSV file R for this (sample, region, q-cutoff)
			csv_file = "{}_{}_{}.csv".format(sample,region,q)
			if "_cleaned_" in amino_freq_csv:
				csv_file = "CLEAN_{}".format(csv_file)
			csv_path = "/tmp/{}".format(csv_file)

			with open(csv_path,"wb") as f_out:
				f_out.write("refseq.aa.pos,coverage\n")
				for row in dataset:
					coverage = 0
					for aa in "ACDEFGHIKLMNPQRSTVWY*":
						coverage += int(row[aa])
					f_out.write("{},{}\n".format(row['refseq.aa.pos'],coverage))

			# Call the R script on the temp csv, then move the png into place + remove the temp csv
			png_path = csv_path.replace(".csv", ".png")
			command = "/usr/bin/env Rscript coverage_plot.R {} {} {} {} {}".format(csv_path, sample, region, q, png_path)
			os.system(command)
			os.remove(csv_path)
			png_destination = "{}/{}".format(images_folder.rstrip("/"),os.path.basename(png_path))
			shutil.move(png_path, png_destination)
