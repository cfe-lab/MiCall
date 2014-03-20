#!/usr/bin/env Rscript

min_coverage = 1000
syntax_error = "Correct syntax: ./coverage_plot.R <input CSV> <sample> <region> <q-cutoff> <file_out>"
args<-commandArgs(TRUE)

if (length(args) != 5) { stop(syntax_error) }
file_path = args[1]
sample = args[2]
region = args[3]
q_cutoff = args[4]
file_out = args[5]

data <- read.csv(file=file_path,head=TRUE,sep=",")
png(file=file_out, width=500, height=500, type="cairo")
plot(NA, xlim=c(1,length(data$refseq.aa.pos)), ylim=c(1,1000000), axes=FALSE, ann=FALSE, xaxs="r", log="y")
abline(h = min_coverage, lty=2)
main_title = paste(sample, "\n[", region, "] q>=", q_cutoff, sep="")
title(main=main_title, xlab="Coordinates", ylab="Coverage", font.lab = 1.4, cex.lab=1.4, cex.main=1.4)
axis(1)
axis(2)
box()
lines(data$coverage + 1)		# Add 1 to allow plotting 0 values on log scale
garbage <- dev.off()
