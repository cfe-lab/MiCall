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

# activate plot device
png(file=file_out, width=500, height=500, type="cairo")

# set up plot region
par(mar=c(5,6,4,2))
plot(NA, xlim=c(1,max(data$refseq.aa.pos)), ylim=c(100,1000000), axes=FALSE, ann=FALSE, xaxs="r", log="y")


main_title = paste(sample, "\n[", region, "] q>=", q_cutoff, sep="")
title(main=main_title, xlab="HXB2 amino acid position", font.lab = 1.4, cex.lab=1.4, cex.main=1.4)
title(ylab='Coverage', cex.lab=1.5, line=4)

#lines(data$coverage + 1)		# Add 1 to allow plotting 0 values on log scale
x <- data$refseq.aa.pos
polygon(x = c(min(x), x, max(x), min(x)), 
	y = c(0.01, data$coverage, 0.01, 0.01),
	col='lightblue')

# draw minimum coverage line
abline(h = min_coverage, lty=2)

axis(1)
axis(2, at=c(1E2, 1E3, 1E4, 1E5), labels=c('100', '1000', '10,000', '100,000'), las=2)
box()

garbage <- dev.off()
