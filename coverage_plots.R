#!/usr/bin/env Rscript

min.coverage = 1000
syntax_error = "Correct syntax: ./coverage_plot.R <input CSV> <output path>"
args <- commandArgs(TRUE)

if (length(args) != 2) { stop(syntax_error) }
input_csv = args[1]
out_path = args[2]

if (!file.exists(out_path)) {
	stop('output path does not exist')
}

# check for trailing directory separator
l <- nchar(out_path)
if (substr(out_path, l, l) != '/') {
	out_path <- paste(out_path, '/', sep='')
}

data <- read.csv(file=input_csv, header=TRUE, sep=',')

#sample,region,q-cutoff,query.aa.pos,refseq.aa.pos,A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*
#50955ARPT-HCV-49537A-INT-PR-RT_S16,HCV1A-H77-core,0,0,1,0,0,0,0,0,0,0,3,0,0,1542,0,0,0,0,0,0,0,0,0,0

alphabet <- c('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y')  # X. represents stop codon '*'

data$coverage <- apply(data[ , which(is.element(names(data), alphabet))], 1, sum)

coverage <- split(data[,which(is.element(names(data), c('refseq.aa.pos', 'q.cutoff', 'coverage')))], f=list(data$region, data$sample), drop=TRUE)

for (i in 1:length(coverage)) {
	label <- names(coverage)[i]
	tokens <- strsplit(label, split='\\.')[[1]]
	region <- tokens[1]
	sample <- tokens[2]
	
	df <- coverage[[i]]
	filename <- paste(sample, region, 'png', sep='.')
	
	# set up plot
	png(file=paste(out_path, filename, sep=''), width=400, height=300, type='cairo')
	par(family='sans', cex=1, mar=c(5,5,1,1))
	plot(NA, xlim=c(1,max(df$refseq.aa.pos)), ylim=c(1,200000), axes=FALSE, ann=FALSE, xaxs="r", log="y")

	#main_title = paste(sample, "\n[", region, "] q>=", qcut, sep="")
	title(xlab="Reference coordinates", font.lab = 1.4, cex.lab=1.4, cex.main=1.4)
	#title(ylab='Coverage', cex.lab=1.5, line=4)

	cutoffs <- sort(as.integer(levels(as.factor(df$q.cutoff))))
	for (j in 1:length(cutoffs)) {
		q.cut <- cutoffs[j]
		lines(x = df$refseq.aa.pos[df$q.cutoff == q.cut], 
			y = df$coverage[df$q.cutoff==q.cut], 
			col=rainbow(length(cutoffs), v=0.8)[j],
			lwd=2
		)
	}

	# draw minimum coverage line
	abline(h = min.coverage, lty=2)
	axis(1)
	axis(2, at=c(1E1, 1E2, 1E3, 1E4, 1E5), labels=c('10', '100', '1000', '10,000', '100,000'), las=2)
	box()
	
	legend(x=0, y=1, legend=paste('q=', cutoffs, sep=''), col=rainbow(length(cutoffs), v=0.8), lty=1, lwd=2, bty='n', yjust=0)
	text(x=max(df$refseq.aa.pos)/2, y=190000, label=paste(sample, region), cex=0.7, col='grey30', adj=c(0.5, 0.5))
	
	garbage <- dev.off()
}
