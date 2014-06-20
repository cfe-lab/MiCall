#!/usr/bin/env Rscript

# parameters
min.coverage <- 1000


syntax.error = "Correct syntax: ./coverage_plot.R <input CSV> <output path>"
args <- commandArgs(TRUE)

if (length(args) != 2) { stop(syntax.error) }
input.csv <- args[1]
out.path <- args[2]

if (!file.exists(out.path)) {
	stop('output path does not exist')
}

# check for trailing directory separator
l <- nchar(out.path)
if (substr(out.path, l, l) != '/') {
	out.path <- paste(out.path, '/', sep='')
}

# load key positions file
key.pos <- read.csv(file='key_positions.csv', header=FALSE)
names(key.pos) <- c('target', 'pos')


# path to CSV file for reporting minimum coverages at key positions
output.csv <- paste(out.path, 'minimum_coverage_at_keys.csv', sep='')
output <- {}

data <- read.csv(file=input.csv, header=TRUE, sep=',')

# the input CSV should look like this:
#sample,region,q-cutoff,query.aa.pos,refseq.aa.pos,A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*
#50955ARPT-HCV-49537A-INT-PR-RT_S16,HCV1A-H77-core,0,0,1,0,0,0,0,0,0,0,3,0,0,1542,0,0,0,0,0,0,0,0,0,0

alphabet <- c('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y')  # X. represents stop codon '*'

data$coverage <- apply(data[ , which(is.element(names(data), alphabet))], 1, sum)

# partition AA frequency table by sample and region
coverage <- split(data[,which(is.element(names(data), c('refseq.aa.pos', 'q.cutoff', 'coverage')))], f=list(data$region, data$sample), drop=TRUE)

# loop through partitions
for (i in 1:length(coverage)) {
	label <- names(coverage)[i]
	tokens <- strsplit(label, split='\\.')[[1]]
	region <- tokens[1]
	sample <- tokens[2]
	
	df <- coverage[[i]]
	filename <- paste(sample, region, 'png', sep='.')
	
	# set up plot
	png(file=paste(out.path, filename, sep=''), width=400, height=300, type='cairo')
	par(family='sans', cex=1, mar=c(5,5,1,1))
	plot(NA, xlim=c(1,max(df$refseq.aa.pos)), ylim=c(1,200000), axes=FALSE, ann=FALSE, xaxs="r", log="y")
	title(xlab="Reference coordinates (AA)", font.lab = 1.4, cex.lab=1.4, cex.main=1.4)

	# indicate key positions for this region
	temp <- key.pos$pos[key.pos$target==region]
	#points(temp, rep(min.coverage, length(temp)), pch=8, cex=1.5)
	for (pos in temp) {
		rect(pos, 500, pos+1, 2e3, col='red', border=NA)
	}

	# draw trend lines for each quality score cutoff
	cutoffs <- sort(as.integer(levels(as.factor(df$q.cutoff))))
	for (j in 1:length(cutoffs)) {
		q.cut <- cutoffs[j]
		df2 <- df[df$q.cutoff == q.cut, ]
		lines(x = df2$refseq.aa.pos, 
			y = df2$coverage, 
			col=rainbow(length(cutoffs), v=0.8)[j],
			lwd=2
		)
		# determine minimum coverage at key positions
		key.coverage <- df2$coverage[is.element(df2$refseq.aa.pos, temp)]
		if (length(key.coverage) > 0) {
			output <- rbind(output, c(sample, region, q.cut, min(key.coverage), temp[which.min(key.coverage)]))
		} else {
			output <- rbind(output, c(sample, region, q.cut, NA, NA))
		}
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


# output minimum coverage at key positions
output <- as.data.frame(output)
names(output) <- c('sample', 'region', 'q.cut', 'min.coverage', 'which.key.pos')
write.csv(output, file=output.csv, quote=FALSE, row.names=FALSE)
