#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
if (length(args) != 2) { 
	stop('Usage: error-metrics.R <input CSV> <output PDF>') 
}
input.csv <- args[1]
output.pdf <- args[2]

# generated from Interop ErrorMetricsOut.bin file with Python
# script parse-interop.py
error <- read.csv(input.csv, header=T)

max.cycle <- max(error$cycle)

# negative cycle values correspond to the second read
error$true.cycle <- ifelse(error$cycle > 0, error$cycle, max.cycle-error$cycle)
by.tile <- split(error, error$tile)
n.tiles <- length(by.tile)

# prepare colouring scheme
hist(log10(error$error)) -> h
n.bins <- length(h$mids)
pal <- rev(rainbow(n.bins, start=0, end=0.35))


# prepare plot region
pdf(file=output.pdf, width=10, height=5)

par(mar=c(5,5,2,5))
plot(NA, xlim=range(error$true.cycle), ylim=c(1, n.tiles), type='n', xlab='Cycle number', ylab='Tile', yaxt='n', cex.lab=1.2)
axis(side=2, at=1:n.tiles, names(by.tile), las=2, cex.axis=0.5)
abline(v=max.cycle)

# perform plot
for (i in 1:n.tiles) {
	tile <- names(by.tile)[i]
	error.rates <- log10(by.tile[[i]]$errorrate)
	max.cycle <- max(by.tile[[i]]$true.cycle, na.rm=T)
	hues <- pal[as.integer(cut(error.rates, breaks=h$breaks))]
	
	# apply functions are supposed to be faster than for-loops...
	sapply(1:max.cycle, function(j) {
		rect(xleft=j-0.5, xright=j+0.5, ybottom=i-0.5, ytop=i+0.5, 
			 col=hues[j], border=NA)
	})
}

# generate legend
y.scaling <- n.tiles / n.bins
max.x <- max(error$true.cycle, na.rm=TRUE)
x.scaling <- max.x / n.tiles


par(xpd=TRUE)
. <- sapply(1:n.bins, function(i) {
	rect(xleft=1.05*max.x, xright=1.05*max.x+x.scaling, ybottom=y.scaling*(i-0.5), ytop=y.scaling*(i+0.5), col=pal[i])
	text(x=1.06*max.x+x.scaling, y=y.scaling*i, label=round(10^h$mid[i], 3), adj=0, cex=0.6)
})
par(xpd=FALSE)

. <- dev.off()  # finished with plotting device, close file
