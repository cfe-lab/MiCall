#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
if (length(args) != 2) {
	stop('Usage: coverage.R <input CSV> <output file prefix>')
}
input.csv <- args[1]
out.prefix <- args[2]

df <- read.csv(input.csv)

# guess if this is a nuc or amino CSV
nucs <- c('A', 'C', 'G', 'T')
aminos <- c('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P',
'Q', 'R', 'S', 'T', 'V', 'W', 'Y', '*')


if (all(is.element(aminos, names(df)))) {
    is.nuc <- FALSE
    df$count <- apply(subset(df, select=aminos), 1, sum)
} else if (all(is.element(nucs, names(df)))) {
    is.nuc <- TRUE
    df$count <- apply(subset(df, select=nucs), 1, sum)
} else {
    stop("Error: input not recognized as a nuc or amino CSV")
}

# generate a coverage plot for each seed and region
temp <- split(df, list(df$region, df$seed), drop=TRUE)


. <- sapply(1:length(temp), function(i) {
    obj <- temp[[i]]
    
    key <- strsplit(names(temp)[i], '\\.')[[1]]
    region <- key[1]
    seed <- key[2]
    
    fn <- paste(out.prefix, seed, region, 'pdf', sep='.')

    # start plot device
    pdf(file=fn, width=6, height=6)
    par(mar=c(5,5,2,2))
    
    if (is.nuc) {
      x <- obj$refseq.nuc.pos
    } else {
      x <- obj$refseq.amino.pos
    }
  
    plot(x, obj$count, type='s',
    xlab='Reference coordinate (nt)', ylab='Coverage')
    title(main=paste(seed, region), cex=0.5, adj=0)

    dev.off()
})


