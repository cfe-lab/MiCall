#!/usr/bin/env Rscript

# Command line expects these arguments:
# 1. input file name (amino frequencies CSV)
# 2. input file, nucleotide frequencies CSV
# 3. output maps archive (tar file where coverage maps will be written)
# 4. output scores (count scores CSV)

# Other dependency: key_positions.csv
# Working files are created in the current directory.


syntax.error = "Correct syntax: ./coverage_plot.R <amino CSV> <nuc CSV> <output maps tar> <output scores CSV>"
args <- commandArgs(TRUE)

if (length(args) != 4) { stop(syntax.error) }
input.csv <- args[1]
nuc.csv <- args[2]
output.maps <- args[3]
output.csv <- args[4]

dir.create('coverage_maps')

# load key positions file
key.pos.ranges <- read.csv(file='key_positions.csv', header=TRUE)
key.pos <- split(key.pos.ranges, f=key.pos.ranges$region)
for (region in names(key.pos)) {
    ranges <- key.pos[[region]]
    positions <- vector()
    for (i in seq(length=length(ranges$ref_pos_start))) {
        ref_pos_start <- ranges$ref_pos_start[i]
        ref_pos_end <- ranges$ref_pos_end[i]
        if (is.na(ref_pos_end)) {
            positions <- append(positions, ref_pos_start)
        }
        else {
            positions <- append(positions, ref_pos_start:ref_pos_end)
        }
    }
    key.pos[[region]] <- positions
}

output.columns <- c(
        'sample',
        'region',
        'q.cut',
        'min.coverage',
        'which.key.pos',
        'off.score',
        'on.score')
output <- matrix(nrow=0, ncol=length(output.columns), dimnames=list(NULL, output.columns))

data <- read.csv(file=input.csv, header=TRUE, sep=',')

# the input CSV should look like this:
#sample,region,q-cutoff,query.aa.pos,refseq.aa.pos,A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*
#50955ARPT-HCV-49537A-INT-PR-RT_S16,HCV1A-H77-core,0,0,1,0,0,0,0,0,0,0,3,0,0,1542,0,0,0,0,0,0,0,0,0,0

alphabet <- c('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y')  # X. represents stop codon '*'

data$coverage <- apply(data[ , which(is.element(names(data), alphabet))], 1, sum)

# zeroes don't show up on a log plot, so we mark no coverage at 0.1
data$coverage <- sapply(data$coverage, function(x) max(x, 0.1))

# partition AA frequency table by sample and region
coverage <- split(data[,which(is.element(names(data), c('refseq.aa.pos', 'q.cutoff', 'coverage')))], f=list(data$region, data$sample), drop=TRUE)

# loop through partitions
for (i in seq(length=length(coverage))) {
    label <- names(coverage)[i]
    tokens <- strsplit(label, split='\\.')[[1]]
    region <- tokens[1]
    sample <- tokens[2]
    
    df <- coverage[[i]]
    filename <- file.path('coverage_maps', paste(sample, region, 'png', sep='.'))
    
    # set up plot
    png(file=filename, width=400, height=300, type='cairo')
    par(family='sans', cex=1, mar=c(5,5,1,1))
    plot(NA, xlim=c(1,max(df$refseq.aa.pos)), ylim=c(1,200000), axes=FALSE, ann=FALSE, xaxs="r", log="y")
    title(xlab="Reference coordinates (AA)", font.lab = 1.4, cex.lab=1.4, cex.main=1.4)
    
    # indicate key positions for this region
    temp <- key.pos[[region]]
    #points(temp, rep(min.coverage, length(temp)), pch=8, cex=1.5)
    for (pos in temp) {
        rect(pos-0.5, 500, pos+0.5, 2e3, col='red', border=NA)
    }
    
    if (length(temp) == 0) {
        temp = seq(max(df$refseq.aa.pos))
    }
    if (substring(region, 1, 4) == 'HLA-') {
        good.coverage <- 100
        coverage.levels <- c(-Inf, 0, 10, good.coverage, Inf)
    }
    else {
        good.coverage <- 1000
        coverage.levels <- c(-Inf, 0, 100, good.coverage, Inf)
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
            min.coverage <- min(key.coverage)
            max.coverage <- max(df2$coverage)
            off.score <- as.character(cut(
                    max.coverage,
                    c(-Inf, 0, 10, 100, Inf),
                    labels=c(0, -1, -2, -3)))
            on.score <- as.character(cut(
                    min.coverage,
                    coverage.levels,
                    labels=c(1, 2, 3, 4)))
            output <- rbind(output, c(
                            sample, 
                            region, 
                            q.cut, 
                            min.coverage, 
                            temp[which.min(key.coverage)],
                            off.score,
                            on.score))
        } else {
            # TODO: Is this useful? Should we only output entries with data?
            output <- rbind(output, c(sample, region, q.cut, NA, NA, 0, 1))
        }
    }
    
    # draw required coverage line
    abline(h = good.coverage, lty=2)
    
    axis(1)
    axis(2, at=c(1E1, 1E2, 1E3, 1E4, 1E5), labels=c('10', '100', '1000', '10,000', '100,000'), las=2)
    box()
    
    legend(x=0, y=1, legend=paste('q=', cutoffs, sep=''), col=rainbow(length(cutoffs), v=0.8), lty=1, lwd=2, bty='n', yjust=0)
    text(x=max(df$refseq.aa.pos)/2, y=190000, label=paste(sample, region), cex=0.7, col='grey30', adj=c(0.5, 0.5))
    
    garbage <- dev.off()
}


# now handle nucleotide coverage for HLA-B (3,337 bp)
region <- 'HLA-B'
data <- read.csv(file=nuc.csv, header=TRUE, sep=',')
data <- data[data$region == 'HLA-B', ]
if (length(data$query.nuc.pos)) {
    alphabet <- c('A', 'C', 'G', 'T')
    sample <- data$sample[1]
    
    data$coverage <- apply(data[ , which(is.element(names(data), alphabet))], 1, sum)
    data$coverage <- sapply(data$coverage, function(x) max(x, 0.1))
    
    good.coverage <- 100
    coverage.levels <- c(-Inf, 0, 10, good.coverage, Inf)
            
    filename <- file.path('coverage_maps', paste(sample, region, 'png', sep='.'))
    
    # set up plot
    png(file=filename, width=400, height=300, type='cairo')
    par(family='sans', cex=1, mar=c(5,5,1,1))
    plot(NA, xlim=c(1,max(data$query.nuc.pos)), ylim=c(1,200000), axes=FALSE, ann=FALSE, xaxs="r", log="y")
    title(xlab="Reference coordinates (bp)", font.lab = 1.4, cex.lab=1.4, cex.main=1.4)
    
    cutoffs <- sort(as.integer(levels(as.factor(data$q.cutoff))))
    for (j in 1:length(cutoffs)) {
        q.cut <- cutoffs[j]
        df2 <- data[data$q.cutoff == q.cut, ]
        lines(x = df2$query.nuc.pos, 
                y = df2$coverage, 
                col=rainbow(length(cutoffs), v=0.8)[j],
                lwd=2
        )
    }
    
    abline(h = good.coverage, lty=2)
    
    axis(1)
    axis(2, at=c(1E1, 1E2, 1E3, 1E4, 1E5), labels=c('10', '100', '1000', '10,000', '100,000'), las=2)
    box()
        
    legend(x=0, y=1, legend=paste('q=', cutoffs, sep=''), col=rainbow(length(cutoffs), v=0.8), lty=1, lwd=2, bty='n', yjust=0)
    text(x=max(data$query.nuc.pos)/2, y=190000, label=paste(sample, region), cex=0.7, col='grey30', adj=c(0.5, 0.5))
        
    garbage <- dev.off()
}

# package outputs
tar(output.maps, 'coverage_maps')


# output minimum coverage at key positions
write.csv(as.data.frame(output), file=output.csv, quote=FALSE, row.names=FALSE)
