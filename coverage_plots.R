#!/usr/bin/env Rscript

# Command line expects these arguments:
# 1. input file name (amino frequencies CSV)
# 2. output maps archive (tar file where coverage maps will be written)
# 3. output scores (count scores CSV)

# Other dependency: projects.json
# Working files are created in the current directory.

library(jsonlite)

syntax.error = "Correct syntax: ./coverage_plots.R <amino CSV> <output maps tar> <output scores CSV>"
args <- commandArgs(TRUE)

if (length(args) != 3) {
    stop(syntax.error)
}
input.csv <- args[1]
output.maps <- args[2]
scores.csv <- args[3]

record.score <- function(
        scores,
        data,
        region.key.pos,
        coverage.levels,
        sample, 
        project.name,
        region,
        seed,
        q.cut) {
    # Add a row of scores to the scores matrix.
    #
    # scores - a matrix that holds data to write out to the scores file.
    # data - A data frame with the following columns:
    #   refseq.aa.pos - the position in the amino acid reference sequence
    #   query.nuc.pos - the position in the nucleotide reference sequence
    #   coverage - the coverage at that position by combining all reads
    #   If the amino acid position column is empty, the nucleotide column
    #   will be used instead.
    # region.key.pos - a vector of key positions within the region. If it is
    #   empty, treat all positions as key.
    # coverage.levels - a vector of four coverage levels that are thresholds
    #   for the four different scores in on-target regions.
    # sample - the sample name
    # project.name - the project name
    # region - the coordinate region name
    # seed - the seed region name
    # q.cut - the quality cutoff used for this data
    # Returns scores with a new row added.
    all.coverage <- data$coverage
    if (length(region.key.pos) == 0) {
        key.coverage <- all.coverage
        min.pos <- which.min(key.coverage)
    }
    else {
        if (length(data$refseq.aa.pos) != 0) {
            positions <- data$refseq.aa.pos
        }
        else {
            positions <- data$query.nuc.pos
        }
        key.coverage <- all.coverage[is.element(
                        positions,
                        region.key.pos)]
        min.pos <- region.key.pos[which.min(key.coverage)]
    }
    
    if (length(key.coverage) > 0) {
        min.coverage <- min(key.coverage)
        max.coverage <- max(all.coverage)
        off.score <- as.character(cut(
                max.coverage,
                c(-Inf, 0, 10, 100, Inf),
                labels=c(0, -1, -2, -3)))
        on.score <- as.character(cut(
                min.coverage,
                coverage.levels,
                labels=c(1, 2, 3, 4)))
        scores <- rbind(scores, c(
                        sample,
                        project.name,
                        region,
                        seed,
                        q.cut, 
                        min.coverage, 
                        min.pos,
                        off.score,
                        on.score))
    } else {
        scores <- rbind(
                scores,
                c(sample, project.name, region, q.cut, NA, NA, 0, 1))
    }
}

prepare.plot <- function(
        xlim,
        x.label,
        good.coverage,
        sample,
        project.name,
        region,
        region.key.pos) {
    # Draw the frame of a plot, along with the key positions.
    #
    # xlim - the maximum value of the x axis
    # x.label - the label text for the x axis
    # good.coverage - the minimum coverage that is considered good
    # sample - sample name
    # project.name - project name that defines the key positions
    # region - region name
    # region.key.pos - vector holding all key positions in the region
    filename <- file.path('coverage_maps', paste(
                    sample,
                    project.name,
                    region,
                    'png',
                    sep='.'))
    
    # set up plot
    png(file=filename, width=400, height=300, type='cairo')
    par(family='sans', cex=1, mar=c(5,5,1,1))
    plot(NA, xlim=c(1,xlim), ylim=c(1,2000000), axes=FALSE, ann=FALSE, xaxs="r", log="y")
    title(xlab=x.label, font.lab = 1.4, cex.lab=1.4, cex.main=1.4)
    
    # indicate key positions for this region
    for (pos in region.key.pos) {
        rect(
                pos-0.505, good.coverage/2,
                pos+0.505, good.coverage*2,
                col='grey30',
                border=NA)
    }
}

dir.create('coverage_maps')

# load project configuration
projects.config <- fromJSON('projects.json')
projects <- projects.config$projects

coverage.levels <- data.frame()
key.positions <- data.frame()  # { pos, project, region.index }
region.indexes <- data.frame() # { seed, coord, region.index }
region.index <- 0
for (project.name in names(projects)) {
    project <- projects[[project.name]]
    region.count <- nrow(project$regions)
    if (is.null(region.count)) {
        region.count <- 0
    }
    for (i in seq_len(region.count)) {
        region <- project$regions[i,]
        seed_regions <- region$seed_region_names[[1]]
        for (seed.index in seq_along(seed_regions)) {
            new.coverage.levels <- data.frame(
                    project=project.name,
                    seed=seed_regions[[seed.index]],
                    coord=region$coordinate_region,
                    red=region$min_coverage1,
                    yellow=region$min_coverage2,
                    green=region$min_coverage3)
            coverage.levels <- rbind(coverage.levels, new.coverage.levels)
            
            new.region.indexes <- data.frame(
                    coord=region$coordinate_region,
                    seed=seed_regions[[seed.index]],
                    region.index=region.index)
            region.indexes <- rbind(region.indexes, new.region.indexes)
        }
        region.pairs <- region$key_positions[[1]]
        if (length(region.pairs) != 0) {
            region.pairs$end_pos[is.na(region.pairs$end_pos)] <- (
                        region.pairs$start_pos[is.na(region.pairs$end_pos)])
        }
        else {
            ref <- projects.config$regions[[region$coordinate_region]]$reference
            region.pairs = data.frame(start_pos=1, end_pos=sum(nchar(ref)))
        }
        for (j in seq_along(region.pairs$start_pos)) {
            new.positions <- data.frame(
                    pos=region.pairs$start_pos[[j]]:region.pairs$end_pos[[j]],
                    project=project.name,
                    region.index=region.index)
            key.positions <- rbind(key.positions, new.positions)
        }
        region.index <- region.index + 1
    }
}

scores.columns <- c(
        'sample',
        'project',
        'region',
        'seed',
        'q.cut',
        'min.coverage',
        'which.key.pos',
        'off.score',
        'on.score')
scores <- matrix(nrow=0, ncol=length(scores.columns), dimnames=list(NULL, scores.columns))

data <- read.csv(file=input.csv, header=TRUE, sep=',')

# the input CSV should look like this:
#sample,seed,region,q-cutoff,query.aa.pos,refseq.aa.pos,A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*
#50955ARPT-HCV-49537A-INT-PR-RT_S16,HCV1A-H77-core-seed,HCV1A-H77-core,0,0,1,0,0,0,0,0,0,0,3,0,0,1542,0,0,0,0,0,0,0,0,0,0

alphabet <- c('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y')  # X. represents stop codon '*'

data$coverage <- apply(data[ , which(is.element(names(data), alphabet))], 1, sum)

# partition AA frequency table by sample and region
coverage <- split(data[,which(is.element(names(data), c('refseq.aa.pos', 'q.cutoff', 'coverage')))], f=list(data$seed, data$region, data$sample), drop=TRUE)

# loop through partitions
for (i in seq_along(coverage)) {
    label <- names(coverage)[i]
    tokens <- strsplit(label, split='\\.')[[1]]
    seed <- tokens[1]
    region <- tokens[2]
    sample <- tokens[3]
    
    df <- coverage[[i]]
    xlim <- max(df$refseq.aa.pos)
    x.label <- "Reference coordinates (AA)"
    related.regions <- region.indexes[
            region.indexes$coord == region & region.indexes$seed == seed,]
    related.positions <- key.positions[
            key.positions$region.index %in% related.regions$region.index,]
    project.positions <- split(
            related.positions$pos,
            f=list(related.positions$project),
            drop=TRUE)
    
    for (project.name in names(project.positions)) {
        region.coverage.levels <- coverage.levels[
                coverage.levels$project == project.name & coverage.levels$seed == seed & coverage.levels$coord == region,]
        score.breaks <- c(
                -Inf,
                region.coverage.levels$red,
                region.coverage.levels$yellow,
                region.coverage.levels$green,
                Inf)
        prepare.plot(
                xlim,
                x.label,
                region.coverage.levels$green,
                sample,
                project.name,
                region,
                project.positions[[project.name]])
        
        # draw trend lines for each quality score cutoff
        cutoffs <- sort(as.integer(levels(as.factor(df$q.cutoff))))
        for (j in 1:length(cutoffs)) {
            q.cut <- cutoffs[j]
            df2 <- df[df$q.cutoff == q.cut, ]
            
            # zeroes don't show up on a log plot, so we plot no coverage at 0.1
            lines(x = df2$refseq.aa.pos, 
                    y = sapply(df2$coverage, function(x) max(x, 0.1)), 
                    col=rainbow(length(cutoffs), v=0.8)[j],
                    lwd=2
            )
            
            # determine minimum coverage at key positions
            scores <- record.score(
                    scores,
                    df2,
                    project.positions[[project.name]],
                    score.breaks,
                    sample,
                    project.name,
                    region,
                    seed,
                    q.cut)
        }
        
        # draw required coverage line
        abline(h = region.coverage.levels$green, lty=2)
        
        axis(1)
        axis(2, at=c(1E1, 1E2, 1E3, 1E4, 1E5, 1E6), labels=c('10', '100', '1000', '10,000', '100,000', '1,000,000'), las=2)
        box()
        
        text(x=max(df$refseq.aa.pos)/2, y=1900000, label=paste(sample, region), cex=0.7, col='grey30', adj=c(0.5, 0.5))
        
        garbage <- dev.off()
    }
}


# package outputs
tar(output.maps, 'coverage_maps')


# output minimum coverage at key positions
write.csv(as.data.frame(scores), file=scores.csv, quote=FALSE, row.names=FALSE)
