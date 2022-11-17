library(tidyverse)
library(devtools)
library(data.table)
library(GenomicRanges)
library(here)

load('../../../LUCAS/rlucas/data/lucas_5mb.rda')
bins <- unique(GRanges(bins5mb %>% select(chr, start, end, arm, bin)))


#--------------------------------------------------------#
get_arm_ranges <- function(assembly){

	mySession <- rtracklayer::browserSession()
	GenomeInfoDb::genome(mySession) <- assembly

	gaps <- rtracklayer::getTable(rtracklayer::ucscTableQuery(mySession, track="cytoBand"))
	gaps$arm <- substr(gaps$name, 1, 1)
	gaps$id <- apply(gaps[,c('chrom', 'arm')], 1, function(x) gsub('chr', '', paste0(x[1], x[2])))

	arm.ranges <- plyr::ddply(gaps, plyr::.(id), plyr::summarize, chrom = unique(chrom), start = min(chromStart), end = max(chromEnd))
	clean.arm.ranges <- subset(arm.ranges,! id %in% c('13p', '14p', '15p', '21p', '22p', 'Xp', 'Xq', 'Yp', 'Yq'))
	clean.arm.ranges <- clean.arm.ranges[,c('chrom', 'start', 'end', 'id')]
	colnames(clean.arm.ranges) <- c('seqnames', 'start', 'end', 'arm')
	return(sort(GRanges(clean.arm.ranges)))
}


find_covered_intervals <- function(intervals, data, min.coverage){
	o <- findOverlaps(intervals, data)
	ow <- width(pintersect(intervals[queryHits(o)], data[subjectHits(o)]))
	base <- intervals[queryHits(o)]
	base$covered.fraction <- ow / width(intervals[queryHits(o)])
	total.feature.coverage <- unlist(lapply(split(base, f = base$bin), FUN = function(x) sum(x$covered.fraction)))
	pos <- intervals[intervals$bin %in% names(which(total.feature.coverage >= min.coverage))]
	return(pos)
}

find_overlapping_genes <- function(ranges, genes, within = 0){
	if (within == 1){
		hits <- unique(as.character(subsetByOverlaps(genes, ranges, type = 'within')$gene_name))
	}
	if (within == 0){
		hits <- unique(as.character(subsetByOverlaps(genes, ranges, type = 'any')$gene_name))
	}
	return(paste(hits, collapse = ','))
}


vectorize <- function(bins, intervals){
  pos.ids <- find_covered_intervals(bins, intervals, 0.9)$bin
  out <- rep(0, length(bins))
  out[bins$bin %in% pos.ids] <- 1
  names(out) <- bins$bin
  return(out)
}
#--------------------------------------------------------#
# determine the cn profile of liver tumors
# threshold to call level changes is based on Davoli et al., Science, supp table S1
date()
liver.thresh <- 0.30

liver <- fread("../data/TCGA/LIHC.txt")
liver$gain <- ifelse(liver$Segment_Mean > liver.thresh, 1, 0)
liver$loss <- ifelse(liver$Segment_Mean < (-1) * liver.thresh, 1, 0)

liver <- split(liver, f = liver$Sample)
liver <- lapply(liver, function(x) {y = x
                                  y$Sample = NULL
                                  y <- subset(y, Chromosome %in% seq(1,22))
                                  y$Chromosome <- sapply(y$Chromosome, function(x) paste0('chr',x))
                                  return(y)})
liver <- lapply(liver, GRanges)

liver.losses <- lapply(liver, function(x) vectorize(bins, subset(x, loss == 1)))
liver.losses <- do.call(rbind, liver.losses)
aliquot <- sapply(rownames(liver.losses), function(x) strsplit(as.character(x), split = '-')[[1]][4])
liver.losses <- liver.losses[aliquot %in% c('01A', '01B', '02A', '02B'),]
liver.loss <- apply(liver.losses, 2, mean)

# 372 liver cases
liver.gains <- lapply(liver, function(x) vectorize(bins, subset(x, gain == 1)))
liver.gains <- do.call(rbind, liver.gains)
aliquot <- sapply(rownames(liver.gains), function(x) strsplit(as.character(x), split = '-')[[1]][4])
liver.gains <- liver.gains[aliquot %in% c('01A', '01B', '02A', '02B'),]
liver.gain <- apply(liver.gains, 2, mean)
#--------------------------------------------------------#

bins$`liver Loss` = liver.loss
bins$`liver Gain` = liver.gain


bins <- data.frame(bins)
bins$pos = apply(bins[,c('start', 'end')], 1, mean)
bins[,c('start', 'end', 'width', 'strand')] = NULL
b = melt(bins, id.vars = c('seqnames', 'pos', 'arm' ,'bin'))
b$arm <- factor(b$arm, levels = unique(bins$arm))
b$disease <- sapply(b$variable, function(x) strsplit(as.character(x), split = '[.]')[[1]][1])
b$change <- sapply(b$variable, function(x) strsplit(as.character(x), split = '[.]')[[1]][2])

b[which(b$change == 'Loss'), 'value'] = b[which(b$change == 'Loss'), 'value'] * (-1)
#--------------------------------------------------------#

fig.data <- b
saveRDS(fig.data, file = '../data/TCGA/fig2c_p2_data.rds')
