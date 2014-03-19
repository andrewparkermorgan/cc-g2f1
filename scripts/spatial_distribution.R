#############################################
## spatial_distribution.R
#
# purpose:	analyze spatial distribution of recombination in CC G2:F1 and other crosses
# author:		APM 1 Jan 2013
#	updated:	APM 24 Sep 2013
#
#############################################

## R libraries
library(GenomicRanges)
library(ggplot2)
library(scales)
library(plyr)

load("seqlengths.mm10.Rdata")

## function to clean up chromosome names and add decoration to enable later plots
rescale.chroms <- function(df, seqlengths, SHORT.CHROMS, ...) {
	
	df$seqnames <- with(df, factor(as.character(seqnames), levels = names(seqlengths) ) )
	df$chr <- with(df, chr <- factor( gsub("chr","", as.character(seqnames)), levels = SHORT.CHROMS[1:20] ) )
	df$ypos <- rnorm(nrow(df))
	df$xpos.start <- with(df, max(seqlengths) - (seqlengths[ as.character(seqnames) ] - start) )
	df$xpos.end <- with(df, max(seqlengths) - (seqlengths[ as.character(seqnames) ] - end) )
	df$norm.start <- with(df, start/seqlengths[ as.character(seqnames) ])
	df$norm.end <- with(df, end/seqlengths[ as.character(seqnames) ])
	
	return(df)

}

## some local constants
SHORT.CHROMS <- c(1:19,"X","Y","M")
CHROMS <- paste("chr", SHORT.CHROMS, sep = "")
SEX.COLORS <- setNames( c("#ff00ff","#0080ff"), c("female","male") )

## pointers to data files (assuming all supplementary files were unzipped into current working directory)
events.file <- "FileS1.csv"
ix.events.file <- "FileS2.csv"
cold.file <- "FileS4.csv"
dsb.file <- "FileS5.csv"
segdup.file <- "FileS6.txt"

## double-strand breaks from Brick et al. 2012
dsb.df <- read.table(dsb.file, header = TRUE, sep = ",")
dsb.df <- rescale.chroms(dsb.df, seqlengths, SHORT.CHROMS)
dsb.gr <- with(dsb.df, GRanges(seqnames = Rle(seqnames), seqlengths = seqlengths,
															 ranges = IRanges(start = start, end = end),
															 strain = strain, intens = strength))

## segmental duplication (SD) density computed with Gepard
segdup.df <- read.table(segdup.file, header = FALSE)
names(segdup.df) <- c("seqnames","start","end","duprate")
segdup.df$seqnames <- factor( names(seqlengths)[ as.integer(segdup.df$seqnames) ], levels = names(seqlengths) )
segdup.gr <- with(segdup.df, GRanges(seqnames = Rle(seqnames), seqlengths = seqlengths,
																		 ranges = IRanges(start = start, end = end),
																		 duprate = duprate))

## cold regions from G2F1
cold.df <- read.table(cold.file, header = TRUE, sep = ",")
cold.df$seqnames <- factor(as.character(cold.df$seqnames), levels = CHROMS, ordered = TRUE)
cold.df$gc <- as.numeric(as.character(cold.df$gc))
cold.gr <- with(cold.df, GRanges(seqnames = Rle(seqnames), seqlengths = seqlengths,
																 ranges = IRanges(start = start, end = end),
																 sd = sd, count.genes = count.genes))

## recombination events in G2F1
events.df <- read.table(events.file, header = TRUE, sep = ",")
events.df$seqnames <- factor( as.character(events.df$seqnames) , levels = CHROMS)
events.df <- transform( events.df, key = paste(seqnames, ":", start, "-", end, sep = "") )
events.df <- rescale.chroms(events.df, seqlengths, SHORT.CHROMS)

## load intercross recombination events (see Supplementary Text)
ix.events.df <- read.table(ix.events.file, header = TRUE, sep = ",")

## make augmented events table -- includes both G2 and intercross events, after scaling chromosomes to unit length
keep.cols <- c("seqnames","start","end","meiosis.sex","id","norm.start","norm.end")
aug.ev <- rescale.chroms(ix.events.df, seqlengths, SHORT.CHROMS)[ ,keep.cols ]
aug.ev$type <- "P" # all of these are paternal meioses
aug.ev$group <- "intercross meioses"
aug.ev <- rbind( aug.ev,
								 transform(subset(events.df, type %in% c("M","P")),
								 					group = "G2 meioses")[ ,c("seqnames","start","end","meiosis.sex","id","type","group","norm.start","norm.end") ] )

## count how many recombs per chromosome
aug.ev <- ddply(aug.ev, .(group, type, seqnames, id), transform, n.recomb = length(start),
								.progress = "text")
aug.ev <- transform(aug.ev, recomb.group = ifelse(n.recomb > 1, "mulitple-recombinant", "single-recombinant"))
aug.ev$type <- factor(aug.ev$type)
aug.ev$group <- factor(aug.ev$group)
aug.ev$recomb.group <- factor(aug.ev$recomb.group)


#############################
## First-pass analysis of spatial distribution of recombination at chromosome scale

## Figure 4
# raw G2F1 recombination data (left panel)
p1 <- ggplot(events.df) +
	geom_point(aes(x = (xpos.start+xpos.end)/(2*1e6), y = ypos, colour = meiosis.sex), size = 0.45) +
	# geom_segment(aes(x = start/(1e6), xend = end/(1e6), y = y.pos, yend = y.pos ), colour = SEX.COLORS[sex]) +
	facet_grid(chr ~ meiosis.sex, drop = FALSE) +
	scale_y_continuous(breaks = NULL) +
	scale_x_continuous(labels = NULL) +
	scale_colour_manual(limits = levels(events.df$meiosis.sex), values = SEX.COLORS) +
	guides(colour = FALSE) +
	xlab("genomic coordinate") + ylab("chromosome") +
	theme_bw()
plot(p1)
# kernel density estimates on recombination data (right panel)
p2 <- ggplot(events.df) +
	# geom_line(data = dsb.df, aes(x = (norm.start+norm.end)/2), stat = "density", colour = "grey", lty = "dashed", lwd = 1) +
	geom_line(aes(x = (norm.start+norm.end)/2, colour = meiosis.sex), stat = "density", lwd = 1) +
	# facet_grid(meiosis.sex ~ .) +
	scale_colour_manual(limits = levels(events.df$meiosis.sex), values = SEX.COLORS) +
	guides(colour = FALSE) +
	xlab("\nnormalised chromosomal coordinate") +
	theme_bw()
plot(p2)

## Figure 5
ggplot(aug.ev) +
	# geom_line(data = dsb.tmp, aes(x = (norm.start+norm.end)/2), stat = "density", colour = "grey", lwd = 1) +
	geom_line(aes(x = (norm.start+norm.end)/2, colour = meiosis.sex), stat = "density", lwd = 1, kernel = "rectangular") +
	facet_grid(group + meiosis.sex ~ recomb.group) +
	scale_colour_manual(limits = levels(events.df$meiosis.sex), values = SEX.COLORS) +
	guides(colour = FALSE) +
	xlab("\nnormalised chromosomal coordinate") +
	theme_bw()


#############################
## Check for enrichment for male or female recombination in bins across genome

## events.to.bins()
# assign events to bins across genome
# <d> = dataframe with columns 'seqnames','pos' and a column with cumulative event count
# <size> = bin size in bp
# <seqlenths> = NAMED vector of chromosome lengths
# <col> = name of column containing cumulative event count
events.to.bins <- function(d, size, seqlengths, col = "cumsum") {
	right.end <- seqlengths[ as.character(d$seqnames[1]) ]
	bins <- seq(1, right.end, size)
	d$bin <- findInterval(d$pos, bins)
	d$max.event <- max(d[,col])
	return(d)
}

## test.sex.enrichment()
# test for enrichment via 1df chi-squared
# <df> = dataframe with columns 'events','max.event','meiosis.sex'
test.sex.enrichment <- function(df, ...) {
	
	events.by.sex <- with(df, tapply(events, meiosis.sex, sum))
	marginal.by.sex <- with(df, tapply(max.event, meiosis.sex, max))
	# print(events.by.sex)
	# print( names(events.by.sex)[ which.max(events.by.sex) ] )
	if( any(sum(is.na(events.by.sex)), sum(is.na(marginal.by.sex)), !sum(events.by.sex)) ) {
		# return(data.frame(p.value = NA, high.sex = NA))
	}
	else {
		return( data.frame(p.value = chisq.test(events.by.sex, p = marginal.by.sex/sum(marginal.by.sex))$p.value,
											 high.sex = names(events.by.sex)[ which.max(events.by.sex) ]) )
	}
	
}

## first convert events to cumulative events
g2f1 <- ddply(events.df, .(seqnames, meiosis.sex), integrate.events)

## Figure S5
## compute chi-squared p-value on enrichment for one sex or another in 1Mbp windows
# conditional on per-chromosome difference in recomb rates
bin.size <- 1e6
maps <- list(G2F1 = g2f1)
for (i in 1:length(maps)) {
	tmp <- ddply(na.omit(subset(maps[[i]], !(seqnames %in% c("chrX")))),
							 .(seqnames, meiosis.sex),
							 events.to.bins, size = bin.size, seqlengths = seqlengths)
	tmp.by.bin <- ddply(tmp,
											.(seqnames, meiosis.sex, bin),
											summarize, max.event = max(max.event), events = diff(range(cumsum)))
	tmp.by.sex <- ddply(tmp.by.bin,
											.(seqnames, bin),
											test.sex.enrichment)
	tmp.by.sex$log.p <- with( tmp.by.sex, ifelse(high.sex == "female", log10(p.value), -1*log10(p.value)) )
	
	p <- ggplot(tmp.by.sex) +
		geom_segment(aes(x = bin*bin.size/1e6, xend = bin*bin.size/1e6, y = 0, yend = log.p, colour = high.sex)) +
		facet_wrap(~ seqnames, ncol = 4) +
		scale_colour_manual(limits = levels(events.df$meiosis.sex), values = SEX.COLORS) +
		scale_y_continuous(labels = abs) +
		guides(colour = guide_legend("enriched sex")) +
		xlab("coordinate (Mbp)") + ylab("-log10( p-value )") + ggtitle(paste("test for local heterochiasmy in", names(maps)[i],"map"))
	plot(p)
	
}

#############################
## Recombination rate vs segmental duplication (SD) content

## use GenomicRanges for overlap calculations, which is lazy but not too slow
events.gr <- with(events.df, GRanges(seqnames = Rle(seqnames), seqlengths = seqlengths,
																		 ranges = IRanges(start = start, end = end),
																		 proximal.strain = proximal.strain, distal.strain = distal.strain,
																		 id = id, funnel = funnel, meiosis.sex = meiosis.sex, type = type))
events.per.bin <- countOverlaps( segdup.gr, events.gr )

## scale event counts to cM/Mbp units
values(segdup.gr)$cMperMbp <- (100*(events.per.bin/(nlevels(events.df$funnel)*8)) / (7/8)) / (width(segdup.gr)/1e6)

## assign intervals to bins according to SD content
sd.df <- as.data.frame(segdup.gr)
sd.df <- transform( sd.df, dup.bins = cut(sd.df$duprate, breaks = c(seq(0, 1, 0.1), Inf), right = FALSE) )
intervals.per.bin <- as.vector(table(sd.df$dup.bins))
sd.df <- transform( sd.df, nintervals = intervals.per.bin[ as.integer(dup.bins) ] )

## Figure S5
ggplot(sd.df) +
	geom_boxplot(aes(y = cMperMbp, x = dup.bins, fill = nintervals)) +
	xlab("\nproportion of sequence covered by a segmental duplication larger than 20 kbp") + ylab("local recombination rate (cM per Mbp)\n") +
	guides(fill = guide_colourbar("number of intervals")) + scale_fill_gradient(trans = "log10") +
	theme_bw()