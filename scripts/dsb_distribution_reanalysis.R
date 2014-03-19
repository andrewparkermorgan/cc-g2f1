#############################################
## dsb_distribution_reanalysis.R
#
# purpose:	compare spatial distribution of DSBs (Smagulova et al. 2011; Brick et al. 2012) with that of G2:F1 recombination events
# author:		APM 22 Dec 2013
#
#############################################
## to set up environment, run <./spatial_distribution.R>

library(quantreg)

## is distribution of recombination events uniform? (by sex, # recombs per chromosome)
alpha <- 0.05
tmp <- ddply(aug.ev, .(seqnames, meiosis.sex, recomb.group, group), function(d) {
	#d <- subset(d, norm.start > 0.2 & norm.end < 0.8)
	mids <- (d$start + d$end)/2
	data.frame(p.value = ks.test(mids, punif, min = min(mids), max = max(mids))$p.value)
})
tmp <- transform(tmp, p.adj = p.adjust(p.value, method = "BH"))
xtabs(~(p.adj < 0.05)+meiosis.sex+recomb.group, subset(tmp, group == "G2 meioses" & seqnames != "chrX"))
ggplot(subset(tmp, group == "G2 meioses" & seqnames != "chrX")) +
	geom_histogram(aes(x = p.adj, y = ..ncount.., fill = recomb.group),
								 position = "dodge") +
	facet_grid(meiosis.sex ~ .)

## is distribution of DSBs different from that of recomb events in males?
tmp <- ddply(dsb.df, .(chr,seqnames), function(d) {
	recomb.mids <- with( subset(events.df, meiosis.sex == "male" & seqnames == d$seqnames[1]), (start+end)/2 )
	p.sim <- replicate(1000, { 
		dsb.mids <- sample( (d$start + d$end)/2, length(recomb.mids), replace = TRUE, prob = d$strength )
		ks.test(recomb.mids, dsb.mids)$p.value
	})
	data.frame(p.value = median(p.sim, na.rm = TRUE))
}, .progress = "text")
tmp <- transform(tmp, p.adj = p.adjust(p.value, method = "BH"))
xtabs(~(p.adj < 0.05), subset(tmp, seqnames != "chrX"))
ggplot(subset(dsb.df, !(seqnames %in% c("chrX","chrY","chrM")))) +
	stat_ecdf(aes(x = (norm.start + norm.end)/2, colour = strain)) +
	facet_wrap(~ chr, ncol = 4) +
	theme_bw()

## recombination events ECDFs by chromosome
ggplot(subset(events.df, chr != "X")) +
	stat_ecdf(aes(x = (norm.start+norm.end)/2, colour = meiosis.sex)) +
	facet_wrap(~ chr, ncol = 4) +
	xlab("\nnormalized coordinate") + ylab("cumulative probability\n") +
	scale_colour_manual("sex", limits = levels(events.df$meiosis.sex), values = SEX.COLORS) +
	theme_bw()

## recombination events ECDFs by chromosome plus DSBs
pvlab <- function(p) {
	logp <- -log10(p)
	round(logp, 1)
}
ggplot(subset(events.df, chr != "X")) +
	stat_ecdf(aes(x = (norm.start+norm.end)/2, colour = meiosis.sex)) +
	stat_ecdf(data = subset(dsb.df, chr != "X"), aes(x = (norm.start+norm.end)/2), colour = "grey80") +
	geom_text(data = subset(tmp, chr != "X"), aes(x = 0.5, y = 0.1, label = pvlab(p.adj))) +
	facet_wrap(~ chr, ncol = 4) +
	xlab("\nnormalized coordinate") + ylab("cumulative probability\n") +
	scale_colour_manual("sex", limits = levels(events.df$meiosis.sex), values = SEX.COLORS) +
	theme_bw()


## use quantile regression to demonstrate that male event distribution is shifted at every quantile
tmp <- ddply(subset(events.df, chr != "X"), .(chr), function(d) {
	cat("\n\n", as.character(d$seqnames[1]), "\n")
	X <- model.matrix(~ 0 + meiosis.sex, d)
	X.alt <- model.matrix(~ meiosis.sex, d)
	y <- with(d, (start+end)/2)
	# m0 <- rq.fit.br(X, y, tau = 0.8, alpha = 0.05, ci = TRUE)
	ldply(seq(0.2, 0.8, 0.05), function(tau) {
		m1 <- rq.fit.br(X.alt, y, tau = tau, alpha = 0.05, ci = TRUE)
		data.frame(tau = tau, sex = rownames(m1$coefficients), m1$coefficients)[2,]
	})
})
# colnames(tmp) <- make.names(colnames(tmp))
# ggplot(tmp) +
# 	geom_pointrange(aes(x = chr, y = coefficients, ymin = lower.bd, ymax = upper.bd, colour = sex)) +
# 	coord_flip() + theme_bw()
ggplot(tmp) +
	geom_hline(yintercept = 0, col = "grey70", lty = "dashed") +
	geom_pointrange(aes(x = tau, y = coefficients/1e6, ymin = lower.bd/1e6, ymax = upper.bd/1e6)) +
	facet_wrap(~chr, ncol = 4) +
	xlab("\nquantile of chromosomal position") + ylab("male-female delta (Mbp)\n") +
	theme_bw() + coord_flip()

df2gr <- function(df, seqnames = "seqnames", start = "start", end = "end", seqlengths = NULL, ...) {
	meta.cols <- setdiff(colnames(df), c(seqnames, start, end))
	gr <- GRanges( seqnames = Rle(df[ ,seqnames ]),
								 ranges = IRanges(start = df[ ,start ], end = df[ ,end ]) )
	
	values(gr) <- df[ ,meta.cols ]
	if (!missing(seqlengths)) {
		seqlengths(gr) <- seqlengths
	}
	
	return(gr)
}

## frequency distribution of DSBs in 500kbp windows
window <- 500e3
step <- 250e3
bins <- ldply(CHROMS, function(x) {
	len <- seqlengths[x]
	starts <- seq(1, len, step)
	ends <- starts + window - 1
	if (max(ends) > len) {
		ends[ which.max(ends) ] <- len
	}
	data.frame(seqnames = x, start = starts, end = ends)
})

## label bins as being "in a cold region" if >50% of the bin is in a cold region
bins.gr <- df2gr(bins)
dsb.gr <- df2gr(dsb.df)
olap <- findOverlaps(dsb.gr, bins.gr)
cold.olap <- countOverlaps(bins.gr, cold.gr, minoverlap = 250e3)
values(bins.gr)$is.cold <- (cold.olap > 0)
dsbs.per.bin <- ddply(as.data.frame(olap), .(subjectHits), function(d) {
	data.frame(reads = sum(values(dsb.gr)$strength[ d$queryHits ], na.rm = TRUE))
}, .progress = "text")
dsbs.per.bin$is.cold <- values(bins.gr)$is.cold[ dsbs.per.bin$subjectHits ]
dsbs.per.bin$coldness <- factor( ifelse(dsbs.per.bin$is.cold, "cold region","background") )

## what is SD rate in these bins?
sd.olap <- findOverlaps(segdup.gr, bins.gr)
sd.per.bin <- ddply(as.data.frame(sd.olap), .(subjectHits), function(d) {
	data.frame(sd = mean(values(segdup.gr)$duprate[ d$queryHits ], na.rm = TRUE))
}, .progress = "text")
dsbs.per.bin <- merge(dsbs.per.bin, sd.per.bin)
dsbs.per.bin$is.sd <- (dsbs.per.bin$sd > 0.5)

## plot frequency distribution of DSBs in cold regions versus not
ggplot(dsbs.per.bin) +
	geom_histogram(aes(x = log10(reads), y = ..ncount.., fill = coldness), position = "dodge") +
	# facet_grid(. ~ is.sd) +
	xlab("\nlog10( DSB tag reads )") + ylab("relative frequency\n") +
	scale_fill_brewer("", palette = "Paired") + scale_y_continuous(labels = function(x) paste(100*x, "%")) +
	theme_bw()

## rank-sum test for difference in location: more or less DSBs in cold regions?
wilcox.test(log10(reads) ~ is.cold, dsbs.per.bin)