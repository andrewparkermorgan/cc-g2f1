setwd("~/Dropbox/G2F1/Andrew/")

library(GenomicRanges)
library(ggplot2)

# load("./peaks.bs.Rdata")
cold.df <- read.table("~/Dropbox/G2F1/Andrew/cold.regions_oldmm9.csv", header = TRUE, sep = ",")
cold.df$seqnames <- factor( as.character(cold.df$seqnames), levels = paste("chr", c(1:19,"X"), sep = "") )
cold.df$Chr <- factor(as.numeric(cold.df$seqnames), levels = as.character(1:20))

cold.gr <- with(cold.df, GRanges(seqnames = Rle(Chr),
																 ranges = IRanges(start = start, end = end)))

bs <- daply(peaks.bs, .(.id), function(df) {
	gr <- with(df, GRanges(seqnames = Rle(Chr),
												 ranges = IRanges(start = Pos.start, end = Pos.end),
												 E = E))
	ifelse(countOverlaps(cold.gr, gr) > 0, 1, 0)
}, .progress = "text")

cold.df$support <- apply(bs, 2, sum)/nrow(bs)

cr <- cold.df
ggplot(cr) + 
	stat_ecdf(aes(x = support), colour = "#0099ff") + 
	xlab("\nbootstrap support level") + ylab("cumulative proportion\n") +
	theme_bw()
qqmath(~ support, data = cr, distribution = qunif,
			 type = "b", pch = 19, aspect = "iso",
			 xlab = "theoretical quantiles (uniform distribution)", ylab = "bootstrap support quantiles",
			 panel = function(...) {
			 	panel.qqmathline(..., col = "grey70", lty = "dashed")
			 	panel.qqmath(..., col = "#0099ff")
			 })