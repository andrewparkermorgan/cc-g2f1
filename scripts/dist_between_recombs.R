## dist_between_recombs.R
##	Author: APM
##	Date: 9 January 2014
##	Purpose: make histogram of distance between adjacent recombination events in the G2F1

library(plyr)
library(ggplot2)

events <- read.csv("./../data/AllEvents6.csv")
events$mid <- with(events, (End + Start)/2)

widths <- ddply(events, .(Chr), function(d) {
	d <- d[ order(d$mid), ]
	to.next <- as.double(diff(d$mid))
	zeros <- which(to.next == as.double(0))
	fakes <- sapply(zeros, function(i) { runif(1, min = 0, max = ) })
	data.frame(to.next = diff(d$mid))
})

logify <- function(x, ...) {
	paste(10^x, "kbp")
}

rng <- -2:9
ggplot(widths) +
	geom_bar(aes(x = log10(to.next/1e3), y = ..ncount..), fill = "grey40") +
	scale_x_continuous(breaks = rng, labels = logify(rng)) +
	xlab("\ndistance between (midpoints of) adjacent recombination events") +
	ylab("proportion\n") +
	theme_bw()