setwd("~/Dropbox/G2F1/Andrew/")

library(ggplot2)

SEX.COLORS <- setNames( c("#ff00ff","#0080ff"), c("female","male") )
bsmaps <- read.table("./boostrap.maplength.txt", header = TRUE)
# bsmaps$generation <- relevel(bsmaps$generation, "G2")
bsmaps$xpos <- as.numeric(bsmaps$generation) + ifelse(bsmaps$sex == "male", -0.02, 0.02)

ggplot(bsmaps, aes(x = xpos, y = len, ymin = lo, ymax = hi, colour = sex, group = sex)) +
	geom_line() +
	# geom_errorbar(width = 0.1) +
	geom_linerange() +
	geom_point(shape = 21, size = 3, fill = "white") +
	xlab("") + ylab("cumulative autosomal map length (cM)\n") +
	scale_colour_manual(limits = names(SEX.COLORS), values = SEX.COLORS) +
	scale_x_continuous(breaks = range(as.numeric(bsmaps$generation)), labels = levels(bsmaps$generation), limits = c(0.8, 2.2)) +
	guides(colour = FALSE) +
	theme_bw() +
	theme(panel.background = element_blank(),
				# panel.border = element_blank(),
				panel.grid = element_blank()
				# axis.line.y = element_line(colour = "black"),
				# axis.ticks.x = element_blank()
				)

plot.new(  )