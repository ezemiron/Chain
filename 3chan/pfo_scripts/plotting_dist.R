#!/usr/bin/env Rscript
##
library(ggplot2)

PATH <- "/home/eze/Documents/papers/chrom_marks/results/"
distn <- list.files(PATH, pattern="distn-summary.csv")
columns <- c("Marker", "class", "AvNorm")
if (length(distn) >= 1) {
	DIST <- read.csv(paste0(PATH, distn[1]))
	
	
	DIST <- cbind(DIST[,columns], 
		rep(sub("distn-summary.csv", "", distn[1]), nrow(DIST)))
		colnames(DIST)[ncol(DIST)] <- "dataset"

for (i in distn[2:length(distn)]) {
	print(i)
	DISTprep <- read.csv(paste0(PATH, i))
	DISTprep <- cbind(DISTprep[,columns],
		rep(sub("distn-summary.csv", "", i), nrow(DISTprep)))
	colnames(DISTprep)[ncol(DISTprep)] <- "dataset"
	DIST <- rbind(DIST, DISTprep)
	}
}

DIST <- data.frame(DIST)
DISTo <- DIST[-grep("crop", DIST$dataset),]

pdf(file="/home/eze/Documents/papers/chrom_marks/results/distribution1.pdf")
ggplot(DISTo, aes(x=X, y=AvNorm, color=dataset, fill=dataset)) + 
	# geom_area(alpha=0.5) +
	geom_ribbon(aes(ymin=0, ymax=AvNorm, x=class), alpha=0.5) +
	theme(panel.background=element_blank(), legend.position=c(0.7, 0.9),
	axis.text.x=element_text(size=14), axis.title=element_text(size=14), 
	legend.text=element_text(size=14), legend.title=element_blank()) +
	xlab("chromatin compaction class") + 
	ylab("proportion of nucleus in compaction class") + 
	scale_color_manual(values=c("blue", "deepskyblue", "olivedrab3","goldenrod1")) +
	scale_fill_manual(values=c("blue", "deepskyblue", "olivedrab3", "goldenrod1"))
dev.off()
