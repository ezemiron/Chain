#! /usr/bin/env RScript
##

library(ggplot2)

## pick up all pf the distribution files that have been created
setwd("/home/eze/Documents/papers/chrom_marks/results")
all=list.files(recursive=TRUE)
allDIS=all[grep("distn.csv", all)]

## make a dataframe to store distribution info
df <- data.frame(cell=NA, condition=NA, marker=NA, file=NA, volume=NA)

for (i in 1:length(allDIS)) {
	dis <- read.csv(allDIS[i])
	df[i,"cell"] <- strsplit(allDIS[i], "/")[[1]][1]
	df[i,"condition"] <-  strsplit(allDIS[i], "/")[[1]][2]
	df[i,"marker"] <-  strsplit(allDIS[i], "/")[[1]][3]
	df[i, "file"] <- strsplit(allDIS[i], "/")[[1]][4]
	df[i,"volume"] <- sum(dis$classnum)
}
df$crop <- rep("uncropped", nrow(df))
df[grep("crop", df$marker),"crop"] <- "cropped"

pdf("/home/eze/Documents/papers/chrom_marks/results/nuclearvolume2.pdf")
ggplot(df) + 
geom_point(aes(x=paste0(cell, "\n", condition, " ", marker), y=volume, color=crop), size=2) +
xlab("") + ylab("nuclear mask volume") + 
theme(panel.background=element_blank(), legend.key=element_blank(),
 legend.title=element_blank(), axis.text.x=element_text(size=14, angle=90), axis.ticks.x=element_blank(), 
axis.text.y=element_text(size=14),  axis.title=element_text(size=15), legend.text=element_text(size=14))
dev.off()
