#!/usr/bin/env Rscript
##
## Copyright (C) 2016 Ezequiel Miron <eze.miron@bioch.ox.ac.uk>
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.



dirchosen <- setwd("home/eze/Documents/papers/chrom_marks/results/")
setwd(dirchosen)

distributions <- list.files(pattern="distn.csv", recursive=TRUE)
fullcompile <- read.csv(distributions[1])
fullcompile$celltype = strsplit(distributions[1], "/")[[1]][1]
fullcompile$timecourse = strsplit(distributions[1], "/")[[1]][2]
fullcompile$marker = strsplit(distributions[1],"/")[[1]][3]

for( i in 2:length(distributions)) {
	comp <-  read.csv(distributions[i])
	comp$celltype = strsplit(distributions[i], "/")[[1]][1]
	comp$timecourse = strsplit(distributions[i], "/")[[1]][2] 
	comp$marker = strsplit(distributions[i],"/")[[1]][3]
	fullcompile <- rbind(fullcompile, comp)
}

count <- table(fullcompile[which(fullcompile$class==1), "timecourse"])
names(count) <- sapply(strsplit(names(count), "_"), "[[",1)

write.csv(fullcompile, file="/home/eze/Documents/papers/chrom_marks/results/fullcompile.csv")


ClassMean <- c(unlist(by(fullcompile$classnumn, paste0(fullcompile$class,"_",fullcompile$timecourse), mean, na.rm=TRUE)))
ClassSd <- c(unlist(by(fullcompile$classnumn, paste0(fullcompile$class,"_",fullcompile$timecourse), sd, na.rm=TRUE)))
if(all(names(ClassMean) == names(ClassSd))) {
	Class <- data.frame(Mean=ClassMean, Sd=ClassSd)
	Class$Compaction = sapply(strsplit(rownames(Class), "_"), "[[",1)
	Class$TrpTxt = sapply(strsplit(rownames(Class), "_"), "[[",2)
}
Class$n <- rep(0, nrow(Class))
for(i in 1:nrow(Class)) {
	Class[i, "n"] <- count[Class[i, "TrpTxt"]]
}
Class$upperCI <- Class$Mean + 1.96*Class$Sd/sqrt(Class$n)
Class$lowerCI <- Class$Mean - 1.96*Class$Sd/sqrt(Class$n)
TEXTSIZE=12; TITLESIZE=14
p <- ggplot(Class) + 
geom_bar(aes(x=Compaction, fill=TrpTxt, y=Mean), stat="identity", position=position_dodge(0.9))+
geom_errorbar(aes(x=Compaction, fill=TrpTxt, ymin=lowerCI, ymax=upperCI), width=0.7, position=position_dodge(0.9)) + 
theme(panel.background=element_blank(), legend.position=c(0.8, 0.8),
	axis.text=element_text(size=TEXTSIZE), axis.title=element_text(size=TITLESIZE), legend.text=element_text(size=TEXTSIZE),
	legend.title=element_text(size=TITLESIZE))+
labs(fill="triptolide \ntreatment",
	x="Chromatin Compaction Class",
	y="Proportion of Nucleus" ) + 
scale_fill_manual(values=c("cadetblue1", "cyan3", "cyan4"), breaks=c("0min", "10min", "30min"), labels=c("0 min", "10 min", "30 min"))

pdf("TrpCompaction.pdf")
p
dev.off()


fullcompile$cellID <- rep(seq(from=1, to=nrow(fullcompile)/7,), each =7)
volumes <- c(unlist(by(fullcompile$classnum, fullcompile$cellID, sum)))
fullcompile$volumes <- rep(0, nrow(fullcompile))
for(i in unique(fullcompile$cellID)) {
	fullcompile[which(fullcompile$cellID == i), "volumes"] <- volumes[i]
}
fullcompile1 <- fullcompile[which(fullcompile$class==1),]
volmean <- c(unlist(by(fullcompile1$classnum, fullcompile1$timecourse, mean, na.rm=TRUE)))
volsd <- c(unlist(by(fullcompile1$classnum, fullcompile1$timecourse, sd, na.rm=TRUE)))
vols <- data.frame(Mean=volmean, Sd=volsd)
vols$timecourse <- sapply(strsplit(rownames(vols), "_"), "[[",1)
vols$n <- count
vols$upperCI <- as.numeric(vols$Mean + 1.96*vols$Sd/sqrt(vols$n))
vols$lowerCI <- as.numeric(vols$Mean - 1.96*vols$Sd/sqrt(vols$n))

p<- ggplot(vols) +
theme(panel.background=element_blank(),
        axis.text=element_text(size=TEXTSIZE), 
	axis.title=element_text(size=TITLESIZE), 
	legend.text=element_text(size=TEXTSIZE),
        legend.title=element_text(size=TITLESIZE))+
geom_bar(aes(x=timecourse, y=Mean, fill=timecourse), position=position_dodge(0.9), stat="identity") + 
#geom_errorbar(aes(x=timecourse, ymin=lowerCI, ymax=upperCI), position=position_dodge(0.9), width=0.7)+
geom_errorbar(aes(x=timecourse, fill=timecourse, ymin=lowerCI, ymax=upperCI), width=0.7, position=position_dodge(0.9)) +
scale_fill_manual(values=c("cadetblue1", "cyan3", "cyan4"), breaks=c("0min", "10min", "30min"), labels=c("0 min", "10 min", "30 min"))



pdf("volumes2.pdf")
p
dev.off()
