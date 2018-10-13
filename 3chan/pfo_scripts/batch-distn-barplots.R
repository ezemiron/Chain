library(ggplot2)
setwd("/home/eze/Documents/papers/chrom_marks/results/distribution/")

markerlist =c("aHP1","H3K9me3","H3K9me2","H3K27me3","Macro-H2A","H4K5ac","H4K20me1","H3K36me2","H3K36me3","H3K4me2","H3K4me3","RNAPIIs2p","Scc1","Smc3","CTCF","SAF-A")
funclist =c("Inact","Inact”,”Inact","Inact","Mid","Mid","Mid","Act","Act","Act","Act","Act", "Struc","Struc","Struc","Struc")
markerLIST <- c("H3K4me3", "RNAPIIs2p", "Smc3")

markerfuncdf<-cbind(markerlist,funclist)
markerfuncdf<- as.data.frame(markerfuncdf)
FuncID <- c("Inact","Mid","Act","Struc")
conditionID <- c("Aux0h", "Aux2h", "ES","MS","LS")

files <- list.files(pattern="distn-compile.csv")
for (i in 1:length(files)) {
	temp <- read.csv(files[i])
	cell <- strsplit(files[i],"_")[[1]][1]
	temp <- cbind(rep(cell, nrow(temp)),temp) 
	colnames(temp)[1] <- "cell"
	assign(paste0("datamap",i),  temp)
	}

distnbind <- rbind(datamap1, datamap2, datamap3, datamap4)
distnbind[which(distnbind$Log2AvNorm=="-Inf"),"Log2AvNorm"] <- 0
pd <- position_dodge(0.9)

for (marker in unique(distnbind$Marker)){
  inddistnbind <- distnbind[distnbind$Marker %in% marker, ]
p <-ggplot(inddistnbind, aes(x=class, y=Log2AvNorm, fill=paste0(cell, " ", condition)))+
  geom_bar(stat = "identity",position = pd)+
  geom_errorbar(aes(ymin=Log2AvNorm-lowLogErr, ymax=Log2AvNorm+uprLogErr), width=.1,position = pd) +
  geom_hline(yintercept=0)+
  geom_vline(xintercept = 1.5)+
  labs(
    x = "Chromatin intensity class",
    y = "Log2-fold enrichment",
    fill = "cell type and treatment"
  )+
  ggtitle(marker)+
  #scale_y_continuous(limits=c(-80, 140),breaks=seq(-80,140,20)) +
  scale_x_continuous(breaks=seq(1,7,1))+
  scale_fill_manual(values=c("cadetblue1", "cadetblue4", "darkolivegreen3", "green4"))+
  theme(axis.line = element_line(colour = "black"),
        panel.background = element_blank())

savename <- paste0("/home/eze/Documents/papers/chrom_marks/results/distribution/", marker,"_condensed-bar.pdf")
pdf(savename)
p
dev.off()
}

for (marker in unique(distnbind$Marker)){
  inddistnbind <- distnbind[distnbind$Marker %in% marker, ]
	inddistnbind <- inddistnbind[which(inddistnbind$cell == "Parental"), ]
p <-ggplot(inddistnbind, aes(x=class, y=Log2AvNorm, fill=paste0(condition)))+
  geom_bar(stat = "identity",position = pd)+
  geom_errorbar(aes(ymin=Log2AvNorm-lowLogErr, ymax=Log2AvNorm+uprLogErr), width=.1,position = pd) +
  geom_hline(yintercept=0)+
  geom_vline(xintercept = 1.5)+
  labs(
    x = "Chromatin intensity class",
    y = "Log2-fold enrichment",
    fill = "treatment"
  )+
  ggtitle(paste0(marker, " in parental cells"))+
  #scale_y_continuous(limits=c(-80, 140),breaks=seq(-80,140,20)) +
  scale_x_continuous(breaks=seq(1,7,1))+
  #scale_fill_manual(values=c("cadetblue1", "cadetblue4", "darkolivegreen3", "green4"))+
  theme(axis.line = element_line(colour = "black"),
        panel.background = element_blank())

savename <- paste0("/home/eze/Documents/papers/chrom_marks/results/distribution/", marker,"_AUXeffect_parentalOnly-condensed-bar.pdf")
pdf(savename)
p
dev.off()
}


#################################
########## FOCI NUMBER ##########
#################################
by(distnbind$spotnum, distnbind$cell, sum)

distnbind$dataset <- paste0(distnbind$cell," ", distnbind$condition)
FOCI <- data.frame("marker"=0, "dataset"=0, "avg"=0, "count"=0)
for (m in unique(distnbind$Marker)) {
	foci <- data.frame(marker=rep(m, length(unique(distnbind$dataset))),
			   dataset=unique(distnbind$dataset), 
                           avg=rep(0, length(unique(distnbind$dataset))), 
                           count=rep(0, length(unique(distnbind$dataset))))
	for(i in 1:nrow(foci)){
		foci[i,"count"] <- sum(distnbind[which(distnbind$dataset==foci[i, "dataset"] & 
							distnbind$Marker==m) ,"spotnum"])
		foci[i, "avg"] <- foci[i, "count"]/distnbind[which(distnbind$dataset==foci[i, "dataset"] & 
								distnbind$Marker==m)[1], "numofcells"] 
	}
	FOCI <- rbind(FOCI, foci)
}
FOCI <- FOCI[-which(FOCI$marker==0),]


fig <- ggplot(FOCI)+
        geom_bar(aes(x=marker, y=avg, fill=dataset), stat="identity", position=pd) +
#geom_errorbar(aes(ymin=Log2AvNorm-lowLogErr, ymax=Log2AvNorm+uprLogErr), width=.1,position = pd) +
   #geom_hline(yintercept=0)+
   labs(
    y = "average foci per cell",
    fill = "cell type and treatment"
  )+
  ggtitle("Foci Counts") + 
  # scale_y_continuous(limits=c(-80, 140),breaks=seq(-80,140,20)) +
  # scale_x_continuous(breaks=seq(1,7,1))+
  scale_fill_manual(values=c("cadetblue1", "cadetblue4", "darkolivegreen3", "green4"))+
  theme(axis.line = element_blank(),
        panel.background = element_blank())

pdf("/home/eze/Documents/papers/chrom_marks/results/distribution/foci_counts.pdf")
fig
dev.off()
