library(ggplot2)
compile <- read.csv("/home/eze/Documents/papers/chrom_marks/results/distribution/full-compiled_June6.csv")
compile$dataset <- paste0(compile$cell, compile$condition)
print(compile$dataset) 
compile <- compile[-grep("mClover_pos", compile$marker), ]
compile$marker_f <- factor(compile$marker, levels=c("Smc3", "RNAPIIs2p", "H3K4me3", "H3K27me3", "H3K9me2"))

#COLORS=c("olivedrab1","darkolivegreen3", "cadetblue1", "cadetblue4")
#COLORS=c("palegreen1", "yellowgreen", "cadetblue1", "cadetblue4")
COLORS=c("olivedrab","olivedrab3","darkslategray","darkslategray4")
LABELS=c("Parental Aux (-)", "Parental Aux (+)", "Scc1mAID Aux (-)", "Scc1mAID Aux (+)")
BREAKS=c("HCT116ParentalAux0h", "HCT116ParentalAux2h", "HCT116Scc1mAIDAux0h", "HCT116Scc1mAIDAux2h")
TITLESIZE=18
TEXTSIZE=14


####################
## CELLS/DATASET ###
####################
datasets <- unique(compile$dataset)
count <- rep(0, length(datasets))
names(count) <- datasets

for(i in datasets) {
	count[i] <- nrow(compile[which(compile$class==1 & compile$dataset == i), ])
	}

count <- (table(compile[which(compile$class==1), "dataset"]))


####################
## NUCLEAR VOLUME ##
####################
compile$cellID <- as.numeric(factor(paste0(compile$cell, compile$condition, compile$marker, compile$cellno)))

volume <- rep(0, max(compile$cellID))
names(volume) <- unique(compile$cellID)
compile$volume <- rep(0, nrow(compile))

for(i in 1:max(compile$cellID) ){
	volume[i] <- sum(compile[which(compile$cellID==i), "classnum"])
	compile[which(compile$cellID==i), "volume"] <- volume[i]
	}

sum(volume)==sum(compile$classnum)
sum(volume)==sum(compile[seq(from=1, to=nrow(compile), by=7), "volume"])

per_cell <- compile[which(compile$class==1), ]

per_condit <- data.frame(dataset=datasets, vol_avg=rep(0, length(datasets)), 
	vol_sd=rep(0, length(datasets)), vol_upperCI=rep(0, length(datasets)),
	vol_lowerCI=rep(0, length(datasets)), n=rep(0, length(datasets)) ) 
rownames(per_condit) <- per_condit$dataset

for(i in per_condit$dataset) {
	per_condit[i, "vol_avg"] <- mean(per_cell[which(per_cell$dataset == per_condit[i, "dataset"]), "volume"])
	per_condit[i, "vol_sd"] <- sd(per_cell[which(per_cell$dataset == per_condit[i, "dataset"]), "volume"])
	per_condit[i, "vol_upperCI"] <- per_condit[i, "vol_avg"] + 1.96 * per_condit[i, "vol_sd"] / sqrt(count[i])
	per_condit[i, "vol_lowerCI"] <- per_condit[i, "vol_avg"] - 1.96 * per_condit[i, "vol_sd"] / sqrt(count[i])
	}

p <- ggplot(per_condit, aes(x= dataset, y=vol_avg, fill=dataset)) +
	#xlab("cell type and treatment") + 
	xlab("")+
	ylab("nuclear mask volume (voxels)") +
	theme(panel.background=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +   
	scale_fill_manual(values=COLORS, breaks=BREAKS, labels=LABELS)+
	geom_bar(position=position_dodge(0.9), stat="identity")+ 
	geom_errorbar(position=position_dodge(0.9), width=0.7, 
		aes(ymin=vol_lowerCI, ymax=vol_upperCI)) 


pdf("/home/eze/Documents/papers/chrom_marks/results/plots/nuclearvolume_june6.pdf")
p
dev.off()

####################
## N. FOCI #########
####################
test<- unlist(by(compile$classnum, compile$cellID, sum))
test==volume

foci <- unlist(by(compile$dist1t, compile$cellID, sum))
compile$total_foci <- rep(foci, each=7)

tester = 49
sum(compile[which(compile$cellID==tester), "dist1t"]) == foci[tester]
all(compile[which(compile$cellID==tester), "total_foci"]== foci[tester])

pdf("/home/eze/Documents/papers/chrom_marks/results/plots/foci_nonnormalized_june9.pdf")
ggplot(compile[which(compile$class=="1"),]) + ylab("total foci/image") + 
	geom_point(aes(x=marker_f, y=total_foci, color=dataset)) + 
	theme(panel.background=element_blank(), legend.key=element_blank(), axis.ticks.x=element_blank()) +
        scale_color_manual(values=COLORS, breaks=BREAKS, labels=LABELS)
dev.off()
	

compile$foci_norm <- compile$total_foci/compile$volume

per_cell <- compile[which(compile$class==1),]
foci_comp <- data.frame(mean= c(by(per_cell$foci_norm, paste0(per_cell$dataset,"*", per_cell$marker), mean)), 
			std= c(by(per_cell$foci_norm, paste0(per_cell$dataset,"*", per_cell$marker), sd)),
			n=c(by(per_cell$foci_norm, paste0(per_cell$dataset,"*", per_cell$marker), length)))

foci_comp$upperCI <- foci_comp$mean + 1.96*foci_comp$std/sqrt(foci_comp$n)
foci_comp$lowerCI <- foci_comp$mean - 1.96*foci_comp$std/sqrt(foci_comp$n)
foci_comp$marker <- unlist(strsplit(rownames(foci_comp), split="\\*"))[seq(from=2, to=nrow(foci_comp)*2,by=2)]
foci_comp$dataset <- unlist(strsplit(rownames(foci_comp), split="\\*"))[seq(from=1, to=nrow(foci_comp)*2, by=2)]
foci_comp$marker_f <- factor(foci_comp$marker, levels=c("Smc3", "RNAPIIs2p", "H3K4me3", "H3K27me3", "H3K9me2"))

p <- ggplot(foci_comp, aes(x=marker_f, y=mean, fill=dataset)) +
	geom_bar(stat="identity", position=position_dodge(0.9) ) + 
	geom_errorbar(aes(ymin=lowerCI, ymax=upperCI), position=position_dodge(0.9), width=0.7) + 
	labs(x= "IF mark",
        	y="Foci number/Nuclear voxel",
		fill="Cell type and treatment") +
	scale_fill_manual(values=COLORS, breaks=BREAKS, labels=LABELS)+ 
        theme(panel.background=element_blank())

pdf("/home/eze/Documents/papers/chrom_marks/results/plots/foci_june6.pdf")
p
dev.off()

p <- ggplot(foci_comp[grep("Aux2h", foci_comp$dataset),], aes(x=marker_f, y=mean, fill=dataset)) +
        geom_bar(stat="identity", position=position_dodge(0.9) ) +
        geom_errorbar(aes(ymin=lowerCI, ymax=upperCI), position=position_dodge(0.9), width=0.7) +
        labs(x= "IF mark",
                y="Foci number/Nuclear voxel",
                fill="Cell type and treatment") +
        scale_fill_manual(values=COLORS[c(2,4)], breaks=BREAKS, labels=LABELS)+
        theme(panel.background=element_blank())

pdf("/home/eze/Documents/papers/chrom_marks/results/plots/foci_AUXpos_june6.pdf")
p
dev.off()

marker <- "Smc3"
marker <- "H3K4me3"
marker <- "H3K27me3"
marker <- "RNAPIIs2p"
marker <- "H3K9me2"

savename <- paste0("/home/eze/Documents/papers/chrom_marks/results/plots/foci_auxpos_", marker,".pdf")
foci_plot <- compile[which(compile$marker == marker),]
foci_plot <- foci_plot[grep("Aux2h", foci_plot$dataset), ]
pdf(savename)
ggplot(foci_plot[which(foci_plot$class=="1"),]) + ylab("total foci / image") + xlab("nuclear volume") + 
	ggtitle(marker) + 
	geom_point(aes(x=volume, y=total_foci, color=dataset), size=2.5) +
	theme(panel.background=element_blank(), legend.key=element_blank(), axis.ticks.x=element_blank(),
	axis.title=element_text(size=TITLESIZE), axis.text=element_text(size=TEXTSIZE), 
	legend.text=element_text(size=TEXTSIZE), legend.title=element_text(size=TITLESIZE), 
	legend.position="none") +
        scale_color_manual(values=COLORS[c(2,4)], breaks=BREAKS, labels=LABELS) 
dev.off()


####################
## CHROM COMPACT ###
####################
compact_mean <- c(unlist(by(compile$classnumn, paste0(compile$dataset, "_", compile$class), mean)))
compact_sd <- c(unlist(by(compile$classnumn, paste0(compile$dataset, "_", compile$class), sd)))

if( all(names(compact_mean) == names(compact_sd))) {
	compaction <- data.frame(data=names(compact_mean), mean=compact_mean, sd=compact_sd)
	compaction$class <- sapply(strsplit(as.character(compaction$data), "_"), "[[", 2)
	compaction$dataset <- sapply(strsplit(as.character(compaction$data), "_"), "[[", 1)
	compaction$n <- rep(0, nrow(compaction))
	for(i in names(count)){
		compaction[which(compaction$dataset==i), "n"] <- count[i]
	}
	compaction$upperCI <- compaction$mean + 1.96 * compaction$sd / sqrt(compaction$n)
	compaction$lowerCI <- compaction$mean - 1.96 * compaction$sd / sqrt(compaction$n)
	
}
p <- ggplot(compaction, aes(x=class, y=mean, fill=dataset)) +
	labs(x= "chromatin compaction class", 
		y= "nuclear proportion", 
		fill="Cell type and treatment") +
        theme(panel.background=element_blank()) +
	scale_fill_manual(values=COLORS, labels=LABELS, breaks=BREAKS)+
        geom_bar(position=position_dodge(0.9), stat="identity")+
        geom_errorbar(position=position_dodge(0.9), width=0.7,
                aes(ymin=lowerCI, ymax=upperCI))

pdf("/home/eze/Documents/papers/chrom_marks/results/plots/compaction_plot_june6.pdf")
p
dev.off()

####################
## EFFECT OF AUX ###
####################
comp <- data.frame(class=unique(compaction$class), 
Parental=rep(0,length(unique(compaction$class))), Scc1=rep(0, length(unique(compaction$class))))
rownames(comp) <- comp$class

for ( i in comp$class ){
	keep <- compaction[which(keep$class==i),]
	comp[i, "Parental"] <- 1- (keep[grep("ParentalAux0h", keep$data), "mean"] / keep[grep("ParentalAux2h", keep$data), "mean"])
	comp[i, "Scc1"] <- 1- (keep[grep("Scc1mAIDAux0h", keep$data), "mean"] / keep[grep("Scc1mAIDAux2h", keep$data), "mean"])
	}


library(reshape); COMP <- melt(comp)
pdf("/home/eze/Documents/papers/chrom_marks/results/plots/chromatincompaction_auxineffect.pdf")
ggplot(COMP, aes(x=class, y=value, fill=variable)) +
	geom_bar(stat="identity", position=position_dodge(0.9)) + 
	labs( fill="cell type",
		x ="chromatin compaction class",
		y = "proportion in class Aux(+)/Aux(-)")  +
	theme(panel.background=element_blank(), legend.key=element_blank(), 
		axis.text=element_text(size=TEXTSIZE), axis.title=element_text(size=TITLESIZE),
		legend.text=element_text(size=TEXTSIZE), legend.title=element_text(size=TITLESIZE) )
dev.off()


####################
### SINGLE BATCH ###
####################
batch <- compile[which(compile$marker=="H3K9me2"), ]
compact_mean <- c(unlist(by(batch$classnumn, paste0(batch$dataset, "_", batch$class), mean)))
compact_sd <- c(unlist(by(batch$classnumn, paste0(batch$dataset, "_", batch$class), sd)))

if( all(names(compact_mean) == names(compact_sd))) {
        compaction <- data.frame(data=names(compact_mean), mean=compact_mean, sd=compact_sd)
        compaction$class <- sapply(strsplit(as.character(compaction$data), "_"), "[[", 2)
        compaction$dataset <- sapply(strsplit(as.character(compaction$data), "_"), "[[", 1)
        compaction$n <- rep(0, nrow(compaction))
	count_batch <- (table(batch[which(batch$class==1), "dataset"]))
        for(i in names(count_batch)){
                compaction[which(compaction$dataset==i), "n"] <- count_batch[i]
        }
        compaction$upperCI <- compaction$mean + 1.96 * compaction$sd / sqrt(compaction$n)
        compaction$lowerCI <- compaction$mean - 1.96 * compaction$sd / sqrt(compaction$n)

}

p <- ggplot(compaction, aes(x=class, y=mean, fill=dataset)) +
        labs(x= "chromatin compaction class",
                y= "nuclear proportion",
                fill="Cell type and treatment") +
	ggtitle("chromatin compaction for slides processed in batch IF")+ 
        theme(panel.background=element_blank()) +
        scale_fill_manual(values=COLORS, labels=LABELS, breaks=BREAKS)+
        geom_bar(position=position_dodge(0.9), stat="identity")+
        geom_errorbar(position=position_dodge(0.9), width=0.7,
                aes(ymin=lowerCI, ymax=upperCI))

pdf("/home/eze/Documents/papers/chrom_marks/results/plots/compaction_batch_plot_june6.pdf")
p
dev.off()



####################
## AUX vs NO-AUX ###
####################
compactAux_mean <- c(unlist(by(compile$classnumn, paste0(compile$condition, "_", compile$class), mean)))
compactAux_sd <- c(unlist(by(compile$classnumn, paste0(compile$condition, "_", compile$class), sd)))

if( all(names(compactAux_mean) == names(compactAux_sd))) {
        compactionA <- data.frame(data=names(compactAux_mean), mean=compactAux_mean, sd=compactAux_sd)
        compactionA$class <- sapply(strsplit(as.character(compactionA$data), "_"), "[[", 2)
        compactionA$dataset <- sapply(strsplit(as.character(compactionA$data), "_"), "[[", 1)
}

p <- ggplot(compactionA, aes(x=class, y=mean, fill=dataset)) +
        labs(x= "chromatin compaction class",
                y= "nuclear proportion",
                fill="Cell treatment") +
        theme(panel.background=element_blank(), legend.position=c(0.8,0.85), 
		axis) +
        ggtitle("Effect of Auxin-treatment, both cell types") +
        scale_fill_manual(values=c("red1", "red4"), labels=c("Aux (+)", "Aux (-)"), breaks=c("Aux2h", "Aux0h"))+
        geom_bar(position=position_dodge(0.9), stat="identity")+
        geom_errorbar(position=position_dodge(0.9), width=0.7,
                aes(ymin=mean-sd, ymax=mean+sd))

pdf("/home/eze/Documents/papers/chrom_marks/results/plots/compaction_AuxvsNoAux2.pdf")
p
dev.off()


####################
### Aux(+), P v S ##
####################
p <- ggplot(compaction[grep("Aux2h", compaction$data),], aes(x=class, y=mean, fill=dataset)) +
        labs(x= "chromatin compaction class",
                y= "nuclear proportion",
                fill="Cell type and treatment") +
        theme(panel.background=element_blank(), legend.position=c(0.8,0.85)) +
        scale_fill_manual(values=COLORS[c(2,4)], labels=LABELS, breaks=BREAKS)+
        geom_bar(position=position_dodge(0.9), stat="identity")+
        geom_errorbar(position=position_dodge(0.9), width=0.7,
                aes(ymin=lowerCI, ymax=upperCI))

pdf("/home/eze/Documents/papers/chrom_marks/results/plots/compaction_Auxpos.pdf")
p
dev.off()


####################
### AUX+ comp'rsn ##
####################
compile2 <-compile[which(compile$dataset %in% c("HCT116Scc1mAIDAux2h", "HCT116ParentalAux2h")),]
compactSvP_mean <- c(unlist(by(compile2$classnumn, paste0(compile2$cell, "_", compile2$class), mean)))
compactSvP_sd <- c(unlist(by(compile2$classnumn, paste0(compile2$cell, "_", compile2$class), sd)))

if( all(names(compactSvP_mean) == names(compactSvP_sd))) {
        compactionS <- data.frame(data=names(compactSvP_mean), mean=compactSvP_mean, sd=compactSvP_sd)
        compactionS$class <- sapply(strsplit(as.character(compactionS$data), "_"), "[[", 2)
        compactionS$dataset <- sapply(strsplit(as.character(compactionS$data), "_"), "[[", 1)
	}

p <- ggplot(compactionS, aes(x=class, y=mean, fill=dataset)) +
        labs(x= "chromatin compaction class",
                y= "nuclear proportion",
                fill="cell type, auxin treated") +
        theme(panel.background=element_blank()) +
	ggtitle("Comparison of Auxin-treated cells") + 
        scale_fill_manual(values=COLORS[c(2,4)])+
        geom_bar(position=position_dodge(0.9), stat="identity")+
        geom_errorbar(position=position_dodge(0.9), width=1,
                aes(ymin=mean-sd, ymax=mean+sd))

pdf("/home/eze/Documents/papers/chrom_marks/results/plots/compaction_AuxTxt.pdf")
p
dev.off()



####################
## IF DISTRIBUT ####
####################
focidist_mean <- c(unlist(by(compile$dist1n, paste0(compile$dataset, "_", compile$class, "_", compile$marker), mean)))
focidist_sd <- c(unlist(by(compile$dist1n, paste0(compile$dataset, "_", compile$class, "_", compile$marker), sd)))
focidist_AvLogNorm <- c(unlist(by(compile$lognorm1, paste0(compile$dataset, "_", compile$class, "_", compile$marker), mean)))


if( all(names(focidist_mean) == names(focidist_sd))) {
        focidist <- data.frame(data=names(focidist_mean), AvLogNorm=focidist_AvLogNorm, mean=focidist_mean, sd=focidist_sd)
        focidist$class <- sapply(strsplit(as.character(focidist$data), "_"), "[[", 2)
        focidist$dataset <- sapply(strsplit(as.character(focidist$data), "_"), "[[", 1)
	focidist$marker <- sapply(strsplit(as.character(focidist$data), "_"), "[[", 3)
        focidist$n <- rep(0, nrow(focidist))
        for(i in names(count)){
                focidist[which(focidist$dataset==i), "n"] <- count[i]
        }
        focidist$upperCI <- focidist$mean + 1.96 * focidist$sd / sqrt(focidist$n)
        focidist$lowerCI <- focidist$mean - 1.96 * focidist$sd / sqrt(focidist$n)

}
focidist$uprLogErr <- sqrt((focidist$AvLogNorm-log2(focidist$upperCI))^2)
focidist$lowLogErr <- sqrt((focidist$AvLogNorm-log2(focidist$lowerCI))^2) 

for(marker in unique(focidist$marker)) {
	savename <- paste0("/home/eze/Documents/papers/chrom_marks/results/plots/focidist_", marker, ".pdf")
	pdf(savename)
	ggplot(focidist[which(focidist$marker == marker),], 
		aes(x=class, y=AvLogNorm, fill=dataset)) +
        	labs(x= "chromatin compaction class",
        	        y= "foci distribution",
			title = marker,
        	        fill="Cell type and treatment") +
        	theme(panel.background=element_blank()) +
        	scale_fill_manual(values=COLORS, breaks=BREAKS, labels=LABELS)+
        	geom_bar(position=position_dodge(0.9), stat="identity")+
        	geom_errorbar(position=position_dodge(0.9), width=0.7,
        	        aes(ymin=AvLogNorm - lowLogErr, ymax=AvLogNorm + uprLogErr))
	
	dev.off()

	savename2 <- paste0("/home/eze/Documents/papers/chrom_marks/results/plots/focidist_Auxpos_", marker, ".pdf")
	pdf(savename2)
	auxpos <- focidist[grep("Aux2h", focidist$dataset),]
	ggplot(auxpos[which(auxpos$marker == marker),],
        	aes(x=class, y=mean, fill=dataset)) +
        	labs(x= "chromatin compaction class",
                	y= "foci distribution",
                	title = marker,
                	fill="Cell type and treatment") +
        	theme(panel.background=element_blank(),legend.position=c(0.75, 0.9),
			legend.title=element_text(size=TITLESIZE), legend.text=element_text(size=TEXTSIZE), 
			axis.title=element_text(size=TITLESIZE), axis.text=element_text(size=TEXTSIZE)) +
       		scale_fill_manual(values=COLORS[c(2,4)], breaks=BREAKS, labels=LABELS)+
        	geom_bar(position=position_dodge(0.9), stat="identity")+
        	geom_errorbar(position=position_dodge(0.9), width=0.7,
                	aes(ymin=lowerCI, ymax=upperCI))

	dev.off()
}





#### STATS: 'effect size'??
#### change in each cell line w.r.t aux: (parental 0aux - parental 2aux) vs (scc1maid 0aux - scc1maid 2aux)
#### what has a bigger effect size: auxin (in population of Scc1maid cell line) or cell type (in population of auxin-treated cells)
