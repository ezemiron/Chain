library(ggplot2)
compile <- read.csv("/home/eze/Documents/papers/chrom_marks/results/EP201712_Glutaraldehyde/Csvs/full_compiled.csv")
compile$dataset <- paste0(compile$cell, compile$condition)
#compile <- compile[-grep("mclover_pos", compile$marker), ]
compile$marker_f <- factor(compile$marker, levels=c("Form_distn", "Glut_distn", "Glut_EtOH_distn"))

#COLORS=c("olivedrab1","darkolivegreen3", "cadetblue1", "cadetblue4")
#COLORS=c("palegreen1", "yellowgreen", "cadetblue1", "cadetblue4")
COLORS=c("olivedrab","olivedrab3","darkslategray")
LABELS=c("C127NucBlue", "C127NucBlue", "C127NucBlue")
BREAKS=c("C127NucBlueC127_Form", "C127NucBlueC127_Glut", "C127NucBlueC127_Glut_EtOH")
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

pdf("/home/eze/Documents/papers/chrom_marks/results/compaction_plot.pdf")
p
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
