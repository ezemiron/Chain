dsn <- list.files("./C127/" ,pattern="distn.csv", recursive=TRUE)
compile <- read.csv(paste0("./C127/", dsn[1]))
compile$file <- rep(dsn[1], nrow(compile)) 
for (i in 2:length(dsn)){
	temp <- read.csv(paste0("./C127/",dsn[i]))
	temp$file <- rep(dsn[i], nrow(temp))
	compile <- rbind(compile, temp)
}
compile$timecourse <- sapply(strsplit(compile$file, "/"), "[[", 1)
compile$marker <- sapply(strsplit(compile$file, "/"), "[[", 2)

 compile$dye <- rep("DAPI", nrow(compile))
 compile[grep("EU", compile$marker), "dye"]<-"SYTOX"



compile <- read.csv("/home/eze/Documents/papers/chrom_marks/results/fullcompile.csv")
COLORS=c("olivedrab","olivedrab3","darkslategray","darkslategray4")
TITLESIZE=18
TEXTSIZE=14

vox=c(0.041, 0.041, 0.125) 


####################
## CELLS/DATASET ###
####################
compile$timecourse <- gsub("^10min$", "20min", compile$timecourse) 
compile$timecourse <- gsub("^30min$", "40min", compile$timecourse) 
compile$timecourse <- gsub("^60min$", "70min", compile$timecourse) 
count <- table(compile[which(compile$class==1), "timecourse"])


####################
## NUCLEAR VOLUME ##
####################
compile$cellID <- rep(seq(from=1, to=nrow(compile)/7), each=7)
sum(count)==max(compile$cellID)


volume <- rep(0, max(compile$cellID))
names(volume) <- unique(compile$cellID)
compile$volume <- rep(0, nrow(compile))

for(i in 1:max(compile$cellID) ){
        volume[i] <- sum(compile[which(compile$cellID==i), "classnum"])
        compile[which(compile$cellID==i), "volume"] <- volume[i]
        }
all(volume==compile[seq(from=1, to=nrow(compile), by=7), "volume"])
compile$volume_micron = compile$volume * vox[1] * vox[2] * vox[3]
compile_norep <-compile[c(compile$marker %in% c("EU", "H3K27me3")), ]
per_cell <- compile_norep[which(compile_norep$class==1), ]


NUCLEARVOL=data.frame(
	MeanVol_micron =c(unlist(by(per_cell$volume_micron, per_cell$timecourse, mean))),
	SdVol_micron = c(unlist(by(per_cell$volume_micron, per_cell$timecourse, sd))),
	n=count)
NUCLEARVOL$CIup = NUCLEARVOL$MeanVol_micron + 1.96*NUCLEARVOL$SdVol_micron/sqrt(NUCLEARVOL$n.Freq)
NUCLEARVOL$CIlo = NUCLEARVOL$MeanVol_micron - 1.96*NUCLEARVOL$SdVol_micron/sqrt(NUCLEARVOL$n.Freq)


pdf("nucvol.pdf")
ggplot(NUCLEARVOL[]) + 
geom_bar(aes(x=n.Var1, y=MeanVol_micron), stat="identity") + 
geom_errorbar(aes(x=n.Var1, ymin=CIlo, ymax=CIup), width=0.6) + 
xlab("Triptolide treatment (min)")
dev.off()


####################
## N. FOCI #########
####################
foci <- unlist(by(compile$dist1t, compile$cellID, sum))
compile$total_foci <- rep(foci, each=7)

tester = 49
sum(compile[which(compile$cellID==tester), "dist1t"]) == foci[tester]

pdf("total-foci.pdf")
ggplot(compile[which(compile$class=="1"),]) + ylab("total foci/image") +
        geom_boxplot(aes(x=marker, y=total_foci, color=timecourse)) +
        theme(panel.background=element_blank(), legend.key=element_blank(), axis.ticks.x=element_blank()) 
dev.off()

compile$foci_norm <- compile$total_foci/compile$volume_micron

per_cell <- compile[which(compile$class==1),]
foci_comp <- data.frame(mean= c(by(per_cell$foci_norm, paste0(per_cell$timecourse,"*", per_cell$marker), mean)),
                        std= c(by(per_cell$foci_norm, paste0(per_cell$timecourse,"*", per_cell$marker), sd)),
                        n=c(by(per_cell$foci_norm, paste0(per_cell$timecourse,"*", per_cell$marker), length)))

foci_comp$upperCI <- foci_comp$mean + 1.96*foci_comp$std/sqrt(foci_comp$n)
foci_comp$lowerCI <- foci_comp$mean - 1.96*foci_comp$std/sqrt(foci_comp$n)
foci_comp$marker <- unlist(strsplit(rownames(foci_comp), split="\\*"))[seq(from=2, to=nrow(foci_comp)*2,by=2)]
foci_comp$dataset <- unlist(strsplit(rownames(foci_comp), split="\\*"))[seq(from=1, to=nrow(foci_comp)*2, by=2)]

p <- ggplot(foci_comp, aes(x=marker, y=mean, fill=dataset)) +
        geom_bar(stat="identity", position=position_dodge(0.9) ) +
        geom_errorbar(aes(ymin=lowerCI, ymax=upperCI), position=position_dodge(0.9), width=0.7) +
        labs(x= "IF mark",
                y="Foci number/cubic micron",
                fill="Cell type and treatment") +
        theme(panel.background=element_blank())

pdf("foci_june9.pdf")
p
dev.off()


####################
## CHROM COMPACT ###
####################
compact_mean <- c(unlist(by(compile$classnumn, paste0(compile$timecourse, "mmm", compile$class), mean)))
compact_sd <- c(unlist(by(compile$classnumn, paste0(compile$timecourse, "mmm", compile$class), sd)))

if( all(names(compact_mean) == names(compact_sd))) {
        compaction <- data.frame(data=names(compact_mean), mean=compact_mean, sd=compact_sd)
        compaction$class <- sapply(strsplit(as.character(compaction$data), "mmm"), "[[", 2)
        compaction$dataset <- sapply(strsplit(as.character(compaction$data), "mmm"), "[[", 1)
        compaction$n <- rep(0, nrow(compaction))
        for(i in names(count)){
                compaction[which(compaction$dataset==i), "n"] <- count[i]
        }
        compaction$upperCI <- compaction$mean + 1.96 * compaction$sd / sqrt(compaction$n)
        compaction$lowerCI <- compaction$mean - 1.96 * compaction$sd / sqrt(compaction$n)
       
 p <- ggplot(compaction, aes(x=class, y=mean, fill=dataset)) +
        labs(x= "chromatin compaction class",
                y= "nuclear proportion",
                fill="triptolide treatment\nlength") +
        theme(panel.background=element_blank()) +
         geom_bar(position=position_dodge(0.9), stat="identity")+
        geom_errorbar(position=position_dodge(0.9), width=0.7,
                aes(ymin=lowerCI, ymax=upperCI))
 
 compaction$dye <- rep("SYTOX", nrow(compaction))
 compaction[grep("split", compaction$dataset), "dye"]<-"DAPI"

  p <- ggplot(compaction[which(compaction$dye=="SYTOX"),], aes(x=class, y=mean, fill=dataset)) +
        labs(x= "chromatin compaction class",
                y= "nuclear proportion",
                fill="triptolide treatment\nlength") +
        theme(panel.background=element_blank()) +
         geom_bar(position=position_dodge(0.9), stat="identity")+
        geom_errorbar(position=position_dodge(0.9), width=0.7,
                aes(ymin=lowerCI, ymax=upperCI))+
                scale_fill_brewer(palette="YlOrRd")
   pdf("SYTOXcompaction.pdf") 
   p
   dev.off()   
 
 
####################
## IF DISTRIBUT ####
####################
count_marker <- table(paste0(compile[which(compile$class==1), "timecourse"], ">", compile[which(compile$class==1), "marker"]))

focidist_mean <- c(unlist(by(compile$dist1n, paste0(compile$timecourse, ">", compile$class, ">", compile$marker), mean)))
focidist_sd <- c(unlist(by(compile$dist1n, paste0(compile$timecourse, ">", compile$class, ">", compile$marker), sd)))
focidist_AvLogNorm <- c(unlist(by(compile$lognorm1, paste0(compile$timecourse, ">", compile$class, ">", compile$marker), mean)))


if( all(names(focidist_mean) == names(focidist_sd))) {
        focidist <- data.frame(data=names(focidist_mean), AvLogNorm=focidist_AvLogNorm, mean=focidist_mean, sd=focidist_sd)
        focidist$class <- sapply(strsplit(as.character(focidist$data), ">"), "[[", 2)
        focidist$dataset <- sapply(strsplit(as.character(focidist$data), ">"), "[[", 1)
        focidist$marker <- sapply(strsplit(as.character(focidist$data), ">"), "[[", 3)
        focidist$n <- rep(0, nrow(focidist))
        for(i in names(count_marker)){
                focidist[which(paste0(focidist$dataset, ">", focidist$marker)==i), "n"] <- count_marker[i]
        }
        focidist$upperCI <- focidist$mean + 1.96 * focidist$sd / sqrt(focidist$n)
        focidist$lowerCI <- focidist$mean - 1.96 * focidist$sd / sqrt(focidist$n)
focidist$uprLogErr <- sqrt((focidist$AvLogNorm-log2(focidist$upperCI))^2)
focidist$lowLogErr <- sqrt((focidist$AvLogNorm-log2(focidist$lowerCI))^2)

marker<- "H3K27me3"
   savename <- paste0("focidist_", marker, ".pdf")
        pdf(savename)
        ggplot(focidist[which(focidist$marker == marker),],
                aes(x=class, y=AvLogNorm, fill=dataset)) +
                labs(x= "chromatin compaction class",
                        y= "log2 fold enrichment",
                        title = marker,
                        fill="triptolide treatment") +
                theme(panel.background=element_blank()) +
              geom_bar(position=position_dodge(0.9), stat="identity")+
              scale_fill_brewer(palette="YlOrRd")+
                geom_errorbar(position=position_dodge(0.9), width=0.7,
                        aes(ymin=AvLogNorm - lowLogErr, ymax=AvLogNorm + uprLogErr))
                        dev.off()
 savename2 <- paste0("focidist_nonlog_", marker, ".pdf")
        pdf(savename2)
        ggplot(focidist[which(focidist$marker == marker),],
                aes(x=class, y=mean, fill=dataset)) +
                labs(x= "chromatin compaction class",
                        y= "foci distribution",
                        title = marker,
                        fill="triptolide treatment") +
                theme(panel.background=element_blank()) +
              scale_fill_brewer(palette="YlOrRd")+
                geom_bar(position=position_dodge(0.9), stat="identity")+
                geom_errorbar(position=position_dodge(0.9), width=0.7,
                        aes(ymin=lowerCI, ymax=upperCI))

        dev.off()

