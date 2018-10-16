#install.packages("ggplot2")
library(ggplot2)
library(reshape2)

setwd("/Users/eze/Dropbox (Schermelleh_Lab)/Schermelleh/EMPaper/graphs-for-paper/C127/2017.01-results/distribution/")
#all markers
markerlist =c("aHP1","H4K20me3","H3K9me3","H3K9me2","H3K27me3","H4K20me1","Macro-H2A","H4K5ac","H3K36me2","H3K4me2","H3K36me3","H3K4me3","Smc3","Scc1","CTCF","SAF-A","30min-BrUTP","RNAP-S2P","hnRNP")

#funclist =c("Inact","Inact","Inact","Inact","Mid","Mid","Mid","Act","Act","Act","Act","Act", "Struc","Struc","Struc","Struc")
pd <- position_dodge(0.9)

markerfuncdf<-cbind(markerlist,funclist)
markerfuncdf<- as.data.frame(markerfuncdf)
FuncID <- c("Inact","Mid","Act","Struc")
conditionID <- c("G1", "NaN3", "Hyper-os", "ES","MS","LS")

###############
#to make and save barplots automatically, one bar plot per marker,
#with G1, NaN3 and Hyper-os only
setwd("/Users/ezejm/Dropbox (Schermelleh_Lab)/Schermelleh/EMPaper/graphs-for-paper/C127/2018.08-Final_eroded-stiched")

datamap<- read.csv(
  "/Users/ezejm/Dropbox (Schermelleh_Lab)/Schermelleh/EMPaper/graphs-for-paper/C127/2018.08-Final_eroded-stiched/G1_distn-compile_volnormE_Final-cropped.csv")
datamap2<- read.csv(
  "/Users/ezejm/Dropbox (Schermelleh_Lab)/Schermelleh/EMPaper/graphs-for-paper/C127/2018.08-Final_eroded-stiched/NaN3_distn-compile_volnormE_Final.csv")
datamap3<- read.csv(
  "/Users/ezejm/Dropbox (Schermelleh_Lab)/Schermelleh/EMPaper/graphs-for-paper/C127/2018.08-Final_eroded-stiched/Hyper-os_distn-compile_volnormE_Final.csv")

distnbind <- rbind(datamap,datamap2,datamap3)



for (marker in markerlist){
  inddistnbind <- distnbind[distnbind$Marker %in% marker, ]
p <-ggplot(inddistnbind, aes(x=class, y=Log2AvNorm, fill=condition))+
  geom_bar(stat = "identity",position = pd)+
  geom_errorbar(aes(ymin=Log2AvNorm-lowLogErr, ymax=Log2AvNorm+uprLogErr), width=.1,position = pd) +
  geom_hline(yintercept=0)+
  geom_vline(xintercept = 1.5)+
  labs(
    x = "Chromatin intensity class",
    y = "Log2-fold enrichment"
  )+
  ggtitle(marker)+
  scale_fill_grey(start = 0.5, end = 0.1)+
  #scale_y_continuous(limits=c(-6, 3),breaks=seq(-6,3,1)) +
  #scale_fill_manual(values=c("#F0E68C"))+
  scale_x_continuous(breaks=seq(1,7,1))+
  theme(axis.line = element_line(colour = "black"),
        panel.background = element_blank())
p

savename <- paste0(marker,"_trip-bar.pdf")
ggsave(savename, plot = p)
}


#######
#to make and save barplots automatically, one bar plot per marker,
#with G1, 6hrAux and 16hrDox-6hrAux Scc1 ablation

#all markers
markerlist =c("H4K20me3","H3K9me3","H3K27me3","H3K4me3","H3K4me2","H3K36me3")

setwd("/Users/ezejm/Dropbox (Schermelleh_Lab)/Schermelleh/EMPaper/graphs-for-paper/HCT116/Roel/")

#datamap<- read.csv(
  #"/Users/ezejm/Dropbox (Schermelleh_Lab)/Schermelleh/EMPaper/graphs-for-paper/C127/2018.08-Final_eroded-stiched/G1_distn-compile_volnormE_Final-cropped.csv")
datamap2<- read.csv(
  "/Users/ezejm/Dropbox (Schermelleh_Lab)/Schermelleh/EMPaper/graphs-for-paper/HCT116/Roel/6haux_distn-compile.csv")
datamap3<- read.csv(
  "/Users/ezejm/Dropbox (Schermelleh_Lab)/Schermelleh/EMPaper/graphs-for-paper/HCT116/Roel/16hdox-6haux_distn-compile.csv")

#datamap$Stain = NULL
datamap2$classnum = NULL
datamap2$spotnum = NULL
datamap3$classnum = NULL
datamap3$spotnum = NULL

#distnbind <- rbind(datamap,datamap2,datamap3)
distnbind <- rbind(datamap2,datamap3)


for (marker in markerlist){
  inddistnbind <- distnbind[distnbind$Marker %in% marker, ]
  p <-ggplot(inddistnbind, aes(x=class, y=Log2AvNorm, fill=condition))+
    geom_bar(stat = "identity",position = pd)+
    geom_errorbar(aes(ymin=Log2AvNorm-lowLogErr, ymax=Log2AvNorm+uprLogErr), width=.1,position = pd) +
    geom_hline(yintercept=0)+
    geom_vline(xintercept = 1.5)+
    labs(
      x = "Chromatin intensity class",
      y = "Log2-fold enrichment"
    )+
    ggtitle(marker)+
    scale_fill_grey(start = 0.5, end = 0.1)+
    #scale_y_continuous(limits=c(-6, 3),breaks=seq(-6,3,1)) +
    #scale_fill_manual(values=c("#F0E68C"))+
    scale_x_continuous(breaks=seq(1,7,1))+
    theme(axis.line = element_line(colour = "black"),
          panel.background = element_blank())
  p
  
  savename <- paste0(marker,"_AID-bar.pdf")
  ggsave(savename, plot = p)
}







##################
#to make and save barplots automatically, one bar plot per marker,
#with G1, ES , MS and LS only
datamap4<- read.csv(
  "/Users/eze/Dropbox (Schermelleh_Lab)/Schermelleh/EMPaper/graphs-for-paper/C127/2017.01-results/distribution/ES_distn-compile.csv")
datamap5<- read.csv(
  "/Users/eze/Dropbox (Schermelleh_Lab)/Schermelleh/EMPaper/graphs-for-paper/C127/2017.01-results/distribution/MS_distn-compile.csv")
datamap6<- read.csv(
  "/Users/eze/Dropbox (Schermelleh_Lab)/Schermelleh/EMPaper/graphs-for-paper/C127/2017.01-results/distribution/LS_distn-compile.csv")


distnbind2 <- rbind(datamap,datamap4,datamap5,datamap6)


pd <- position_dodge(0.9)

for (marker in markerlist){
  inddistnbind2 <- distnbind2[distnbind2$Marker %in% marker, ]
  p <-ggplot(inddistnbind2, aes(x=class, y=Log2AvNorm, fill=condition))+
    geom_bar(stat = "identity",position = pd)+
    geom_errorbar(aes(ymin=Log2AvNorm-lowLogErr, ymax=Log2AvNorm+uprLogErr), width=.1,position = pd) +
    geom_hline(yintercept=0)+
    geom_vline(xintercept = 1.5)+
    labs(
      x = "Chromatin intensity class",
      y = "Log2-fold enrichment"
    )+
    ggtitle(marker)+
    scale_fill_manual(values=c("#F8766D", "#D24F44", "#A62E25","#7B1C15"))+
    #scale_y_continuous(limits=c(-80, 140),breaks=seq(-80,140,20)) +
    scale_x_continuous(breaks=seq(1,7,1))+
    theme(axis.line = element_line(colour = "black"),
          panel.background = element_blank())
  
  savename <- paste0(marker,"_sphase-bar.png")
  ggsave(savename, plot = p)
}


###############
#to make and save barplots, one bar plot per collection of markers
#with either G1, NaN3 and Hyper-os, ES, MS, LS
datamap<- read.csv(
  "/Users/eze/Dropbox (Schermelleh_Lab)/Schermelleh/EMPaper/graphs-for-paper/C127/2017.01-results/distribution/G1_distn-compile.csv")
datamap2<- read.csv(
  "/Users/eze/Dropbox (Schermelleh_Lab)/Schermelleh/EMPaper/graphs-for-paper/C127/2017.01-results/distribution/NaN3_distn-compile.csv")
datamap3<- read.csv(
  "/Users/eze/Dropbox (Schermelleh_Lab)/Schermelleh/EMPaper/graphs-for-paper/C127/2017.01-results/distribution/Hyper-os_distn-compile.csv")
datamap4<- read.csv(
  "/Users/eze/Dropbox (Schermelleh_Lab)/Schermelleh/EMPaper/graphs-for-paper/C127/2017.01-results/distribution/ES_distn-compile.csv")
datamap5<- read.csv(
  "/Users/eze/Dropbox (Schermelleh_Lab)/Schermelleh/EMPaper/graphs-for-paper/C127/2017.01-results/distribution/MS_distn-compile.csv")
datamap6<- read.csv(
  "/Users/eze/Dropbox (Schermelleh_Lab)/Schermelleh/EMPaper/graphs-for-paper/C127/2017.01-results/distribution/LS_distn-compile.csv")

distnbind <- rbind(datamap,datamap2,datamap3,datamap4,datamap5,datamap6)


pd <- position_dodge(0.9)

for (marker in markerlist){
  inddistnbind <- distnbind[distnbind$Marker %in% marker, ]

  for (condition in conditionID){
    inddistnbind2 <- inddistnbind[inddistnbind$condition %in% condition, ]
  p <-ggplot(inddistnbind2, aes(x=class, y=Log2AvNorm))+
    geom_bar(stat = "identity",position = pd)+
    geom_errorbar(aes(ymin=Log2AvNorm-lowLogErr, ymax=Log2AvNorm+uprLogErr), width=.1,position = pd) +
    geom_hline(yintercept=0)+
    geom_vline(xintercept = 1.5)+
    labs(
      x = "Chromatin intensity class",
      y = "Log2-fold enrichment"
    )+
    ggtitle(paste(marker,condition))+
    #scale_y_continuous(limits=c(-80, 140),breaks=seq(-80,140,20)) +
    scale_x_continuous(breaks=seq(1,7,1))+
    theme(axis.line = element_line(colour = "black"),
          panel.background = element_blank())
  
   savename <- paste0(marker,"_",condition,"_single-bar.png")
   ggsave(savename, plot = p)
  }
}

###############
#to make and save barplots, one bar plot per GROUPED collection of markers
#with either G1, NaN3 and Hyper-os, ES, MS, LS
datamap<- read.csv(
  "/Users/eze/Dropbox (Schermelleh_Lab)/Schermelleh/EMPaper/graphs-for-paper/C127/2017.01-results/distribution/G1_distn-compile.csv")
datamap2<- read.csv(
  "/Users/eze/Dropbox (Schermelleh_Lab)/Schermelleh/EMPaper/graphs-for-paper/C127/2017.01-results/distribution/NaN3_distn-compile.csv")
datamap3<- read.csv(
  "/Users/eze/Dropbox (Schermelleh_Lab)/Schermelleh/EMPaper/graphs-for-paper/C127/2017.01-results/distribution/Hyper-os_distn-compile.csv")
datamap4<- read.csv(
  "/Users/eze/Dropbox (Schermelleh_Lab)/Schermelleh/EMPaper/graphs-for-paper/C127/2017.01-results/distribution/ES_distn-compile.csv")
datamap5<- read.csv(
  "/Users/eze/Dropbox (Schermelleh_Lab)/Schermelleh/EMPaper/graphs-for-paper/C127/2017.01-results/distribution/MS_distn-compile.csv")
datamap6<- read.csv(
  "/Users/eze/Dropbox (Schermelleh_Lab)/Schermelleh/EMPaper/graphs-for-paper/C127/2017.01-results/distribution/LS_distn-compile.csv")

distnbind <- rbind(datamap,datamap2,datamap3,datamap4,datamap5,datamap6)


pd <- position_dodge(0.9)

for (func in FuncID){
  markerfunc <- markerfuncdf[markerfuncdf$funclist %in% func, ]
  markerfuncvect<-as.vector(markerfunc$markerlist)
  inddistnbind <- distnbind[distnbind$Marker %in% markerfuncvect, ]

    for (condition in conditionID){
    inddistnbind2 <- inddistnbind[inddistnbind$condition %in% condition, ]
    #inddistnbind2 <- inddistnbind[inddistnbind$marker %in% markerfunc, ]
    
    p <-ggplot(inddistnbind2, aes(x=class, y=Log2AvNorm, fill=Marker))+
      geom_bar(stat = "identity",position = pd)+
      geom_errorbar(aes(ymin=Log2AvNorm-lowLogErr, ymax=Log2AvNorm+uprLogErr), width=.1,position = pd) +
      geom_hline(yintercept=0)+
      geom_vline(xintercept = 1.5)+
      labs(
        x = "Chromatin intensity class",
        y = "Log2-fold enrichment"
      )+
      ggtitle(paste(func,condition))+
      scale_fill_grey(start = 0.2, end = 0.8)+
      #scale_y_continuous(limits=c(-80, 140),breaks=seq(-80,140,20)) +
      scale_x_continuous(breaks=seq(1,7,1))+
      theme(axis.line = element_line(colour = "black"),
            panel.background = element_blank())
    
    savename <- paste0(func,"_",condition,"_bar.png")
    ggsave(savename, plot = p)
  }  
    
}







###############
#to make and save barplots automatically for IMR90, one bar plot per marker,
#with C127 G1, IMR90 G1 and IMR90 G0

setwd("/Users/eze/Dropbox (Schermelleh_Lab)/Schermelleh/EMPaper/graphs-for-paper/IMR90/2017.01-results/distribution/")

datamap<- read.csv(
  "/Users/eze/Dropbox (Schermelleh_Lab)/Schermelleh/EMPaper/graphs-for-paper/IMR90/2017.01-results/distribution/G1_distn-compile.csv")
datamap2<- read.csv(
  "/Users/eze/Dropbox (Schermelleh_Lab)/Schermelleh/EMPaper/graphs-for-paper/IMR90/2017.01-results/distribution/G0_distn-compile.csv")

distnbind <- rbind(datamap,datamap2)


pd <- position_dodge(0.9)

for (marker in markerlist){
  inddistnbind <- distnbind[distnbind$Marker %in% marker, ]
  p <-ggplot(inddistnbind, aes(x=class, y=Log2AvNorm, fill=condition))+
    geom_bar(stat = "identity",position = pd)+
    geom_errorbar(aes(ymin=Log2AvNorm-lowLogErr, ymax=Log2AvNorm+uprLogErr), width=.1,position = pd) +
    geom_hline(yintercept=0)+
    geom_vline(xintercept = 1.5)+
    labs(
      x = "Chromatin intensity class",
      y = "Log2-fold enrichment"
    )+
    ggtitle(marker)+
    #scale_y_continuous(limits=c(-80, 140),breaks=seq(-80,140,20)) +
    scale_x_continuous(breaks=seq(1,7,1))+
    theme(axis.line = element_line(colour = "black"),
          panel.background = element_blank())
  
  savename <- paste0(marker,"_senescent-bar.png")
  ggsave(savename, plot = p)
}


###############
#to make and save barplots automatically for HCT116, one bar plot per marker,
#with parental or Scc1mAID

setwd("/Users/eze/Dropbox (Schermelleh_Lab)/Schermelleh/EMPaper/graphs-for-paper/HCT116/distributions/")
COLORS=c("olivedrab","olivedrab3","darkslategray","darkslategray4")
datamap<- read.csv(
  "/Users/eze/Dropbox (Schermelleh_Lab)/Schermelleh/EMPaper/graphs-for-paper/HCT116/HCT116_distn-compiled.csv")

distnbind <- datamap

pd <- position_dodge(0.9)

for (marker in markerlist){
  inddistnbind <- distnbind[distnbind$Marker %in% marker, ]
  p <-ggplot(inddistnbind, aes(x=class, y=Log2AvNorm, fill=condition))+
    geom_bar(stat = "identity",position = pd)+
    geom_errorbar(aes(ymin=Log2AvNorm-lowLogErr, ymax=Log2AvNorm+uprLogErr), width=.1,position = pd) +
    geom_hline(yintercept=0)+
    geom_vline(xintercept = 1.5)+
    labs(
      x = "Chromatin intensity class",
      y = "Log2-fold enrichment"
    )+
    ggtitle(marker)+
    #scale_y_continuous(limits=c(-80, 140),breaks=seq(-80,140,20)) +
    scale_x_continuous(breaks=seq(1,7,1))+
    scale_fill_manual(values=COLORS[c(2,4)])+
    theme(axis.line = element_line(colour = "black"),
          panel.background = element_blank())
  
  savename <- paste0(marker,"_Auxin-bar.png")
  ggsave(savename, plot = p)
}




###############
#to make and save a barplot for WF controls
markerlist =c("SIR","PWF-SIR","PWF")
distnbind <- read.csv("WF-distribution-summary.csv")
conditionID <- c("WF-ctrl")
condition <- conditionID
inddistnbind2 <- distnbind

p <-ggplot(inddistnbind2, aes(x=class, y=lognorm1, fill=Marker))+
  geom_bar(stat = "identity",position = pd)+
  #geom_errorbar(aes(ymin=Log2AvNorm-lowLogErr, ymax=Log2AvNorm+uprLogErr), width=.1,position = pd) +
  geom_hline(yintercept=0)+
  geom_vline(xintercept = 1.5)+
  labs(
    x = "Chromatin intensity class",
    y = "Log2-fold enrichment",
    fill = "Resolution"
  )+
  #ggtitle(paste(marker,condition))+
  #scale_y_continuous(limits=c(-80, 140),breaks=seq(-80,140,20)) +
  scale_x_continuous(breaks=seq(1,7,1))+
  scale_fill_grey(start = 0.8, end = 0.2)+
  theme(axis.line = element_line(colour = "black"),
        panel.background = element_blank())
p
savename <- paste0(condition,"_bar.png")
ggsave(savename, plot = p)
 