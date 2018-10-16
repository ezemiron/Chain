library(plyr) 

setwd("/Users/eze/Dropbox (Schermelleh_Lab)/Schermelleh/EMPaper/graphs-for-paper/C127/2017.01-results/distribution/")
G1_distn-compile

markerlist =c("aHP1","H3K9me3","H3K9me2","H3K27me3","Macro-H2A","H4K5ac","H4K20me1","H3K36me2","H3K36me3","H3K4me2","H3K4me3","RNAP-S2P","Scc1","Smc3","CTCF","SAF-A")

datamap<- read.csv(
  "/Users/eze/Dropbox (Schermelleh_Lab)/Schermelleh/EMPaper/graphs-for-paper/C127/2017.01-results/distribution/G1_distn-compile.csv")
datamap2<- read.csv(
  "/Users/eze/Dropbox (Schermelleh_Lab)/Schermelleh/EMPaper/graphs-for-paper/C127/2017.01-results/distribution/NaN3_distn-compile.csv")
datamap3<- read.csv(
  "/Users/eze/Dropbox (Schermelleh_Lab)/Schermelleh/EMPaper/graphs-for-paper/C127/2017.01-results/distribution/Hyper-os_distn-compile.csv")

#####
###to plot table of sum of spot number for each marker in a given condition
Mnum=ddply(datamap, c("Marker"),summarise, spotnumber= sum(spotnum)) 


write.csv(Mnum,file="C127-G1_spotnum.csv")
