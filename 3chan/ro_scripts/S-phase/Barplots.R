readline(prompt = "This script should be run from the directory containing the distn-compile.csv file")

#install.packages("ggplot2")
library(ggplot2)

dirchosen = getwd()
csv_file = "G1_distn-compile.csv"
csv_file1 = "ES_EdU_distn-compile.csv"
csv_file2 = "MS_EdU_distn-compile.csv"
csv_file3 = "LS_EdU_distn-compile.csv" 
#topPath1 = paste(dirchosen, csv_file sep = "/")

setwd(dirchosen)
markerlist = c("30min-BrUTP","H3K9me3","H3K9me2","H4K20me3","H3K27me3","H3K4me3","H4K20me3", "RNAP-S2P", "SCC1", "CTCF", "H3K36me3", "hnRNP")
funclist =c("Inact","Inact","Inact","Inact","Inact","Inact","Inact","Inact","Inact","Mid","Mid","Mid","Mid","Mid","Mid","Mid","Mid","MID","Act","Act","Act","Act","Act","Act","Act","Act","Act","Struc","Struc","Struc","Struc","Struc","Struc","Struc","Struc","Struc")

markerfuncdf<-cbind(markerlist,funclist)
markerfuncdf<- as.data.frame(markerfuncdf)
FuncID <- c("Inact", "Mid", "Act", "Struc")
conditionID <- c("15min EdU")

###############
#to make and save barplots automatically, one bar plot per marker,
#with G1, NaN3 and Hyper-os only
datamap<- read.csv(csv_file)
datamap1<- read.csv(csv_file1)
datamap2<- read.csv(csv_file2)
datamap3<- read.csv(csv_file3)
#datamap2<- read.csv(topPath2)
#  ("/Users/eze/Desktop/Dropbox (Schermelleh_Lab)/graphs-for-paper/C127/2017.01-results/distribution/NaN3_distn-compile.csv")
#datamap3<- read.csv()
#  ("/Users/eze/Desktop/Dropbox (Schermelleh_Lab)/graphs-for-paper/C127/2017.01-results/distribution/Hyper-os_distn-compile.csv")

# distnbind <- rbind(datamap)
distnbind <- rbind(datamap)

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
  scale_fill_manual(breaks=c("ES","MS","LS"),values=c("#000000","#969696","#636363"))+
  scale_y_continuous(expand = c(0,0), limits=c(-10,2)) +
  scale_x_continuous(breaks=seq(1,7,1))+
  theme(axis.line = element_line(colour = "black"),
        panel.background = element_blank(),legend.position="none")
 
savename <- paste0(marker,"_condensed-bar.pdf")
ggsave(savename, plot = p, device = "pdf")
}
