readline(prompt = "This script should be run from the directory containing the distn-compile.csv file")

#install.packages("ggplot2")
library(ggplot2)

dirchosen = getwd()
csv_file1 = "ES/ES_EdU_spot-class-compile.csv"
csv_file2 = "MS/MS_EdU_spot-class-compile.csv"
csv_file3 = "LS/LS_EdU_spot-class-compile.csv" 
#topPath1 = paste(dirchosen, csv_file, sep = "/")

setwd(dirchosen)
markerlist = c("H4K20me3","H3K27me3","H3K4me3","H4K20me3", "RNAP-S2P", "SCC1", "CTCF", "H3K36me3", "hnRNP")
funclist =c("Inact","Inact","Inact","Inact","Inact","Inact","Inact","Inact","Inact","Mid","Mid","Mid","Mid","Mid","Mid","Mid","Mid","MID","Act","Act","Act","Act","Act","Act","Act","Act","Act","Struc","Struc","Struc","Struc","Struc","Struc","Struc","Struc","Struc"))
markerfuncdf<-cbind(markerlist,funclist)
markerfuncdf<- as.data.frame(markerfuncdf)
FuncID <- c("Inact", "Mid", "Act", "Struc")
conditionID <- c("2haux")

###############
#to make and save barplots automatically, one bar plot per marker,
#with G1, NaN3 and Hyper-os only
datamap1<- read.csv(csv_file1)
datamap2<- read.csv(csv_file2)
datamap3<- read.csv(csv_file3)
#datamap2<- read.csv(topPath2)
#  ("/Users/eze/Desktop/Dropbox (Schermelleh_Lab)/graphs-for-paper/C127/2017.01-results/distribution/NaN3_distn-compile.csv")
#datamap3<- read.csv()
#  ("/Users/eze/Desktop/Dropbox (Schermelleh_Lab)/graphs-for-paper/C127/2017.01-results/distribution/Hyper-os_distn-compile.csv")

distnbind <- rbind(datamap1,datamap2,datamap3)

pd <- position_dodge(0.9)

for (marker in markerlist){
  inddistnbind <- distnbind[distnbind$Marker %in% marker, ]
p <-ggplot(inddistnbind, aes(x=Marker, y=spotnum_vol_avg_um, fill=condition))+
  geom_bar(stat = "identity",position = pd)+
  geom_errorbar(aes(ymin=spotnum_vol_avg_um-spotnum_vol_SD_um, ymax=spotnum_vol_avg_um-spotnum_vol_SD_um), width=.1,position = pd) +
  labs(
    x = marker,
    y = "Marker number"
  )+
  ggtitle(marker)+
  scale_y_continuous(expand = c(0,0), limits=c(0,30)) +
  theme(axis.line = element_line(colour = "black"),
        panel.background = element_blank())

savename <- paste0(marker,"_EdU-foci-number.pdf")
ggsave(savename, plot = p, device = "pdf")
}
