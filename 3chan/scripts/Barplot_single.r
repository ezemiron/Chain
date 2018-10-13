#install.packages("ggplot2")
library(ggplot2)

dirchosen=getwd()
setwd(dirchosen)
csv_file1 = "DAPI_distn-compile.csv"
topPath1 = paste(dirchosen, csv_file1, sep = "/")

markerlist =c("s1_hnRNPC-488_DAPI", "s3_H4K20me3-488_DAPI", "s5_H3K27me3-488_DAPI")
funclist =c("Inact","Mid","Act")

markerfuncdf<-cbind(markerlist,funclist)
markerfuncdf<- as.data.frame(markerfuncdf)
FuncID <- c("Inact","Mid","Act")
conditionID <- c("DAPI")

datamap<- read.csv(topPath1)

pd <- position_dodge(0.9)

for (marker in markerlist){
  inddatamap <- datamap[datamap$Marker %in% marker, ]
p <-ggplot(datamap, aes(x=class, y=Log2AvNorm, fill=condition))+
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

savename <- paste0(marker,"_condensed-bar.png")
ggsave(savename, plot = p)
}



