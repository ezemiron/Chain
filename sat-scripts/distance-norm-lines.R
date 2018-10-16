library(ggplot2)

setwd("/Users/eze/Desktop/Dropbox (Schermelleh_Lab)/graphs-for-paper/C127/2017.01-results/d2b-norm/")


d2bnorm<- read.csv(
  "/Users/eze/Desktop/Dropbox (Schermelleh_Lab)/graphs-for-paper/C127/2017.01-results/d2b-norm/G1_d2bnorm-compile.csv")
d2bnorm2<- read.csv(
  "/Users/eze/Desktop/Dropbox (Schermelleh_Lab)/graphs-for-paper/C127/2017.01-results/d2b-norm/NaN3_d2bnorm-compile.csv")
d2bnorm3<- read.csv(
  "/Users/eze/Desktop/Dropbox (Schermelleh_Lab)/graphs-for-paper/C127/2017.01-results/d2b-norm/Hyper-os_d2bnorm-compile.csv")

d2bnormbind <- rbind(d2bnorm,d2bnorm2,d2bnorm3)

markerlist =c("aHP1","H3K9me3","H3K9me2","H3K27me3","Macro-H2A","H4K5ac","H4K20me1","H3K36me2","H3K36me3","H3K4me2","H3K4me3","RNAP-S2P","Scc1","Smc3","CTCF","SAF-A")

#to plot averaged then logged results with shaded error bands:

for (marker in markerlist){
  indd2bnormbind <- d2bnormbind[d2bnormbind$Marker %in% marker, ]
p<-ggplot(data=indd2bnormbind, aes(x=stddist, y=Log2AvNorm, colour=Condition, fill=Condition )) + 
  geom_line()+
  geom_ribbon(aes(ymin=lowLogErr, ymax=uprLogErr), linetype=2, alpha=0.3)+
  geom_hline(yintercept=0)+
  geom_vline(xintercept=0)+
  labs(
    x = "Distance from IC border (nm)",
    y = "Log2-fold enrichment"
  )+
  ggtitle(marker)+
  scale_x_continuous(limits=c(-400, 400),breaks=seq(-400,400,50),expand = c(0, 0))+
  theme(axis.line = element_line(colour = "black"),
        panel.background = element_blank())

  
savename <- paste0(marker,"-d2bnorm-line.png")
ggsave(savename, plot = p)
}
