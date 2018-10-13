library(ggplot2)

d2bdata <- read.csv("MS/MS_mask_md2b-compile.csv")
d2bdata2 <- read.csv("MS/MS_dif_md2b-compile.csv")
d2bdata3 <- read.csv("MS/MS_EdU_md2b-compile.csv")

markerlist =c("hnRNP","SCC1","CTCF","RNAP-S2P","H3K36me3","H3K4me3","H3K27me3","H4K20me3")
##markerlist = rev(markerlist)

d2bbind <- rbind(d2bdata,d2bdata2,d2bdata3)
d2bbind <- d2bdata


#new script for compiled absolute discrete values:
#only one set at a time, need to change d2bdata, ggtitle and also errorbar and point colours to match:
ggplot(d2bdata, aes(x=Marker, y=mu.d2b.nm))+
geom_errorbar(aes(ymin=(mu.d2b.nm-CI95), ymax=(mu.d2b.nm+CI95)), width=.1) +
  geom_point()+
  geom_hline(yintercept=0)+
  labs(
    x = "Markers",
    y = "nm to IC surface"
  )+
  ggtitle("C127 MS")+
  scale_y_continuous(limits=c(-80, 140),breaks=seq(-80,140,20)) +
  scale_x_discrete(limits=markerlist)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.line = element_line(colour = "black"),
        panel.background = element_blank())
