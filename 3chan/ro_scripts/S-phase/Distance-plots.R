library(ggplot2)

#d2bdata <- read.csv("ES/ES_dif_md2b-compile.csv")
#d2bdata2 <- read.csv("MS/MS_dif_md2b-compile.csv")
#d2bdata3 <- read.csv("LS/LS_dif_md2b-compile.csv")
d2bdata4 <- read.csv("G1_md2b-compile.csv")

markerlist =c("30min-BrUTP","hnRNP","SCC1","CTCF","RNAP-S2P","H3K36me3","H3K4me3","H3K27me3","H3K9me2","H3K9me3","H4K20me3")
##markerlist = rev(markerlist)

G1col <- "#F8766D"
NaN3col <- "#00BA38"
Hyperoscol <-"#619CFF"
Color4 <- "#FF4F00"

#d2bbind <- rbind(d2bdata4,d2bdata,d2bdata2,d2bdata3)
d2bbind <- d2bdata4

#all conditioins together:
pd <- position_dodge(0.5)

ggplot(d2bbind, aes(x=Marker, y=mu.d2b.nm, colour=Condition, fill=Condition))+
  geom_hline(yintercept=0)+
  geom_errorbar(aes(ymin=(mu.d2b.nm-CI95), ymax=(mu.d2b.nm+CI95)), width=.1, position=pd) +
  geom_point(position=pd)+
    labs(
    x = "Markers",
    y = "nm to IC surface"
  )+
  ggtitle("EdU-mask all conditions")+
  scale_color_manual(values=c("#000000"))+
  scale_fill_manual(values=c("#000000"))+
  scale_y_continuous(limits=c(-80, 140),breaks=seq(-80,140,20)) +
  scale_x_discrete(limits=markerlist)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.line = element_line(colour = "black"),
        panel.background = element_blank())
