library(ggplot2)

d2bdata <- read.csv("G1_md2b-compile.csv")

markerlist =c("hnRNP","H3K4me3", "H3K27me3","H4K20me3")
markerlist = rev(markerlist)

G1col <- "#F8766D"

d2bbind <- rbind(d2bdata)

#new script for compiled absolute discrete values:
#only one set at a time, need to change d2bdata, ggtitle and also errorbar and point colours to match:
ggplot(d2bdata, aes(x=Marker, y=mu.d2b.nm))+
geom_errorbar(aes(ymin=(mu.d2b.nm-CI95), ymax=(mu.d2b.nm+CI95)), width=.1, colour= G1col) +
  geom_point( colour= G1col)+
  geom_hline(yintercept=0)+
  labs(
    x = "Markers",
    y = "nm to IC surface"
  )+
  ggtitle("C127 G1")+
  scale_y_continuous(limits=c(-80, 140),breaks=seq(-80,140,20)) +
  scale_x_discrete(limits=markerlist)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.line = element_line(colour = "black"),
        panel.background = element_blank())

#all conditioins together:
pd <- position_dodge(0.5)

ggplot(d2bbind, aes(x=Marker, y=mu.d2b.nm, colour=Condition))+
  geom_hline(yintercept=0)+
  geom_errorbar(aes(ymin=(mu.d2b.nm-CI95), ymax=(mu.d2b.nm+CI95)), width=.1, position=pd) +
  geom_point(position=pd)+
    labs(
    x = "Markers",
    y = "nm to IC surface"
  )+
  ggtitle("C127 all conditions")+
  scale_y_continuous(limits=c(-80, 140),breaks=seq(-80,140,20)) +
  scale_x_discrete(limits=markerlist)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.line = element_line(colour = "black"),
        panel.background = element_blank())

