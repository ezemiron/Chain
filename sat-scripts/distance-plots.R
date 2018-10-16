library(ggplot2)
library(RColorBrewer)

d2bdata <- read.csv("/Users/ezejm/Dropbox (Schermelleh_Lab)/Schermelleh/EMPaper/graphs-for-paper/C127/2018.08-Final_eroded-stiched/G1_md2b-compile_Final.csv")
d2bdata2 <- read.csv("/Users/ezejm/Dropbox (Schermelleh_Lab)/Schermelleh/EMPaper/graphs-for-paper/C127/2018.08-Final_eroded-stiched/NaN3_md2b-compile_Final.csv")
d2bdata3 <- read.csv("/Users/ezejm/Dropbox (Schermelleh_Lab)/Schermelleh/EMPaper/graphs-for-paper/C127/2018.08-Final_eroded-stiched/Hyper-os_md2b-compile_Final.csv")
d2bdata4 <- read.csv("/Users/ezejm/Dropbox (Schermelleh_Lab)/Schermelleh/EMPaper/graphs-for-paper/C127/2018.08-Final_eroded-stiched/Trip/Trip_md2b-compile.csv")

markerlist =c("aHP1","H4K20me3","H3K9me3","H3K9me2","H3K27me3","H4K20me1","Macro-H2A","H4K5ac","H3K36me2","H3K4me2","H3K36me3","H3K4me3","Smc3","Scc1","CTCF","SAF-A","30min-BrUTP","RNAP-S2P","hnRNP")
markerlist = rev(markerlist)

#new markerlist for plotting smaller subset of markers:
markerlist =c("aHP1","H3K9me3","H3K9me2","H3K27me3","Macro-H2A","H3K4me2","H3K36me3","H3K4me3","Smc3","Scc1","CTCF","SAF-A","RNAP-S2P")
markerlist = rev(markerlist)


G1col <- "#F8766D"
NaN3col <- "#00BA38"
Hyperoscol<-"#619CFF"

d2bbind <- rbind(d2bdata,d2bdata2,d2bdata3)



#new script for compiled absolute discrete values:
#only one set at a time, need to change d2bdata, ggtitle and also errorbar and point colours to match:
ggplot(d2bdata, aes(x=Marker, y=mu.d2b.nm))+
geom_errorbar(aes(ymin=(mu.d2b.nm-CI95), ymax=(mu.d2b.nm+CI95)), width=.1, colour= G1col) +
  geom_point( colour= G1col)+
  geom_hline(yintercept=0)+
  geom_vline(xintercept = 4.5)+
  geom_vline(xintercept = 7.5)+
  geom_vline(xintercept = 13.5)+
  geom_vline(xintercept = 16.5)+
  labs(
    x = "Markers",
    y = "nm to IC surface"
  )+
  ggtitle("C127 G1")+
  scale_y_continuous(limits=c(-50, 180),breaks=seq(-60,120,20)) +
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
  geom_hline(yintercept=0)+
  geom_vline(xintercept = 2.5)+
  geom_vline(xintercept = 5.5)+
  geom_vline(xintercept = 9.5)+
  geom_vline(xintercept = 11.5)+
    labs(
    x = "Markers",
    y = "nm to IC surface"
  )+
#  scale_color_brewer(palette="Dark2")+
  ggtitle("C127 all conditions")+
  scale_color_grey(start = 0.5, end = 0.1)+
  scale_y_continuous(limits=c(-80, 200),breaks=seq(-80,140,20)) +
  scale_x_discrete(limits=markerlist)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.line = element_line(colour = "black"),
        panel.background = element_blank())














#prior scripts for continuous distances.
fewmarker <- subset(d2bdata2,d2bdata2$Marker %in% c("H3K27me3" , "H3K4me3"))


ggplot(d2bdata, aes(d2bdata$nm2bound, fill= d2bdata$Marker))+
  geom_histogram(binwidth = 30)+
  scale_x_continuous(name="d2b", limits=c(-300, 300),breaks=seq(-300,300,20)) +
#  scale_y_continuous(name="Stopping distance", limits=c(0, 150))#geom_density()

ggplot(d2bdata, aes(d2bdata$Marker, d2bdata$nm2bound, fill=d2bdata$Marker))+
  geom_boxplot()+#+ coord_flip()
labs(
  x = "Markers",
  y = "nm to surface",
  color = "Markers"
)+
  theme_bw()+
  scale_x_discrete(limits=markerlist)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

ggplot(d2bdata, aes(d2bdata$nm2bound, colour= d2bdata$Marker))+
  geom_freqpoly(binwidth = 10)+
labs(
  x = "nm to surface",
  y = "count",
  color = "Markers"
)+
  xlim(-300,300)+
  theme_bw()

#d2bnames <- names(d2bdata)
d2bcomp <- rbind(d2bdata,d2bdata2,d2bdata3)

ggplot(d2bcomp, aes(d2bcomp$Marker, d2bcomp$nm2bound, fill=d2bcomp$condition))+
  geom_boxplot()+# coord_flip()
  labs(
    x = "Markers",
    y = "nm to surface",
    fill = "Chromatin"
  )+
  scale_x_discrete(limits=markerlist)+
  theme_bw()+
  ylim(0,250)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))




fewmarker <- subset(d2bdata2,d2bdata2$Marker %in% c("H3K27me3" , "H3K4me3"))


ggplot(fewmarker, aes(fewmarker$nm2bound, fill= fewmarker$Marker))+
  geom_histogram(binwidth = 30)+
  scale_x_continuous(name="d2b", limits=c(-300, 300),breaks=seq(-300,300,20)) 
