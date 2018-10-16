library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(dplyr)


SAdata <- read.csv("/Users/ezejm/Dropbox (Schermelleh_Lab)/Schermelleh/EMPaper/graphs-for-paper/C127/2018.08-Final_eroded-stiched/C127_G1_SA2vol_raw_Final.csv")
SAdata2 <- read.csv("/Users/ezejm/Dropbox (Schermelleh_Lab)/Schermelleh/EMPaper/graphs-for-paper/C127/2018.08-Final_eroded-stiched/C127_NaN3_SA2vol_raw_Final.csv")
SAdata3 <- read.csv("/Users/ezejm/Dropbox (Schermelleh_Lab)/Schermelleh/EMPaper/graphs-for-paper/C127/2018.08-Final_eroded-stiched/C127_Hyper-os_SA2vol_raw_Final.csv")
SAdata4 <- read.csv("/Users/ezejm/Dropbox (Schermelleh_Lab)/Schermelleh/EMPaper/graphs-for-paper/C127/2018.08-Final_eroded-stiched/Trip/C127_Trip_SA2vol_raw.csv")


SAdata$condition <- "G1"
SAdata2$condition <- "NaN3"
SAdata3$condition <- "Hyper-os"
SAdata4$condition <- "Trip"

SAcomp <- rbind(SAdata,SAdata2,SAdata3)



#####
#plot violin SA 2 vol:
ggplot(SAcomp, aes(SAcomp$condition, SAcomp$SA2vol, fill=SAcomp$condition))+
  #geom_boxplot()+# coord_flip()
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75))+
  labs(
    x = "Chromatin conditions",
    y = "SA:Vol",
    fill = "Chromatin"
  )+
  theme_bw()+
  scale_x_discrete(limits=c("G1","NaN3","Hyper-os"))+
  #ylim(0,1)+
  scale_fill_manual(breaks=c("G1","NaN3","Hyper-os"),values=c("#F8766D", "#619CFF","#00BA38"))
  #theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=1))

#plot single condition
ggplot(SAdata, aes(SAdata$condition, SAdata$SA2vol))+
  geom_violin(fill='#F8766D',draw_quantiles = c(0.25, 0.5, 0.75))+
  labs(
    x = "Chromatin condition",
    y = "Chromatin SA:Volume"#,
    #  fill = "Chromatin"
  )+
  theme(axis.line = element_line(colour = "black"),panel.background = element_blank())+
  scale_x_discrete(limits=c("G1"))+
  #ylim(0,1)+
  scale_fill_manual(breaks=c("G1"),values=c("#F8766D", "#619CFF","#00BA38"))
#theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=1))



#dot plot of mean with 95CI
gdata <- group_by(SAcomp, condition)
sdata <- summarise(gdata, mean=mean(SA2vol), sd=sd(SA2vol), num=n())
df <- transform(sdata, CI95 =1.96*sd/sqrt(num))
df <- transform(df, lowErr = mean-CI95)
df <- transform(df, uprErr = mean+CI95)

pd <- position_dodge(0.5)

ggplot(df, aes(x=condition, y=mean, colour=condition))+
  geom_errorbar(aes(ymin=(mean-CI95), ymax=(mean+CI95)), width=.1) +
  geom_point()+geom_line()+
  #geom_hline(yintercept=0)+
  #geom_vline(xintercept = 4.5)+
  #geom_vline(xintercept = 7.5)+
  #geom_vline(xintercept = 13.5)+
  #geom_vline(xintercept = 16.5)+
  labs(
    x = "Condition",
    y = "Surface area to volume ratio"
  )+
  #ggtitle("C127 G1")+
  #scale_color_grey(start = 0.5, end = 0.1)+
  scale_x_discrete(limits=c("G1","Trip"))+
  scale_y_continuous(limits=c(0.28, 0.38),breaks=seq(0.28,0.38,0.02)) +
  scale_color_manual(breaks=c("G1","Trip"),values=c("#7f7f7f","#191919"))+
  #scale_x_continuous(breaks=seq(1,7,1))+
  theme(axis.line = element_line(colour = "black"),
        panel.background = element_blank())



#####
#plot chrom2vol
ggplot(SAcomp, aes(SAcomp$condition, SAcomp$chrom2vol, fill=SAcomp$condition))+
  #geom_boxplot()+# coord_flip()
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75))+
  labs(
    x = "Chromatin conditions",
    y = "chromatin:nucleus",
    fill = "Chromatin"
  )+
  theme_bw()+
  scale_x_discrete(limits=c("G1","NaN3","Hyper-os"))+
  #ylim(0,1)+
  scale_fill_manual(breaks=c("G1","NaN3","Hyper-os"),values=c("#F8766D", "#619CFF","#00BA38"))
#theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=1))

#plot single condition
ggplot(SAdata, aes(SAdata$condition, SAdata$chrom2vol))+
  #geom_boxplot(fill='#F8766D')+# coord_flip()
  geom_violin(fill='#F8766D',draw_quantiles = c(0.25, 0.5, 0.75))+
  labs(
    x = "Chromatin condition",
    y = "chromatin:nucleus"#,
    #  fill = "Chromatin"
  )+
  theme(axis.line = element_line(colour = "black"),panel.background = element_blank())+
  scale_x_discrete(limits=c("G1"))+
  #ylim(0,1)+
  scale_fill_manual(breaks=c("G1"),values=c("#F8766D", "#619CFF","#00BA38"))
#theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=1))


#####
#plot whole nuclear volume:

voxnm3 <- 41*41*125
voxum3 <- voxnm3/1000000000


ggplot(SAcomp, aes(SAcomp$condition, SAcomp$nucvol, fill=SAcomp$condition))+
  geom_boxplot()+# coord_flip()
  labs(
    x = "Chromatin conditions",
    y = "Nuclear volume",
    fill = "Chromatin"
  )+
  theme_bw()+
  scale_x_discrete(limits=c("G1","NaN3","Hyper-os"))+
  #ylim(0,1)+
  scale_fill_manual(breaks=c("G1","NaN3","Hyper-os"),values=c("#F8766D", "#619CFF","#00BA38"))
#theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=1))


  #plot single condition:

SAdatavox <- data.frame(SAdata$nucvol*voxum3, SAdata$condition)
names(SAdatavox) <- c("nucvol","condition")

ggplot(SAdatavox, aes(SAdatavox$condition, SAdatavox$nucvol))+
  geom_violin(fill='#F8766D',draw_quantiles = c(0.25, 0.5, 0.75))+# coord_flip()
  labs(
    x = "Chromatin condition",
    y = "Nuclear volume (um^3)"#,
    #fill = "Chromatin"
  )+
  theme(axis.line = element_line(colour = "black"),panel.background = element_blank())+
    #scale_x_discrete(limits=c("G1"))+
  #ylim(0,1)+
  scale_fill_manual(breaks=c("G1"),values=c("#F8766D", "#619CFF","#00BA38"))
#theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=1))
#export as 390 by 358 aspect ratio to get square.

median(SAdatavox$nucvol)
mean(SAdatavox$nucvol)

#############

#plot absolute changes in each chromatin density class
datamap<- read.csv(
  "/Users/ezejm/Dropbox (Schermelleh_Lab)/Schermelleh/EMPaper/graphs-for-paper/C127/2018.08-Final_eroded-stiched/G1_distn-compile_volnormE_Final.csv")
datamap2<- read.csv(
  "/Users/ezejm/Dropbox (Schermelleh_Lab)/Schermelleh/EMPaper/graphs-for-paper/C127/2018.08-Final_eroded-stiched/NaN3_distn-compile_volnormE_Final.csv")
datamap3<- read.csv(
  "/Users/ezejm/Dropbox (Schermelleh_Lab)/Schermelleh/EMPaper/graphs-for-paper/C127/2018.08-Final_eroded-stiched/Hyper-os_distn-compile_volnormE_Final.csv")
datamap4 <- read.csv(
  "/Users/ezejm/Dropbox (Schermelleh_Lab)/Schermelleh/EMPaper/graphs-for-paper/C127/2018.08-Final_eroded-stiched/Trip/Trip_distn-compile.csv")



datamapcomp <- rbind(datamap,datamap2, datamap3)

ggplot(datamapcomp, aes(factor(datamapcomp$class), datamapcomp$classnum, fill=datamapcomp$condition))+
  geom_boxplot()+# coord_flip()
  labs(
    x = "Chromatin class",
    y = "Absolute class volume",
    fill = "Chromatin"
  )+
  #theme_bw()+
  theme(axis.line = element_line(colour = "black"),
        panel.background = element_blank())
  #scale_x_discrete(limits=c("G1","NaN3","Hyper-os"))+
  #ylim(0,1)+
  scale_fill_manual(breaks=c("G1","NaN3","Hyper-os"),values=c("#F8766D", "#00BA38","#619CFF"))
#theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=1))

  #to plot single condition
  ggplot(datamap, aes(factor(datamap$class), datamap$classnum))+
    geom_boxplot(fill='#F8766D')+# coord_flip()
    labs(
      x = "Chromatin class",
      y = "Absolute class volume (voxels)",
      fill = "Chromatin"
    )+
    #theme_bw()+
    theme(axis.line = element_line(colour = "black"),
          panel.background = element_blank())
  #scale_x_discrete(limits=c("G1","NaN3","Hyper-os"))+
  #ylim(0,1)+
  scale_fill_manual(breaks=c("G1","NaN3","Hyper-os"),values=c("#F8766D", "#00BA38","#619CFF"))
  
  
  
  #dot plot of mean with 95CI
  gdata <- group_by(datamapcomp, condition, class)
  sdata <- summarise(gdata, mean=mean(classnumn), sd=sd(classnumn), num=n())
  df <- transform(sdata, CI95 =1.96*sd/sqrt(num))
  df <- transform(df, lowErr = mean-CI95)
  df <- transform(df, uprErr = mean+CI95)
  
  pd <- position_dodge(0.5)
  
  ggplot(df, aes(x=class, y=mean, colour=condition))+
    geom_errorbar(aes(ymin=(mean-CI95), ymax=(mean+CI95)), width=.1) +
    geom_point()+geom_line()+
    #geom_hline(yintercept=0)+
    #geom_vline(xintercept = 4.5)+
    #geom_vline(xintercept = 7.5)+
    #geom_vline(xintercept = 13.5)+
    #geom_vline(xintercept = 16.5)+
    labs(
      x = "Class",
      y = "Class volume / nuclear volume"
     )+
    ggtitle("C127 G1")+
    scale_color_grey(start = 0.5, end = 0.1)+
    scale_y_continuous(limits=c(0, 0.4),breaks=seq(0,0.4,0.1)) +
  scale_x_continuous(breaks=seq(1,7,1))+
    theme(axis.line = element_line(colour = "black"),
         panel.background = element_blank())

#######
#plot changes in each chromatin density class normalised to whole nuclear volume for different conditions:

  datamapnucvol <-  sum(datamap$classnum)/length(unique(datamap$Marker))
  datamap2nucvol <-  sum(datamap2$classnum)/length(unique(datamap2$Marker))
  datamap3nucvol <-  sum(datamap3$classnum)/length(unique(datamap3$Marker))
  datamap4nucvol <-  sum(datamap4$classnum)/length(unique(datamap4$Marker))

datamapa <- transform(datamap, classvolnorm = datamap$classnum/datamapnucvol)
datamap2a <- transform(datamap2, classvolnorm = datamap2$classnum/datamap2nucvol)
datamap3a <- transform(datamap3, classvolnorm = datamap3$classnum/datamap3nucvol)
datamap4a <- transform(datamap4, classvolnorm = datamap4$classnum/datamap4nucvol)

datamapcompa <- rbind(datamapa,datamap2a, datamap3a)

ggplot(datamapcompa, aes(factor(datamapcompa$class), datamapcompa$classvolnorm, fill=datamapcompa$condition))+
  geom_boxplot()+# coord_flip()
  labs(
    x = "Chromatin class",
    y = "Volume : Nuclear volume",
    fill = "Chromatin"
  )+
  #theme_bw()+
  theme(axis.line = element_line(colour = "black"),
        panel.background = element_blank())
#scale_x_discrete(limits=c("G1","NaN3","Hyper-os"))+
#ylim(0,1)+
scale_fill_manual(breaks=c("G1","NaN3","Hyper-os"),values=c("#F8766D", "#00BA38","#619CFF"))
#theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=1))
#

  #to plot single condition
ggplot(datamapa, aes(factor(datamapa$class), datamapa$classvolnorm))+
  geom_boxplot(fill='#F8766D')+# coord_flip()
  labs(
    x = "Chromatin class",
    y = "Chromatin Volume : Nuclear volume",
    fill = "Chromatin"
  )+
  #theme_bw()+
  theme(axis.line = element_line(colour = "black"),
        panel.background = element_blank())
#scale_x_discrete(limits=c("G1","NaN3","Hyper-os"))+
#ylim(0,1)+
scale_fill_manual(breaks=c("G1","NaN3","Hyper-os"),values=c("#F8766D", "#00BA38","#619CFF"))
#theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=1))



#######
#plot changes in each chromatin density class normalised to whole nuclear volume for Sphase:
datamap<- read.csv(
  "/Users/eze/Dropbox (Schermelleh_Lab)/Schermelleh/EMPaper/graphs-for-paper/C127/2017.01-results/distribution/G1_distn-compile.csv")
datamap4<- read.csv(
  "/Users/eze/Dropbox (Schermelleh_Lab)/Schermelleh/EMPaper/graphs-for-paper/C127/2017.01-results/distribution/ES_distn-compile.csv")
datamap5<- read.csv(
  "/Users/eze/Dropbox (Schermelleh_Lab)/Schermelleh/EMPaper/graphs-for-paper/C127/2017.01-results/distribution/MS_distn-compile.csv")
datamap6<- read.csv(
  "/Users/eze/Dropbox (Schermelleh_Lab)/Schermelleh/EMPaper/graphs-for-paper/C127/2017.01-results/distribution/LS_distn-compile.csv")

datamapnucvol <-  sum(datamap$classnum)/length(unique(datamap$Marker))
datamap4nucvol <-  sum(datamap4$classnum)/length(unique(datamap4$Marker))
datamap5nucvol <-  sum(datamap5$classnum)/length(unique(datamap5$Marker))
datamap6nucvol <-  sum(datamap6$classnum)/length(unique(datamap6$Marker))

datamapa <- transform(datamap, classvolnorm = datamap$classnum/datamapnucvol)
datamap4a <- transform(datamap4, classvolnorm = datamap4$classnum/datamap4nucvol)
datamap5a <- transform(datamap5, classvolnorm = datamap5$classnum/datamap5nucvol)
datamap6a <- transform(datamap6, classvolnorm = datamap6$classnum/datamap6nucvol)

#choose to put G1 datamapa or not:
#datamapcompa <- rbind(datamapa,datamap4a, datamap5a, datamap6a)
datamapcompa <- rbind(datamapa, datamap4a, datamap5a, datamap6a)

ggplot(datamapcompa, aes(factor(datamapcompa$class), datamapcompa$classvolnorm, fill=datamapcompa$condition))+
  geom_boxplot()+# coord_flip()
  labs(
    x = "Chromatin class",
    y = "Volume : Nuclear volume",
    fill = "Chromatin"
  )+
  scale_fill_manual(values=c("#F8766D", "#D24F44", "#A62E25","#7B1C15"))+
  #scale_y_continuous(limits=c(-80, 140),breaks=seq(-80,140,20)) +
  scale_x_discrete(breaks=seq(1,7,1))+
  theme(axis.line = element_line(colour = "black"),
        panel.background = element_blank())
#scale_x_discrete(limits=c("G1","NaN3","Hyper-os"))+
#ylim(0,1)+
#scale_fill_manual(breaks=c("G1","NaN3","Hyper-os"),values=c("#F8766D", "#00BA38","#619CFF"))
#theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=1))

################
#To plot remodeller size to marker as well as peak density class
library("ggrepel")
set.seed(42)

markerlist =c("H4K20me3","H3K9me3","H3K9me2","H3K27me3","H4K5ac","H4K20me1","H3K36me2","H3K4me2","H3K36me3","H3K4me3","RNAP-S2P")

  remod <- read.csv("/Users/ezejm/Dropbox (Schermelleh_Lab)/Schermelleh/EMPaper/graphs-for-paper/C127/2018.08-Final_eroded-stiched/remod-size-table_Final.csv")

p <- ggplot(remod, aes(y=remod$Marker, x=remod$Size..kDa.))+
  geom_point()+
  geom_text_repel(aes(label=remod$Complex, color=factor(remod$Density.enrichment..peak.class.)))+
  scale_y_discrete(limits=markerlist)+
  labs(
    x = "Size of Complex (kDa)",
    y = "Marker",
    color = "Peak\nClass"
  )+
  #scale_color_brewer(palette="Dark2")+
  scale_color_hue(l=65, c=100)+
  theme(axis.line = element_line(colour = "black"),
        panel.background = element_blank())
  #geom_hline(yintercept = 4.5)+
  #geom_hline(yintercept = 6.5)
  
p

  ########
#To plot WF vs SIR class volumes
datamapcomp <- read.csv("WF-distribution-summary.csv")
voxnm3 <- 41*41*125
voxum3 <- voxnm3/1000000000
pd <- position_dodge(0.9)


p <- ggplot(datamapcomp, aes(x=class, y=classnumn, fill=Marker))+
  geom_bar(stat = "identity",position = pd)+# coord_flip()
  labs(
    x = "Chromatin class",
    y = "Relative class volume",
    fill = "Resolution"
  )+
  scale_fill_grey(start = 0.8, end = 0.2)+
  scale_x_continuous(breaks=seq(1,7,1))+
    theme(axis.line = element_line(colour = "black"),
        panel.background = element_blank())
#scale_x_discrete(limits=c("G1","NaN3","Hyper-os"))+
#ylim(0,1)+
#scale_fill_manual(values=c("#F8766D", "#00BA38","#619CFF"))
#theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=1))

p
savename <- paste0(condition,"_classnumn-bar.png")
ggsave(savename, plot = p)
