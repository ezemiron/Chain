library(ggplot2)
##plot changes in each chromatin density class normalised to whole nuclear volume for Sphase:
##datamap <- read.csv("G1/G1_distn-compile_S-phase.csv")
##datamap4 <- read.csv("ES/ES_distn-compile_S-phase.csv")
##datamap5 <- read.csv("MS/MS_distn-compile_S-phase.csv")
##datamap6 <- read.csv("LS/LS_distn-compile_S-phase.csv")

##datamapnucvol <-  sum(datamap$classnum)/length(unique(datamap$Marker))
##datamap4nucvol <-  sum(datamap4$classnum)/length(unique(datamap4$Marker))
##datamap5nucvol <-  sum(datamap5$classnum)/length(unique(datamap5$Marker))
##datamap6nucvol <-  sum(datamap6$classnum)/length(unique(datamap6$Marker))

##datamapa <- transform(datamap, classvolnorm = datamap$classnum/datamapnucvol)
##datamap4a <- transform(datamap4, classvolnorm = datamap4$classnum/datamapnucvol)
##datamap5a <- transform(datamap5, classvolnorm = datamap5$classnum/datamapnucvol)
##datamap6a <- transform(datamap6, classvolnorm = datamap6$classnum/datamapnucvol)

##datamapcompa <- rbind(datamap4a, datamap5a, datamap6a)

datamap<- read.csv("ES/ES_distn-compile.csv")
##datamap2<- read.csv("18/18-distn-compile.csv")
##datamap3<- read.csv("LS/LS_distn-compile_S-phase.csv")

datamapcompa <- rbind(datamap)

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
