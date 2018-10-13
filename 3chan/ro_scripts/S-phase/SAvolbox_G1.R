library(ggplot2)

#SAdata3 <- read.csv("C127_MS_SA2vol_raw.csv")
SAdata <- read.csv("C127_G1_SA2vol_raw.csv")
#SAdata2 <- read.csv("C127_ES_SA2vol_raw.csv")
#SAdata4 <- read.csv("C127_LS_SA2vol_raw.csv")


SAdata$condition <- "G1"
#SAdata2$condition <- "ES"
#SAdata3$condition <- "MS"
#SAdata4$condition <- "LS"

#SAcomp <- rbind(SAdata,SAdata2,SAdata3,SAdata4)
SAcomp <- SAdata


###plot SA 2 vol:
ggplot(SAcomp, aes(SAcomp$condition, SAcomp$SA2vol, fill=SAcomp$condition))+
  geom_boxplot()+# coord_flip()
  labs(
    x = "Chromatin conditions",
    y = "SA:Vol",
    fill = "Chromatin"
  )+
  theme(axis.line = element_line(colour = "black"),panel.background = element_blank(),legend.position="none")+
  scale_x_discrete(limits=c("G1"))+
  #ylim(0,1)+
  scale_fill_manual(breaks=c("G1"),values=c("#000000"))
  #theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=1))


###plot chrom2vol
ggplot(SAcomp, aes(SAcomp$condition, SAcomp$chrom2vol, fill=SAcomp$condition))+
  geom_boxplot()+# coord_flip()
  labs(
    x = "Chromatin conditions",
    y = "chromatin:nucleus",
    fill = "Chromatin"
  )+
  theme(axis.line = element_line(colour = "black"),panel.background = element_blank(),legend.position="none")+
  scale_x_discrete(limits=c("G1"))+
  #ylim(0,1)+
  scale_fill_manual(breaks=c("G1"),values=c("#000000"))



###plot whole nuclear volume:
ggplot(SAcomp, aes(SAcomp$condition, SAcomp$nucvol, fill=SAcomp$condition))+
  geom_boxplot()+# coord_flip()
  labs(
    x = "Chromatin conditions",
    y = "Nuclear volume",
    fill = "Chromatin"
  )+
  theme(axis.line = element_line(colour = "black"),panel.background = element_blank(),legend.position="none")+
  #theme_bw()+
  scale_x_discrete(limits=c("G1"))+
  #ylim(0,1)+
  scale_fill_manual(breaks=c("G1"),values=c("#000000"))+
  theme(axis.line = element_line(colour = "black"),
           panel.background = element_blank(),legend.position = "none")

#theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=1))



#############

#plot absolute changes in each chromatin density class
datamap<- read.csv("G1_distn-compile.csv")
#datamap2<- read.csv("ES_mask_distn-compile.csv")
#datamap3<- read.csv("MS_mask_distn-compile.csv")
#datamap4<- read.csv("LS_mask_distn-compile.csv")

#datamapcomp <- rbind(datamap,datamap2,datamap3,datamap4)
datamapcomp <-datamap

ggplot(datamapcomp, aes(factor(datamapcomp$class), datamapcomp$classnum, fill=datamapcomp$condition))+
  geom_boxplot()+# coord_flip()
  labs(
    x = "Chromatin class",
    y = "Absolute class volume",
    fill = "Chromatin"
  )+
  theme(axis.line = element_line(colour = "black"),panel.background = element_blank(),legend.position="none")+
  #theme_bw()+
  #scale_x_discrete(limits=c("G1"))+
  #ylim(0,1)+
  scale_fill_manual(breaks=c("G1"),values=c("#000000"))


#theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=1))



#######
#plot changes in each chromatin density class normalised to whole nuclear volume for different conditions:

  datamapnucvol <-  sum(datamap$classnum)/length(unique(datamap$Marker))
#  datamap2nucvol <-  sum(datamap2$classnum)/length(unique(datamap2$Marker))
#  datamap3nucvol <-  sum(datamap3$classnum)/length(unique(datamap3$Marker))
#  datamap4nucvol <-  sum(datamap4$classnum)/length(unique(datamap4$Marker))  


datamapa <- transform(datamap, classvolnorm = datamap$classnum/datamapnucvol)
#datamap2a <- transform(datamap2, classvolnorm = datamap2$classnum/datamap2nucvol)
#datamap3a <- transform(datamap3, classvolnorm = datamap3$classnum/datamap3nucvol)
#datamap4a <- transform(datamap4, classvolnorm = datamap4$classnum/datamap4nucvol)

#datamapcompa <- rbind(datamapa,datamap2a,datamap3a,datamap4a)
datamapcompa <- datamapa

ggplot(datamapcompa, aes(factor(datamapcompa$class), datamapcompa$classvolnorm, fill=datamapcompa$condition))+
  geom_boxplot()+# coord_flip()
  labs(
    x = "Chromatin class",
    y = "Volume : Nuclear volume",
    fill = "Chromatin"
  )+
  #theme_bw()+
  theme(axis.line = element_line(colour = "black"),panel.background = element_blank(),legend.position="none")+
#scale_x_discrete(limits=c("G1","ES","MS","LS"))+
#ylim(0,1)+
scale_fill_manual(breaks=c("G1"),values=c("#000000"))

#theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=1))
#

######
#plot changes in each chromatin network size:
library(plyr)

d2bnorm<- read.csv("G1_network-compile.csv")
#d2bnorm2<- read.csv("ES_mask_network-compile.csv")
#d2bnorm3<- read.csv("MS_mask_network-compile.csv")
#d2bnorm4<- read.csv("LS_mask_network-compile.csv")

#netbind <- rbind(d2bnorm,d2bnorm2,d2bnorm3,d2bnorm4)
netbind <- d2bnorm

netbindply <- ddply(netbind, c ("Condition","stddist"), summarise, meanNet = mean(AvNetwork), MnegCI95 = mean(negCI95), MposCI95 = mean(posCI95))

    p<-ggplot(data=netbindply, aes(x=stddist, y=meanNet, colour=Condition, fill=Condition )) + 
    geom_line()+
    geom_ribbon(aes(ymin= MnegCI95, ymax=MposCI95), linetype=2, alpha=0.3)+
    #geom_hline(yintercept=0.789)+
    #geom_vline(xintercept=82.314)+
    labs(
      x = "Half-network width (nm)",
      y = "Frequency"
    )+
    #ggtitle(marker)+
    scale_x_continuous(limits=c(0, 1000),breaks=seq(0,1000,50),expand = c(0, 0))+
    scale_y_continuous(limits=c(0, 1),breaks=seq(0,1,0.1),expand = c(0, 0))+
theme(axis.line = element_line(colour = "black"),panel.background = element_blank(),legend.position="none")
    
     closest<-function(xv,sv){
     xv[which(abs(xv-sv)==min(abs(xv-sv)))] }
     
    empty <- NA
    condlist <- c("G1")
       for (cond in condlist){
       condtable<- netbindply[netbindply$Condition %in% cond, ]
       hist <- (condtable$stddist+1)*condtable$meanNet
        #hist <- hist -1
       index <- mean(hist)

       halfnet <- sum(hist)/sum (condtable$meanNet)
       empty <- cbind(empty,halfnet)  
       }
    
    empty <- empty[,2:length(empty)]
    meandists <- cbind(condlist,as.vector(empty))
    
    empty2 <- NA
    empty2 <- cbind(empty2,empty2,empty2)
    for (cond in condlist){
      condtable<- netbindply[netbindply$Condition %in% cond, ]
      meandist <- as.numeric(meandists[which(cond==meandists[,1]),2])
      truemeandist <- closest(condtable$stddist,meandist)
                             
      posfreqNet <- condtable[truemeandist,4]
      negfreqNet <- condtable[truemeandist,5]
      negdist <- (2+sqrt(2))*which(closest(condtable$MnegCI95,negfreqNet)==condtable$MnegCI95)
      posdist <- (2+sqrt(2))*which(closest(condtable$MposCI95,posfreqNet)==condtable$MposCI95)
      truemeandist <- (2+sqrt(2))*truemeandist
      
      prep <- c(truemeandist,negdist,posdist)
      empty2<- rbind(empty2,prep)
    }
    
    empty2 <- as.data.frame(empty2[2:dim(empty2)[1],])
    rownames(empty2)<- c()
    colnames(empty2) <- c("meandist","negdist","posdist")
    Network <- cbind(empty2,condlist)
    
    write.csv(Network,file="C127-G1-network-abs.csv")
    
    ggplot(Network, aes(x=condlist, y=meandist, colour=condlist))+
      geom_point()+
      scale_y_continuous(limits=c(200, 350),breaks=seq(200,350,25)) +
      geom_errorbar(aes(ymin=negdist, ymax=posdist, width=.1))+ 
    labs(
      x = "Chromatin conditions",
      y = "Distance (nm)",
      fill = "Chromatin"
    )+
    theme(axis.line = element_line(colour = "black"),panel.background = element_blank(),legend.position="none")+
      scale_x_discrete(limits=condlist)+
      scale_color_manual(breaks=c("G1"),values=c("#000000"))
    
    
    
