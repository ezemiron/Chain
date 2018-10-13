library(ggplot2)

SAdata <- read.csv("C127_G1_SA2vol_raw.csv")
##SAdata2 <- read.csv("/Documents/papers/chrom_marks/results/C127/ES/C127_ES_SA2vol_raw.csv")
##SAdata3 <- read.csv("/Users/eze/Desktop/Dropbox (Schermelleh_Lab)/Schermelleh_/EMPaper/graphs-for-paper/C127/2017.01-results/SA2vol/C127_Hyper-os_SA2vol_raw.csv")


SAdata$condition <- "G1"
##SAdata3$condition <- "Hyper-os"
##SAdata2$condition <- "NaN3"

##SAcomp <- rbind(SAdata,SAdata2,SAdata3)
SAcomp <- rbind(SAdata)


###plot SA 2 vol:
ggplot(SAcomp, aes(SAcomp$condition, SAcomp$SA2vol, fill=SAcomp$condition))+
  geom_boxplot()+# coord_flip()
  labs(
    x = "Chromatin conditions",
    y = "SA:Vol",
    fill = "Chromatin"
  )+
  theme_bw()+
  scale_x_discrete(limits=c("G1"))+
  #ylim(0,1)+
  scale_fill_manual(breaks=c("G1"),values=c("#F8766D"))
  #theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=1))


###plot chrom2vol
ggplot(SAcomp, aes(SAcomp$condition, SAcomp$chrom2vol, fill=SAcomp$condition))+
  geom_boxplot()+# coord_flip()
  labs(
    x = "Chromatin conditions",
    y = "chromatin:nucleus",
    fill = "Chromatin"
  )+
  theme_bw()+
  scale_x_discrete(limits=c("G1"))+
  #ylim(0,1)+
  scale_fill_manual(breaks=c("G1"),values=c("#F8766D"))
#theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=1))



###plot whole nuclear volume:
ggplot(SAcomp, aes(SAcomp$condition, SAcomp$nucvol, fill=SAcomp$condition))+
  geom_boxplot()+# coord_flip()
  labs(
    x = "Chromatin conditions",
    y = "Nuclear volume",
    fill = "Chromatin"
  )+
  theme_bw()+
  scale_x_discrete(limits=c("G1"))+
  #ylim(0,1)+
  scale_fill_manual(breaks=c("G1"),values=c("#F8766D"))
#theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=1))



#############

#plot absolute changes in each chromatin density class
datamap<- read.csv("G1_distn-compile.csv")

datamapcomp <- rbind(datamap)

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
  scale_fill_manual(breaks=c("G1"),values=c("#F8766D"))
#theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=1))



#######
#plot changes in each chromatin density class normalised to whole nuclear volume for different conditions:

  datamapnucvol <-  sum(datamap$classnum)/length(unique(datamap$Marker))

  
datamapa <- transform(datamap, classvolnorm = datamap$classnum/datamapnucvol)


datamapcompa <- rbind(datamapa)

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
scale_fill_manual(breaks=c("G1"),values=c("#F8766D"))
#theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=1))
#



#######
#plot changes in each chromatin density class normalised to whole nuclear volume for Sphase:
datamap<- read.csv("G1_distn-compile.csv")

datamapnucvol <-  sum(datamap$classnum)/length(unique(datamap$Marker))


datamapa <- transform(datamap, classvolnorm = datamap$classnum/datamapnucvol)


datamapcompa <- rbind(datamapa)

ggplot(datamapcompa, aes(factor(datamapcompa$class), datamapcompa$classvolnorm, fill=datamapcompa$condition))+
  geom_boxplot()+# coord_flip()
  labs(
    x = "Chromatin class",
    y = "Volume : Nuclear volume",
    fill = "Chromatin"
  )+
  scale_fill_manual(values=c("#F8766D"))+
  #scale_y_continuous(limits=c(-80, 140),breaks=seq(-80,140,20)) +
  scale_x_discrete(breaks=seq(1,7,1))+
  theme(axis.line = element_line(colour = "black"),
        panel.background = element_blank())
#scale_x_discrete(limits=c("G1","NaN3","Hyper-os"))+
#ylim(0,1)+
#scale_fill_manual(breaks=c("G1","NaN3","Hyper-os"),values=c("#F8766D", "#00BA38","#619CFF"))
#theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=1))







#######
#plot changes in each chromatin network size:
library(plyr)

d2bnorm<- read.csv("G1_network-compile.csv")

netbind <- rbind(d2bnorm)

netbindply <- ddply(netbind, c ("Condition"), summarise, meanNet = mean(AvNetwork), MnegCI95 = mean(negCI95), MposCI95 = mean(posCI95))

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
    scale_x_continuous(limits=c(0, 400),breaks=seq(0,400,50),expand = c(0, 0))+
    scale_y_continuous(limits=c(0, 1),breaks=seq(0,1,0.1),expand = c(0, 0))+
    theme(axis.line = element_line(colour = "black"),
          panel.background = element_blank())

    
    
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
      negdist <- 2*which(closest(condtable$MnegCI95,negfreqNet)==condtable$MnegCI95)
      posdist <- 2*which(closest(condtable$MposCI95,posfreqNet)==condtable$MposCI95)
      truemeandist <- 2* truemeandist
      
      prep <- c(truemeandist,negdist,posdist)
      empty2<- rbind(empty2,prep)
    }
    
    empty2 <- as.data.frame(empty2[2:dim(empty2)[1],])
    rownames(empty2)<- c()
    colnames(empty2) <- c("meandist","negdist","posdist")
    Network <- cbind(empty2,condlist)
    
    write.csv(Network,file="C127-network-abs.csv")
    
    ggplot(Network, aes(x=condlist, y=meandist, colour=condlist))+
      geom_point()+
      scale_y_continuous(limits=c(100, 250),breaks=seq(100,250,25)) +
      geom_errorbar(aes(ymin=negdist, ymax=posdist, width=.1))+ 
    labs(
      x = "Chromatin conditions",
      y = "Distance (nm)",
      fill = "Chromatin"
    )+
    theme(axis.line = element_line(colour = "black"),panel.background = element_blank())+
      scale_x_discrete(limits=condlist)+
      scale_color_manual(breaks=c("G1"),values=c("#F8766D""))
    
