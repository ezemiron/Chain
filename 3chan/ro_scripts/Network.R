#plot changes in each chromatin network size:
library(plyr)
library(ggplot2)

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


 
   
