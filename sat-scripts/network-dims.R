library(ggplot2)
library(plyr)

#set wd eg for C127:
setwd("/Users/ezejm/Dropbox (Schermelleh_Lab)/Schermelleh/EMPaper/graphs-for-paper/C127/2017.01-results/SA2vol/")

#######
#plot changes in each chromatin network size:


d2bnorm<- read.csv(
  "/Users/ezejm/Dropbox (Schermelleh_Lab)/Schermelleh/EMPaper/graphs-for-paper/C127/2017.01-results/SA2vol/G1_network-compile.csv")
d2bnorm2<- read.csv(
  "/Users/ezejm/Dropbox (Schermelleh_Lab)/Schermelleh/EMPaper/graphs-for-paper/C127/2017.01-results/SA2vol/NaN3_network-compile.csv")
d2bnorm3<- read.csv(
  "/Users/ezejm/Dropbox (Schermelleh_Lab)/Schermelleh/EMPaper/graphs-for-paper/C127/2017.01-results/SA2vol/Hyper-os_network-compile.csv")
d2bnorm4 <- read.csv(
  "/Users/ezejm/Dropbox (Schermelleh_Lab)/Schermelleh/EMPaper/graphs-for-paper/C127/2018.08-Final_eroded-stiched/Trip/Trip_network-compile.csv")

netbind <- rbind(d2bnorm,d2bnorm2,d2bnorm3)

netbindply <- ddply(netbind, c ("Condition","stddist"), summarise, meanNet = mean(AvNetwork), MnegCI95 = mean(negCI95), MposCI95 = mean(posCI95))

p<-ggplot(data=netbindply, aes(x=stddist, y=meanNet, colour=Condition, fill=Condition )) + 
  geom_line()+
  geom_ribbon(aes(ymin= MnegCI95, ymax=MposCI95), linetype=2, alpha=0.3)+
  #geom_hline(yintercept=0.789)+
  #geom_vline(xintercept=82.314)+
  labs(
    x = "Distance to border (nm)",
    y = "Frequency"
  )+
  #ggtitle(marker)+
  scale_x_continuous(limits=c(0, 400),breaks=seq(0,400,50),expand = c(0, 0))+
  scale_y_continuous(limits=c(0, 1),breaks=seq(0,1,0.1),expand = c(0, 0))+
  theme(axis.line = element_line(colour = "black"),
        panel.background = element_blank())

savename <- "C127-rndm-d2b-line.png"
ggsave(savename, plot = p)

  #to plot single condition
netbindply2 <- netbindply[netbindply$Condition %in% "G1", ]

p<-ggplot(data=netbindply2, aes(x=stddist, y=meanNet, colour=Condition, fill=Condition )) + 
  geom_line()+
  geom_ribbon(aes(ymin= MnegCI95, ymax=MposCI95), linetype=2, alpha=0.3)+
  #geom_hline(yintercept=0.789)+
  #if radius was calculated then the delta distance can be plotted as vertical line:
  #geom_vline(xintercept=empty[1])+
  labs(
    x = "Distance to border (nm)",
    y = "Frequency"
  )+
  #ggtitle(marker)+
  scale_x_continuous(limits=c(0, 400),breaks=seq(0,400,50),expand = c(0, 0))+
  scale_y_continuous(limits=c(0, 1),breaks=seq(0,1,0.1),expand = c(0, 0))+
  theme(axis.line = element_line(colour = "black"),
        panel.background = element_blank())

######
#For working out network radius
empty <- NA
#condlist <- c("G1","NaN3","Hyper-os")
condlist <- c("G1","Trip")
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

#define scale factor by modelling tubular network
ratio <- 2+sqrt(2)

#define closest function to figure out error bands
closest<-function(xv,sv){
  xv[which(abs(xv-sv)==min(abs(xv-sv)))] }


for (cond in condlist){
  condtable<- netbindply[netbindply$Condition %in% cond, ]
  meandist <- as.numeric(meandists[which(cond==meandists[,1]),2])
  truemeandist <- closest(condtable$stddist,meandist)
  
  posfreqNet <- condtable[truemeandist,4]
  negfreqNet <- condtable[truemeandist,5]
  negdist <- ratio*which(closest(condtable$MnegCI95,negfreqNet)==condtable$MnegCI95)
  posdist <- ratio*which(closest(condtable$MposCI95,posfreqNet)==condtable$MposCI95)
  truemeandist <- ratio* truemeandist
  
  prep <- c(truemeandist,negdist,posdist)
  empty2<- rbind(empty2,prep)
}

empty2 <- as.data.frame(empty2[2:dim(empty2)[1],])
rownames(empty2)<- c()
colnames(empty2) <- c("meandist","negdist","posdist")
Network <- cbind(empty2,condlist)

write.csv(Network,file="C127-network-rad-abs.csv")

######
#to plot network radius:

Network<- read.csv(
  "/Users/ezejm/Dropbox (Schermelleh_Lab)/Schermelleh/EMPaper/graphs-for-paper/C127/2018.08-Final_eroded-stiched/C127-network-rad-abs_Final.csv")

p <- ggplot(Network, aes(x=condlist, y=meandist, colour=condlist))+
  geom_point()+
  scale_y_continuous(limits=c(100, 500),breaks=seq(100,500,50)) +
  geom_errorbar(aes(ymin=negdist, ymax=posdist, width=.1))+ 
  labs(
    x = "Chromatin conditions",
    y = "Network radius (nm)",
    fill = "Chromatin"
  )+
  theme(axis.line = element_line(colour = "black"),panel.background = element_blank())+
  scale_x_discrete(limits=condlist)+
  #scale_colour_grey(breaks=c("G1","NaN3","Hyper-os"),start = 0.5, end = 0.1)
  scale_color_manual(breaks=c("G1","NaN3","Hyper-os"),values=c("#7f7f7f","#191919","#5e5e5e"))
p
savename <- "C127-net-rad-abs.png"
ggsave(savename, plot = p)

  #to plot single condition
Networka <- Network[Network$condlist %in% "G1", ]

p <- ggplot(Networka, aes(x=condlist, y=meandist, colour=condlist))+
  geom_point()+
  scale_y_continuous(limits=c(100, 500),breaks=seq(100,500,50)) +
  geom_errorbar(aes(ymin=negdist, ymax=posdist, width=.1))+ 
  labs(
    x = "Chromatin condition",
    y = "Network radius (nm)",
    fill = "Chromatin"
  )+
  theme(axis.line = element_line(colour = "black"),panel.background = element_blank())+
  scale_x_discrete(limits="G1")+
  scale_color_manual(breaks=c("G1","NaN3","Hyper-os"),values=c("#F8766D", "#619CFF","#00BA38"))
  #export as 300 by 358 ratio

######
#to plot network diameter:
Network2 <- empty2*2
Network2 <- cbind(Network2,condlist)

write.csv(Network2,file="C127-network-diam-abs.csv")

p<-ggplot(Network2, aes(x=condlist, y=meandist, colour=condlist))+
  geom_point()+
  scale_y_continuous(limits=c(200, 1000),breaks=seq(200,1000,100)) +
  geom_errorbar(aes(ymin=negdist, ymax=posdist, width=.1))+ 
  labs(
    x = "Chromatin conditions",
    y = "Network diameter (nm)",
    fill = "Chromatin"
  )+
  theme(axis.line = element_line(colour = "black"),panel.background = element_blank())+
  scale_x_discrete(limits=condlist)+
  scale_color_manual(breaks=c("G1","NaN3","Hyper-os"),values=c("#F8766D", "#619CFF","#00BA38"))

savename <- "C127-net-diam-abs.png"
ggsave(savename, plot = p)

  #to plot single condition
Network2a <- Network2[Network2$condlist %in% "G1", ]

p<-ggplot(Network2a, aes(x=condlist, y=meandist, colour=condlist))+
  geom_point()+
  scale_y_continuous(limits=c(200, 1000),breaks=seq(200,1000,100)) +
  geom_errorbar(aes(ymin=negdist, ymax=posdist, width=.1))+ 
  labs(
    x = "Chromatin condition",
    y = "Network diameter (nm)",
    fill = "Chromatin"
  )+
  theme(axis.line = element_line(colour = "black"),panel.background = element_blank())+
  scale_x_discrete(limits="G1")+
  scale_color_manual(breaks=c("G1","NaN3","Hyper-os"),values=c("#F8766D", "#619CFF","#00BA38"))
