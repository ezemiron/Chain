library(ggplot2)
#plot absolute changes in each chromatin density class
datamap<- read.csv("ES_distn-compile.csv")
condition<-("ES")
savename<-paste0(condition,"_distribution.pdf")

datamapcomp <- rbind(datamap)

p <- ggplot(datamapcomp, aes(factor(datamapcomp$class), datamapcomp$classnumn, fill=datamapcomp$condition))+
  geom_boxplot()+# coord_flip()
  labs(
    x = "Chromatin class",
    y = "Chromatin : Nuclear volume",
    fill = "Chromatin"
  )+
  #scale_fill_manual(values=c("#F8766D", "#D24F44", "#A62E25"))+
  #scale_y_continuous(limits=c(-80, 140),breaks=seq(-80,140,20)) +
  scale_x_discrete(breaks=seq(1,7,1))+
  #theme_bw()+
  theme(axis.line = element_line(colour = "black"),
        panel.background = element_blank())
#scale_x_discrete(limits=c("G1","NaN3","Hyper-os"))+
#ylim(0,1)+
#scale_fill_manual(breaks=c("G1","NaN3","Hyper-os"),values=c("#F8766D", "#00BA38","#619CFF"))
#theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=1))

ggsave(savename, device="pdf", plot = p)
