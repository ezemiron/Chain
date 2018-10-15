#install.packages("ggplot2")
library(ggplot2)
library(RColorBrewer)
library(scales)

setwd("/Users/eze/Dropbox (Schermelleh_Lab)/Schermelleh/EMPaper/graphs-for-paper/C127/2017.01-results/distribution/")

#for IMR90
setwd("/Users/eze/Dropbox (Schermelleh_Lab)/Schermelleh/EMPaper/graphs-for-paper/IMR90/2017.01-results/distribution/")

# hmcols10<-colorRampPalette(c('#000000','#000000','#000000','#000000','#000000','#2F2C62','#2F2C62','#2F2C62','#2F2C62','#2F2C62','#2F2C62','#2F2C62','#2F2C62','#2F2C62','#2F2C62', '#42399B', '#4A52A7', '#59AFEA', '#7BCEB8','#FFFF88', '#F5952D', '#E93131', '#D70131', '#D70131', '#D70131', '#D70131', '#D70131', '#D70131', '#D70131', '#D70131', '#D70131', '#D70131', '#D70131', '#D70131', '#D70131', '#D70131','#D70131','#D70131'))(256)
# hmcols10 = colorRampPalette(brewer.pal(11,"RdYlBu"))(100)
hmcolsA = colorRampPalette(brewer.pal(9,"YlOrRd"))(100)
hmcolsB = colorRampPalette(brewer.pal(9,"Blues"))(100)
hmcolsA <- rev(hmcolsA)
hmcols10 <- c(hmcolsA, hmcolsB)

b <- c(seq(from=-3, to=2, by =1 ))
#to set the markers manually you have to make a list in reverse order, the first will appear at the bottom of the geom_tile:
#unique(datamap$Marker)
markerlist =c("aHP1","H3K9me3","H3K9me2","H3K27me3","Macro-H2A","H4K5ac","H4K20me1","H3K36me2","H3K36me3","H3K4me2","H3K4me3","RNAP-S2P","Scc1","Smc3","CTCF","SAF-A")

#marker list for IMR90:
markerlist =c("H3K9me3","H3K9me2","H3K27me3","H4K5ac","H3K36me2","H3K36me3","H3K4me2","H3K4me3","RNAP-S2P","Smc3","SAF-A")


############
#to plot parts of markerlist heatmap

datamap<- read.csv(
  "/Users/eze/Dropbox (Schermelleh_Lab)/Schermelleh/EMPaper/graphs-for-paper/C127/2017.01-results/distribution/G1_distn-compile.csv")
datamap2<- read.csv(
  "/Users/eze/Dropbox (Schermelleh_Lab)/Schermelleh/EMPaper/graphs-for-paper/C127/2017.01-results/distribution/NaN3_distn-compile.csv")
datamap3<- read.csv(
  "/Users/eze/Dropbox (Schermelleh_Lab)/Schermelleh/EMPaper/graphs-for-paper/C127/2017.01-results/distribution/Hyper-os_distn-compile.csv")
datamap4<- read.csv(
  "/Users/eze/Dropbox (Schermelleh_Lab)/Schermelleh/EMPaper/graphs-for-paper/C127/2017.01-results/distribution/ES_distn-compile.csv")
datamap5<- read.csv(
  "/Users/eze/Dropbox (Schermelleh_Lab)/Schermelleh/EMPaper/graphs-for-paper/C127/2017.01-results/distribution/MS_distn-compile.csv")
datamap6<- read.csv(
  "/Users/eze/Dropbox (Schermelleh_Lab)/Schermelleh/EMPaper/graphs-for-paper/C127/2017.01-results/distribution/LS_distn-compile.csv")

##for all
#distnbind <- rbind(datamap,datamap2,datamap3,datamap4,datamap5,datamap6)

##for S phase:
distnbind <- rbind(datamap6,datamap5,datamap4,datamap)
filedescribe <- "Sphase"

#for condensation:
distnbind <- rbind(datamap3,datamap2,datamap)
filedescribe <- "Compact"

condnames<- as.vector((unique(distnbind$condition)))
lineh <- length(condnames)-0.5
#min(datamap$Log2AvNorm)
#max(datamap$Log2AvNorm)

for (marker in markerlist){
  inddistnbind <- distnbind[distnbind$Marker %in% marker, ]

p <- ggplot(inddistnbind,aes(x = class, y = condition, fill = Log2AvNorm)) +
  geom_tile()+
  geom_vline(xintercept = 1.5)+
  geom_hline(yintercept = lineh)+
  scale_fill_gradientn(limits =c(-6,6), colours = hmcols10, breaks=b, name="Log2")+
  scale_x_discrete(limits=c(1,2,3,4,5,6,7),expand = c(0, 0))+
  #scale_x_discrete(expand = c(0, 0))+
  scale_y_discrete(limits=condnames,expand = c(0, 0))+
  #scale_y_discrete(expand = c(0, 0))+
  labs(
    x = "Chromatin intensity class",
    y = "Condition"
  )+
  coord_equal()+
  ggtitle(marker)+
  theme_bw()
  
savename <- paste0(marker,"_",filedescribe,"-heatmap.png")
ggsave(savename, plot = p)
}

############
#to plot whole markerlist heatmap
datamap<- read.csv(
  "/Users/eze/Dropbox (Schermelleh_Lab)/Schermelleh/EMPaper/graphs-for-paper/C127/2017.01-results/distribution/G1_distn-compile.csv")
datamap<- read.csv(
  "/Users/eze/Dropbox (Schermelleh_Lab)/Schermelleh/EMPaper/graphs-for-paper/IMR90/2017.01-results/distribution/G0_distn-compile.csv")

distnbind <- datamap

# correcting extreme values for better LUT display
distnbind["Log2AvNorm"] <- lapply(distnbind["Log2AvNorm"], function(x) ifelse(x<(-3.3), (-3.3), x))
distnbind["Log2AvNorm"] <- lapply(distnbind["Log2AvNorm"], function(x) ifelse(x>1, 1.1, x))


p <- ggplot(distnbind,aes(x = class, y = Marker, fill = Log2AvNorm)) +
  geom_tile()+
  facet_grid(. ~ Stain)+
  geom_vline(xintercept = 1.5)+
  geom_hline(yintercept = 4.5)+
  geom_hline(yintercept = 8.5)+
  geom_hline(yintercept = 12.5)+
  scale_fill_gradientn(values=rescale(c(1.3,0,-3.3)), 
                       limits =c(-3.3,1.3), 
                       colours = hmcols10, 
                       breaks=b, name="Log2")+
  scale_x_discrete(limits=c(1,2,3,4,5,6,7),expand = c(0, 0))+
  scale_y_discrete(limits=markerlist,expand = c(0, 0))+
  coord_equal()+
  ggtitle("C127-G1")+
  theme_bw()+
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank())

p

ggsave("C127-G1-unstitch_LUT.png", plot = p)



#for data stiched by different chromatin labels:
setwd("/Users/eze/Documents/Paper work/")

b <- c(seq(from=-3, to=2, by =1 ))
#all markers
markerlist =c("aHP1","H3K9me3","H3K9me2","H3K27me3","H4K20me1","Macro-H2A","H4K5ac","H3K36me2","H3K4me2","H3K36me3","H3K4me3","RNAP-S2P","Scc1","Smc3","CTCF","SAF-A")
markerlist =c("H4K20me3","aHP1","H3K9me3","H3K9me2","H3K27me3","H4K20me1","Macro-H2A","H4K5ac","H3K36me2","H3K4me2","H3K36me3","H3K4me3","RNAP-S2P","Scc1","Smc3","CTCF","SAF-A","hnRNP")

datamap<- read.csv(
  "/Users/eze/Documents/Paper work/G1_distn-compile_volnormE.csv")
datamap2<- read.csv(
  "/Users/eze/Documents/Paper work/NaN3_distn-compile_volnormE.csv")
datamap3<- read.csv(
  "/Users/eze/Documents/Paper work/Hyper-os_distn-compile_volnormE.csv")
datamap4<- read.csv(
  "/Users/eze/Documents/Paper work/ES_distn-compile_volnormE.csv")
datamap5<- read.csv(
  "/Users/eze/Documents/Paper work/MS_distn-compile_volnormE.csv")
datamap6<- read.csv(
  "/Users/eze/Documents/Paper work/LS_distn-compile_volnormE.csv")

##for S phase:
distnbind <- rbind(datamap,datamap4,datamap5,datamap6)

#for condensation:
distnbind <- rbind(datamap,datamap2,datamap3)


# correcting extreme values for better LUT display
distnbind["Log2AvNorm"] <- lapply(distnbind["Log2AvNorm"], function(x) ifelse(x<(-3.3), (-3.3), x))
distnbind["Log2AvNorm"] <- lapply(distnbind["Log2AvNorm"], function(x) ifelse(x>1, 1.1, x))

#change facet_grid to Stain or condition depending on
#if you are plotting differences in DAPI vs SYTOX to unstitch the data in datamap
#or if plotting the distnbind for S phase or condensation in distnbind
p <- ggplot(distnbind,aes(x = class, y = Marker, fill = Log2AvNorm)) +
  geom_tile()+
  facet_grid(. ~ condition)+
  geom_vline(xintercept = 1.5)+
  geom_hline(yintercept = 4.5)+
  geom_hline(yintercept = 8.5)+
  geom_hline(yintercept = 12.5)+
  scale_fill_gradientn(values=rescale(c(1.3,0,-3.3)), 
                       limits =c(-3.3,1.3), 
                       colours = hmcols10, 
                       breaks=b, name="Log2")+
  scale_x_discrete(limits=c(1,2,3,4,5,6,7),expand = c(0, 0))+
  scale_y_discrete(limits=markerlist,expand = c(0, 0))+
  coord_equal()+
  ggtitle("C127-Compaction")+
  theme_bw()+
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank())

p

  #save as stich or unstich for DAPI vs SYTOX comparison
ggsave("C127-comp-stitch_LUT.png", plot = p)

#save with conditions as facets
ggsave("C127-G1-unstitch_LUT.png", plot = p)



