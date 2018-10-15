# Script to test LUTs for heatmaps from simulated data:

library(ggplot2)
library(scales)
library(RColorBrewer)

#colorRampPalette(brewer.pal(n_palette, "palette_name"))(n_plot)
#where palette_name is the name of the preset palette
#and n_palette and n_plot is the number of colours wanted from the palette and into the plot, respectively
hmcols3 = colorRampPalette(brewer.pal(11,"RdYlBu"))(100)


b <- c(seq(from=-2, to=2, by =1 ))

pp <- function (n,r=4) {
  x <- seq(-r*pi, r*pi, len=n)
  df <- expand.grid(x=x, y=x)
  df$r <- sqrt(df$x^2 + df$y^2)
  df$z <- (cos(df$r^10)*exp(-df$r/50))*6
  df
}

# ggplot(pp(15),aes(x=x,y=y))+geom_tile(aes(fill=z))+
# scale_fill_gradientn(values=c(1,0), limits =c(-2,2), colours = hmcols3, breaks=b, name="Log2")

M <- pp(15)
M["z"] <- lapply(M["z"], function(x) ifelse(x<(-2), (-2.5), x))
M["z"] <- lapply(M["z"], function(x) ifelse(x>2, 2.5, x))

ggplot(M,aes(x=x,y=y))+geom_tile(aes(fill=z))+
  scale_fill_gradientn(values=c(1,0), limits =c(-2.5,2.5), colours = hmcols3, breaks=b, name="Log2")






#######
# Trial with real data

setwd("/Users/ezejm/Dropbox (Schermelleh_Lab)/Schermelleh/EMPaper/graphs-for-paper/C127/2018.08-Final_eroded-stiched")

#all markers
markerlist =c("aHP1","H4K20me3","H3K9me3","H3K9me2","H3K27me3","H4K20me1","Macro-H2A","H4K5ac","H3K36me2","H3K4me2","H3K36me3","H3K4me3","Smc3","Scc1","CTCF","SAF-A","30min-BrUTP","RNAP-S2P","hnRNP")

datamap<- read.csv(
  "/Users/ezejm/Dropbox (Schermelleh_Lab)/Schermelleh/EMPaper/graphs-for-paper/C127/2018.08-Final_eroded-stiched/G1_distn-compile_volnormE_Final.csv")
datamap2<- read.csv(
  "/Users/ezejm/Dropbox (Schermelleh_Lab)/Schermelleh/EMPaper/graphs-for-paper/C127/2018.08-Final_eroded-stiched/NaN3_distn-compile_volnormE_Final.csv")
datamap3<- read.csv(
  "/Users/ezejm/Dropbox (Schermelleh_Lab)/Schermelleh/EMPaper/graphs-for-paper/C127/2018.08-Final_eroded-stiched/Hyper-os_distn-compile_volnormE_Final.csv")

#for condensation:
distnbind <- rbind(datamap,datamap2,datamap3)


#######
#first solution:
hmcols10 = colorRampPalette(brewer.pal(11,"RdYlBu"))(100)

b <- c(seq(from=-3, to=2, by =1 ))
# correcting extreme values for better LUT display
distnbind["Log2AvNorm"] <- lapply(distnbind["Log2AvNorm"], function(x) ifelse(x<(-3.5), (-4), x))
distnbind["Log2AvNorm"] <- lapply(distnbind["Log2AvNorm"], function(x) ifelse(x>1.5, 2, x))

ggplot(distnbind,aes(x = class, y = Marker, fill = Log2AvNorm)) +
  geom_tile()+
  facet_grid(. ~ condition)+
  geom_vline(xintercept = 1.5)+
  geom_hline(yintercept = 3.5)+
  geom_hline(yintercept = 7.5)+
  geom_hline(yintercept = 11.5)+
  scale_fill_gradientn(values=rescale(c(2,0,-4)), 
                       limits =c(-4,2), 
                       colours = c("red","white","blue"), 
                       breaks=b, name="Log2")+
  scale_x_discrete(limits=c(1,2,3,4,5,6,7),expand = c(0, 0))+
  #scale_x_discrete(expand = c(0, 0))+
  scale_y_discrete(limits=markerlist,expand = c(0, 0))+
  #scale_y_discrete(expand = c(0, 0))+
  coord_equal()+
  ggtitle("C127-Compaction")+
  theme_bw()



#######
#Second solution:
hmcolsA = colorRampPalette(brewer.pal(9,"YlOrRd"))(100)
hmcolsB = colorRampPalette(brewer.pal(9,"Blues"))(100)
hmcolsA <- rev(hmcolsA)
hmcols10 <- c(hmcolsA, hmcolsB)
#hmcols10[1]  is "#800026" dark red
#hmcols10[200]  is "#08306B" dark blue

b <- c(seq(from=-3, to=2, by =1 ))
# correcting extreme values for better LUT display
distnbind["Log2AvNorm"] <- lapply(distnbind["Log2AvNorm"], function(x) ifelse(x<(-3.3), (-3.3), x))
distnbind["Log2AvNorm"] <- lapply(distnbind["Log2AvNorm"], function(x) ifelse(x>1, 1.1, x))

ggplot(distnbind,aes(x = class, y = Marker, fill = Log2AvNorm)) +
  geom_tile()+
  facet_grid(. ~ condition)+
  geom_vline(xintercept = 1.5)+
  geom_hline(yintercept = 3.5)+
  geom_hline(yintercept = 7.5)+
  geom_hline(yintercept = 11.5)+
  scale_fill_gradientn(values=rescale(c(1.3,0,-3.3)), 
                       limits =c(-3.3,1.3), 
                       colours = hmcols10, 
                       breaks=b, name="Log2")+
  scale_x_discrete(limits=c(1,2,3,4,5,6,7),expand = c(0, 0))+
  #scale_x_discrete(expand = c(0, 0))+
  scale_y_discrete(limits=markerlist,expand = c(0, 0))+
  #scale_y_discrete(expand = c(0, 0))+
  coord_equal()+
  ggtitle("C127-Compaction")+
  theme_bw()


############## 
#third solution:

#Reds:
hmcolsA = colorRampPalette(brewer.pal(9,"YlOrRd"))(100)
hmcolsA <- rev(hmcolsA)
#Blues:
hmcolsB = colorRampPalette(brewer.pal(9,"Blues"))(100)

#hmcols10[200]  is "#08306B" dark blue needs to extend to "#000000" black
hmcolsC<-colorRampPalette(c('#000000','#08306B'))(50)
hmcolsC <- rev(hmcolsC)
#hmcols10[1]  is "#800026" dark red needs to extend to "#FF00FF" Magenta
hmcolsD<-colorRampPalette(c('#800026','#800026'))(50)
hmcolsD <- rev(hmcolsD)
#Stitch LUT:
hmcols10 <- c(hmcolsD, hmcolsA, hmcolsB,hmcolsC)

distnbind <- rbind(datamap,datamap2,datamap3)
b <- c(seq(from=-4, to=2, by =1 ))
# correcting extreme values for better LUT display
distnbind["Log2AvNorm"] <- lapply(distnbind["Log2AvNorm"], function(x) ifelse(x<(-4.8), (-4.8), x))
distnbind["Log2AvNorm"] <- lapply(distnbind["Log2AvNorm"], function(x) ifelse(x>1.6, 1.6, x))

ggplot(distnbind,aes(x = class, y = Marker, fill = Log2AvNorm)) +
  geom_tile()+
  facet_grid(. ~ condition)+
  geom_vline(xintercept = 1.5)+
  geom_hline(yintercept = 3.5)+
  geom_hline(yintercept = 6.5)+
  geom_hline(yintercept = 12.5)+
  geom_hline(yintercept = 15.5)+
  scale_fill_gradientn(values=rescale(c(2.3,0,-4.8)), 
                       limits =c(-4.8,2.3), 
                       colours = hmcols10, 
                       breaks=b, name="Log2")+
  scale_x_discrete(limits=c(1,2,3,4,5,6,7),expand = c(0, 0))+
  #scale_x_discrete(expand = c(0, 0))+
  scale_y_discrete(limits=markerlist,expand = c(0, 0))+
  #scale_y_discrete(expand = c(0, 0))+
  coord_equal()+
  ggtitle("C127-Compaction")+
  theme_bw()



###########
#new markerlist for plotting smaller subset of markers:

distnbind<- read.csv(
  "/Users/ezejm/Dropbox (Schermelleh_Lab)/Schermelleh/EMPaper/graphs-for-paper/C127/2018.08-Final_eroded-stiched/Trip/Trip_distn-compile.csv")

markerlist =c("H4K20me3","H3K9me3","H3K27me3","H3K4me3","H3K36me3","H3K4me2","CTCF","Scc1","RNAP-S2P","RNAP","hnRNP")


ggplot(distnbind,aes(x = class, y = Marker, fill = Log2AvNorm)) +
  geom_tile()+
  facet_grid(. ~ condition)+
  geom_vline(xintercept = 1.5)+
  geom_hline(yintercept = 2.5)+
  geom_hline(yintercept = 3.5)+
  geom_hline(yintercept = 6.5)+
  geom_hline(yintercept = 8.5)+
  #geom_hline(yintercept = 15.5)+
  scale_fill_gradientn(values=rescale(c(2.3,0,-4.8)), 
                       limits =c(-4.8,2.3), 
                       colours = hmcols10, 
                       breaks=b, name="Log2")+
  scale_x_discrete(limits=c(1,2,3,4,5,6,7),expand = c(0, 0))+
  #scale_x_discrete(expand = c(0, 0))+
  scale_y_discrete(limits=markerlist,expand = c(0, 0))+
  #scale_y_discrete(expand = c(0, 0))+
  coord_equal()+
  ggtitle("C127")+
  theme_bw()
