library(tidyr)
library(gplots)

cell <- "C127"
dirchosen <- paste0("/Users/eze/Desktop/Dropbox (Schermelleh_Lab)/heatmaps/",cell)
setwd(dirchosen)

#hmcol <- colorRampPalette(brewer.pal(9, "Blues"))(255)
#hmcol<-colorRampPalette(c('#000000','#000000','#000000','#000000','#000000','#2F2C62','#2F2C62','#2F2C62','#2F2C62','#2F2C62','#2F2C62','#2F2C62','#2F2C62','#2F2C62','#2F2C62', '#42399B', '#4A52A7', '#59AFEA', '#7BCEB8','#FFFF88', '#F5952D', '#E93131', '#D70131', '#D70131', '#D70131', '#D70131', '#D70131', '#D70131', '#D70131', '#D70131', '#D70131', '#D70131', '#D70131', '#D70131', '#D70131', '#D70131','#D70131','#D70131'))(256)
hmcol<-colorRampPalette(c('#D70131','#D70131','#D70131', '#D70131', '#D70131', '#D70131', '#D70131', '#D70131', '#D70131', '#D70131', '#D70131', '#D70131', '#D70131', '#D70131','#D70131','#D70131','#E93131','#F5952D','#FFFF88','#7BCEB8','#59AFEA','#4A52A7','#42399B','#2F2C62','#2F2C62','#2F2C62','#2F2C62','#2F2C62','#2F2C62','#2F2C62','#2F2C62','#2F2C62','#2F2C62','#000000','#000000','#000000','#000000','#000000'))(256)

#datamap<- read.csv("G1_distn-compile.csv")
datamap<- read.csv("G1_distn-compile.csv", colClasses = c("NULL","NULL", NA, NA,NA, "NULL", "NULL", "NULL", "NULL", "NULL"))
datamap<- tidyr::spread(datamap,class,Log2AvNorm)
#datamapb <- cbind(datamap$Marker, datamap$class, datamap$Log2AvNorm)

#datamapmat <- read.csv("G1_class-dist.csv")
datamapmat<- datamap[,2:length(datamap)]
datamapmat <- as.matrix(datamapmat)
rownames(datamapmat) <- datamap$Marker
colnames(datamapmat) <- c(1,2,3,4,5,6,7)

datamapmat <- as.matrix(datamapmat)
datamapmat <- as.numeric(datamapmat)

#parameters need tweaking to make pretty:
heatmap.2(datamapmat2, Colv = FALSE, Rowv = FALSE, trace="none", col = rev(hmcol), margin=c(13, 13))
heatmap.2(datamapmat2, Colv = FALSE, trace="none", col = rev(hmcol), margin=c(13, 13))
#heatmap.2(datamapmat2, Rowv = FALSE, trace="none", col = rev(hmcol), margin=c(13, 13))
