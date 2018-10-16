#install.packages("ggplot2")
library(ggplot2)

#Uncomment one of the following colour schemes
# hmcols<-colorRampPalette(c("blue","white","red"))(256)
# hmcols<-colorRampPalette(c('#2F2C62', '#42399B', '#4A52A7', '#59AFEA', '#7BCEB8', '#A7DA64','#EFF121', '#F5952D', '#E93131', '#D70131', '#D70131'))(256)
# hmcols<-colorRampPalette(c('#2F2C62','#2F2C62','#2F2C62','#2F2C62', '#42399B', '#4A52A7', '#59AFEA', '#7BCEB8', '#A7DA64','#EFF121', '#F5952D', '#E93131', '#D70131', '#D70131', '#D70131', '#D70131'))(256)
# hmcols<-colorRampPalette(c('#2F2C62', '#42399B', '#4A52A7', '#59AFEA', '#7BCEB8', "white",'#EFF121', '#F5952D', '#E93131', '#D70131', '#D70131'))(256)
# hmcols<-colorRampPalette(c('#2F2C62', '#42399B', '#4A52A7', '#59AFEA', '#7BCEB8', '#A7DA64','#FFFF88', '#F5952D', '#E93131', '#D70131', '#D70131'))(256)
# hmcols<-colorRampPalette(c('#2F2C62','#2F2C62','#2F2C62','#2F2C62','#2F2C62','#2F2C62', '#42399B', '#4A52A7', '#59AFEA', '#7BCEB8', '#A7DA64','#FFFF88', '#F5952D', '#E93131', '#D70131', '#D70131', '#D70131', '#D70131', '#D70131', '#D70131', '#D70131'))(256)
# hmcols<-colorRampPalette(c('#2F2C62','#2F2C62','#2F2C62','#2F2C62','#2F2C62','#2F2C62', '#42399B', '#4A52A7', '#59AFEA', '#7BCEB8','#FFFF88', '#F5952D', '#E93131', '#D70131', '#D70131', '#D70131', '#D70131', '#D70131', '#D70131', '#D70131'))(256)
#hmcols<-colorRampPalette(c('#2F2C62','#2F2C62','#2F2C62','#2F2C62','#2F2C62','#2F2C62','#2F2C62','#2F2C62','#2F2C62','#2F2C62','#2F2C62', '#42399B', '#4A52A7', '#59AFEA', '#7BCEB8','#FFFF88', '#F5952D', '#E93131', '#D70131', '#D70131', '#D70131', '#D70131', '#D70131', '#D70131', '#D70131', '#D70131', '#D70131', '#D70131', '#D70131', '#D70131'))(256)
#hmcols<-colorRampPalette(c('#000000','#000000','#2F2C62','#2F2C62','#2F2C62','#2F2C62','#2F2C62','#2F2C62','#2F2C62','#2F2C62','#2F2C62','#2F2C62', '#42399B', '#4A52A7', '#59AFEA', '#7BCEB8','#FFFF88', '#F5952D', '#E93131', '#D70131', '#D70131', '#D70131', '#D70131', '#D70131', '#D70131', '#D70131', '#D70131', '#D70131', '#D70131', '#D70131', '#D70131', '#D70131'))(256)
hmcols<-colorRampPalette(c('#000000','#000000','#000000','#000000','#000000','#2F2C62','#2F2C62','#2F2C62','#2F2C62','#2F2C62','#2F2C62','#2F2C62','#2F2C62','#2F2C62','#2F2C62', '#42399B', '#4A52A7', '#59AFEA', '#7BCEB8','#FFFF88', '#F5952D', '#E93131', '#D70131', '#D70131', '#D70131', '#D70131', '#D70131', '#D70131', '#D70131', '#D70131', '#D70131', '#D70131', '#D70131', '#D70131', '#D70131', '#D70131','#D70131','#D70131'))(256)
b <- c(seq(from=-8, to=8, by =2 ))

cell <- "IMR90"
dirchosen <- paste0("/Users/eze/Desktop/Dropbox (Schermelleh_Lab)/heatmaps/",cell)
setwd(dirchosen)

allfiles<- list.files()
files <- grep("compile.csv", allfiles)
compfiles <- allfiles[files]

for (compfile in compfiles){ 
condition <-sub("_distn-compile.csv", "",compfile)
datamap<- read.csv(compfile)

# min(datamap$Log2AvNorm)
# max(datamap$Log2AvNorm)


graphtitle <- paste(cell,condition,sep="-")
p <- ggplot(datamap,aes(x = class, y = Marker, fill = Log2AvNorm)) +
  geom_tile()+
  scale_fill_gradientn(limits =c(-6,6), colours = hmcols, breaks=b)+
  scale_x_discrete(expand = c(0, 0))+
  scale_y_discrete(expand = c(0, 0))+
  coord_equal()+
  ggtitle(graphtitle)+
  theme_bw()

savename <- paste0(graphtitle,"_heatmap.png")
ggsave(savename, plot = p)


}


