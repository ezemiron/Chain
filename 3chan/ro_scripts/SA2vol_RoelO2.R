#!/usr/bin/env Rscript
##
## Copyright (C) 2016 Ezequiel Miron <eze.miron@bioch.ox.ac.uk>
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.

library(tiff)

dirchosen  <- readline(prompt = "This script should be run from the results directory as the working directory. \n  ie: \"~/Documents/papers/chrom_marks/results\".\nIf already in results directory hit enter\n")

dirchosen <- getwd()

celltype = "C127"
condition = "Trip"
topPATH = paste(celltype, condition, sep = "/")

dir_names <- list.dirs(topPATH, recursive = TRUE)
dirchosen2 <- paste(dirchosen, dir_names[1], sep = "/")
dir_names <- dir_names[2:length(dir_names)]

#var_names <- list.dirs("C127/sub2", recursive = FALSE, full.names=FALSE)


df <- as.data.frame(NA)
df <- cbind(df,df,df)
names(df)<- NA

for (dir_name in dir_names){

  newdir <- paste(dirchosen, dir_name, sep ="/")
  setwd(newdir)

  #dir_num <- which(dir_names == dir_name)
  #var_name <- var_names[dir_num]


  filesandpaths <- list.files(pattern= "*_bound.tif", full.names=TRUE, recursive = FALSE);
  if (length(filesandpaths) > 0){

    chromclass <- list.files(pattern= "*spot-class-count.csv", full.names=TRUE, recursive = FALSE);

    if ((length(chromclass ==1)&(file.info(chromclass)$size> 0))){
      chromclass <- read.csv(chromclass[1])

      if (max(chromclass$cell)==length(filesandpaths)){

        for (fileandpath in filesandpaths){
          if (file.info(fileandpath)$size > 0){

          ##find volume of classes 2 to 7:
          cellnum <- which(fileandpath == filesandpaths)
          cellLogic = chromclass$cell == cellnum
          chromcell <- chromclass [cellLogic,]
          chromvol <- sum(chromcell$classnum[2:7])
          nucvol <- sum(chromcell$classnum[1:7])
	  
          ##find surface area of corresponding cell:
          bound <- readTIFF(fileandpath, all= TRUE,convert=TRUE)
          bound <- as.data.frame(bound)
          boundsum <- sum(bound)

          dfA <- as.data.frame(cbind(nucvol, chromvol, boundsum))
          titles <- names(dfA)
          names(dfA) <- NA

          df <- rbind(df, dfA)
          } #else {boundary file: fileandpath is empty}
        }

      } #else {number of cells in var_name spot-class-count does not match numner of boundary files}

    } #else {no var_name spot-class-count found, or empty file}

  } #else { no var_name boundaries found}

}

df <- df[2:nrow(df),]
names(df) <- titles
df <- transform(df,  SA2vol= ((df$boundsum)/(df$chromvol)))
df <- transform(df,  chrom2vol= ((df$chromvol)/(df$nucvol)))

setwd(dirchosen2)
savename  <- paste0(celltype,"_",condition,"_SA2vol_raw.csv")
write.csv(df,savename);



Statsum <- function(x){
  meanX <- mean(x)
  stDev <- sd(x)
  lCI <- mean(x)+(-1.96)*sd(x)/sqrt(length(x))
  uCI <- mean(x)+(1.96)*sd(x)/sqrt(length(x))
  list(meanX,stDev,lCI,uCI)
}

dfStat <- as.data.frame(Statsum(df$SA2vol))
names(dfStat) <- c("mean-SA2vol","sdDev","lowCI95","uprCI95")
dfStat <- cbind(celltype, condition, dfStat)
savename2  <- paste0(celltype,"_",condition,"_SA2vol_stat.csv")
write.csv(dfStat,savename2);

dfStat <- as.data.frame(Statsum(df$chrom2vol))
names(dfStat) <- c("mean-chrom2vol","sdDev","lowCI95","uprCI95")
dfStat <- cbind(celltype, condition, dfStat)
savename2  <- paste0(celltype,"_",condition,"_chrom2vol_stat.csv")
write.csv(dfStat,savename2);
