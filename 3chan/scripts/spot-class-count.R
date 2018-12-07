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



dirchosen <- readline(prompt = "This script should be run two directories above the marks")

dirchosen <- getwd()



Cell = "Fena"
Condition = "G2"

topPATH = paste(Cell, Condition, sep = "/")

dir_names <- list.dirs(topPATH, recursive = TRUE)
#dir_names <- dir_names[3:length(dir_names)]

var_names <- list.dirs(topPATH, recursive = TRUE, full.names=FALSE)
#var_names <- var_names[3:length(var_names)]

for (dir_name in dir_names){

    newdir <- paste(dirchosen, dir_name, sep ="/")
    setwd(newdir)
    dir_num <- which(dir_names == dir_name)
    var_name <- var_names[dir_num]

    df <- data.frame(NA)
    
    names(df) <- NA
    df <- cbind(df,df,df,df,df,df)
    
    allfilesandpaths <- list.files(pattern= var_name, full.names=TRUE);
    if (length(allfilesandpaths) > 0){
    
      filesandpaths  <- grep("*-subdistn.csv",allfilesandpaths, value = TRUE)
      if (length(filesandpaths) > 0){
    
        for (fileandpath in filesandpaths){
            if (length(fileandpath) > 0){
              csv <- read.csv(fileandpath, colClasses=c("NULL",NA,NA,NA,NA,NA,"NULL","NULL"))
              cell_num <- which(filesandpaths == fileandpath)
              csv$cell <- cell_num
              names(csv) <- rep(NA, ncol(csv))
              df <- rbind(df,csv)
            }
        }
        df <- df[2:nrow(df),]
        names(df) <- c("class","spotnum","spotnumn","classnum","classnumn","cell")
#    df <- transform(df, AvNorm = rowMeans(df[,-1], na.rm = TRUE))
#    df <- transform(df, Log2AvNorm = log2(df$AvNorm))

        savename  <- paste0(var_name,"_spot-class-sub-count.csv")
        write.csv(df,savename);

        df2 <-aggregate(df[,2:3:4:5], list(class =df$class), mean)
        df2 <-transform(df2, numofcells = rep(max(df$cell),dim(df2)[1]))
      
        savename2  <- paste0(var_name,"_spot-class-sub-summary.csv")
        write.csv(df2,savename2);
        }
    }
}

