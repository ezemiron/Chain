#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
##
## Copyright (C) 2016 Ezequiel Miron <eze.miron@bioch.ox.ac.uk>
## modified by Phoebe Oldach 2017
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



dirchosen  <- readline(prompt = "This script should be run from the results directory as the working directory. \n  ie: \"~/Documents/papers/chrom_marks/results\".\nIf already in results directory hit enter\n")

dirchosen <- getwd()

if (length(args)==0) {
	celltype="HCT116Scc1mAID"
	condition="Aux2h"
	} else {
	celltype = args[1]
	condition = args[2]
}

topPATH = paste(celltype, condition, sep = "/")

dir_names <- list.dirs(topPATH, recursive = TRUE)
dir_names <- dir_names[2:length(dir_names)]

var_names <- list.dirs(topPATH, recursive = TRUE, full.names=FALSE)
var_names <- var_names[2:length(var_names)]

for (dir_name in dir_names){

    newdir <- paste(dirchosen, dir_name, sep ="/")
    setwd(newdir)
    dir_num <- which(dir_names == dir_name)
    var_name <- var_names[dir_num]
    Marker <- rep(var_name,7)
    df  <- as.data.frame(Marker)
    df  <- transform(df, class = c(1:7))

    allfilesandpaths <- list.files(full.names=TRUE);
    if (length(allfilesandpaths) > 0){
    
      filesandpaths  <- grep("*distn.csv",allfilesandpaths, value = TRUE)
      
    if (length(filesandpaths) > 0){
    
      for (fileandpath in filesandpaths){

          csv <- read.csv(fileandpath, colClasses=c("NULL","NULL",NA,"NULL","NULL","NULL", NA,"NULL")) 
	  df <- cbind(df,csv)
	 }

	if (length(filesandpaths) > 2) {
      df <- transform(df, AvNorm = rowMeans(df[,grep("norm1", colnames(df))], na.rm = TRUE))
      #log base 2 of the average:
      df <- transform(df, Log2AvNorm = log2(df$AvNorm))
    #standard deviation of the average:
      df <- transform(df, StDv = apply (df[,grep("norm1", colnames(df))], 1, sd, na.rm =TRUE))
    #95% confidence interval:
      df <- transform(df, CI95 = apply(df[, grep("norm1", colnames(df))], 1, function(x){1.96*sd(x)/sqrt(length(x))}))
    #lower and upper 95% confidence intervals:
      df <- transform(df, negCI95 = apply(df[,grep("norm", colnames(df))], 1, function(x){mean(x)+(-1.96)*sd(x)/sqrt(length(x))}))
      df <- transform(df, posCI95 = apply(df[,grep("norm", colnames(df))], 1, function(x){mean(x)+c(1.96)*sd(x)/sqrt(length(x))}))
    
#lower and upper log errors for the log average from the 95CIs:
      df <- transform(df, lowLogErr = sqrt((df$Log2AvNorm-log2(df$negCI95))^2))
      df <- transform(df, uprLogErr = sqrt((df$Log2AvNorm-log2(df$posCI95))^2))

  df <- transform(df, AvCount = rowMeans(df[,grep("dist", colnames(df))], na.rm = TRUE))
  df <- transform(df, StdDvCount = apply(df[,grep("dist", colnames(df))], 1, sd, na.rm = TRUE))
    } else { df <- cbind(df, AvNorm=df$norm1) }
    
      savename  <- paste0("/home/eze/Documents/papers/chrom_marks/results/csvs/", celltype, condition, Marker[1],"_distn-summary.csv")
      write.csv(df, file =savename);
    }
  }

}
