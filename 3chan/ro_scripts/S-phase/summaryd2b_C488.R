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



dirchosen  <- readline(prompt = "This script should be run from the results directory as the working directory. \n  ie: \"~/Documents/papers/chrom_marks/results\".\nIf already in results directory hit enter\n")

dirchosen <- getwd()

##celltype =  "Cardiomyocytes"
##condition = "37C"
celltype = "HCT116-SCC1mAID"
condition = "16hdox-6haux"
topPATH = paste(celltype, condition, sep = "/")

dir_names <- list.dirs(topPATH, recursive = TRUE)
#dir_names <- list.dirs("C127/sub2")
dir_names <- dir_names[2:length(dir_names)]

var_names <- list.dirs(topPATH, recursive = TRUE, full.names=FALSE)
#var_names <- list.dirs("C127/sub2", full.names=FALSE)
var_names <- var_names[2:length(var_names)]

for (dir_name in dir_names){

  newdir <- paste(dirchosen, dir_name, sep ="/")
  setwd(newdir)
  dir_num <- which(dir_names == dir_name)
  var_name <- var_names[dir_num]

  stddist <- seq(-400,400)
  df <- as.data.frame(stddist)
  
  df2a <- as.data.frame(NA)
  df2b <- df2a
  df3a <- df2a
  df3b <- df2a
  df4 <- df
#  names(df2) <- c(NA,NA)
  allfilesandpaths <- list.files()
  
##  allfilesandpaths <- list.files(pattern= var_name, full.names=TRUE);

  if (length(allfilesandpaths) > 0){
    
###########################
#for normalised d2b only:
    filesandpaths  <- grep("*C488-d2bnorm.csv",allfilesandpaths, value = TRUE)
    if (length(filesandpaths) > 0){

    for (fileandpath in filesandpaths){
        if (length(fileandpath) > 0){
          
          d2b <- read.csv(fileandpath, header=TRUE,colClasses=c("NULL","NULL","NULL","NULL",NA,"NULL"))
#          cell_num <- which(filesandpaths == fileandpath)
#          d2b$cell <- cell_num
#          names(d2b) <- rep(NA, ncol(d2b))
          df <- cbind(df,d2b)
        }
    }

    df <- transform(df, AvNorm = rowMeans(df[,2:length(df)], na.rm = TRUE))
    #log base 2 of the average:
    df <- transform(df, Log2AvNorm = log2(df$AvNorm))
    #standard deviation of the average:
    df <- transform(df, StDv = apply (df[,2:(length(df)-2)], 1, sd, na.rm =TRUE))
    #95% confidence interval:
    df <- transform(df, CI95 = apply(df[,2:(length(df)-3)], 1, function(x){1.96*sd(x)/sqrt(length(x))}))
    #lower and upper 95% confidence intervals:
    df <- transform(df, negCI95 = apply(df[,2:(length(df)-4)], 1, function(x){mean(x)+(-1.96)*sd(x)/sqrt(length(x))}))
    df <- transform(df, posCI95 = apply(df[,2:(length(df)-5)], 1, function(x){mean(x)+c(1.96)*sd(x)/sqrt(length(x))}))
    
#lower and upper log errors for the log average from the 95CIs:
    df <- transform(df, lowLogErr = df$Log2AvNorm - sqrt((df$Log2AvNorm-log2(df$negCI95))^2))
    df <- transform(df, uprLogErr = df$Log2AvNorm + sqrt((df$Log2AvNorm-log2(df$posCI95))^2))
    
    ##to plot averaged then logged results with shaded error bands:
    #library(ggplot2)
    #p<-ggplot(data=df, aes(x=df$stddist, y=df$Log2AvNorm)) + geom_line()
    #p<-p+geom_ribbon(aes(ymin=df$lowLogErr, ymax=df$uprLogErr), linetype=2, alpha=0.1)
    #p

#p<-ggplot(data=data, aes(x=interval, y=OR, colour=Drug)) + geom_point() + geom_line()
#p<-p+geom_ribbon(aes(ymin=data$lower, ymax=data$upper), linetype=2, alpha=0.1)

    savename  <- paste0(var_name,"_C488_d2bnorm-summary.csv")
    write.csv(df,savename);
    }

###########################
#for absolute marker and random d2b only:

filesandpaths2  <- grep("_C488-d2bfit.csv",allfilesandpaths, value = TRUE)
    if (length(filesandpaths2) > 0){
    numcells <- 0

    for (fileandpath2 in filesandpaths2){
        if (length(fileandpath2) > 0){
          
          fitdf <- read.csv(fileandpath2, header=TRUE)
#          cell_num <- which(filesandpaths == fileandpath)
#          d2b$cell <- cell_num

          df2a <- rbind(df2a, fitdf$mfitted[1])
          df2b <- rbind(df2b, fitdf$mfitted[2])
          
          df3a <- rbind(df3a, fitdf$rfitted[1])
          df3b <- rbind(df3b, fitdf$rfitted[2])
          numcells <- numcells + 1
        }
    }

    Markermu <- mean(df2a[2:nrow(df2a),1])
    Randommu <- mean(df3a[2:nrow(df3a),1])
#to get the average standard deviation you have to: square them to get the variance for each, then get their average (ie, sum them and divide by the number of inputs) and finally square root them. This is NOT the same as just taking the mean of all the standard deviations!
    Markersigma <- sqrt(sum(df2b[2:nrow(df2b),1]^2)/(nrow(df2b)-1))
    Randomsigma <- sqrt(sum(df3b[2:nrow(df3b),1]^2)/(nrow(df3b)-1))

    Mabsdf <- cbind(Markermu,Markersigma,numcells)
    rownames(Mabsdf) <- var_name

    Rabsdf <- cbind(Randommu,Randomsigma,numcells)
    rownames(Rabsdf) <- var_name
    
    savename2  <- paste0(var_name,"_C488_abs-md2b_summary.csv")
    write.csv(Mabsdf,savename2);
    savename3  <- paste0(var_name,"_C488_abs-rd2b_summary.csv")
    write.csv(Rabsdf,savename3);
    }



###########################
#for POSITVE fitted random d2b only, metric for chromatin network width:

    filesandpaths  <- grep("*_C488-d2bnorm.csv",allfilesandpaths, value = TRUE)
    if (length(filesandpaths) > 0){

    for (fileandpath in filesandpaths){
        if (length(fileandpath) > 0){
          
          d2b <- read.csv(fileandpath, header=TRUE,colClasses=c("NULL","NULL","NULL",NA,"NULL","NULL"))
#          cell_num <- which(filesandpaths == fileandpath)
#          d2b$cell <- cell_num
#          names(d2b) <- rep(NA, ncol(d2b))
          df4 <- cbind(df4,d2b)
        }
    }

    #remove negative distances:
    df4ind <- df4$stddist > -1
    df4 <- df4[df4ind,]


    df4 <- transform(df4, AvNetwork = rowMeans(df4[,2:length(df4)], na.rm = TRUE))

    df4 <- transform(df4, StDv = apply (df4[,2:(length(df4)-1)], 1, sd, na.rm =TRUE))
    #95% confidence interval:
    df4 <- transform(df4, CI95 = apply(df4[,2:(length(df4)-2)], 1, function(x){1.96*sd(x)/sqrt(length(x))}))
    #lower and upper 95% confidence intervals:
    df4 <- transform(df4, negCI95 = apply(df4[,2:(length(df4)-3)], 1, function(x){mean(x)+(-1.96)*sd(x)/sqrt(length(x))}))
    df4 <- transform(df4, posCI95 = apply(df4[,2:(length(df4)-4)], 1, function(x){mean(x)+c(1.96)*sd(x)/sqrt(length(x))}))
    

    ##to plot averaged then logged results with shaded error bands:
    #library(ggplot2)
    #p<-ggplot(data=df, aes(x=df$stddist, y=df$Log2AvNorm)) + geom_line()
    #p<-p+geom_ribbon(aes(ymin=df$lowLogErr, ymax=df$uprLogErr), linetype=2, alpha=0.1)
    #p

#p<-ggplot(data=data, aes(x=interval, y=OR, colour=Drug)) + geom_point() + geom_line()
#p<-p+geom_ribbon(aes(ymin=data$lower, ymax=data$upper), linetype=2, alpha=0.1)

    savename  <- paste0(var_name,"_C488_network-summary.csv")
    write.csv(df4,savename);
    }







  }
}

