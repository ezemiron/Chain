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


#This script needs to have a prior setwd() to that of the RESULTS.
# that is so that ./"celltype" is posible

dirchosen  <- readline(prompt = "This script should be run from the results directory as the working directory. \n  ie: \"~/Documents/papers/chrom_marks/results\".\nIf already in results directory hit enter\n")

dirchosen <- getwd()

celltype = "C127"


topdir_names <- list.dirs(celltype, recursive = FALSE)
#topdir_names <- topdir_names[2:length(topdir_names)]

con_names <- list.dirs(celltype, recursive = FALSE,full.names=FALSE)
#con_names <- con_names[2:length(con_names)]


for (topdir_name in topdir_names){

  newdir <- paste(dirchosen, topdir_name, sep ="/")
  setwd(newdir)
  topdir_num <- which(topdir_names == topdir_name)
  con_name <- con_names[topdir_num]

  dir_names <- list.dirs(recursive = FALSE)
  var_names <- list.dirs(recursive = FALSE,full.names=FALSE)

  df <- data.frame(NA)
  names(df) <- NA
  dfA1 <- cbind(df,df,df,df)

  dirchosen2 <- getwd()


#######
#Process all d2b markers first:

  for (dir_name in dir_names){

    newdir2 <- paste(dirchosen2, dir_name, sep ="/")
    setwd(newdir2)
    dir_num <- which(dir_names == dir_name)
    var_name <- var_names[dir_num]

    allfilesandpaths <- list.files(pattern= var_name, full.names=TRUE);

    if (length(allfilesandpaths) > 0){
      d2bFsum  <- grep("*abs-md2b_summary.csv",allfilesandpaths, value = TRUE)

      if (length(d2bFsum) > 0){

        dfB1 <- read.csv(d2bFsum)

#        Marker <- rep(var_name,nrow(dfB1))
#        Marker <- as.data.frame(Marker)
#        dfC1 <- cbind(Marker, dfB1)
        names(dfB1) <- rep(NA, ncol(dfB1))
        dfA1 <- rbind(dfA1, dfB1)
      }
      else{
      warn <- paste("\n",con_name, var_name, "not found.","\n", "Will not add",var_name, "(md2b-abs) to", con_name ,"compile list.\n Run summaryd2b.R with", con_name,"/",var_name,".")
      warning(warn)
      }
    }
  }
condition <- rep(con_name, nrow(dfA1))
condition <- as.data.frame(condition)
dfA1 <- cbind(condition, dfA1)

#rename
  dfA1 <- dfA1[2:nrow(dfA1),]
  names(dfA1) <- c("Condition","Marker", "mu-d2b-nm", "sigma-d2b-nm","numcells")
  dfA1 <- transform(dfA1, CI95 = (1.96*dfA1[4])/sqrt(dfA1$numcells))
  names(dfA1) <- c("Condition","Marker", "mu-d2b-nm", "sigma-d2b-nm","numcells","CI95")


setwd(dirchosen2)
savename  <- paste0(con_name,"_md2b-compile.csv")
write.csv(dfA1,savename);




#######
##Repeat process d2b norm only:

setwd(dirchosen)

  dfA1 <- cbind(df,df,df,df,df)
#  
  for (dir_name in dir_names){

    newdir2 <- paste(dirchosen2, dir_name, sep ="/")
    setwd(newdir2)
    dir_num <- which(dir_names == dir_name)
    var_name <- var_names[dir_num]

    allfilesandpaths <- list.files(pattern= var_name, full.names=TRUE);

    if (length(allfilesandpaths) > 0){
      d2bFsum  <- grep("*d2bnorm-summary.csv",allfilesandpaths, value = TRUE)

      if (length(d2bFsum) > 0){

        dfB1 <- read.csv(d2bFsum)
        dfB1 <- dfB1[c("stddist", "Log2AvNorm", "lowLogErr", "uprLogErr")]
        Marker <- rep(var_name,nrow(dfB1))
        Marker <- as.data.frame(Marker)
        dfB1 <- cbind(Marker, dfB1)
        names(dfB1) <- rep(NA, ncol(dfB1))
        dfA1 <- rbind(dfA1, dfB1)
      }
      else{
      warn <- paste("\n",con_name, var_name, "not found.","\n", "Will not add",var_name, "(d2bnorm) to", con_name ,"compile list.\n Run summaryd2b.R with", con_name,"/",var_name,".")
      warning(warn)
      }
    }
  }
condition <- rep(con_name, nrow(dfA1))
condition <- as.data.frame(condition)
dfA1 <- cbind(condition, dfA1)

#rename
  dfA1 <- dfA1[2:nrow(dfA1),]
  names(dfA1) <- c("Condition","Marker", "stddist", "Log2AvNorm","lowLogErr", "uprLogErr")


setwd(dirchosen2)
savename  <- paste0(con_name,"_d2bnorm-compile.csv")
write.csv(dfA1,savename);




#######
##Repeat process d2b norm only:

setwd(dirchosen)

  dfA1 <- cbind(df,df,df,df,df)
#  
  for (dir_name in dir_names){

    newdir2 <- paste(dirchosen2, dir_name, sep ="/")
    setwd(newdir2)
    dir_num <- which(dir_names == dir_name)
    var_name <- var_names[dir_num]

    allfilesandpaths <- list.files(pattern= var_name, full.names=TRUE);

    if (length(allfilesandpaths) > 0){
      d2bFsum  <- grep("*_network-summary.csv",allfilesandpaths, value = TRUE)

      if (length(d2bFsum) > 0){

        dfB1 <- read.csv(d2bFsum)
        dfB1 <- dfB1[c("stddist", "AvNetwork", "negCI95", "posCI95")]
        Marker <- rep(var_name,nrow(dfB1))
        Marker <- as.data.frame(Marker)
        dfB1 <- cbind(Marker, dfB1)
        names(dfB1) <- rep(NA, ncol(dfB1))
        dfA1 <- rbind(dfA1, dfB1)
      }
      else{
      warn <- paste("\n",con_name, var_name, "not found.","\n", "Will not add",var_name, "(network) to", con_name ,"compile list.\n Run summaryd2b.R with", con_name,"/",var_name,".")
      warning(warn)
      }
    }
  }
condition <- rep(con_name, nrow(dfA1))
condition <- as.data.frame(condition)
dfA1 <- cbind(condition, dfA1)

#rename
  dfA1 <- dfA1[2:nrow(dfA1),]
  names(dfA1) <- c("Condition","Marker", "stddist", "AvNetwork", "negCI95", "posCI95")


setwd(dirchosen2)
savename  <- paste0(con_name,"_network-compile.csv")
write.csv(dfA1,savename);


}
