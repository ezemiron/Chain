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
    dfA <- cbind(df,df,df)
    dirchosen2 <- getwd()

    for (dir_name in dir_names){

      newdir2 <- paste(dirchosen2, dir_name, sep ="/")
      setwd(newdir2)
      dir_num <- which(dir_names == dir_name)
      var_name <- var_names[dir_num]


      allfilesandpaths <- list.files(pattern= var_name, full.names=TRUE);
      
      if (length(allfilesandpaths) > 0){
        distnsum  <- grep("*_distn-summary.csv",allfilesandpaths, value = TRUE)

        if (length(distnsum) > 0){
          spotsum <- sub ("distn", "spot-num", distnsum)
          
          if (length (which (allfilesandpaths == spotsum))> 0){
           # condition <- con_name
           # condition <- as.data.frame(condition)

#            dfB <- read.csv(distnsum)
#            dfB <- dfB[, c("Marker","class","Log2AvNorm","lowLogErr","uprLogErr")]

            dfC <- read.csv(spotsum)
	   # dfC <- read.csv(spotsum, colClasses=c("NULL","NULL",NA,NA))
            dfC <- dfC[, c("marker","cell","spotnum")]

	    #dfC <- cbind(con_name,dfC)

            names(dfC) <- rep(NA, ncol(dfC))
            dfA <- rbind(dfA, dfC)
          }
          else{
          warn <- paste("\n",con_name, spotsum , "not found.","\n", "Will not add",distnsum, "to compile list.\n Run spot-class-count.R with", con_name,"/",var_name,".")
          warning(warn)
          }
        }
      }
    }
  setwd(dirchosen2)


#rename
        dfA <- dfA[2:nrow(dfA),]
	names(dfA) <- c("Marker", "cell","spotnum")

	dfZ <- dfA
	
	dfB <- data.frame(
	aggregate(dfA$spotnum, list(Marker =dfA$Marker), sum),
	aggregate(dfA$spotnum, list(Marker =dfA$Marker), mean),
	aggregate(dfZ$spotnum, list(Marker =dfZ$Marker), sd),
	aggregate(dfZ$cell, list(Marker =dfZ$Marker), max)
	)


	dfB <- cbind(con_name,dfB)
	dfB <- dfB[,c(1:3,5,7,9)]
        names(dfB) <- c("condition","Marker","TSpotnum","AvSpotnum","sdSpotnum" ,"cellnum")

        savename  <- paste0(con_name,"_spot-num_compile.csv")
        write.csv(dfB,savename);
  }



