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

celltype = "Fena"

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
    dfA <- cbind(df,df,df,df,df,df,df,df,df,df,df)
    dirchosen2 <- getwd()

    for (dir_name in dir_names){

      newdir2 <- paste(dirchosen2, dir_name, sep ="/")
      setwd(newdir2)
      dir_num <- which(dir_names == dir_name)
      var_name <- var_names[dir_num]


      allfilesandpaths <- list.files(pattern= var_name, full.names=TRUE);
      
      if (length(allfilesandpaths) > 0){
        distnsum  <- grep("*_subdistn-chrom-summary.csv",allfilesandpaths, value = TRUE)

        if (length(distnsum) > 0){
          spotsum <- sub ("subdistn-chrom-summary", "spot-class-sub-summary", distnsum)
          
          if (length (which (allfilesandpaths == spotsum))> 0){
            condition <- rep(con_name, 7)
            condition <- as.data.frame(condition)

            dfB <- read.csv(distnsum)
            dfB <- dfB[, c("Marker","class","StDv","negCI95","posCI95")]

            dfC <- read.csv(spotsum)
            dfC <- dfC[, c("spotnum","spotnumn","classnum","classnumn","numofcells")]

            dfB <- cbind(condition, dfB, dfC)
            names(dfB) <- rep(NA, ncol(dfB))
            dfA <- rbind(dfA, dfB)
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
        names(dfA) <- c("condition","Marker", "class", "StDv","negCI95","posCI95","spotnum","spotnumn","classnum","classnumn","numofcells")

        savename  <- paste0(con_name,"_distn-compile-subchrom_volnorm.csv")
        write.csv(dfA,savename);
  }
