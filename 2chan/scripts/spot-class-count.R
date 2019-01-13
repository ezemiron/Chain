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

#This script needs to have a prior setwd() to that of the Makefile.
# that is so that ./results is posible
dirchosen  <- readline(prompt = "Enter path to directory of Makefile so that ./results is possible. \n  ie: \"~/Documents/papers/chrom_marks\".\nIf al
ready in Makefile directory hit enter\n")

if (! dirchosen == ""){
setwd(dirchosen)
} else {
dirchosen <- getwd()
}

#dirchosen <- "/home/eze/Documents/papers/chrom_marks/results"

cell="C127"
condition="numtest"

dir = paste0(cell,"/", condition, "/")
dir_names <- list.dirs(dir, recursive = FALSE)
#dir_names <- dir_names[2:length(dir_names)]

var_names <- list.dirs(dir, recursive = FALSE, full.names=FALSE)
#var_names <- var_names[2:length(var_names)]

for (dir_name in dir_names){

    newdir <- paste(dirchosen, dir_name, sep ="/")
    setwd(newdir)
    dir_num <- which(dir_names == dir_name)
    var_name <- var_names[dir_num]

    df <- data.frame(NA)
    
    names(df) <- NA
    df <- cbind(df,df,df,df)
    
    allfilesandpaths <- list.files(pattern= var_name, full.names=TRUE);
    if (length(allfilesandpaths) > 0){
    
      filesandpaths  <- grep("*distn.csv",allfilesandpaths, value = TRUE)
      if (length(filesandpaths) > 0){
    
        for (fileandpath in filesandpaths){
            if (length(fileandpath) > 0){
              csv <- read.csv(fileandpath, colClasses=c("NULL",NA,NA,"NULL",NA,"NULL", "NULL","NULL"))
              cell_num <- which(filesandpaths == fileandpath)
              csv$cell <- cell_num
              names(csv) <- rep(NA, ncol(csv))
              df <- rbind(df,csv)
            }
        }
        df <- df[2:nrow(df),]
        names(df) <- c("class", "spotnum","classnum","cell")
#    df <- transform(df, AvNorm = rowMeans(df[,-1], na.rm = TRUE))
#    df <- transform(df, Log2AvNorm = log2(df$AvNorm))

#aggregate for spot-class-count:
        savename  <- paste0(var_name,"_spot-class-count.csv")
        write.csv(df,savename);

#aggregate for spot-sum-summary:
	df3 <- data.frame(var_name, aggregate(df[,2], list(cell =df$cell), sum))

#	df3 <-aggregate(df[,2], list(cell =df$cell), sum)
	names(df3) <- c("marker","cell","spotnum")

#	df3 <- data.frame(var_name, sum(df$spotnum), max(df$cell))
#	names(df3) <- c("marker", "spotnum", "cellnum")

        savename3  <- paste0(var_name,"_spot-num-summary.csv")
        write.csv(df3,savename3);

#aggregate for spot-class-summary:
        df2 <-aggregate(df[,2:3], list(class =df$class), mean)
        df2 <-transform(df2, numofcells = rep(max(df$cell),dim(df2)[1]))
      
        savename2  <- paste0(var_name,"_spot-class-summary.csv")
        write.csv(df2,savename2);






        }
    }
}

