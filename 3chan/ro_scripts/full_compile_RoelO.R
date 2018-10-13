#!/usr/bin/env Rscript

# This program is free software: you can redistribute it and/or modify
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


args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
        celltype="C127"
        } else {
        celltype = args[1]
}

setwd("/home/eze/Documents/papers/chrom_marks/results/EP20171221-2_Gluteraldehyde2/DAPI/Csvs")
FILES <- list.files(pattern="distn-full.csv", recursive=TRUE)

for (i in 1:length(FILES)) {
	file <- FILES[i]
	print(file)
	if(i ==1) {
		compile <- read.csv(file)
		Marker <-  strsplit(file, "-") [[1]] [1]
                Cell <- "C127"
                Condition <- "DAPI"
		compile <- cbind(
			File      =rep(file, nrow(compile)),
			cell      =rep(Cell, nrow(compile)),
			condition =rep(Condition, nrow(compile)),
			marker    =rep(Marker, nrow(compile)),
			compile) 
	} else {
	  csv <- read.csv(file)
                Marker <- strsplit(file, "-") [[1]] [1]
                Cell <- "C127" 
                Condition <- "DAPI"
                csv <- cbind(
			File = rep(file, nrow(csv)),
			cell=rep(Cell, nrow(csv)),
                        condition=rep(Condition, nrow(csv)),
                        marker=rep(Marker, nrow(csv)),
                        csv) 
		compile <- rbind(compile, csv)

	}

}
compile <- compile[, -grep("^X", colnames(compile))]     

write.csv(compile, "/home/eze/Documents/papers/chrom_marks/results/EP20171221-2_Gluteraldehyde2/DAPI/Csvs/full_compiled.csv")
