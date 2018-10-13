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


main <-function (argv)
{
  oldw <- getOption("warn") # to suppress warnings.
  options(warn = -1)
  library (EBImage, warn.conflicts = FALSE, verbose = FALSE)
  library (bioimagetools, warn.conflicts = FALSE, verbose = FALSE)  

  MCNRmask_fpath <- argv[1]
  csv_fpath <- argv[2]
  ## TODO figure out a way to read (and write) the resolution from the segmented image
  if (length (argv) < 2)
    stop ('usage error: main_example requires at least 2 arguments')

  xlength <- as.numeric (41) #set to 41nm in Oxford OMX V3 (should be 39.5 in Munich)
  ylength <- as.numeric (41) #set to 41nm in Oxford OMX V3 (should be 39.5 in Munich)
  zlength <- as.numeric (125) #set to 125nm in Oxford OMX V3 (and Munich)


  ## read segmented image
  MCNRmask <- readTIF (MCNRmask_fpath) #This opens the seg 16bit tiff with all channels correclty

  ## read the csvs:
  csv  <- read.csv (csv_fpath, header=TRUE)

  mask <- MCNRmask
  mask = mask==1

  csvr <- round (data.frame (csv[1]/xlength, csv[2]/ylength, csv[3]/zlength))

  ind1 <- ((dim(mask)[1] * dim(mask)[2] * (csvr$z - 1)) + ((csvr$x -1) * dim(mask)[1]) + csvr$y)

  imind1 <- mask[ind1]

  csv2<-csv[imind1,]


  ##write(compile1, file=stdout())
  #write.csv(csv2, file ="")
  write.csv(csv2, "outputcsv.csv")
  options(warn = oldw)
  return (0)
}

if (identical (environment (), globalenv ()))
  quit (status=main (commandArgs (trailingOnly = TRUE)))

