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

  centroids.fpath <- argv[1]
  segimg.fpath <- argv[2]
  mask <- argv [3]
   
  if (length (argv) < 3)
    stop ('usage error: main_example requires at least 3 arguments')

  MCNRmask <- readTIF (mask)
  MCNRmask <- MCNRmask*255

if(length(MCNRmask)>0){
    img <- readTIF (segimg.fpath)
    csv1 <- read.csv(centroids.fpath, header=TRUE)
    csv1 = csv1[2:7]

    MCNRlogic = (1==MCNRmask)   

    vox  <- c(41, 41, 125)
    x.microns <- (dim(MCNRmask)[2]*vox[2])/1000
    y.microns <- (dim(MCNRmask)[1]*vox[1])/1000
    z.microns <- (dim(MCNRmask)[3]*vox[3])/1000
    img <- round(img*7)
    
    csv1r <- round (data.frame (csv1[1]/vox[2], csv1[2]/vox[1], csv1[3]/vox[3]))

  ## for 3D linear index: ind <- (m.n)(z-1)+m(c-1)+r
  ## remember that
  ## m = length of columns, determined by the 1st dim of img
  ## n = length of rows, determined by the 2nd dim of img
  ## c = column number, determined by the x coordinate
  ## r = row number, determined by the y coordinate
  ## z = optical plane number, determined by z coordinate

  ind1 <- ((dim(MCNRmask)[1] * dim(MCNRmask)[2] * (csv1r$z - 1)) + ((csv1r$x -1) * dim(MCNRmask)[1]) + csv1r$y)

#linear index on MCNRlogical mask gets logical vector
  MCNRind1 <- MCNRlogic[ind1]
 
#use logical vector to remove centroids outside MCNR
  points1<-csv1[MCNRind1,]
    
  write.csv(points1, file ="")
##  outpath <- gsub("inv.csv","invvol.csv",centroidsinv.fpath)
##  write.csv(points2, file = outpath) 

  }
  else{
    print ("MCNR mask and segmented image are of different dimensions.")
  }
  options(warn = oldw)
  return (0)
}

if (identical (environment (), globalenv ()))
  quit (status=main (commandArgs (trailingOnly = TRUE)))
