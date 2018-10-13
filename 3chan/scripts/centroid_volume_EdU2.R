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
  MCNRmask.fpath <- argv[2]
  segimg.fpath <- argv[3]
  dif.mask <- argv[4]
  edu.mask <- argv [5]
  edu.fpath = gsub("mask-centroids.csv","EdU-centroidsMCN.csv",centroids.fpath)
  dif.fpath = gsub("mask-centroids.csv","dif-centroidsMCN.csv",centroids.fpath)
  mask.out = gsub("mask-centroids.csv","masktrial.csv",centroids.fpath)
  
  if (length (argv) < 5)
    stop ('usage error: main_example requires at least 3 arguments')


  MCNRmask <- readTIF (MCNRmask.fpath)
  MCNRmask <- MCNRmask[,,2,]

if(length(MCNRmask)>0){
    img <- readTIF (segimg.fpath)
    csv1 <- read.csv(centroids.fpath, header=TRUE)

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

  csv2r<-csv1r[MCNRind1,]

#second trimmed linear index
  ind2 <- ((dim(MCNRmask)[1] * dim(MCNRmask)[2] * (csv2r$z - 1)) + ((csv2r$x -1) * dim(MCNRmask)[1]) + csv2r$y)
 

 #use logical vector to remove centroids outside MCNR
  points1<-csv1[MCNRind1,]
 
  #indexing the image using the second trimmed linear index
  classnum <- img[ind2]
  compile <- data.frame(points1,classnum)
 
  write.csv(compile, file ="")
##  outpath <- gsub("inv.csv","invvol.csv",centroidsinv.fpath)
##  write.csv(points2, file = outpath) 

    csv1 <- compile

 EdUmask <- readTIF (edu.mask)
 EdUmask <- EdUmask*255
 EdUmasklogic = (1==EdUmask)
 difmask <- readTIF (dif.mask)
 difmask <- difmask*255
 difmasklogic = (1==EdUmask)
   
 csv1r <- round (data.frame (csv1[1]/vox[2], csv1[2]/vox[1], csv1[3]/vox[3]))

 ind4 <- ((dim(EdUmask)[1] * dim(EdUmask)[2] * (csv1r$z - 1)) + ((csv1r$x -1) * dim(EdUmask)[1]) + csv1r$y)
 EdUind1 <- EdUmasklogic[ind4]
 csv2r<-csv1r[EdUind1,]
## points1<-csv1[EdUind1,]
 write.csv(points1, file =edu.fpath)

 ind3 <- ((dim(difmask)[1] * dim(difmask)[2] * (csv1r$z - 1)) + ((csv1r$x -1) * dim(difmask)[1]) + csv1r$y)
 difind1 <- difmasklogic[ind3]
 csv2r<-csv1r[difind1,]
 points1<-csv1[difind1,]
 write.csv(csv2r, file =dif.fpath)
 

  }
  else{
    print ("MCNR mask and segmented image are of different dimensions.")
  }
  options(warn = oldw)
  return (0)
}

if (identical (environment (), globalenv ()))
  quit (status=main (commandArgs (trailingOnly = TRUE)))
