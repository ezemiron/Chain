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
  oldw <- getOption("warn") # to suppress warnings. (If the file was writen in Fiji instead of IJ then there are some unknown fields with tags 50838 and 50839 which cause a warning when read in R.)
  options(warn = -1)

library (EBImage, warn.conflicts = FALSE, verbose = FALSE)
library (bioimagetools, warn.conflicts = FALSE, verbose = FALSE)


  centroids.fpath <- argv[1]

  seg.fpath <- argv[2]

  MCNRmask.fpath <- argv[3]

  if (length (argv) < 3){
    stop ('usage error: main_example requires at least 3 arguments')
  }
  
  classA <- as.numeric(1)
  classB <- NULL

  img.classes <- readTIF (seg.fpath);
  MCNRmask <- readTIF (MCNRmask.fpath)
  MCNRmask <- MCNRmask[,,1,]
#  MCNRlogic <- as.logical(round(MCNRmask))
  
  if (length(img.classes)==length(MCNRmask)){

    csv1 <- read.csv(centroids.fpath, header=TRUE);
    MCNRlogic = MCNRmask==1
    
    vox  <- c(41, 41, 125)
    x.microns <- (dim(MCNRmask)[2]*vox[2])/1000
    y.microns <- (dim(MCNRmask)[1]*vox[1])/1000
    z.microns <- (dim(MCNRmask)[3]*vox[3])/1000
    img.classes <- round(img.classes*7)
    
    csv1r <- round (data.frame (csv1[1]/vox[2], csv1[2]/vox[1], csv1[3]/vox[3]))
    ind1 <- ((dim(MCNRmask)[1] * dim(MCNRmask)[2] * (csv1r$z - 1)) + ((csv1r$x -1) * dim(MCNRmask)[1]) + csv1r$y)

MCNRind1 <- MCNRlogic[ind1]

points<-csv1[MCNRind1,]/1000


  #reshape from xyz to yxz for use in distance2border (ReadTIF reads yxz)
  #points <- data.frame(points[2],points[1],points[3])
  

  d2b <- distance2border(points,img.classes,x.microns,y.microns,z.microns,class1=classA,class2=classB, silent=TRUE)

  #currently, negative values in d2b represent points outside classA, 
  #but it should be the other way round according to 
  #the documentation, so fixed it: 
  d2b <- d2b*(-1)


  ##write(compile1, file=stdout())
  write.csv(d2b, file ="")

  }
  else{
    print ("MCNR mask and segmented image are of different dimensions.")
  }
  options(warn = oldw)
  
  return (0)
}

if (identical (environment (), globalenv ()))
  quit (status=main (commandArgs (trailingOnly = TRUE)))
