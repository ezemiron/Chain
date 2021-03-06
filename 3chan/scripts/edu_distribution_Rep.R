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
  MCNRmask.fpath <- argv[3]
  
  
  if (length (argv) < 3)
    stop ('usage error: main_example requires at least 3 arguments')


  MCNRmask <- readTIF (MCNRmask.fpath)

##comment out next two lines after finishing manual controls as these are not used in automated run. and uncomment the last line.  Also swap the linear index lines 
# MCNRmask <-round(MCNRmask[,,1,seq(1, as.numeric(attributes(MCNRmask)$images), 2)])
#  MCNRmask <-round(aperm(MCNRmask[,,1,seq(1, as.numeric(attributes(MCNRmask)$images), 2)], c(2,1,3)))

# MCNRmask <- MCNRmask[,,1,]

## read segmented image
  img <- readTIF (segimg.fpath) #This opens the seg 16bit tiff with all channels correctly
  if (length(img)==length(MCNRmask)){

    csv1 <- read.csv(centroids.fpath, header=TRUE);
    csv1 <- csv1[,1:3]
    MCNRlogic = (1==MCNRmask)
    
    vox  <- c(41, 41, 125)
print(dim(MCNRmask)[2])
    x.microns <- (dim(MCNRmask)[2]*vox[2])/1000
    y.microns <- (dim(MCNRmask)[1]*vox[1])/1000
    z.microns <- (dim(MCNRmask)[3]*vox[3])/1000
    img <- round(img*7)
    
##comment out next three lines after finishing manual controls as these are not used in automated run  and uncoment the final line
#  csv1names <- colnames(csv1)
#  csv1r <- round (data.frame (csv1[2]/vox[2], csv1[1]/vox[1], csv1[3]/vox[3]))
#  colnames(csv1r) <- csv1names

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
  csv2r<-csv1r[MCNRind1,]
  
  #second trimmed linear index
  ind2 <- ((dim(MCNRmask)[1] * dim(MCNRmask)[2] * (csv2r$z - 1)) + ((csv2r$x -1) * dim(MCNRmask)[1]) + csv2r$y)


#  points<-csv1[MCNRind1,]/1000

  #indexing the image using the second trimmed linear index
  imind1 <- img[ind2]



  ## for normalisation find how many pixels in each class,ie:
  ## count the size of bins, and remove background:
    classnum <- as.numeric (table (img))
    if (0 == min (img))
      classnum <- classnum[2:length(classnum)]

    classnumn <- classnum / sum (classnum)

    ## count how many spots fell in each bin
    dist1t <- tabulate (imind1, nbins=length (classnum))

    ## Removes spots found in background if any:
    if (0 == min (imind1)){
      bgnum <- as.numeric (table (imind1)[1])
      tnum <- sum (table (imind1))
    }

    dist1n <- dist1t / sum (dist1t)

    norm1 <- dist1n / classnumn
    lognorm1 <- log2 (norm1)
    compile1 <- data.frame (class=1:length(classnum), dist1t, dist1n, classnum, classnumn, norm1, lognorm1)

    ##write(compile1, file=stdout())
    write.csv(compile1, file ="")
  }
  else{
    print ("MCNR mask and segmented image are of different dimensions.")
  }
  options(warn = oldw)
  return (0)
}

if (identical (environment (), globalenv ()))
  quit (status=main (commandArgs (trailingOnly = TRUE)))
