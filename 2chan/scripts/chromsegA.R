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

  dapi_channel_n <- as.numeric (argv[1])
  mask_fpath <- argv[2]
  img_fpath <- argv[3]
  maska_fpath <- argv[4]
  out_fpath <- argv[5]

  if (length (argv) < 4)
    stop ('usage error: m5in_example requires at least 5 arguments')


  
  img <- readTIF (img_fpath) #This opens the 16bit tiff with all channels correclty

  ## FIXME use dapi_n_channel to index chromatin
  chromatin <- img[ , ,dapi_channel_n, ] #finds the chromatin channel

  mask <- readTIF (mask_fpath) # opens the mask
  mask <- mask*255
  mask = mask==1
  ## this should no longer happen given that both of these are taken from
  ## the same file and not manually selected
  if (! identical (dim (mask), dim (chromatin)))
    stop ('usage error: Chromatin landscape and MASK are of different sizes')

  ## The segmentation:
  seg <- segment (chromatin, 7, 0.1, 1/3, mask=mask, maxit=20, varfixed=TRUE,
                  inforce.nclust=TRUE, start="equal",silent=TRUE)

  ## Rescaling so imageJ lut can read the output
  segscale = seg$class / length (seg$mu)

  ## save as an 8 bit tif
  writeTIF (segscale, out_fpath)

  options(warn = oldw)
  
  return (0)
}

if (identical (environment (), globalenv ()))
  quit (status=main (commandArgs (trailingOnly = TRUE)))
