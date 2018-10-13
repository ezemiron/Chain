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

#eg:
#fit_fpath <- "fittable"

  fit_fpath <- argv[1]


  if (length (argv) != 1)
    stop ('usage error: Requires 1 argument')

  fitdf <- read.csv(fit_fpath, header=TRUE);

#get fitted mu and sigma:
  mfitted <- fitdf$mfitted
  rfitted <- fitdf$rfitted


  #generate a set of standard distances and find what the normal on these are:
  stddist <- seq(-400,400)
  markerfit <- dnorm(stddist, mfitted[1], mfitted[2])
  randomfit <- dnorm(stddist, rfitted[1], rfitted[2])
  #normalise the fit to unity:
  markerfit <- markerfit/max(markerfit)
  randomfit <- randomfit/max(randomfit)
  ##To check them graphically:
  #plot (markerfit~stddist)
  #plot (randomfit~stddist)



  d2bnorm <- markerfit/randomfit
  logd2bnorm <- log2(d2bnorm)

  normtable <- cbind(stddist, markerfit, randomfit, d2bnorm, logd2bnorm)
  ##To check them graphically:
  #plot (d2bnorm~stddist)
  #plot (logd2bnorm~stddist)

  ##write to standard output:
  write.csv(normtable, file ="")
#write.csv(normtable, file = "EM16-12-A_C127_S4_1514_H3K4me3-594_Sytox_G1_02_SIR_EAL_THR_d2bnorm.csv")

  options(warn = oldw)
  
  return (0)
}

if (identical (environment (), globalenv ()))
  quit (status=main (commandArgs (trailingOnly = TRUE)))


