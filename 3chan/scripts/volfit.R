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



# Compute mu and sigma for fitted gaussian:
fitting <- function(marker) {
#  rmarker <- round(marker)
  tabmarker <- table(marker)

  x <- as.numeric(rownames(tabmarker))
  #tab <- as.data.frame(t(tabract))
  #x <- as.vector(tab[2])
  rownames(tabmarker) <- c()
  p <- as.vector(tabmarker/max(tabmarker))
  prob <- cumsum(p);
  prob <- prob / prob[length(prob)]

# Compute sum of squared residuals to a fit
  f <- function(q) {
    res <- pnorm(x, q[1], q[2]) - prob
    sum(res * res)
  }

  #estimates:
  mu <- mean (marker[,1])
  stdev <-sd (marker[,1])

  # Find the least squares fit using function f defined previously
  coeff <-(fit <- nlm(f, c(mu, stdev)))$estimate

  ## To plot the cumulative fit
  #plot(x, prob) + curve(pnorm(x, coeff[1], coeff[2]), add=TRUE)
  ##to check the sum of residuals:
  #sumresid <- sum ((pnorm(x, coeff[1],coeff[2])-prob)^2)
  return(coeff)
}


main <-function (argv)
{
  oldw <- getOption("warn") # to suppress warnings. (If the file was writen in Fiji instead of IJ then there are some unknown fields with tags 50838 and 50839 which cause a warning when read in R.)
  options(warn = -1)

#eg:
#markerd2b_fpath <- "EM16-12-A_C127_S4_1514_H3K4me3-594_Sytox_G1_02_SIR_EAL_THR_md2b.csv"
#randomd2b_fpath <- "EM16-12-A_C127_S4_1514_H3K4me3-594_Sytox_G1_02_SIR_EAL_THR_rd2b.csv"

#marker <- rbind(active,active2)
  markerd2b_fpath <- argv[1]
  randomd2b_fpath <- argv[2]
  

  if (length (argv) != 2)
    stop ('usage error: Requires at least 2 arguments')

  markerd2b <- read.csv(markerd2b_fpath, header=FALSE);
  randomd2b <- read.csv(randomd2b_fpath, header=FALSE);

#get fitted mu and sigma:
  mfitted <- fitting(markerd2b)
  rfitted <- fitting(randomd2b)

  fittable <- cbind(mfitted,rfitted)
  rownames(fittable) <- c("mu","sigma")
  ##write to standard output:
  write.csv(fittable, file ="")
#write.csv(normtable, file = "EM16-12-A_C127_S4_1514_H3K4me3-594_Sytox_G1_02_SIR_EAL_THR_d2bnorm.csv")

  options(warn = oldw)
  
  return (0)
}

if (identical (environment (), globalenv ()))
  quit (status=main (commandArgs (trailingOnly = TRUE)))


