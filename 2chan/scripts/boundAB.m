#! /usr/bin/env octave
##
## Copyright (C) 2018 Ezequiel Miron <eze.miron@bioch.ox.ac.uk>
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, see <http://www.gnu.org/licenses/>.

pkg load image;

function [rv] = main (argv)

  if (numel (argv) != 3)
    error ("usage: bound_FPATH maskA/B_FPATH boundA/B_FPATH");
  endif
 
  bound_fpath = argv{1};
  mask_fpath = argv{2};
  out_fpath = argv{3};

  ## We are assuming that the image has 1 channels.  Maybe we shouldn't?
  #n_img = numel (imfinfo (in_fpath));
  #n_img = 1;
  bound = imread (bound_fpath, "Index", "all");
  mask = imread (mask_fpath, "Index", "all");


##puts them together:
bound_mask = bound - mask;
#bound_mask_clean = or(bound_mask,bound_mask3);



  ##we are having to write this logical mask as an 8bit due to R unable to open "1bit" in subsequent steps
  imwrite (uint8(bound_mask), out_fpath);

  rv = 0;
  return;
endfunction

main (argv ());
