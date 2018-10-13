
#! /usr/bin/env octave
##
## Copyright (C) 2016-2018 Ezequiel Miron <eze.miron@bioch.ox.ac.uk>
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
  pkg load statistics;

function [rv] = main (argv)

  if (numel (argv) != 3)
  error ("script requires 3 arguments");
  endif

  nucleus_mask = argv{1};
  EdU_mask = argv{2};
  out_fpath = argv{3};

  nucleus = imread (nucleus_mask, "index","all");
  edu = imread (EdU_mask, "index","all");

  mask = nucleus - edu;

  imwrite (EdU_mask, out_fpath);

    rv = 0;
  return;
endfucntion
 
main (argv ());

