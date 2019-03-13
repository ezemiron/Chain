#! /usr/bin/env octave
##
## Copyright (C) 2014-2016 Carnë Draug <carandraug+dev@gmail.com>
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

#old:
function edu_mask_clean = get_edu_mask_old (edu)
  se_2d = strel ("disk", 1, 0);
  edu_mask = im2bw (imdilate (edu, se_3d), graythresh (edu(:)));

  edu_mask = im2bw (edu, graythresh (edu(:)));
  edu_mask = reshape (edu_mask, size (edu));
  edu_mask = imclose (edu_mask, se_2d);
  
  edu_mask_clean = edu_mask;
  #edu_mask_clean = bwareaopen (edu_mask, 10000, ones (3, 3));
endfunction



#normal:
function edu_mask_clean = get_edu_mask (edu)
  se_2d = strel ("disk", 1, 0);
  se_3d = strel ("arbitrary", repmat (getnhood (se_2d), [1 1 3]));

  edu_mask = im2bw (edu, graythresh (edu(:)));
  edu_mask = reshape (edu_mask, size (edu));
  edu_mask = imclose (edu_mask, se_2d);
  
  #edu_mask_clean = edu_mask;

#size exclusion for background spots
  edu_mask_med_large = bwareaopen (edu_mask, 50, ones (3, 3));

#size exclusion for large spots
  edu_mask_large = bwareaopen (edu_mask_med_large, 150, ones (3, 3));
  edu_mask_clean = minus(edu_mask_med_large,edu_mask_large);
endfunction



#loose:
function edu_mask_clean = get_edu_mask_loose (edu)
  se_2d = strel ("disk", 1, 0);
  se_3d = strel ("arbitrary", repmat (getnhood (se_2d), [1 1 3]));

  edu_mask = im2bw (edu, graythresh (edu(:)));
  #edu_mask = im2bw (imdilate (edu, se_3d), graythresh (edu(:)));
  edu_mask = reshape (edu_mask, size (edu));
  edu_mask = imclose (edu_mask, se_2d);
  
  #edu_mask_clean = edu_mask;

#size exclusion for background spots
  edu_mask_med_large = bwareaopen (edu_mask, 50, ones (3, 3));

#size exclusion for large spots
  edu_mask_large = bwareaopen (edu_mask_med_large, 150, ones (3, 3));
  edu_mask_clean = minus(edu_mask_med_large,edu_mask_large);

#loosen mask
  edu_mask_clean = imdilate (edu_mask_clean, se_3d);
  edu_mask_clean = imdilate (edu_mask_clean, se_3d);
  edu_mask_clean = imdilate (edu_mask_clean, se_3d);
  edu_mask_clean = imdilate (edu_mask_clean, se_3d);
  edu_mask_clean = imdilate (edu_mask_clean, se_3d);
  edu_mask_clean = imerode (edu_mask_clean, se_3d);
  edu_mask_clean = imerode (edu_mask_clean, se_3d);
  edu_mask_clean = imerode (edu_mask_clean, se_3d);
  edu_mask_clean = imerode (edu_mask_clean, se_3d);

#fill holes
edu_mask_clean = imfill (edu_mask_clean, "holes");
#edu_mask_clean = bwfill (edu_mask_clean, "holes", 8);
endfunction






function [rv] = main (argv)

  if (numel (argv) != 3)
    error ("usage: nucleus_mask edu_CHANNEL IN_FPATH MASK_FPATH");
  endif

  ## TODO: Could use bioformats to default to the channel closest to edu?
  ##          Only if we use the dv files. The tif files saved by ImageJ
  ##          that we got Yolanda had lost that information (instead they
  ##          simply attached a blue LUT to the channel).
  edu_channel = str2double (argv{1});
  in_fpath = argv{2};
  mask_fpath = argv{3};

  if (isnan (edu_channel))
    error ("nucleus_mask: invalid edu_CHANNEL number");
  elseif (fix (edu_channel) != edu_channel)
    error ("nucleus_mask: edu_CHANNEL must be an integer");
  endif

  ## We are assuming that the image has 2 channels.  Maybe we shouldn't?
  n_img = numel (imfinfo (in_fpath));
  edu = imread (in_fpath, "Index", edu_channel:3:n_img);

#choose "normal" or loose mask
  mask = get_edu_mask_old (edu);
  #mask = get_edu_mask (edu);
  #mask = get_edu_mask_loose (edu);

  ##we are having to write this logical mask as an 8bit due to R unable to open "1bit" in subsequent steps
  imwrite (uint8(mask), mask_fpath);

  rv = 0;
  return;
endfunction

main (argv ());
