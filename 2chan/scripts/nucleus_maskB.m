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

function dapi_mask_clean = get_dapi_mask (dapi)
  se_2d = strel ("disk", 10, 0);
  se_2d2 = strel ("disk", 1, 0);
  se_2d3 = strel ("disk", 3, 0);
  se_3d = strel ("arbitrary", repmat (getnhood (se_2d), [1 1 3]));
  se_3d2 = strel ("arbitrary", repmat (getnhood (se_2d2), [1 1 3]));
  se_3d3 = strel ("arbitrary", repmat (getnhood (se_2d3), [1 1 3]));

  dapi_mask = im2bw (imdilate (dapi, se_3d), graythresh (dapi(:)));
  dapi_mask = bwfill (dapi_mask, "holes", 8);
  dapi_mask = reshape (dapi_mask, size (dapi));
  dapi_mask = imerode (dapi_mask, se_2d);
  dapi_mask = imclose (dapi_mask, se_2d);
  dapi_mask = bwareaopen (dapi_mask, 10000, ones (3, 3));

#gets lamina segmented:
dapi_mask2 = dapi_mask;
dapi_mask2 = imerode (dapi_mask2, se_3d);
dapi_mask2 = imerode (dapi_mask2, se_3d);
dapi_mask2 = imerode (dapi_mask2, se_3d);
dapi_mask2 = imerode (dapi_mask2, se_3d);
dapi_mask2 = imerode (dapi_mask2, se_3d);
dapi_mask2 = imerode (dapi_mask2, se_3d);
dapi_mask2 = imerode (dapi_mask2, se_3d);
dapi_mask2 = imerode (dapi_mask2, se_3d);
dapi_mask2 = imerode (dapi_mask2, se_3d);
dapi_mask2 = imdilate (dapi_mask2, se_2d);
dapi_mask2 = imdilate (dapi_mask2, se_2d);
dapi_mask2 = imdilate (dapi_mask2, se_2d);
dapi_mask2 = imdilate (dapi_mask2, se_2d);
dapi_mask2 = imdilate (dapi_mask2, se_2d);
dapi_mask2 = imdilate (dapi_mask2, se_2d);
dapi_mask2 = imdilate (dapi_mask2, se_2d);
dapi_mask2 = imdilate (dapi_mask2, se_2d);


#segments chromocenters:
mx = max(unique (dapi));
mxthof8 = mx*0.4;
mxof16b = im2double (mxthof8);
dapi_mask3 = im2bw (dapi(:), mxof16b);

dapi_mask3 = bwfill (dapi_mask3, "holes", 8);
dapi_mask3 = reshape (dapi_mask3, size (dapi));
dapi_mask3 = imdilate (dapi_mask3, se_3d2);
dapi_mask3 = imdilate (dapi_mask3, se_3d2);
dapi_mask3 = imerode (dapi_mask3, se_2d2);
dapi_mask3 = imerode (dapi_mask3, se_2d2);
dapi_mask3 = imerode (dapi_mask3, se_2d2);
dapi_mask3 = imerode (dapi_mask3, se_2d2);
dapi_mask3 = imerode (dapi_mask3, se_2d2);
dapi_mask3 = imerode (dapi_mask3, se_2d2);
dapi_mask3 = imdilate (dapi_mask3, se_3d2);
dapi_mask3 = imdilate (dapi_mask3, se_3d2);
dapi_mask3 = imdilate (dapi_mask3, se_3d2);
dapi_mask3 = imdilate (dapi_mask3, se_3d3);


##puts them together:
dapi_mask = dapi_mask - dapi_mask2;
dapi_mask_clean = or(dapi_mask,dapi_mask3);


endfunction

function [rv] = main (argv)

  if (numel (argv) != 3)
    error ("usage: nucleus_mask DAPI_CHANNEL IN_FPATH MASK_FPATH");
  endif
 
    dapi_channel = str2double (argv{1});
  in_fpath = argv{2};
  mask_fpath = argv{3};

  if (isnan (dapi_channel))
    error ("nucleus_mask: invalid DAPI_CHANNEL number");
  elseif (fix (dapi_channel) != dapi_channel)
    error ("nucleus_mask: DAPI_CHANNEL must be an integer");
  endif

  ## We are assuming that the image has 2 channels.  Maybe we shouldn't?
  n_img = numel (imfinfo (in_fpath));
  dapi = imread (in_fpath, "Index", dapi_channel:2:n_img);

  mask = get_dapi_mask (dapi);

  ##we are having to write this logical mask as an 8bit due to R unable to open "1bit" in subsequent steps
  imwrite (uint8(mask), mask_fpath);

  rv = 0;
  return;
endfunction

main (argv ());
