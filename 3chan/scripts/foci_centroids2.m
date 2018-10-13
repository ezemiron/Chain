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
pkg load bioformats;
javaMethod ("enableLogging", "loci.common.DebugTools", "ERROR");

function [centroids, r_area, meanintensity] = get_centroids (im, mask)
  if (! isinteger (im))
    ## we could just use abs() but if data is floating point, we may
    ## be sure we are not messing up something.
    error ("data is not integer class");
  endif

  bw = im2bw (im, graythresh (im(mask)));
  im(! bw | ! mask) = 0;
  imc = imcomplement (im);
  w_lines = ! watershed (imc);
  im(w_lines) = 0;
  bw = bwareaopen (logical (im), 10);

  props = regionprops (bw, im, {"Area", "WeightedCentroid","MeanIntensity"});
  centroids = cell2mat ({props(:).WeightedCentroid}(:));
  r_area = [props(:).Area];
  meanintensity = [props(:).MeanIntensity];
endfunction

# function [voxel_sizes] = get_pixel_sizes (img_fpath);
#  reader = bfGetReader (img_fpath); # FIXME: this sometimes segfaults
#  metadata_store = reader.getMetadataStore ();

#  micron = java_get ("ome.units.UNITS", "MICROM");
#  voxel_sizes = [
#    metadata_store.getPixelsPhysicalSizeX(0).value(micron)
#    metadata_store.getPixelsPhysicalSizeY(0).value(micron)
#    metadata_store.getPixelsPhysicalSizeZ(0).value(micron)
#  ]';
# endfunction

function [rv] = main (argv)

  if (numel (argv) != 4)
    error ("usage: foci_centroids FOCI_CHANNEL MASK_FPATH IMG_FPATH num_of_channels");
  endif

  ## TODO: Could use bioformats to default to the channel closest to DAPI?
  ##          Only if we use the dv files. The tif files saved by ImageJ
  ##          that we got Yolanda had lost that information (instead they
  ##          simply attached a blue LUT to the channel).
  foci_channel = str2num (argv{1});
  #foci_channel = argv{1};
  mask_fpath = argv{2};
  img_fpath = argv{3};
  num_channels = str2num (argv{4});
  
  if (isnan (foci_channel))
    error ("foci_centroids: invalid FOCI_CHANNEL number");
  elseif (fix (foci_channel) != foci_channel)
    error ("foci_centroids: FOCI_CHANNEL must be an integer");
  endif

  dapi_mask = imread (mask_fpath, "index", "all");
  dapi_mask = logical(dapi_mask);  
#  dapi_mask2 = imcomplement(dapi_mask);

## Untick when you want tot save the inverted mask as a inv-mask.tif
#  invmask_out_fpath = strrep(mask_fpath, 'mask.tif', 'inv-mask.tif');
#  imwrite (uint8(dapi_mask2), invmask_out_fpath);



  n_img = numel (imfinfo (img_fpath));
  foci = imread (img_fpath, "Index", foci_channel:num_channels:n_img);

  if (! size_equal (dapi_mask, foci))
    error ("foci_centroids: MASK and IMG are images of different sizes");
  endif

  [centroids, r_area, meanintensity] = get_centroids (foci, dapi_mask);
#  voxel_sizes = get_pixel_sizes (img_fpath);

## Multiply the cooridinates with voxel sizes
  voxel_sizes = [
   41
   41
   125
  ]';
  centroids(:,3) = []; # Remove coordinates from singleton dimension 3
  centroids .*= voxel_sizes;

  
## Multiply the number of voxels witht the volume of the a voxel (41*41*125)
  r_area .*= 210125;
  r_volume = [r_area(:)];
  meanintensity = [meanintensity(:)];

   printf ("x, y, z, volume (nm3), meanIntensity\n")
   printf ("%f,%f,%f,%f,%f\n", [centroids,r_volume,meanintensity]')

#  [centroids2, r_area2, meanintensity2] = get_centroids (foci, dapi_mask2);
#  centroids2(:,3) = [];
#  centroids2 .*= voxel_sizes;

#  r_area2 .*= 210125;
#  r_volume2 = [r_area2(:)];
#  meanintensity2 = [meanintensity2(:)];

#  csv_fpath = strrep(mask_fpath, 'mask.tif', 'centroids-inv.csv');
#  fid= fopen(csv_fpath, "w");
#  fprintf (fid, "x, y, z, volume (nm3), meanintensity\n");
#  fclose(fid);
#  csvwrite(csv_fpath, [centroids2 r_volume2(:) meanintensity2(:)], "-append");

  rv = 0;
  return;
endfunction

main (argv ());
