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


#the d2b (distance to boundary) function:

function [d2boundcomp] = d2b (cen1sing, cen2sing) #centroids2 are those of the boundary
  cen1sing = cen1sing(2:length(cen1sing),:);
  mag = floor(log10 (length(cen1sing)));
  if (mag < 2)
    div = 1; #do it in a single chunk
    cen1div = floor(length (cen1sing)/div);
  endif
  if (mag == 2)
    div = 10; # do it in 10 chunks
    cen1div = floor(length (cen1sing)/div);
  endif
  if (mag > 2)
    cen1div = 100; #for large number of spots always have divisions of 100 elements max and work out number of divisions accordingly:
    div = floor(length(cen1sing)/cen1div); 
  endif
#break d2b in chunks to avoid running out of Octave index elements or memory:
  i = 1;
  comp = zeros(1,1);
  while (i < (div+1))
    cen1singa = cen1sing((((i*cen1div)-cen1div)+1):(i*cen1div),:);
    d2bound = min (pdist2 (cen2sing, cen1singa, "sqeuclidean"));
    d2bound = permute(d2bound, [2,1]);
    comp = cat(1,comp,d2bound);
    i +=1;
  endwhile
#calculate any remaining coordinates not factorised to last chunk, if any:
  dif = length(cen1sing)-(length(cen1div)*div);
  if (dif > 0)
    j = i-1;
    cen1singa = cen1sing(((cen1div*j)+1):length(cen1sing),:);
    d2bound = min (pdist2 (cen2sing, cen1singa, "sqeuclidean"));
    d2bound = permute(d2bound, [2,1]);
    comp = cat(1,comp,d2bound);
  endif
#tidy:
  #comp = comp (2:length(comp),:);
  d2boundcomp = sqrt (comp);
endfunction




#the main script function:

function [rv] = main (argv)

  if (numel (argv) != 7)
    error ("script requires 7 arguments");
  endif
  
  refclass = str2double (argv{1});
  SEG_fpath = argv{2};
  cen_fpath = argv{3};
  bound_fpath = argv{4};
  mask_fpath = argv{5};
  d2bmarker_fpath = argv{6};
  d2brandom_fpath = argv{7};
  vox = [20,20,20];

#read:

  SEG = imread (SEG_fpath, "index", "all");

  cen1 = csvread (cen_fpath); 

  bound = imread (bound_fpath, "index", "all");

#prep:
  refclass = floor((refclass/(length (unique(SEG))-1))*max (unique(SEG)));

  cen1 = [cen1(:,2),cen1(:,1),cen1(:,3)];
  #to get cen1 to the closest voxel and ignore sub voxel positions:
#  cen1 = round(cen1./vox).*vox;
  cen1 = cen1(2:length(cen1),:);

  SEGlogic = (SEG>refclass);
  #SEGlogic = logical(SEGlogic);

  bound_dim = size(bound);
#drop singleton 3rd dimension:
  bound_dim(:,3) = [];
# or could do:
#bound_dim = size(squeeze (bound));

#Turning bound into yxz cen2:
ind = find(bound); % linear index
[row,col,pag] = ind2sub(bound_dim ,ind); % sub indices
cen2 = [row,col,pag];
cen2 = cen2.*vox;


#finding distances using d2b function defined earlier:
  cen1sing = single(cen1);
  cen2sing = single(cen2);

  md2bound = d2b(cen1sing, cen2sing);

#Filter centroids into those within and without reference class:
  cen1vox = cen1./vox;
  #if cen1 is not rounded before then cen1vox needs rounding now:
  cen1vox = round(cen1vox);

  cen1ind = sub2ind( bound_dim,cen1vox(:,1), cen1vox(:,2), cen1vox(:,3));
  ind_mask = SEGlogic(cen1ind);
  #d2binternal = md2bound (ind_mask);
  md2bound(!ind_mask) *= -1;
  
  ind2_mask = bound(cen1ind);
  halfind2 = randi([0 1],length (ind2_mask),1);
  ind3_mask = and(ind2_mask,halfind2);
  md2bound(ind3_mask) *= -1;

#make a random distribution of centroids at 1/200th concentration:
#A = randi([0 1], m,n,o) where m, n and o are dimensions
mask = imread (mask_fpath, "index", "all");
mask = squeeze(mask);
randm = randi([-198 1], bound_dim(1),bound_dim(2),bound_dim(3));
randmlogic = randm>0;
maskedrandm = randmlogic.*mask;
#Turning random spots into yxz cen3:
randmind = find(maskedrandm); % linear index
[row,col,pag] = ind2sub(bound_dim ,randmind); % sub indices
cen3 = [row,col,pag];
cen3 = cen3.*vox;
cen3sing = single(cen3);
#get sub voxel "resolution":
subvox1 = randi([(floor(vox(1)/2)*-1) floor(vox(1)/2)], length(cen3),1);
subvox2 = randi([(floor(vox(2)/2)*-1) floor(vox(2)/2)], length(cen3),1);
subvox3 = randi([(floor(vox(3)/2)*-1) floor(vox(3)/2)], length(cen3),1);
randmextra = [subvox1, subvox2, subvox3];
cen4 = cen3.+randmextra;
cen4sing = single(cen4);

#find distances:
  rd2bound = d2b(cen4sing, cen2sing);
  randmind_mask = SEGlogic(randmind);
#  randmind_mask = randmind_mask (2:length(randmind_mask));
  rd2bound(!randmind_mask) *= -1;

  randmind2_mask = bound(randmind);

  halfrandmind2 = randi([0 1],length (randmind2_mask),1);
  randmind3_mask = and(randmind2_mask,halfrandmind2);
  rd2bound(randmind3_mask) *= -1;

#write 2 csvs in one go:
  csvwrite(d2bmarker_fpath, md2bound);
  csvwrite(d2brandom_fpath, rd2bound);
## write to sdt out:
#  if (int_true ==1)
#    d2binternal = d2boundcomp (ind_mask);
#    printf ("d2bInt-nm\n");
#    printf ("%f\n", d2binternal');
#  else
#    d2boundcomp(!ind_mask) *= -1;
#    printf ("d2bFull-nm\n");
#    printf ("%f\n", d2boundcomp');
#  endif
  rv = 0;
  return;
endfunction

main (argv ());
