pkg load image;
pkg load statistics;

function [rv] = main (argv)

  if(numel (argv) !=4)
    error("script requires 4 arguments");
  endif

  centroid_C1_fpath = argv{1};
  centroid_C2_fpath = argv{2};
  NNA_C1_fpath = argv{3};
  NNA_C2_fpath = argv{4};

#read centroid file and, if needed, filter on size of the segmented centroid (column 4)
  cen11 = csvread (centroid_C1_fpath);
##  cen11 = cen11(2:end,:);
##  cen11 = cen11(cen11(:,4) >= 11600000,:);
  cen1 = cen11(:, [2,3,4]);
  cen22 = csvread (centroid_C2_fpath);
  cen22 = cen22(2:end,:);
##  cen22 = cen22(cen22(:,4) >= 13000000,:);
  cen2 = cen22(:, [2,3,4]);


############################################################
############################################################

#finding nearest neighbour (nn_dist) and the index number of the closest centroid (nn_idx) using euclidean measurements
  [nn_dist_1to2, nn_idx_1to2] = min (pdist2 (cen2, cen1, "sqeuclidean"));
    nn_dist_1to2 = sqrt (nn_dist_1to2);

  [nn_dist_2to1, nn_idx_2to1] = min (pdist2 (cen1, cen2, "sqeuclidean"));
    nn_dist_2to1 = sqrt (nn_dist_2to1);

#save two files simultaniously (distance C1-C2 and distance C2-C1) 
  fid= fopen(NNA_C1_fpath, "w");
  fprintf (fid, "x, y, z, volume, meanintensity, NN_dist_1to2, NN_idx1to2\n");
  fclose(fid);
  csvwrite(NNA_C1_fpath, [cen11 nn_dist_1to2(:) nn_idx_1to2(:)], "-append");

  fid= fopen(NNA_C2_fpath, "w");
  fprintf (fid, "x, y, z, volume, meanintensity, NN_dist_2to1, NN_idx_2to1\n");
  fclose(fid);
  csvwrite(NNA_C2_fpath, [cen22 nn_dist_2to1(:) nn_idx_2to1(:)], "-append");

  rv = 0;
  return;

endfunction

main (argv ());
