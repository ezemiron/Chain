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
  cen1 = cen11(:, [1,2,3]);
  cen22 = csvread (centroid_C2_fpath);
##  cen22 = cen22(2:end,:);
##  cen22 = cen22(cen22(:,4) >= 13000000,:);
  cen2 = cen22(:, [1,2,3]);


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

############################################################
############################################################

#finding all the distances between centroids of a single channel
  [nn_dist_1to1] = pdist2 (cen1, cen1, "sqeuclidean");
    nn_dist_1to1 = sqrt (nn_dist_1to1);

  [nn_dist_2to2] = pdist2 (cen2, cen2, "sqeuclidean");
    nn_dist_2to2 = sqrt (nn_dist_2to2);

#Sort the matrix ascending and take the second-eleventh measure (fisrt measure is 0 because it finds itself)    
  nn_dist_1to1 = sort(nn_dist_1to1,'ascend');
    nn_dist_1to1 = nn_dist_1to1(2:11,:);
      nn_dist_1to1 = sort(nn_dist_1to1,'descend');
      avg_dist_1to1 = mean(nn_dist_1to1);
        nn_dist_1to1 = rot90(nn_dist_1to1,3);
        avg_dist_1to1 = rot90(avg_dist_1to1,3);

  nn_dist_2to2 = sort(nn_dist_2to2,'ascend');
    nn_dist_2to2 = nn_dist_2to2(2:11,:);
      nn_dist_2to2 = sort(nn_dist_2to2,'descend');
      avg_dist_2to2 = mean(nn_dist_2to2);
        nn_dist_2to2 = rot90(nn_dist_2to2,3);
        avg_dist_2to2 = rot90(avg_dist_2to2,3);

#save two files simultaniously (X,Y,Z coordinate with 10 distances)
  csv_fpath1 = strrep(NNA_C1_fpath, 'C1-C2_NNA.csv', 'C1-C1_NNA.csv');
  fid= fopen(csv_fpath1, "w");
  fprintf (fid, "x,y,z,volume,meanintensity,nn_dist1,nn_dist2,nn_dist3,nn_dist4,nn_dist5,nn_dist6,nn_dist7,nn_dist8,nn_dist9,nn_dist10,AVG_nn_dist\n");
  fclose(fid);
  csvwrite(csv_fpath1, [cen11 nn_dist_1to1 avg_dist_1to1], "-append");

csv_fpath2 = strrep(NNA_C1_fpath, 'C1-C2_NNA.csv', 'C2-C2_NNA.csv');
  fid= fopen(csv_fpath2, "w");
  fprintf (fid, "x,y,z,volume,meanintensity,nn_dist1,nn_dist2,nn_dist3,nn_dist4,nn_dist5,nn_dist6,nn_dist7,nn_dist8,nn_dist9,nn_dist10,AVG_nn_dist\n");
  fclose(fid);
  csvwrite(csv_fpath2, [cen22 nn_dist_2to2 avg_dist_2to2], "-append");

  rv = 0;
  return;

endfunction

main (argv ());
