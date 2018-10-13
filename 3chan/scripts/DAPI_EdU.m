  pkg load image;
  pkg load statistics;

function [rv] = main (argv)

  if (numel (argv) != 3)
  error ("script requires 3 arguments");
  endif

  nucleus_fpath = argv{1};
  edu_fpath = argv{2};
  out_fpath = argv{3};

  nucleus = imread (nucleus_fpath, "index", "all");
  edu = imread (edu_fpath, "index", "all");

  mask = minus(nucleus,edu);
  imwrite(uint8(mask),out_fpath);

    rv = 0;
  return;
endfunction

main (argv ());
