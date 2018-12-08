# Cha*i*N - Cha*i*N analysis of the *in situ* Nucleome

## What?
A tool to dissect the spatial distribution of the mammalian nucleome relative to chromatin from super resolution 3D structured illumination microscopy (3D SIM) data. This tool was developed as part of my doctoral thesis and was used in a related publication.

## Where?
The original versions of Cha*i*N are currently in the Oxford biochemistry Micron servers under `wolf4192/Chain/`
And in github: git@github.com:ezemiron/Chain.git
https://github.com/ezemiron/Chain.git

One version exists for 2 Channel data with chromatin in second channel and spots channel in first
and another version for 3 channel data, with chromatin in third channel and spots channel in first and second and channel.

##### Older versions:
Only for 2 channels:
`wolf4192/backup-chain/eze7134/Documents/papers/chrom_marks/`

### Requirements and set up
First you should check you have all the necessary...
#### R:  	at least R 3.3
	bioimagetools: (source code: https://r-forge.r-project.org/R/?group_id=1143 )
 
	 $ R
	    R version 3.3.3 (2017-03-06) -- "Another Canoe"
	    Copyright (C) 2017 The R Foundation for Statistical Computing
	    Platform: x86_64-pc-linux-gnu (64-bit)
	    [...]
	    > setRepositories(ind=c(1,2))
	    > install.packages("bioimagetools")
	    Installing package into ‘/usr/local/lib/R/site-library’
	    (as ‘lib’ is unspecified)
	    Warning in install.packages("bioimagetools") :
	      'lib = "/usr/local/lib/R/site-library"' is not writable
	    Would you like to use a personal library instead?  (y/n) y
	    Would you like to create a personal library
	    ~/R/x86_64-pc-linux-gnu-library/3.3
	    to install packages into?  (y/n) y
	    --- Please select a CRAN mirror for use in this session ---
	    [...]
	    Selection: 67
	    also installing the dependencies ‘jpeg’, ‘locfit’, ‘fftwtools’, ‘tiff’, ‘EBImage’
	    [...]
	    * DONE (bioimagetools)
	    [...]
	    > library('bioimagetools')
	    Bioimagetools 1.1.3

>**Troubleshooting:**
If R complains when installing bioimagetools it may lack system-wide libraries. These would have to be installed from terminal, probably with admin rights. eg:
libssl-dev: `sudo apt-get install libssl-dev`
libcurl4: `sudo apt-get install libcurl4`
libcurl4-openssl-dev: `sudo apt-get install libcurl4-openssl-dev`
libtiff5-dev: `sudo apt-get install libtiff5-dev`
EBImage:
`if (!requireNamespace("BiocManager", quietly = TRUE))`
    `install.packages("BiocManager")`
`BiocManager::install("EBImage", version = "3.8")`
Or for older versions:
`source("https://bioconductor.org/biocLite.R")`
`biocLite("EBImage")`
####Octave:	at least Octave 4.0.3
	Octave packages:
	Octave-image 2.6.1
	Octave-statistics 1.3.0

####Bioformats:
	$ wget https://downloads.openmicroscopy.org/bio-formats/5.9.2/artifacts/bioformats-octave-5.9.2.tar.gz
	$ octave --eval 'pkg install bioformats-octave-5.9.2.tar.gz'

####Java:
Need to also make a text file in your home: `~/javaclasspath.txt`
with this inside: `/usr/local/share/java/bioformats_package.jar`
To guide bioformats to the location of the bioformats_package.jar located on the machine.
If you are doing this on your own machine the PATH you write inside the javaclasspath file will change.

see: https://docs.openmicroscopy.org/bio-formats/5.8.1/users/octave/index.html


####GraphicsMagick-1.3.23
A library for working with images.

> ####Misc
There is a chance I've forgotten some basic library in this brief desciption. 
Basically if Cha*i*N complains at the start because it requires a library it should tell you what is missing.

## How to Cha*i*N:
Lemme show you. You'll need two types of *"raw"* data.
Files ending in `*SIR_THR_ALN.tif` & `ALN_MCNR-mask.tif`
These are not true raw `.tif` files. The first is the reconstructed raw (SIR) which has been processed by the script `16THR.ijm`. This is an custom ImageJ macro for Fiji which utilises the **SIMcheck** plugging from Micron to threshold the image to the modal value (assumed to be background) and converts the image to a 16 bit tiff.

In the process it also uses the true raw for **SIMcheck** to generate a modulation contrast map, the second file required.

Both files need to be aligned, currently using **Chromagnon**.
 
###Where is the raw data?
The raw data is far too large to keep under version control (each SIM result is around 150MB after conversion to 16bit). At the moment it is assumed that the files exist in the `/data` directory but we will come up with something better later. Maybe using our omero server to serve this files. To meet this assumption currently you just have to copy the data to be processed into the actual Cha*i*n `/data` directory.
The output will be in the Cha*i*n `/results` directory.

>Previously when Cha*i*N was in a different machine than the micron server, data was mounted as read only using `sshfs` from whichever directory it was in, into Cha*i*N's `/data` directory. 
>Can't use `sshfs` because 2125 does not have it installed
Can't used mount (below) because I need to be root.
>This mounting does not currently work:
Mount data as read-only
`$ mount -Br ~/path/to/data ~/Chain/2chan/data`
Mount results
`mount -B ~/path/to/results ~/Chain/2chan/results`


### How to structure and name raw data
The raw data to be processed by Cha*i*N can be organised in whatever architecture is required. The results data will be saved in an architecture mimicking that of the input.
However, for ease of compiling al results `*.csv` tables afterwards it is recommeded to organise the data under the following structure: `/data/cell-type/exp-condition/chrom-marker`. eg: `/data/C127/G1/CTCF` 

Make sure that the name of the directories is the same as the tags included in the file name:
`[experiment-code]_[microscope-slide-number]_[cell-type]_[experimental-condition]_[acquisition-oil]_[wavelength]_[chromatin-marker]_[acquisition-number]_[FILE-TYPE]`
i.e., for the previous example, a file would be: 
`2018-07-03_1_C127_Trip_G1_514_Scc1-488_CTCF-594_04_SIR_THR_ALN.tif`

###The Makefile
The numerous R and Octave scripts that form the basis of image analysis in this pipeline are managed under the **Make** build system. One can still call up any of these scripts individually to process a single file (or number of files a script requires) but this is not necessary.
In the Cha*i*N directory there is a file called `Makefile` with the instructions about what each script requires and what it should output. 
Furthermore the `Makefile` has instructions to look for the files appropriate `.tif` files (as mentioned previously) in the `/data` directory. It will mimick the subdirectory architecture when saving all output files to the `/results` directory.
On top of this, the relationship between the ouptut of some of the scripts as inputs for others downstream is coded in the makefile, such that **Make** can generate all required intermediate files if asked for an ouptut.
If running this analysis on a machine with enough memory and ram to handle paralleled image processing **Make** is also able to parallelize jobs, greatly reducing the time required for processing large numbers of images.
>	$ make command
or
	
>	$ make -j**n** command

>where **n** is the number of parallel jobs requred

Finally **Make** is able to recognize when files have been processed already do calling a function twice does not result in unnecessary processing of the same data that had previously been analysed. Only new data will be processed.


The `Makefile` includes a short guideline of the analyses available which can be called by:
	$ make help

Here are some examples of how to process files:
> Depending on the details of the Makefile the actual file title endings may be slighly differently processed

####Nucleus mask:
	$ make mask


**Input:** Reconstructed, thresholded and aligned data.
>the order of alignment and thresholding has changed due to the alignment software changes, thus files could be `*_SIR_THR_ALN.tif` or `*_SIR_EAL_THR.tif`. This **WILL** affect the ability of the Makefile to find these files so the recipes in the Makefile will have to be changed accordingly (see section below).

**script:** `nucleus_mask.m`

To make a `.tif` mask from the chromatin stain which will be the outline of the nucleus.
**Output:**	`*_mask.tif`

####Replication mask:
	$ make submask
	 

Will process all `EdU` files.
To make a `.tif`mask from the replication regions if required to separate nascent replication from non replicating regions.
Output:	`*EdU-submask.tif`

####Chromatin density segmentation:
	$ make segmented-chromatin
	

**Input:** 	`results/.../*_mask.tif` & `data/.../*_SIR_THR_ALN.tif` or `*_SIR_EAL_THR.tif`
**script:** `chromseg.R`
To make a `.tif`segmented chromatin landscape, binned into arbitrary classes inside the nuclear mask made from `make mask`
**Output:**	`*_SEG.tif`

####Centroids:
	$ make centroids
	

**Input:** 	`results/.../*_mask.tif` & `data/.../*_SIR_THR_ALN.tif` or `*_SIR_EAL_THR.tif`
**script:** `foci_centroids.m`
To make a `.csv` file of `xyz` coordinates of the intensity weighted centroids from IF signals against epigenetic markers.
**Output:**	`*_centroids.csv`

####Density ditributions:
	$ make distributions
	

**Input:** `results/.../*_centroids.csv` & `results/.../*_SEG.tif` & `data/.../*_MCNR-mask.tif`
**script:** `distribution.R`
To take the centroid coordinates from `make centroids`and index them on the segmented landscape from `make segmented-chromatin`, generating a `.csv` distribution of the number of centroids that fall in each density class.
**Output:**	`*_distn.csv`

####Chromatin-IC boundaries:
	$ make boundaries
	

**Input:** `results/.../*_SEG.tif`
**script:** `bound.m`
To take the the segmented landscape from `make segmented-chromatin` and generate a `.tif`binary of the boundary between class 1 and all other classes.
**Output:**	`*_bound.tif`

####Distances to boundaries:
	$ make d2b
	

**Input:** `results/.../*_centroids.csv` & `results/.../*_bound.tif`& `results/.../*_masktif`
**script:** `d2bound.m`
To take the centroid coordinates from `make centroids`and run a nearest neighbour eucledian algorithm to voxels in the segmented boundaires from `make boundaries`, generating a `.csv` distribution of the number of centroids that fall in each density class (*md2b*).
This is also repeated for a random sample of centroids (same number as those for the biological marker) to get a baseline for the random distribution (*rd2b*)
**Output:**	`*_md2b.csv` & `*_rd2b.csv`

>**Troubleshooting:**
The Make workflow will work automatically until this point.
That is to say that typing `make d2b` will result in all files being made which are required to make the d2b files. This includes masks, SEGs, boundaries etcs.
The next two scripts (d2bfit and d2bnorm) are also managed by the Make system. However typing `make d2bfit` or `make d2bnorm` without having all files required will result in an error saying that there is no recipe to make the required output files of these scripts. 
This is related to the dual outputs of the `d2b` script which are the dual inputs of `d2bfit` command.
Thus if the output of `d2bfit` or `d2bnorm` are required, it is adviced to do this in two steps.
1. `make d2b`
2. `make d2bfit` or `make d2bnorm`


####Distance fit to model
	$ make d2bfit
	

**Input:** `results/.../*_md2b.csv` & `results/.../*_rd2b.csv`
**script:** `d2bfit.R`
To use a least squares regression to fit a *mu* and *sigma* for the available distances of the marker and the random sample from `make d2b`.
**Output:**	`*_d2bfit.csv`

####Distances normalized
	$ make d2bnorm
	

**Input:** `results/.../*_d2bfit.csv`
**script:** `d2bnorm.R`
To make normalised distances `make d2b`as log2fold enrichment simulated over +/- 400nm from the *mu* and *sigma* from `make d2bfit`.
**Output:**	`*_d2bnorm.csv`

### Tweaking Makefile (when/if required)
>If you can help it, don't do it.

Currently need to copy all of G1 data from all the trip directories in bioc1090/data which also end in the "*SIR_THR_ALN.tif" and "*MCNR-mask.tif":

find LS2018-07-03_1_C127_Trip/ -type f \( -name "*G1*" -and \( -name "*SIR_THR_ALN.tif" -or -name "*MCNR-mask.tif" \) \) -exec cp {} ~/papers/data/C127/Trip \;

Then need to scp to a machine that can run Fiji to open all and save them as tifs again in imageJ with a OS.ijm macro.

Then scp to copy them back to data in Cha*i*n.

## Final data management and aggregation
After processing with Cha*i*N the `/results` directory (with any created sub directories) will be filled with the output files from the image analysis.

Then to aggregate the distribution and network data employ these scipts from the directory where the makefile is, in this order:

* Summary.R (**input:** `*_distn.csv`, **output:** `*_distn-summary.csv`& `*_distn-full.csv`)
	* spot-class-count.R (**input:** `*_distn.csv`, **output:** `*_spot-class-count.csv`& `*_spot-class-summaryl.csv`)
		* dist_compile.R (**input:**  `*_distn-summary.csv`& `*_spot-class-summary.csv`, **output:** `*_distn-compile.csv`)
			* SA2vol.R (**input:** `*_spot-class-summary.csv` & `*_bound.tif`, **output:** `*_SA2vol_raw.csv` & `*_SA2vol_stat.csv` & `*_chrom2vol_stat.csv`)

To aggregate d2b data:

* summaryd2b.R (**input:** `.csv`, **output:** `.csv`)
	* d2b_compile.R (**input:** `*_network-summary.csv` & `*_abs-md2b_summary.csv` & `*_d2bnorm-summary.csv` **output:** `*_md2b-compile.csv` & `*_d2bnorm-compile.csv` & `*_network-compile.csv`)

For changes in cell type and condition each script needs to be manually adjusted by changing those variables at the beginning of the script.

##Making plots: 
Below are the scripts used to generate the graphs in the paper. These scripts cannot be run from terminal as `Rscript [path-to-script]`.
They should be opened in R (or R studio) to asses them. They are all comprised of multiple semgents each outputting different kind of graphs and requiring input of different data.



####batch-distn-barplots.R
Required packages:

		ggplot2
		reshape2
		
* Automatically produce a barplot for each marker of its enrichment/depletion in each chromatin class ([condition]_distn-compile*.csv).
* Depending on the ggplot used and data uploaded the barplots can be for single condition or multiple conditions in parallel.
* Groups of markers can also be plotted in parallel for a single conditions (to group by function)
* Bar plots of WF vs SIR enrichment/depletion data for controls ("WF-distribution-summary.csv").


####heatmapLUT.R
Required packages:

		ggplot2
		scales
		RColorBrewer
* This will produce the heatmaps with the corrected loops (for extreme values)
* The variety of small ggplots to choose from include all markers for a single condition. or comparing across conditions in different panels etc.
* It requires the type of file with the density distributions for all markers: [condition]_distn-compile*.csv (eg: G1_distn-compile_volnormE_Final.csv)

####distance-plots.R
Required packages:

		ggplot2
		RColorBrewer
* This will produce the dot plot of distance from IC surface for all markers in a single condition or multiple conditions in parallel.
* It requires the type of file with the distances for each marker: [condition]_md2b-compile*. (eg: csv G1_md2b-compile_Final.csv)

####SAvolbox.R
Required packages:
	
		ggplot2
		RColorBrewer
		reshape2
		dplyr
* Box or violin plots of SA:vol ratios distribution (C127_[G1/NaN3/Hyper-os/Trip]_SA2vol_raw_Final.csv).
* Dot plot of mean SA:vol with 95CI (C127_[G1/NaN3/Hyper-os/Trip]_SA2vol_raw_Final.csv).
* Box or violin of chromatin vol:nuclear vol (C127_[G1/NaN3/Hyper-os/Trip]_SA2vol_raw_Final.csv).
* Box or violin of whole nuclear volume (C127_[G1/NaN3/Hyper-os/Trip]_SA2vol_raw_Final.csv).
* Changes of volumes for each density class as absolute values or normalised to nuclear volume for G1 or Sphase ([G1/ES/MS/LS]_distn-compile.csv"). 
* Remodeller size to marker distn (from remod-size-table_Final.csv).
* Bar plots of WF vs SIR chromatin data for controls ("WF-distribution-summary.csv").

####network-dims.R
Required packages:

		ggplot2
		plyr
* To plot changes in each chromatin network size as continuous line (G1_network-compile.csv)
* To work out the network radius from a tubular model and plot the radius as point with 95CI bars or the network diameter

###Useful tools:
General Bash commands for manipulating data in batch:

Use this to **MOVE** all files not tifs or dvs from directory "." to another directory:
`find . -type f ! -name "*.tif" ! -name "*.dv" -exec mv {} ~/Documents/Microscopy/2015/2015.08/EM15-08-B/ \;`

Use this to **COPY** all files in this directory and below which have “*histone mod*” and end in X or Y to somewhere else
`find . -type f \( -name "*H3K4me3*" -and \( -name "*EAL_THR.tif" -or -name "*MCNR-mask.tif" \) \) -exec cp {} /Volumes/wolf4192/papers/chrom_marks... \;`

`find . -type f \( -name "*EAL_THR.tif" -or -name "*_THR_mask.tif" \) -exec cp {} /Volumes/wolf4192/data/OLDSYS \;`
`find . -type f -name "*MCNR-mask.tif" -name "*SIR_EAL_THR.tif" -exec cp {} ~/data/ \;`

finds all files that finish with centroids.csv and are also empty and removes them:
`find . -size 0b -a -type f -name "*centroids.csv" -exec rm -f {} \;`

counts the number of files of a particular type
`find results/ -type f -name *d2b* | wc -l`

>**Troubleshooting:**
Make works with the use of `.dirstamp` files to alert it to the change in input data. This is done on purpose. As if the input data is updated then the result data is no longer relevant to the new dataset so the analysis is redone making all new files from the point of change and downstream files that depend on this.
For a change in a script this will be all the files depending on (and including) the output of the scripts. For a change in the input data this will be ALL files.
If the user called all `make` commands they thought necessary and then removed the data from the `/data` directory to then realise they wanted to add another `make` command they can put the relevant data back in the `/data` directory.
However, this new data will have a date later than the `.dirstamp` of the `/results` directory. Despite the fact that the data is the same, Make believes it is updated and will try to re-run the whole analysis.
In such scenario, to avoid this, find all `.dirstamp` files in the `/results` directory and all directories below and delete them:
`find results/ -type f -name ".dirstamp*" -exec rm {} +`
This will allow Makefile to match the results to the data by the filename only and update a new `.dirstamp` so all the analysis does not need to be repeated.








