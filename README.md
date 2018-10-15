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

####Nucleus mask
	$ make mask


Will process all files...
####Replication mask
	$ make submask
	

Will process all files...
####Chromatin density segmentation
	$ make segmented-chromatin
	

Will process all files...
####Centroids
	$ make centroids
	

Will process all files...
####Density ditributions
	$ make distributions
	

Will process all files...
####Chromatin-IC boundaries
	$ make boundaries
	

Will process all files...
####Distanes to 
	$ make distributions
	

Will process all files...
####Density ditributions
	$ make distributions
	

Will process all files...


### Tweaking Makefile when required
>If you can help it, don't do it.

Currently need to copy all of G1 data from all the trip directories in bioc1090/data which also end in the "*SIR_THR_ALN.tif" and "*MCNR-mask.tif":

find LS2018-07-03_1_C127_Trip/ -type f \( -name "*G1*" -and \( -name "*SIR_THR_ALN.tif" -or -name "*MCNR-mask.tif" \) \) -exec cp {} ~/papers/data/C127/Trip \;

Then need to scp to a machine that can run Fiji to open all and save them as tifs again in imageJ with a OS.ijm macro.

Then scp to copy them back to data in Cha*i*n.

### Final data management
After processing with Cha*i*N the `/results` directory (with any created sub directories) will be filled with the output files from the image analysis.

Do `make distribution` (or `make -j8 distribution` for parallelization)
make d2b
make d2bfit
make d2bnorm

Then to aggregate the distribution and network data:
Summary.R
	spot-class-count.R
		dist_compile.R
		SA2vol.R

To aggregate d2b data:
summaryd2b.R
	d2b_compile.R


###Making plots: 
Below are the scripts used to generate the graphs in the paper. These scripts cannot be run from terminal as Rscript [path-to-script]. 
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










