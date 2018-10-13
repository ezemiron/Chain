#Cha*i*N - Cha*i*N analysis of the *in situ* Nucleome


##What?
A tool to dissect the spatial distribution of the mammalian nucleome relative to chromatin from super resolution 3D structured illumination microscopy (3D SIM) data.

##Where?
The original versions of Cha*i*N are currently in the Oxford biochemistry Micron servers under `wolf4192/Chain/`
And in github: git@github.com:ezemiron/Chain.git
https://github.com/ezemiron/Chain.git

One version exists for 2 Channel data with chromatin in second channel and spots channel in first
and another version for 3 channel data, with chromatin in third channel and spots channel in first and second and channel.

#####Older versions:
`~/backup-chain/eze7134`
or
`~/sandbox/C127`

Test output
`~/Chain-backup/results ...`

Can't use sshfs because 2125 does not have it installed
Can't used mount (below) because I need to be root.
Currently just have to copy the data into the actual Cha*i*n data folder

(this mounting does not currently work:
Mount data as read-only from micron 5 to bioch7134
mount -Br ~/Chain-backup/test2125/data ~/Chain/eze/Documents/papers/chrom_marks/data

Mount results
mount -B ~/Chain-backup/test2125/results ~/Chain/eze/Documents/papers/chrom_marks/results/
)

Currently need to copy all of G1 data from all the trip directories in bioc1090/data which also end in the "*SIR_THR_ALN.tif" and "*MCNR-mask.tif":

find LS2018-07-03_1_C127_Trip/ -type f \( -name "*G1*" -and \( -name "*SIR_THR_ALN.tif" -or -name "*MCNR-mask.tif" \) \) -exec cp {} ~/papers/data/C127/Trip \;

Then need to scp to a machine that can run Fiji to open all and save them as tifs again in imageJ with a OS.ijm macro.

Then scp to copy them back to data in Cha*i*n.


Do make distribution (make -j8 distribution)
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

## How to:
Lemme show you. But first you should check you have all the necessary...
###Requirements

####R:  	at least R 3.3
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
(in a terminal)
	wget https://downloads.openmicroscopy.org/bio-formats/5.9.2/artifacts/bioformats-octave-5.9.2.tar.gz
    	octave --eval 'pkg install bioformats-octave-5.9.2.tar.gz'


####Java:
Need to also make a text file in your home: `~/javaclasspath.txt`
with this inside: `/usr/local/share/java/bioformats_package.jar`
to guide bioformats to the location of the bioformats_package.jar located on the machine.
If you are doing this on your own machine the PATH you write inside the javaclasspath file will change.


####GraphicsMagick-1.3.23
A library for working with images.

####Misc
There is a chance I've forgotten some basic library in this brief desciption. Basically if Cha*i*N complains at the start because it requires a library it should tell you what is missing.


##Graphs: 
Below are the scripts used to generate the graphs in the paper. These scripts cannot be run from terminal as Rscript [path-to-script]. 
They should be opened in R (or R studio) to asses them. They are all comprised of multiple semgents each outputting different kind of graphs and requiring input of different data.



batch-distn-barplots.R
	Required packages:
		ggplot2
		reshape2
	Automatically produce a barplot for each marker of its enrichment/depletion in each chromatin class ([condition]_distn-compile*.csv).
	Depending on the ggplot used and data uploaded the barplots can be for single condition or multiple conditions in parallel.
	Groups of markers can also be plotted in parallel for a single conditions (to group by function)
	Bar plots of WF vs SIR enrichment/depletion data for controls ("WF-distribution-summary.csv").


heatmapLUT.R
	Required packages:
		ggplot2
		scales
		RColorBrewer
	This will produce the heatmaps with the corrected loops (for extreme values)
	The variety of small ggplots to choose from include all markers for a single condition. or comparing across conditions in different panels etc.
	It requires the type of file with the density distributions for all markers: [condition]_distn-compile*.csv (eg: G1_distn-compile_volnormE_Final.csv)

distance-plots.R
	Required packages:
		ggplot2
		RColorBrewer
	This will produce the dot plot of distance from IC surface for all markers in a single condition or multiple conditions in parallel.
	It requires the type of file with the distances for each marker: [condition]_md2b-compile*. (eg: csv G1_md2b-compile_Final.csv)

SAvolbox.R
	Required packages:
		ggplot2
		RColorBrewer
		reshape2
		dplyr
	Box or violin plots of SA:vol ratios distribution (C127_[G1/NaN3/Hyper-os/Trip]_SA2vol_raw_Final.csv).
	Dot plot of mean SA:vol with 95CI (C127_[G1/NaN3/Hyper-os/Trip]_SA2vol_raw_Final.csv).
	Box or violin of chromatin vol:nuclear vol (C127_[G1/NaN3/Hyper-os/Trip]_SA2vol_raw_Final.csv).
	Box or violin of whole nuclear volume (C127_[G1/NaN3/Hyper-os/Trip]_SA2vol_raw_Final.csv).
	Changes of volumes for each density class as absolute values or normalised to nuclear volume for G1 or Sphase ([G1/ES/MS/LS]_distn-compile.csv"). 
	Remodeller size to marker distn (from remod-size-table_Final.csv).
	Bar plots of WF vs SIR chromatin data for controls ("WF-distribution-summary.csv").

network-dims.R
	Required packages:
		ggplot2
		plyr
	To plot changes in each chromatin network size as continuous line (G1_network-compile.csv)
	To work out the network radius from a tubular model and plot the radius as point with 95CI bars or the network diameter










