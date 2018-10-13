## Copyright (C) 2015 David Pinto <david.pinto@bioch.ox.ac.uk>
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

## FORCE_MAKE is used as dependency that is never built to keep a target
## always out of date.  This is used to leave the finding of dependencies
## and deciding for out of date up to latexmk (see latexmk man page)

.PHONY: publication clean help FORCE_MAKE \
        mask centroids segmented_chromatin distribution boundaries d2b

help:
	@echo "Targets:"
	@echo "   publication         - Create PDF of paper for publication"
	@echo
	@echo "   centroids           - CSV files with centroid coordinates"
	@echo "   mask                - Binary masks for nucleus (DAPI or Sytox)"
	@echo "   submask             - Binary masks for replication foci"
	@echo "   segmented_chromatin - Multi-level segmentation of chromatin."
	@echo "   distribution        - CSV file with the total number of centroids found in each class of the segmented nucleus"
	@echo "   boundaries          - TIF of 3D boundaries of specified chromatin class, default class 1"
	@echo "   d2b                 - CSV files with distance from centroids to the surface of interchromatin (IC, class 1) md2b, and distribution of random set of distances for normalisation rd2b."
	@echo "   d2bfit              - CSV file with fitted mu and sigma for the marker and a random sample of available distances."
	@echo "   d2bnorm             - CSV file with normalised distances over random from centroids to the surface of interchromatin and as log2fold enrichment."
	@echo "   clean               - Remove paper and all intermediary files (results)."
	@echo "   NNA                 - Euclidean nearest neighbor analysis on two centroid CSV files."
	@echo "   centroid_volume     - Determine the volume of centroids based on the MCNR-mask
	@echo "   difmask	      - Difference between nuclear and EdU-mask
	@echo
	@echo "Environment variables:"
	@echo "   OCTAVE  - octave interpreter program (default: octave)."
	@echo "   RSCRIPT - program to call R scripts (default: Rscript)."

##
## Environment variables and target alias.
##

OCTAVE ?= octave
RSCRIPT ?= Rscript

publication: sim_a.pdf

clean:
	$(RM) -r results sim_a.pdf dirstamp-file
	ls -1 sim_a* | grep -v tex | xargs $(RM)


##
## The real rules to build the project.
##

%.pdf: %.tex FORCE_MAKE
	latexmk -pdf -silent $<

## Our base images (by now they have already been reconstructed, aligned,
## and thresholded to 16 bit integers).
THR_FILES := $(shell find data/ -type f -name '*ALN_THR.tif')

## The images in this files were acquired incorrectly.  They will have
## chromatin staining on the first channel instead of the second.
THR_CHROMATIN_FIRST_FILES := $(shell find data/ -type f -name 'EM15-10-B*_THR.tif')

##Distinguish between images of 2 or 3 (EdU) channels:
##EDU_ONLY_FILES := $(shell find data/ -type f -name 'PO*EdU*_THR.tif')
#NNA_ONLY_FILES := $(shell find data/ -type f -name '*Cloud*_THR_MCF_ALN.tif')
HCT_ONLY_FILES := $(shell find data/ -type f -name '*G1*THR_ALN.tif')
##4C_EDU_ONLY := $(shell find data/ -type f -name '*4-color-EdU_*THR_ALN.tif')
##THR_EDU_FILES := $(filter-out $(HCT_ONLY_FILES), $(THR_CHROMATIN_FIRST_FILES), $(THR_FILES))
##THR_ONLY_FILES := $(filter-out $(HCT_ONLY_FILES), $(EDU_ONLY_FILES), $(THR_FILES))


## Makefile oddities... We do this so we can use them int pattern rules
## because we can't escape stuff there http://stackoverflow.com/a/19324181/1609556
UNDERSCORE = _


mask: $(patsubst data/%.tif, results/%$(UNDERSCORE)mask.tif, $(THR_ONLY_FILES)) $(patsubst data/%.tif, results/%$(UNDERSCORE)EdU-mask.tif, $(EDU_ONLY_FILES)) $(patsubst data/%.tif, results/%$(UNDERSCORE)mask.tif, $(HCT_ONLY_FILES)) $(patsubst data/%.tif, results/%$(UNDERSCORE)4color-mask.tif, $(4C_EDU_ONLY)) $(patsubst data/%.tif, results/%$(UNDERSCORE)NNA-mask.tif, $(NNA_ONLY_FILES))
## $(patsubst data/%.tif, results/%$(UNDERSCORE)FISH-mask.tif, $(RFISH_ONLY_FILES)) $(patsubst data/%.tif, results/%$(UNDERSCORE)mask.tif, $(HCT_ONLY_FILES))
submask: $(patsubst data/%.tif, results/%$(UNDERSCORE)EdU-submask.tif, $(EDU_ONLY_FILES))
difmask: $(patsubst data/%.tif, results/%$(UNDERSCORE)EdU-difmask.tif, $(EDU_ONLY_FILES)) 
## $(patsubst data/%.tif, results/%$(UNDERSCORE)FISH-submask.tif, $(RFISH_ONLY_FILES))
centroids: $(patsubst data/%.tif, results/%$(UNDERSCORE)centroids.csv, $(THR_ONLY_FILES)) $(patsubst data/%.tif, results/%$(UNDERSCORE)Rep-centroids.csv, $(EDU_ONLY_FILES)) $(patsubst data/%.tif, results/%$(UNDERSCORE)C594-centroids.csv, $(HCT_ONLY_FILES)) $(patsubst data/%.tif, results/%$(UNDERSCORE)C488-centroids.csv, $(HCT_ONLY_FILES)) $(patsubst data/%.tif, results/%$(UNDERSCORE)C1-centroids.csv, $(4C_EDU_ONLY)) $(patsubst data/%.tif, results/%$(UNDERSCORE)C2-centroids.csv, $(4C_EDU_ONLY)) $(patsubst data/%.tif, results/%$(UNDERSCORE)C3-centroids.csv, $(4C_EDU_ONLY)) $(patsubst data/%.tif, results/%$(UNDERSCORE)C4-centroids.csv, $(4C_EDU_ONLY)) $(patsubst data/%.tif, results/%$(UNDERSCORE)C1-centroids.csv, $(NNA_ONLY_FILES)) $(patsubst data/%.tif, results/%$(UNDERSCORE)C2-centroids.csv, $(NNA_ONLY_FILES))
## $(patsubst data/%.tif, results/%$(UNDERSCORE)mask-centroids.csv, $(EDU_ONLY_FILES)) $(patsubst data/%.tif, results/%$(UNDERSCORE)EdU-centroids.csv, $(EDU_ONLY_FILES)) 
## $(patsubst data/%.tif, results/%$(UNDERSCORE)mask-centroids.csv, $(EDU_ONLY_FILES))
## $(patsubst data/%.tif, results/%$(UNDERSCORE)C1-centroids.csv, $(NNA_ONLY_FILES)) $(patsubst data/%.tif, results/%$(UNDERSCORE)C2-centroids.csv, $(NNA_ONLY_FILES)) $(patsubst data/%.tif, results/%$(UNDERSCORE)FISH-centroids.csv, $(RFISH_ONLY_FILES)) $(patsubst data/%.tif, results/%$(UNDERSCORE)C594-centroids.csv, $(HCT_ONLY_FILES)) $(patsubst data/%.tif, results/%$(UNDERSCORE)C488-centroids.csv, $(HCT_ONLY_FILES))
centroid_volume: $(patsubst data/%.tif, results/%$(UNDERSCORE)centroidsvol.csv, $(THR_ONLY_FILES)) $(patsubst data/%.tif, results/%$(UNDERSCORE)mask-centroidsMCN.csv, $(EDU_ONLY_FILES)) $(patsubst data/%.tif, results/%$(UNDERSCORE)C488-centroidsMCN.csv, $(HCT_ONLY_FILES))
centroid_MCN: $(patsubst data/%.tif, results/%$(UNDERSCORE)dif-centroidsMCN.csv, $(EDU_ONLY_FILES)) $(patsubst data/%.tif, results/%$(UNDERSCORE)EdU-centroidsMCN.csv, $(EDU_ONLY_FILES))
## $(patsubst data/%.tif, results/%$(UNDERSCORE)FISH-centroidsvol.csv, $(RFISH_ONLY_FILE)) $(patsubst data/%.tif, results/%$(UNDERSCORE)EdU-centroidsvol.csv, $(EDU_ONLY_FILES))
segmented_chromatin: $(patsubst data/%.tif, results/%$(UNDERSCORE)SEG.tif, $(THR_ONLY_FILES)) $(patsubst data/%.tif, results/%$(UNDERSCORE)EdU-SEG.tif, $(EDU_ONLY_FILES)) $(patsubst data/%.tif, results/%$(UNDERSCORE)FISH-SEG.tif, $(RFISH_ONLY_FILES)) $(patsubst data/%.tif, results/%$(UNDERSCORE)SEG.tif, $(HCT_ONLY_FILES))
boundaries: $(patsubst data/%.tif, results/%$(UNDERSCORE)bound.tif, $(THR_ONLY_FILES)) $(patsubst data/%.tif, results/%$(UNDERSCORE)mask-bound.tif, $(EDU_ONLY_FILES)) $(patsubst data/%.tif, results/%$(UNDERSCORE)bound.tif, $(HCT_ONLY_FILES))
#$(patsubst data/%.tif, results/%$(UNDERSCORE)EdU-bound.tif, $(EDU_ONLY_FILES))
distribution: $(patsubst data/%.tif, results/%$(UNDERSCORE)Rep-distn.csv, $(EDU_ONLY_FILES)) $(patsubst data/%.tif, results/%$(UNDERSCORE)EdU-distn.csv, $(EDU_ONLY_FILES)) $(patsubst data/%.tif, results/%$(UNDERSCORE)dif-distn.csv, $(EDU_ONLY_FILES)) $(patsubst data/%.tif, results/%$(UNDERSCORE)C594-distn.csv, $(HCT_ONLY_FILES)) $(patsubst data/%.tif, results/%$(UNDERSCORE)C488-distn.csv, $(HCT_ONLY_FILES))
#$(patsubst data/%.tif, results/%$(UNDERSCORE)distn.csv, $(THR_ONLY_FILES))
d2b: $(patsubst data/%.tif, results/%$(UNDERSCORE)MCN_md2b.csv, $(THR_ONLY_FILES)) $(patsubst data/%.tif, results/%$(UNDERSCORE)C594-md2b.csv, $(HCT_ONLY_FILES)) $(patsubst data/%.tif, results/%$(UNDERSCORE)C488-md2b.csv, $(HCT_ONLY_FILES)) $(patsubst data/%.tif, results/%$(UNDERSCORE)mask-md2b.csv, $(EDU_ONLY_FILES)) $(patsubst data/%.tif, results/%$(UNDERSCORE)dif-md2b.csv, $(EDU_ONLY_FILES)) $(patsubst data/%.tif, results/%$(UNDERSCORE)dif-md2b.csv, $(EDU_ONLY_FILES)) 
#$(patsubst data/%.tif, results/%$(UNDERSCORE)EdU-md2b.csv, $(EDU_ONLY_FILES))
#$(patsubst data/%.tif, results/%$(UNDERSCORE)rd2b.csv, $(THR_ONLY_FILES)) $(patsubst data/%.tif, results/%$(UNDERSCORE)EdU-rd2b.csv, $(EDU_ONLY_FILES))
#$(patsubst data/%.tif, results/%$(UNDERSCORE)EdU-md2b.csv, $(EDU_ONLY_FILES)) $(patsubst data/%.tif, results/%$(UNDERSCORE)EdU-rd2b.csv, $(EDU_ONLY_FILES))
d2bfit: $(patsubst data/%.tif, results/%$(UNDERSCORE)MCN_d2bfit.csv, $(THR_ONLY_FILES)) $(patsubst data/%.tif, results/%$(UNDERSCORE)mask-d2bfit.csv, $(EDU_ONLY_FILES)) $(patsubst data/%.tif, results/%$(UNDERSCORE)C594-d2bfit.csv, $(HCT_ONLY_FILES)) $(patsubst data/%.tif, results/%$(UNDERSCORE)C488-d2bfit.csv, $(HCT_ONLY_FILES))
d2bnorm: $(patsubst data/%.tif, results/%$(UNDERSCORE)MCN_d2bnorm.csv, $(THR_ONLY_FILES)) $(patsubst data/%.tif, results/%$(UNDERSCORE)mask-d2bnorm.csv, $(EDU_ONLY_FILES)) $(patsubst data/%.tif, results/%$(UNDERSCORE)C594-d2bnorm.csv, $(HCT_ONLY_FILES)) $(patsubst data/%.tif, results/%$(UNDERSCORE)C488-d2bnorm.csv, $(HCT_ONLY_FILES))
NNA: $(patsubst data/%.tif, results/%$(UNDERSCORE)EdUtoMark_NNA.csv, $(EDU_ONLY_FILES)) $(patsubst data/%.tif, results/%$(UNDERSCORE)C594-NNA.csv, $(HCT_ONLY_FILES)) $(patsubst data/%.tif, results/%$(UNDERSCORE)C1-C2_NNA.csv, $(NNA_ONLY_FILES)) 

## FIXME It seems that the dapi/sytox channel is not always the same.  And
##      we can't guess it from filename, or wavelength (e.g., G1/H3K36me2
##      has Sytox on first channel by accident).
##      Figure out a way to use the same rules below for THR_FILES and
##      THR_CHROMATIN_FIRST files.  Set the channel number argument to a
##      variable and delay the expansion?

.SECONDEXPANSION:
##results/%$(UNDERSCORE)mask.tif: data/%.tif | $$(dir $$@).dirstamp
##	$(OCTAVE) scripts/nucleus_mask.m 2 $^ $@

##results/%$(UNDERSCORE)centroids.csv: results/%$(UNDERSCORE)mask.tif data/%.tif
##	$(OCTAVE) scripts/foci_centroids2.m 1 $^ 2 > $@ 

##results/%$(UNDERSCORE)SIR$(UNDERSCORE)ALN$(UNDERSCORE)THR$(UNDERSCORE)centroidsvol.csv: results/%$(UNDERSCORE)SIR$(UNDERSCORE)ALN$(UNDERSCORE)THR$(UNDERSCORE)centroids.csv data/%$(UNDERSCORE)ALN$(UNDERSCORE)MCNR-mask.tif results/%$(UNDERSCORE)SIR$(UNDERSCORE)ALN$(UNDERSCORE)THR$(UNDERSCORE)SEG.tif
##	$(RSCRIPT) scripts/centroid_volume.R $^ > $@

##results/%$(UNDERSCORE)SEG.tif: results/%$(UNDERSCORE)mask.tif data/%.tif
##	$(RSCRIPT) scripts/chromseg.R 2 $^ $@

##results/%$(UNDERSCORE)SIR$(UNDERSCORE)ALN$(UNDERSCORE)THR$(UNDERSCORE)distn.csv: results/%$(UNDERSCORE)SIR$(UNDERSCORE)ALN$(UNDERSCORE)THR$(UNDERSCORE)centroids.csv results/%$(UNDERSCORE)SIR$(UNDERSCORE)ALN$(UNDERSCORE)THR$(UNDERSCORE)SEG.tif data/%$(UNDERSCORE)ALN$(UNDERSCORE)MCNR-mask.tif
##	$(RSCRIPT) scripts/distribution.R $^ > $@

##results/%$(UNDERSCORE)SIR$(UNDERSCORE)ALN$(UNDERSCORE)THR$(UNDERSCORE)bound.tif: results/%$(UNDERSCORE)SIR$(UNDERSCORE)ALN$(UNDERSCORE)THR$(UNDERSCORE)SEG.tif
##	$(OCTAVE) scripts/bound.m 1 $^ $@

##results/%$(UNDERSCORE)SIR$(UNDERSCORE)ALN$(UNDERSCORE)THR$(UNDERSCORE)MCN_md2b.csv: results/%$(UNDERSCORE)SIR$(UNDERSCORE)ALN$(UNDERSCORE)THR$(UNDERSCORE)SEG.tif results/%$(UNDERSCORE)SIR$(UNDERSCORE)ALN$(UNDERSCORE)THR$(UNDERSCORE)centroidsvol.csv results/%$(UNDERSCORE)SIR$(UNDERSCORE)ALN$(UNDERSCORE)THR$(UNDERSCORE)bound.tif results/%$(UNDERSCORE)SIR$(UNDERSCORE)ALN$(UNDERSCORE)THR$(UNDERSCORE)mask.tif
##	$(OCTAVE) scripts/d2bound.m 1 $^ results/$*$(UNDERSCORE)SIR$(UNDERSCORE)ALN$(UNDERSCORE)THR$(UNDERSCORE)MCN_md2b.csv results/$*$(UNDERSCORE)SIR$(UNDERSCORE)ALN$(UNDERSCORE)THR$(UNDERSCORE)MCN_rd2b.csv

##results/%$(UNDERSCORE)SIR$(UNDERSCORE)ALN$(UNDERSCORE)THR$(UNDERSCORE)rd2b.csv: results/%$(UNDERSCORE)SIR$(UNDERSCORE)ALN$(UNDERSCORE)THR$(UNDERSCORE)md2b.csv
#fake recipe, does nothing because rd2b is already made as the final argument and target of md2b

##results/%$(UNDERSCORE)SIR$(UNDERSCORE)ALN$(UNDERSCORE)THR$(UNDERSCORE)MCN_d2bfit.csv: results/%$(UNDERSCORE)SIR$(UNDERSCORE)ALN$(UNDERSCORE)THR$(UNDERSCORE)MCN_md2b.csv results/%$(UNDERSCORE)SIR$(UNDERSCORE)ALN$(UNDERSCORE)THR$(UNDERSCORE)MCN_rd2b.csv
##	$(RSCRIPT) scripts/d2bfit.R $^ > $@

##results/%$(UNDERSCORE)SIR$(UNDERSCORE)ALN$(UNDERSCORE)THR$(UNDERSCORE)MCN_d2bnorm.csv: results/%$(UNDERSCORE)SIR$(UNDERSCORE)ALN$(UNDERSCORE)THR$(UNDERSCORE)MCN_d2bfit.csv
##	$(RSCRIPT) scripts/d2bnorm.R $^ > $@

## EdU-color
##results/%$(UNDERSCORE)4color-mask.tif: data/%.tif | $$(dir $$@).dirstamp
##	$(OCTAVE) scripts/nucleus_mask_NNA.m 1 $^ $@

##results/%$(UNDERSCORE)C1-centroids.csv: results/%$(UNDERSCORE)4color-mask.tif data/%.tif
##	$(OCTAVE) scripts/foci_centroids2.m 1 $^ 4 > $@

##results/%$(UNDERSCORE)C2-centroids.csv: results/%$(UNDERSCORE)4color-mask.tif data/%.tif
##	$(OCTAVE) scripts/foci_centroids2.m 2 $^ 4 > $@

##results/%$(UNDERSCORE)C3-centroids.csv: results/%$(UNDERSCORE)4color-mask.tif data/%.tif
##	$(OCTAVE) scripts/foci_centroids2.m 3 $^ 4 > $@

##results/%$(UNDERSCORE)C4-centroids.csv: results/%$(UNDERSCORE)4color-mask.tif data/%.tif
##	$(OCTAVE) scripts/foci_centroids2.m 4 $^ 4 > $@


#replication data recipes:
##results/%$(UNDERSCORE)EdU-mask.tif: data/%.tif | $$(dir $$@).dirstamp
##	$(OCTAVE) scripts/edu_mask2.m 3 $^ $@

##results/%$(UNDERSCORE)EdU-submask.tif: data/%.tif | $$(dir $$@).dirstamp
##	$(OCTAVE) scripts/edu_mask.m 1 $^ $@

##results/%$(UNDERSCORE)SIR$(UNDERSCORE)ALN$(UNDERSCORE)THR$(UNDERSCORE)EdU-difmask.tif: results/%$(UNDERSCORE)SIR$(UNDERSCORE)ALN$(UNDERSCORE)THR$(UNDERSCORE)EdU-mask.tif results/%$(UNDERSCORE)SIR$(UNDERSCORE)ALN$(UNDERSCORE)THR$(UNDERSCORE)EdU-submask.tif
##	$(OCTAVE) scripts/DAPI_EdU.m $^ $@

##results/%$(UNDERSCORE)Rep-centroids.csv: results/%$(UNDERSCORE)mask.tif data/%.tif
##	$(OCTAVE) scripts/foci_centroids2.m 1 $^ 2 > $@ 

##results/%$(UNDERSCORE)EdU-SEG.tif: results/%$(UNDERSCORE)EdU-mask.tif data/%.tif
##	$(RSCRIPT) scripts/chromseg.R 3 $^ $@

##results/%$(UNDERSCORE)SIR$(UNDERSCORE)ALN$(UNDERSCORE)THR$(UNDERSCORE)mask-distn.csv: results/%$(UNDERSCORE)SIR$(UNDERSCORE)ALN$(UNDERSCORE)THR$(UNDERSCORE)mask-centroidsMCN.csv results/%$(UNDERSCORE)SIR$(UNDERSCORE)ALN$(UNDERSCORE)THR$(UNDERSCORE)EdU-SEG.tif data/%$(UNDERSCORE)ALN$(UNDERSCORE)MCNR-mask.tif
##	$(RSCRIPT) scripts/edu_distribution_mask.R $^ > $@

##results/%$(UNDERSCORE)SIR$(UNDERSCORE)ALN$(UNDERSCORE)THR$(UNDERSCORE)dif-distn.csv: results/%$(UNDERSCORE)SIR$(UNDERSCORE)ALN$(UNDERSCORE)THR$(UNDERSCORE)dif-centroidsMCN.csv results/%$(UNDERSCORE)SIR$(UNDERSCORE)ALN$(UNDERSCORE)THR$(UNDERSCORE)EdU-SEG.tif data/%$(UNDERSCORE)ALN$(UNDERSCORE)MCNR-mask.tif results/%$(UNDERSCORE)SIR$(UNDERSCORE)ALN$(UNDERSCORE)THR$(UNDERSCORE)EdU-difmask.tif
##	$(RSCRIPT) scripts/edu_distribution.R $^ > $@

##results/%$(UNDERSCORE)SIR$(UNDERSCORE)ALN$(UNDERSCORE)THR$(UNDERSCORE)EdU-distn.csv: results/%$(UNDERSCORE)SIR$(UNDERSCORE)ALN$(UNDERSCORE)THR$(UNDERSCORE)EdU-centroidsMCN.csv results/%$(UNDERSCORE)SIR$(UNDERSCORE)ALN$(UNDERSCORE)THR$(UNDERSCORE)EdU-SEG.tif data/%$(UNDERSCORE)ALN$(UNDERSCORE)MCNR-mask.tif results/%$(UNDERSCORE)SIR$(UNDERSCORE)ALN$(UNDERSCORE)THR$(UNDERSCORE)EdU-submask.tif
##	$(RSCRIPT) scripts/edu_distribution.R $^ > $@

##results/%$(UNDERSCORE)SIR$(UNDERSCORE)ALN$(UNDERSCORE)THR$(UNDERSCORE)Rep-distn.csv: results/%$(UNDERSCORE)SIR$(UNDERSCORE)ALN$(UNDERSCORE)THR$(UNDERSCORE)Rep-centroids.csv results/%$(UNDERSCORE)SIR$(UNDERSCORE)ALN$(UNDERSCORE)THR$(UNDERSCORE)SEG.tif results/%$(UNDERSCORE)SIR$(UNDERSCORE)ALN$(UNDERSCORE)THR$(UNDERSCORE)mask.tif
##	$(RSCRIPT) scripts/edu_distribution_Rep.R $^ > $@

##results/%$(UNDERSCORE)SIR$(UNDERSCORE)ALN$(UNDERSCORE)THR$(UNDERSCORE)mask-centroidsMCN.csv: results/%$(UNDERSCORE)SIR$(UNDERSCORE)ALN$(UNDERSCORE)THR$(UNDERSCORE)mask-centroids.csv data/%$(UNDERSCORE)ALN$(UNDERSCORE)MCNR-mask.tif results/%$(UNDERSCORE)SIR$(UNDERSCORE)ALN$(UNDERSCORE)THR$(UNDERSCORE)EdU-SEG.tif
##	$(RSCRIPT) scripts/centroid_volume_EdU.R $^ > $@

##results/%$(UNDERSCORE)SIR$(UNDERSCORE)ALN$(UNDERSCORE)THR$(UNDERSCORE)EdU-centroidsMCN.csv: results/%$(UNDERSCORE)SIR$(UNDERSCORE)ALN$(UNDERSCORE)THR$(UNDERSCORE)mask-centroidsMCN.csv results/%$(UNDERSCORE)SIR$(UNDERSCORE)ALN$(UNDERSCORE)THR$(UNDERSCORE)EdU-SEG.tif results/%$(UNDERSCORE)SIR$(UNDERSCORE)ALN$(UNDERSCORE)THR$(UNDERSCORE)EdU-submask.tif
##	$(RSCRIPT) scripts/MCNcentroids.R $^ > $@

##results/%$(UNDERSCORE)SIR$(UNDERSCORE)ALN$(UNDERSCORE)THR$(UNDERSCORE)dif-centroidsMCN.csv: results/%$(UNDERSCORE)SIR$(UNDERSCORE)ALN$(UNDERSCORE)THR$(UNDERSCORE)mask-centroidsMCN.csv results/%$(UNDERSCORE)SIR$(UNDERSCORE)ALN$(UNDERSCORE)THR$(UNDERSCORE)EdU-SEG.tif results/%$(UNDERSCORE)SIR$(UNDERSCORE)ALN$(UNDERSCORE)THR$(UNDERSCORE)EdU-difmask.tif
##	$(RSCRIPT) scripts/MCNcentroids.R $^ > $@

##results/%$(UNDERSCORE)SIR$(UNDERSCORE)ALN$(UNDERSCORE)THR$(UNDERSCORE)mask-bound.tif: results/%$(UNDERSCORE)SIR$(UNDERSCORE)ALN$(UNDERSCORE)THR$(UNDERSCORE)EdU-SEG.tif
##	$(OCTAVE) scripts/bound.m 1 $^ $@

##results/%$(UNDERSCORE)SIR$(UNDERSCORE)ALN$(UNDERSCORE)THR$(UNDERSCORE)mask-md2b.csv: results/%$(UNDERSCORE)SIR$(UNDERSCORE)ALN$(UNDERSCORE)THR$(UNDERSCORE)EdU-SEG.tif results/%$(UNDERSCORE)SIR$(UNDERSCORE)ALN$(UNDERSCORE)THR$(UNDERSCORE)mask-centroidsMCN.csv results/%$(UNDERSCORE)SIR$(UNDERSCORE)ALN$(UNDERSCORE)THR$(UNDERSCORE)mask-bound.tif results/%$(UNDERSCORE)SIR$(UNDERSCORE)ALN$(UNDERSCORE)THR$(UNDERSCORE)EdU-mask.tif
##	$(OCTAVE) scripts/d2bound.m 1 $^ results/$*$(UNDERSCORE)SIR$(UNDERSCORE)ALN$(UNDERSCORE)THR$(UNDERSCORE)mask-md2b.csv results/$*$(UNDERSCORE)SIR$(UNDERSCORE)ALN$(UNDERSCORE)THR$(UNDERSCORE)mask-rd2b.csv

##results/%$(UNDERSCORE)SIR$(UNDERSCORE)ALN$(UNDERSCORE)THR$(UNDERSCORE)dif-md2b.csv: results/%$(UNDERSCORE)SIR$(UNDERSCORE)ALN$(UNDERSCORE)THR$(UNDERSCORE)EdU-SEG.tif results/%$(UNDERSCORE)SIR$(UNDERSCORE)ALN$(UNDERSCORE)THR$(UNDERSCORE)dif-centroidsMCN.csv results/%$(UNDERSCORE)SIR$(UNDERSCORE)ALN$(UNDERSCORE)THR$(UNDERSCORE)mask-bound.tif results/%$(UNDERSCORE)SIR$(UNDERSCORE)ALN$(UNDERSCORE)THR$(UNDERSCORE)EdU-difmask.tif
##	$(OCTAVE) scripts/d2bound.m 1 $^ results/$*$(UNDERSCORE)SIR$(UNDERSCORE)ALN$(UNDERSCORE)THR$(UNDERSCORE)dif-md2b.csv results/$*$(UNDERSCORE)SIR$(UNDERSCORE)ALN$(UNDERSCORE)THR$(UNDERSCORE)dif-rd2b.csv

##results/%$(UNDERSCORE)SIR$(UNDERSCORE)ALN$(UNDERSCORE)THR$(UNDERSCORE)EdU-md2b.csv: results/%$(UNDERSCORE)SIR$(UNDERSCORE)ALN$(UNDERSCORE)THR$(UNDERSCORE)EdU-SEG.tif results/%$(UNDERSCORE)SIR$(UNDERSCORE)ALN$(UNDERSCORE)THR$(UNDERSCORE)EdU-centroidsMCN.csv results/%$(UNDERSCORE)SIR$(UNDERSCORE)ALN$(UNDERSCORE)THR$(UNDERSCORE)mask-bound.tif results/%$(UNDERSCORE)SIR$(UNDERSCORE)ALN$(UNDERSCORE)THR$(UNDERSCORE)EdU-submask.tif
##	$(OCTAVE) scripts/d2bound.m 1 $^ results/$*$(UNDERSCORE)SIR$(UNDERSCORE)ALN$(UNDERSCORE)THR$(UNDERSCORE)EdU-md2b.csv results/$*$(UNDERSCORE)SIR$(UNDERSCORE)ALN$(UNDERSCORE)THR$(UNDERSCORE)EdU-rd2b.csv

##results/%$(UNDERSCORE)SIR$(UNDERSCORE)ALN$(UNDERSCORE)THR$(UNDERSCORE)EdU-d2bfit.csv: results/%$(UNDERSCORE)SIR$(UNDERSCORE)ALN$(UNDERSCORE)THR$(UNDERSCORE)EdU-md2b.csv results/%$(UNDERSCORE)SIR$(UNDERSCORE)ALN$(UNDERSCORE)THR$(UNDERSCORE)EdU-rd2b.csv
##	$(RSCRIPT) scripts/d2bfit.R $^ > $@

##results/%$(UNDERSCORE)SIR$(UNDERSCORE)ALN$(UNDERSCORE)THR$(UNDERSCORE)dif-d2bfit.csv: results/%$(UNDERSCORE)SIR$(UNDERSCORE)ALN$(UNDERSCORE)THR$(UNDERSCORE)dif-md2b.csv results/%$(UNDERSCORE)SIR$(UNDERSCORE)ALN$(UNDERSCORE)THR$(UNDERSCORE)dif-rd2b.csv
##	$(RSCRIPT) scripts/d2bfit.R $^ > $@

##results/%$(UNDERSCORE)SIR$(UNDERSCORE)ALN$(UNDERSCORE)THR$(UNDERSCORE)mask-d2bfit.csv: results/%$(UNDERSCORE)SIR$(UNDERSCORE)ALN$(UNDERSCORE)THR$(UNDERSCORE)mask-md2b.csv results/%$(UNDERSCORE)SIR$(UNDERSCORE)ALN$(UNDERSCORE)THR$(UNDERSCORE)mask-rd2b.csv
##	$(RSCRIPT) scripts/d2bfit.R $^ > $@

##results/%$(UNDERSCORE)SIR$(UNDERSCORE)ALN$(UNDERSCORE)THR$(UNDERSCORE)EdU-d2bnorm.csv: results/%$(UNDERSCORE)SIR$(UNDERSCORE)ALN$(UNDERSCORE)THR$(UNDERSCORE)EdU-d2bfit.csv
##	$(RSCRIPT) scripts/d2bnorm.R $^ > $@

##results/%$(UNDERSCORE)SIR$(UNDERSCORE)ALN$(UNDERSCORE)THR$(UNDERSCORE)dif-d2bnorm.csv: results/%$(UNDERSCORE)SIR$(UNDERSCORE)ALN$(UNDERSCORE)THR$(UNDERSCORE)dif-d2bfit.csv
##	$(RSCRIPT) scripts/d2bnorm.R $^ > $@

##results/%$(UNDERSCORE)SIR$(UNDERSCORE)ALN$(UNDERSCORE)THR$(UNDERSCORE)mask-d2bnorm.csv: results/%$(UNDERSCORE)SIR$(UNDERSCORE)ALN$(UNDERSCORE)THR$(UNDERSCORE)mask-d2bfit.csv
##	$(RSCRIPT) scripts/d2bnorm.R $^ > $@

##results/%$(UNDERSCORE)SIR$(UNDERSCORE)ALN$(UNDERSCORE)THR$(UNDERSCORE)EdUtoMark_NNA.csv: results/%$(UNDERSCORE)SIR$(UNDERSCORE)ALN$(UNDERSCORE)THR$(UNDERSCORE)Rep-centroids.csv results/%$(UNDERSCORE)SIR$(UNDERSCORE)ALN$(UNDERSCORE)THR$(UNDERSCORE)EdU-centroidsMCN.csv
##	$(OCTAVE) scripts/NNA_S-phase.m $^ results/$*$(UNDERSCORE)SIR$(UNDERSCORE)ALN$(UNDERSCORE)THR$(UNDERSCORE)EdUtoMark_NNA.csv results/$*$(UNDERSCORE)SIR$(UNDERSCORE)ALN$(UNDERSCORE)THR$(UNDERSCORE)MarktoEdU_NNA.csv

## results/%$(UNDERSCORE)SIR$(UNDERSCORE)EAL$(UNDERSCORE)THR$(UNDERSCORE)EdU-rd2b.csv: results/%$(UNDERSCORE)SIR$(UNDERSCORE)EAL$(UNDERSCORE)THR$(UNDERSCORE)EdU-md2b.csv
## fake recipe, does nothing because EdU-rd2b is already made as the final argument and target of EdU-md2b


##HCT116-SCC1-mAID data recipes:
results/%$(UNDERSCORE)mask.tif: data/%.tif | $$(dir $$@).dirstamp
	$(OCTAVE) scripts/edu_mask2.m 3 $^ $@

results/%$(UNDERSCORE)C594-centroids.csv: results/%$(UNDERSCORE)mask.tif data/%.tif
	$(OCTAVE) scripts/foci_centroids2.m 1 $^ 3 > $@

results/%$(UNDERSCORE)C488-centroids.csv: results/%$(UNDERSCORE)mask.tif data/%.tif
	$(OCTAVE) scripts/foci_centroids2.m 2 $^ 3 > $@

results/%$(UNDERSCORE)SEG.tif: results/%$(UNDERSCORE)mask.tif data/%.tif
	$(RSCRIPT) scripts/chromseg.R 3 $^ $@

results/%$(UNDERSCORE)SIR$(UNDERSCORE)THR$(UNDERSCORE)ALN$(UNDERSCORE)C594-distn.csv: results/%$(UNDERSCORE)SIR$(UNDERSCORE)THR$(UNDERSCORE)ALN$(UNDERSCORE)C594-centroids.csv results/%$(UNDERSCORE)SIR$(UNDERSCORE)THR$(UNDERSCORE)ALN$(UNDERSCORE)SEG.tif data/%$(UNDERSCORE)ALN$(UNDERSCORE)MCNR-mask.tif
	$(RSCRIPT) scripts/distribution.R $^ > $@

results/%$(UNDERSCORE)SIR$(UNDERSCORE)THR$(UNDERSCORE)ALN$(UNDERSCORE)C488-distn.csv: results/%$(UNDERSCORE)SIR$(UNDERSCORE)THR$(UNDERSCORE)ALN$(UNDERSCORE)C488-centroids.csv results/%$(UNDERSCORE)SIR$(UNDERSCORE)THR$(UNDERSCORE)ALN$(UNDERSCORE)SEG.tif data/%$(UNDERSCORE)ALN$(UNDERSCORE)MCNR-mask.tif
	$(RSCRIPT) scripts/distribution_HCT116.R $^ > $@

results/%$(UNDERSCORE)SIR$(UNDERSCORE)THR$(UNDERSCORE)ALN$(UNDERSCORE)bound.tif: results/%$(UNDERSCORE)SIR$(UNDERSCORE)THR$(UNDERSCORE)ALN$(UNDERSCORE)SEG.tif
	$(OCTAVE) scripts/bound.m 1 $^ $@

results/%$(UNDERSCORE)SIR$(UNDERSCORE)THR$(UNDERSCORE)ALN$(UNDERSCORE)C488-md2b.csv: results/%$(UNDERSCORE)SIR$(UNDERSCORE)THR$(UNDERSCORE)ALN$(UNDERSCORE)SEG.tif results/%$(UNDERSCORE)SIR$(UNDERSCORE)THR$(UNDERSCORE)ALN$(UNDERSCORE)C488-centroids.csv results/%$(UNDERSCORE)SIR$(UNDERSCORE)THR$(UNDERSCORE)ALN$(UNDERSCORE)bound.tif results/%$(UNDERSCORE)SIR$(UNDERSCORE)THR$(UNDERSCORE)ALN$(UNDERSCORE)mask.tif
	$(OCTAVE) scripts/d2bound-HCT116.m 1 $^ results/$*$(UNDERSCORE)SIR$(UNDERSCORE)THR$(UNDERSCORE)ALN$(UNDERSCORE)C488-md2b.csv results/$*$(UNDERSCORE)SIR$(UNDERSCORE)THR$(UNDERSCORE)ALN$(UNDERSCORE)C488-rd2b.csv

results/%$(UNDERSCORE)SIR$(UNDERSCORE)THR$(UNDERSCORE)ALN$(UNDERSCORE)C594-md2b.csv: results/%$(UNDERSCORE)SIR$(UNDERSCORE)THR$(UNDERSCORE)ALN$(UNDERSCORE)SEG.tif results/%$(UNDERSCORE)SIR$(UNDERSCORE)THR$(UNDERSCORE)ALN$(UNDERSCORE)C594-centroids.csv results/%$(UNDERSCORE)SIR$(UNDERSCORE)THR$(UNDERSCORE)ALN$(UNDERSCORE)bound.tif results/%$(UNDERSCORE)SIR$(UNDERSCORE)THR$(UNDERSCORE)ALN$(UNDERSCORE)mask.tif
	$(OCTAVE) scripts/d2bound-HCT116.m 1 $^ results/$*$(UNDERSCORE)SIR$(UNDERSCORE)THR$(UNDERSCORE)ALN$(UNDERSCORE)C594-md2b.csv results/$*$(UNDERSCORE)SIR$(UNDERSCORE)THR$(UNDERSCORE)ALN$(UNDERSCORE)C594-rd2b.csv


##results/%$(UNDERSCORE)SIR$(UNDERSCORE)THR$(UNDERSCORE)THR$(UNDERSCORE)C488-centroidsMCN.csv: results/%$(UNDERSCORE)SIR$(UNDERSCORE)ALN$(UNDERSCORE)THR$(UNDERSCORE)C488-centroids.csv data/%$(UNDERSCORE)ALN$(UNDERSCORE)MCNR-mask.tif results/%$(UNDERSCORE)SIR$(UNDERSCORE)ALN$(UNDERSCORE)THR$(UNDERSCORE)SEG.tif
##	$(RSCRIPT) scripts/centroid_volume_C488.R $^ > $@

results/%$(UNDERSCORE)SIR$(UNDERSCORE)THR$(UNDERSCORE)ALN$(UNDERSCORE)C594-d2bfit.csv: results/%$(UNDERSCORE)SIR$(UNDERSCORE)THR$(UNDERSCORE)ALN$(UNDERSCORE)C594-md2b.csv results/%$(UNDERSCORE)SIR$(UNDERSCORE)THR$(UNDERSCORE)ALN$(UNDERSCORE)C594-rd2b.csv
	$(RSCRIPT) scripts/d2bfit.R $^ > $@

results/%$(UNDERSCORE)SIR$(UNDERSCORE)THR$(UNDERSCORE)ALN$(UNDERSCORE)C488-d2bfit.csv: results/%$(UNDERSCORE)SIR$(UNDERSCORE)THR$(UNDERSCORE)ALN$(UNDERSCORE)C488-md2b.csv results/%$(UNDERSCORE)SIR$(UNDERSCORE)THR$(UNDERSCORE)ALN$(UNDERSCORE)C488-rd2b.csv
	$(RSCRIPT) scripts/d2bfit.R $^ > $@

results/%$(UNDERSCORE)SIR$(UNDERSCORE)THR$(UNDERSCORE)ALN$(UNDERSCORE)C594-d2bnorm.csv: results/%$(UNDERSCORE)SIR$(UNDERSCORE)THR$(UNDERSCORE)ALN$(UNDERSCORE)C594-d2bfit.csv
	$(RSCRIPT) scripts/d2bnorm.R $^ > $@

results/%$(UNDERSCORE)SIR$(UNDERSCORE)THR$(UNDERSCORE)ALN$(UNDERSCORE)C488-d2bnorm.csv: results/%$(UNDERSCORE)SIR$(UNDERSCORE)THR$(UNDERSCORE)ALN$(UNDERSCORE)C488-d2bfit.csv
	$(RSCRIPT) scripts/d2bnorm.R $^ > $@

##results/%$(UNDERSCORE)SIR$(UNDERSCORE)ALN$(UNDERSCORE)THR$(UNDERSCORE)C594-NNA.csv: results/%$(UNDERSCORE)SIR$(UNDERSCORE)ALN$(UNDERSCORE)THR$(UNDERSCORE)C488-centroidsMCN.csv results/%$(UNDERSCORE)SIR$(UNDERSCORE)ALN$(UNDERSCORE)THR$(UNDERSCORE)C594-centroidsMCN.csv 
##	$(OCTAVE) scripts/NNA_HCT.m $^ results/$*$(UNDERSCORE)SIR$(UNDERSCORE)ALN$(UNDERSCORE)THR$(UNDERSCORE)C488-C594_NNA.csv results/$*$(UNDERSCORE)SIR$(UNDERSCORE)ALN$(UNDERSCORE)THR$(UNDERSCORE)C594-C488_NNA.csv


#Raserfish data recipes:
##results/%$(UNDERSCORE)FISH-mask.tif: data/%.tif | $$(dir $$@).dirstamp
##	$(OCTAVE) scripts/nucleus_mask.m 2 $^ $@

##results/%$(UNDERSCORE)FISH-centroids.csv: results/%$(UNDERSCORE)FISH-mask.tif data/%.tif
##	$(OCTAVE) scripts/foci_centroids2.m 1 $^ 2 > $@

##results/%$(UNDERSCORE)FISH-submask.tif: data/%.tif | $$(dir $$@).dirstamp
##	$(OCTAVE) scripts/FISH_mask2.m 1 $^ $@

##results/%$(UNDERSCORE)SIR$(UNDERSCORE)ALN$(UNDERSCORE)THR$(UNDERSCORE)FISH-centroidsvol.csv: results/%$(UNDERSCORE)SIR$(UNDERSCORE)ALN$(UNDERSCORE)THR$(UNDERSCORE)FISH-centroids.csv data/%$(UNDERSCORE)ALN$(UNDERSCORE)MCNR-mask.tif results/%$(UNDERSCORE)SIR$(UNDERSCORE)ALN$(UNDERSCORE)THR$(UNDERSCORE)SEG.tif
##	$(RSCRIPT) scripts/centroid_volume.R $^ > $@

##results/%$(UNDERSCORE)FISH-SEG.tif: results/%$(UNDERSCORE)FISH-mask.tif data/%.tif
##	$(RSCRIPT) scripts/chromseg.R 2 $^ $@

##results/%$(UNDERSCORE)SIR$(UNDERSCORE)ALN$(UNDERSCORE)THR$(UNDERSCORE)FISH-distn.csv: results/%$(UNDERSCORE)SIR$(UNDERSCORE)ALN$(UNDERSCORE)THR$(UNDERSCORE)FISH-centroids.csv results/%$(UNDERSCORE)SIR$(UNDERSCORE)ALN$(UNDERSCORE)THR$(UNDERSCORE)FISH-SEG.tif data/%$(UNDERSCORE)ALN$(UNDERSCORE)MCNR-mask.tif
##	$(RSCRIPT) scripts/distribution.R $^ > $@

#NNA data recipes:
##results/%$(UNDERSCORE)NNA-mask.tif: data/%.tif | $$(dir $$@).dirstamp
##	$(OCTAVE) scripts/nucleus_mask_NNA.m 2 $^ $@

##results/%$(UNDERSCORE)C1-centroids.csv: results/%$(UNDERSCORE)NNA-mask.tif data/%.tif
##	$(OCTAVE) scripts/foci_centroids2.m 1 $^ 2 > $@

##results/%$(UNDERSCORE)C2-centroids.csv: results/%$(UNDERSCORE)NNA-mask.tif data/%.tif
##	$(OCTAVE) scripts/foci_centroids2.m 2 $^ 2 > $@

##results/%$(UNDERSCORE)THR$(UNDERSCORE)MCF$(UNDERSCORE)ALN$(UNDERSCORE)C1-C2_NNA.csv: results/%$(UNDERSCORE)THR$(UNDERSCORE)MCF$(UNDERSCORE)ALN$(UNDERSCORE)C1-centroids.csv results/%$(UNDERSCORE)THR$(UNDERSCORE)MCF$(UNDERSCORE)ALN$(UNDERSCORE)C2-centroids.csv
##	$(OCTAVE) scripts/NNA.m $^ results/$*$(UNDERSCORE)THR$(UNDERSCORE)MCF$(UNDERSCORE)ALN$(UNDERSCORE)C1-C2_NNA.csv results/$*$(UNDERSCORE)THR$(UNDERSCORE)MCF$(UNDERSCORE)ALN$(UNDERSCORE)C2-C1_NNA.csv

#results/%$(UNDERSCORE)C1-C2_NNA.csv: results/%$(UNDERSCORE)C1-centroids.csv results/%$(UNDERSCORE)C2-centroids.csv
#	$(OCTAVE) scripts/NNA.m $^ results/$*$(UNDERSCORE)C1-C2_NNA.csv results/$*$(UNDERSCORE)C2-C1_NNA.csv


## Generate required parent directories in results/
##

## For this to work, recipes that generate files in results/ from data/
## must have an order-only prerequisite on a `.dirstamp` file on the
## directory of the target. Example:
##
## results/foo/bar/xxx.csv: data/foo/bar/xxx.tif | results/foo/bar/.dirstamp
##     $(OCTAVE) scripts/analysis.m $^ > $@
##
## Why we want to have data and results separated:
##    * Organization  --- we don't want to mix the two things.
##    * Allow us to have data as read-only to avoid accidentally removing data.
##
## Details on this solution:
##    * We use 'mkdir -p' however this is not reliable in non-GNU
##      systems and definitely not reliable over smb mounts.
##    * Can't copy data into results since that's several TB.
##
## See our question and answer at: 
##  http://stackoverflow.com/questions/35995714/writing-a-recursive-make-recipe-with-prerequisite-on-parent-directory

RESULT_DIRS := $(patsubst data/%, results/%, $(shell find data/* -type d -print))
DIRSTAMPS := $(addsuffix /.dirstamp, $(RESULT_DIRS))

##Uncomment both of the "@mkdir..." lines below if a new directory was added to data/
results/.dirstamp:
#	@mkdir -p $(dir $@)
	@touch $@

.SECONDEXPANSION:
$(DIRSTAMPS): $$(dir $$(patsubst %/.dirstamp, %, $$@)).dirstamp
#	@mkdir $(dir $@)
	@touch $@
