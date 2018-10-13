#!/usr/bin/env Rscript

setwd("/home/eze/Documents/papers/chrom_marks/results/")
all_centroid_files <- list.files(pattern="centroids.csv", recursive=TRUE)

#marker <- "H3K27me3" ; downsample <- "Scc1mAID"; SIZE <- 500
marker <- "H3K4me3" ; downsample <- "Parental"; SIZE <- 4000
#marker <- "Smc3" ; downsample <- "none"
#marker <- "RNAPIIs2p"; downsample <- "none"
#marker <- "H3K9me2"; downsample <- "Scc1mAID"; SIZE <- 4750

all_centroid_files <- all_centroid_files[grep("/Aux2h/", all_centroid_files)]
all_centroid_files <- all_centroid_files[grep(marker, all_centroid_files)]
compiled <- data.frame(nfoci=rep(0, length(all_centroid_files)),
		total_mean=rep(0, length(all_centroid_files)),
		total_sd=rep(0, length(all_centroid_files)),
		IQ1=rep(0, length(all_centroid_files)),
		 IQ2=rep(0, length(all_centroid_files)),
		 IQ3=rep(0, length(all_centroid_files)),
		 IQ4=rep(0, length(all_centroid_files)),
		IQ5=rep(0, length(all_centroid_files))
		)
	rownames(compiled) <- all_centroid_files


to_downsample <- all_centroid_files[grep(downsample, all_centroid_files)]
downsampled <- data.frame(nfoci=rep(0, length(to_downsample)),
                total_mean=rep(0, length(to_downsample)),
                total_sd=rep(0, length(to_downsample)),
                IQ1=rep(0, length(to_downsample)),
                IQ2=rep(0, length(to_downsample)),
		IQ3=rep(0, length(to_downsample)),
		IQ4=rep(0, length(to_downsample)),
		IQ5=rep(0, length(to_downsample))
		)
	rownames(downsampled) <- to_downsample

for(centroid_file in all_centroid_files) {
	centroids <- read.csv(centroid_file)
	stack <- range(centroids$z)
	## cut off the top and bottom 1/8th of the z stack 
	## (for a bit less computational work, and because its really noisy)
	lowerTrunc <- stack[1] + (stack[2] -stack[1])/8
	upperTrunc <- stack[2] - (stack[2] - stack[1])/8

	cent <- centroids[which(centroids$z > lowerTrunc & centroids$z < upperTrunc),]
	distance <- matrix(nrow=nrow(cent), ncol=nrow(cent))

	for(i in 1:nrow(cent)) {
		distance[i,i] <- 0
			if(i < nrow(cent)){
				for(j in seq(from=(i+1), to=nrow(cent))) {
						distance[i,j] = distance[j,i] = ( 
						(cent[i,"x"]-cent[j,"x"])^2 + 
						(cent[i,"y"]-cent[j,"y"])^2 + 
						(cent[i, "z"]-cent[j,"z"])^2)^(1/2)
				}
			}
	}
	
	compiled[centroid_file, "nfoci"] <- nrow(cent)
	compiled[centroid_file, "total_mean"] <- mean(distance, na.rm=TRUE)
	compiled[centroid_file, "total_sd"] <- sd(distance, na.rm=TRUE)
	compiled[centroid_file, "IQ1"] <- quantile(distance, na.rm=TRUE)[1]
	compiled[centroid_file, "IQ2"] <- quantile(distance, na.rm=TRUE)[2]
	compiled[centroid_file, "IQ3"] <- quantile(distance, na.rm=TRUE)[3]
	compiled[centroid_file, "IQ4"] <- quantile(distance, na.rm=TRUE)[4]
	compiled[centroid_file, "IQ5"] <- quantile(distance, na.rm=TRUE)[5] 

	if( centroid_file %in% to_downsample ) {
		if(nrow(cent) > SIZE) {	
	
		keep = sample(seq(from=1, to=nrow(cent)), size=SIZE)
		distanceD <- distance[keep, keep]

		downsampled[centroid_file, "nfoci"] <- length(keep)
        	downsampled[centroid_file, "total_mean"] <- mean(distanceD, na.rm=TRUE)
       		downsampled[centroid_file, "total_sd"] <- sd(distanceD, na.rm=TRUE)
        	downsampled[centroid_file, "IQ1"] <- quantile(distanceD, na.rm=TRUE)[1]
        	downsampled[centroid_file, "IQ2"] <- quantile(distanceD, na.rm=TRUE)[2]
        	downsampled[centroid_file, "IQ3"] <- quantile(distanceD, na.rm=TRUE)[3]
        	downsampled[centroid_file, "IQ4"] <- quantile(distanceD, na.rm=TRUE)[4]
        	downsampled[centroid_file, "IQ5"] <- quantile(distanceD, na.rm=TRUE)[5]	
		} else { downsampled[centroid_file, "nfoci" ] <- nrow(cent) - downsample }
	}
}

compiled$celltype = sapply(strsplit(rownames(compiled), "/"), "[", 1)
compiled$condition = sapply(strsplit(rownames(compiled), "/"), "[", 2)
compiled$marker= sapply(strsplit(rownames(compiled), "/"), "[", 3)

savename <- paste0("/home/eze/Documents/papers/chrom_marks/results/", "focidistances_", marker, "_June26.csv")
write.csv(compiled,savename)

downsampled$celltype = sapply(strsplit(rownames(downsampled), "/"), "[", 1)
downsampled$condition = sapply(strsplit(rownames(downsampled), "/"), "[", 2)
downsampled$marker= sapply(strsplit(rownames(downsampled), "/"), "[", 3)

savename2 <- paste0("/home/eze/Documents/papers/chrom_marks/results/", "focidistances_", marker, "downsampled_June26.csv")
write.csv(downsampled,savename2)
