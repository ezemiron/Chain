library (EBImage, warn.conflicts = FALSE, verbose = FALSE)
library (bioimagetools, warn.conflicts = FALSE, verbose = FALSE)

setwd("/home/eze/Documents/papers/chrom_marks/results/")
nuclearVol=list.files(pattern="mask.tif", recursive=TRUE)
nuclearVols=data.frame(file=nuclearVol, vol=rep(0, length(nuclearVol)))

for(i in 1:length(nuclearVol)){
	DAPIm = readTIF(nuclearVol[i])
	nuclearVols[i, "vol"]=sum(DAPIm > 0)
}

write.csv(nuclearVols, "/home/eze/Documents/papers/chrom_marks/results/nuclearmaskVolumes.csv")


## compare nuclear mask pixel volume to SEG pixel volume
nuclearVols$filefix <- gsub("HCT116.*/Aux[02]h/.*/", "", nuclearVol$file)
for(unique(nuclearVols$file)){
	nuclearVols[which(nuclearVols$file==i),"vol"]
}
