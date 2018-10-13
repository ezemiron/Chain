dirchosen <- "/home/eze/Documents/papers/chrom_marks/results"

Cell = "C127"
Condition = "Hyper-os"

topPATH = paste(Cell, Condition, sep = "/")

dir_names <- list.dirs(topPATH, recursive = TRUE)
#dir_names <- dir_names[3:length(dir_names)]

var_names <- list.dirs(topPATH, recursive = TRUE, full.names=FALSE)
#var_names <- var_names[3:length(var_names)]

for (dir_name in dir_names){

    newdir <- paste(dirchosen, dir_name, sep ="/")
    setwd(newdir)
    dir_num <- which(dir_names == dir_name)
    var_name <- var_names[dir_num]

    df <- data.frame(NA)
    
    names(df) <- NA
    df <- cbind(df,df,df,df)
    
    allfilesandpaths <- list.files(pattern= var_name, full.names=TRUE);
    if (length(allfilesandpaths) > 0){
    
      filesandpaths  <- grep("*distn.csv",allfilesandpaths, value = TRUE)
      if (length(filesandpaths) > 0){
    
        for (fileandpath in filesandpaths){
            if (length(fileandpath) > 0){
              csv <- read.csv(fileandpath, colClasses=c("NULL",NA,"NULL",NA,"NULL",NA,"NULL","NULL"))
              cell_num <- which(filesandpaths == fileandpath)
              csv$cell <- cell_num
              names(csv) <- rep(NA, ncol(csv))
              df <- rbind(df,csv)
            }
        }
        df <- df[2:nrow(df),]
        names(df) <- c("class", "spotnumn","classnumn","cell")
#    df <- transform(df, AvNorm = rowMeans(df[,-1], na.rm = TRUE))
#    df <- transform(df, Log2AvNorm = log2(df$AvNorm))

        savename  <- paste0(var_name,"_spot-class-count_volnorm.csv")
        write.csv(df,savename);

        df2 <-aggregate(df[,2:3], list(class =df$class), mean)
        df2 <-transform(df2, numofcells = rep(max(df$cell),dim(df2)[1]))
      
        savename2  <- paste0(var_name,"_spot-class-summary_volnorm.csv")
        write.csv(df2,savename2);
        }
    }
}

