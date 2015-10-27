##
## Main script, adjust filenames, paths etc. here
##

##
# required libraries
##

library(GenomicRanges)
library(Rsamtools)
library(parallel)

##
# all functions needed
##
source("functions.R")

##
# all parameters needed
##

# a one-by-line list of absolute filenames to all bam files containing the aligned reads
bamfilelist_filename <- commandArgs(trailingOnly=TRUE)[1] # "data/bamfilelist.txt"
if(!file.exists(bamfilelist_filename)){
	stop(paste("File does not exist:",bamfilelist_filename))
}

# main directory for output files (inkl. bams, counts, etc) 
outputdir <- "results/"
subfolder <- strsplit(basename(bamfilelist_filename),"\\.")[[1]][1]

# filename of file containing the annotation of the chromatin state ranges, bed_file_style with last colum=chromatin state of bin
chromBins_filename <- "data/all_chromstate_bins.txt"

# filename of file containing hms/cms, sample name, corresponding input name
mapInputSample_filename <- "data/inputSampleMap.txt"

# number of available cores
nrcores <- 10  

# cutoff: the factor for multiplying the standard deviation of the features in the Elasticnet 
sdf <- 2 

# files containing the hm and cm names
hms <- as.character(read.table("data/hmnames.txt",header=FALSE)[,1])
cms <- as.character(read.table("data/cmnames.txt",header=FALSE)[,1])

# create results paths
dir.create(file.path(outputdir))
outputdir <- paste0(outputdir,subfolder,"/")
dir.create(file.path(outputdir))
dir.create(file.path(outputdir,"counts"))

##
# Read and process bam files into count matrix
##

# list of all available bedfiles 
bamfilelist <- as.character(read.table(bamfilelist_filename)[,1])

# read the chromatin segments
bins <- read.table(chromBins_filename, header=TRUE)
chromranges <- GRanges(IRanges(start=as.numeric(bins[,2]), end=as.numeric(bins[,3])), seqname=bins[,1], mcols=bins[,6])

# count the reads overlapping each chromatin segment for each file
count_list <- mclapply(bamfilelist, countsFromBam, ranges=chromranges, mc.cores=nrcores)
names(count_list) <- sapply(bamfilelist, function(x){ tail(strsplit(x,"\\/")[[1]],1) })

##
# combine replicates just by adding up the counts
# note this will also amplify the input and create temporary files
##

mapInputSample <- read.table(mapInputSample_filename, header=F)
mapInputSample <- mapInputSample[apply(mapInputSample[,2:3], 1, function(x){ any(grepl(x[1], bamfilelist)) & any(grepl(x[2], bamfilelist)) }),]
all_comps <- unique(as.character(mapInputSample[,1]))
countfilelist <- c()
for(comp in all_comps){
	# select all replicates
       	select <- mapInputSample[,1]%in%comp 
	# IP:	
	countfilelist <- c(countfilelist, combineCoverage(count_list[as.character(mapInputSample[select,2])], outfile=paste0(outputdir, "counts/", comp, "_bins_ranges_coverage.bed")))
	# Input:
        countfilelist <- c(countfilelist, combineCoverage(count_list[as.character(mapInputSample[select,3])], outfile=paste0(outputdir, "counts/", comp, "_input_bins_ranges_coverage.bed")))
}

##
# read in the data matrix
##

# collect from the temporary files above
all_data <- collectData(countfilelist)

# give them easier names
cutNames <- function(x){ gsub("_bins_ranges_coverage.bed","",tail(strsplit(x,"\\/")[[1]],n=1)) }
names(all_data) <- sapply(names(all_data), cutNames)

# construct a map mapping each IP to Input
allsamples <- grep(paste(c(hms,cms),collapse="|"), names(all_data), value=T)
select <- grepl("input",allsamples)
sampleInputmap <- cbind( allsamples[!select], sapply(allsamples[!select], function(x){allsamples[select&grepl(paste0("^",x,"_"),allsamples)]} ))
colnames(sampleInputmap) <- c("sample","input")
 
##
# normalize
##

all_data_norm <- normalizeData(all_data, sampleInputmap)
save(list=ls(), file=paste0(outputdir,"all_data_norm.RData"))

##
# Get for each chrom state the merged enet and spcn results
##

# run prediction and save results
all_states <- sort(unique(mcols(chromranges)[,1]))
res <- lapply(all_states, predictByChromState, chromstate_ranges=chromranges, data=all_data_norm, hms=hms, cms=cms, sdfactor=sdf, mccores=nrcores)
names(res) <- all_states
save(res, file=paste0(outputdir, "all_data_norm_all_states.RData"))

# write the final network
nets <- lapply(res, "[[", 1)
res_mat <- lapply(names(nets), function(x, net){ mat <- cbind(net[[x]], rep(x, nrow(net[[x]]))); colnames(mat) <- c(colnames(net[[x]]),"STATE"); return(mat) }, net=nets)
write.table(do.call(rbind,res_mat), file=paste0(outputdir, "chromnet.txt"), quote=FALSE, row.names=FALSE)

##
# Make the supplementary tables
##
library(xlsx)
library(gdata)

forsupp_enweights <- array(NA, c(length(all_states), length(cms), length(c(cms,hms))), dimnames=list( paste("state",all_states,sep="_"), cms, c(cms,hms) ))
forsupp_pcor <- array(NA, c(length(all_states), length(c(cms,hms)), length(c(cms,hms))), dimnames=list( paste("state",all_states,sep="_"), c(cms,hms), c(cms,hms) ))

for(s in all_states){
       
        # ENET 
        cv_coeffs_median <- lapply(res[[s]][["enet_res"]][["coeffs"]],function(x){(apply(x,1,median))})
        for(i in c(hms,cms)){
                forsupp_enweights[paste("state",s,sep="_"),names(cv_coeffs_median[[i]]),i] <- cv_coeffs_median[[i]]
        }

        # SPCN
        forsupp_pcor[paste("state",s,sep="_"),row.names(res[[s]][["spcn_res"]][["coeffs"]]), colnames(res[[s]][["spcn_res"]][["coeffs"]])] <- res[[s]][["spcn_res"]][["coeffs"]]

}

writeAllStates <- function( states, restable, outfile){	
	for(i in states){ 
        	if(i%in%states[1]){ write.xlsx(as.data.frame(restable[i,,]), file=outfile, sheetName=i, showNA=F)
        	}else{ write.xlsx(as.data.frame(restable[i,,]), file=outfile, sheetName=i, append=T, showNA=F) }
	}
}

writeAllStates(dimnames(forsupp_pcor)[[1]], forsupp_pcor, paste0(outputdir,"supptable_pcor.xls"))
writeAllStates(dimnames(forsupp_enweights)[[1]], forsupp_enweights, paste0(outputdir,"supptable_enweights.xls"))

