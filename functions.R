##
## All functions used by main.R
##

library(glmnet)

countsFromBam <- function(bfile, ranges){
        
        # current bed_file name:
        bam_file_id <- tail(strsplit(bfile,"\\/")[[1]],1)
        print(paste("Counting for file",bam_file_id))

        # we need to select chromosomes according to bam-header information
        cur_ranges <- ranges[seqnames(ranges)%in%names(scanBamHeader(bfile)[[1]]$targets)]

        # set the parameters for counting only within the current ranges
        param <- ScanBamParam(which=cur_ranges, what=c("pos", "qwidth"), flag = scanBamFlag(isUnmappedQuery=FALSE, isNotPassingQualityControls=FALSE, isDuplicate=FALSE))

        # count and save
        bin_coverage <- scanBam(bfile, param=param)
        ids <- paste(seqnames(cur_ranges),":",start(cur_ranges),"-",end(cur_ranges),sep="")
        bin_coverage <- bin_coverage[match(ids, names(bin_coverage))]

	# write out counts and annotation
        bin_coverage_table <- cbind(as.data.frame(cur_ranges)[,1:5],sapply(bin_coverage,function(x){length(x$pos)}))
       	names(bin_coverage_table) <- c("seqnames", "start", "end", "width", "strand","coverage") 

	return(bin_coverage_table)
}

combineCoverage <- function(counts_sublist, outfile){
        cur_ranges <- counts_sublist[[1]]
        if(length(counts_sublist)>1){
                for(i in 2:length(counts_sublist)){
                        cur_ranges[,"coverage"] <- cur_ranges[,"coverage"] + counts_sublist[[i]][,"coverage"]
                }
        }
        write.table(cur_ranges, file=outfile)
        return(outfile)
}

readData <- function(filename, reorderBy=NULL){

	# read data and separate into ids and coverage for easy access
	d <- read.table(filename, header=T)
	ids <- d[, names(d)%in%c("seqnames", "start", "end", "width", "strand") ]		

	# if order provided, reorder data	
	if(!is.null(reorderBy)){	
		reorder <- match(reorderBy, apply(ids,paste,collapse="_"))
		ids <- ids[reorder,]
		d <- d[reorder]	
	}

	return( list( 	coverage=as.numeric(d[,"coverage"]), 
			ids=ids) )
}


collectData <- function(countfilelist){
	
	# read in first sample to get 'reference'	
	d <- readData(countfilelist[1])

	# build all data matrix
	all_data <- matrix(NA, ncol=length(countfilelist), nrow=length(d$coverage), dimnames=list(NULL,countfilelist) )
	all_data[,1] <- d$coverage 
	for(filename in countfilelist){
	        d <- readData(filename)
		all_data[,filename] <- d$coverage
	}
	all_data <- as.data.frame( cbind(d$ids, all_data) )

	return(all_data)
}

widthNorm <- function(datamat, widths){

	datamat_norm <- datamat
	for(i in 1:ncol(datamat)){
		datamat_norm[,i] <- as.numeric(datamat[,i]) / as.numeric(widths) 
	}
	
	return(datamat_norm)
}


ratioNorm <- function(datamat, inputmat, si_match){

        datamat_norm <- datamat

	# normalize for each pair of sample/input 
        for(i in 1:nrow(si_match)){
		sampleID <- si_match[i,"sample"]

		# data of sample and corresponding input
		data <- datamat[,sampleID]
		input <- inputmat[,si_match[i,"input"]]

		# subselect/remove outliers from estimation of median	
                sub <- input>0 & data>0 & data<=quantile(data,0.99) & input<=quantile(input,0.99)

		# if there are enough points for estimation...
                if(sum(sub)>100){
			
			# estimate median
                        m <- median( data[sub] / input[sub] )

			# normalize with median
                        datamat_norm[,sampleID] <- (data+1) / (input+1)*m
                }else{
               		# not enough points to make a good approximation 
		        datamat_norm[,sampleID] <- NA
                }               
        }
        
	return(datamat_norm)
}

log1Transform <- function(datamat){

	datamat_norm <- datamat

        for (i in 1:ncol(datamat)){ 
       		datamat_norm[,i] <- log( datamat[,i] + 1 ) 
	}

        return(datamat_norm)
}

zscoreTransform <- function(datamat){

        datamat_norm <- datamat

        for (i in 1:ncol(datamat)){
                datamat_norm[,i] <- (datamat[,i] - mean(datamat[,i])) / sd(datamat[,i]) 
        }

        return(datamat_norm)
}

checkNAStatus <- function(datamat){

	select <- !apply(datamat, 2, function(x){ any(is.na(x)) })

        return( datamat[,select] )
}

normalizeData <- function(datamat, sampleInputMatch){

	# subselect by sample/input
	onlysamples <- unlist(sampleInputMatch[,"sample"])
	onlyinput <- unlist(sampleInputMatch[,"input"]) 
	if( any(duplicated(onlysamples)|duplicated(onlyinput)) ){ print("Warning: Sample or Input name duplicated in sampleInputMatch") }
 	if( any(! names(datamat)%in%c("seqnames", "start", "end", "width", "strand", onlysamples, onlyinput)) ){ print("Warning: Some columns in the matrix are not part of either annotation or sampleInputMatch matrix") }	

	# normalize each segment by width
	datamat[,c(onlysamples,onlyinput)] <- widthNorm(datamat[,c(onlysamples,onlyinput)], datamat[,"width"])
	
	# throw out all regions with no signal for most samples (excluding input)
	sub <- apply( datamat[,onlysamples], 1, function(x) { sum(x>0)>1 } )
	datamat <- datamat[sub,]

	# Perform ratio normalization
	datamat[,onlysamples] <- ratioNorm(datamat[,onlysamples], datamat[,onlyinput], sampleInputMatch) 
	datamat <- datamat[,!names(datamat)%in%onlyinput]
	datamat <- checkNAStatus( datamat )

	# Log the data adding a pseudocount of 1
	datamat[,onlysamples] <- log1Transform(datamat[,onlysamples])
	datamat <- checkNAStatus( datamat )

	# Z-score transformation
	datamat[,onlysamples] <- zscoreTransform(datamat[,onlysamples]) 
	datamat <- checkNAStatus( datamat )
	
	return( datamat ) 	
}


predictByChromState <- function(chromstate, data, chromstate_ranges, hms, cms, spcn_path=NULL, alphas=c(seq(0,0.1,0.01),seq(0.2,0.9,0.1)), cv_n=10, sdfactor=1, mccores=1){

	# overlap with chromatin state assignment
        chrom_ranges <- chromstate_ranges[ mcols(chromstate_ranges)[,1]==chromstate ]
        data_ranges <- GRanges(IRanges(start=as.numeric(data[,"start"]), end=as.numeric(data[,"end"])), seqname=data[,"seqnames"])
        subselect <- overlapsAny(data_ranges,chrom_ranges,minoverlap=2)
        sub_data <- data[subselect,]

	# normalize subset
	sub_data[,c(hms,cms)] <- zscoreTransform(sub_data[,c(hms,cms)])  
	sub_data <- checkNAStatus(sub_data)

	if(mccores>=2){
		# run both in parallel
		tmp <- mclapply( c(	"spcnWrapper(as.matrix(sub_data[,c(cms,hms)]), spcn_path)",
				 	"enetWrapper(sub_data, hms, cms, alphas, cv_n, sdfactor, mccores-1)" ),
				 function(x) { eval( parse(text=x) )}, mc.cores=2)
		spcn_res <- tmp[[1]]
		enet_res <- tmp[[2]]
	}else{
	        # run SPCN
        	spcn_res <- spcnWrapper(as.matrix(sub_data[,c(cms,hms)]), spcn_path)
        	# run ENET
        	enet_res <- enetWrapper(sub_data, hms, cms, alphas, cv_n, sdfactor)
	}
	
	# combine the results
	res <- combineAll(enet_res$net, spcn_res$net, hms, cms) 

	return( list( net=res, spcn_res=spcn_res, enet_res=enet_res ) )
}


spcnWrapper <- function(tmpD, spcn_path){

        ## load SPCN script
	if(is.null(spcn_path)){
		source("http://spcn.molgen.mpg.de/code/sparse_pcor.R")
        }else{
		source(spcn_path)
	}
	
        # do the SPCN
        res <- sparse_pcor(tmpD,rank=TRUE,sparse=TRUE)
        spcn_mat <- res$sparse_partial

	# reformat into a table
        nr_vars <- ncol(tmpD)
        spc_net <- matrix(NA,nrow=nr_vars^2,4)
        for(i in 1:nrow(spcn_mat)){
                for(j in 1:ncol(spcn_mat)){
                        spc_net[j+((i-1)*nr_vars),1] = dimnames(spcn_mat)[[1]][i]
                        spc_net[j+((i-1)*nr_vars),2] = dimnames(spcn_mat)[[2]][j]
                        if(j>i){
                                spc_net[j+((i-1)*nr_vars),3] = round(spcn_mat[i,j],4)
                       		spc_net[j+((i-1)*nr_vars),4] = as.numeric(spcn_mat[i,j]>0) 
			}else{
                                spc_net[j+((i-1)*nr_vars),3] = 0
                        }
                }
        }
        spc_net <- spc_net[!is.na(spc_net[,3]),]
        spc_net <- spc_net[spc_net[,3]!="0",]
	
	return( list(net=spc_net, coeffs=res$partial) )
}

makeCVids <- function(nr_ys, cv_n){
	
	ids <- rep(1:cv_n, each=nr_ys%/%cv_n)
        if(length(ids)!=nr_ys){ ids <- c(ids, 1:(nr_ys-length(ids))) }
        ids <- sample(ids,length(ids),replace=F)

	return(ids)
}

findBestParams <- function(Xmat, y, alphas, cv_n){

	# get cross-folds 
	a_ids <- makeCVids(length(y), cv_n=cv_n)

	# remember error and standard deviation for each crossvalidation 
	cvcss_alphas <- rep(NA, length(alphas))
	cvse_alphas <- rep(NA, length(alphas))
       	cv_models <- vector("list", length(alphas)) 

	# run glmnet for each alpha                
	for(i in 1:length(alphas)){

		#find best lambda given alpha
		model <- cv.glmnet(Xmat, y, family="gaussian", alpha=alphas[i], foldid=a_ids)
		cvcss_alphas[i] <- model$cvm[which(model$lambda==model$lambda.1se)]
		cvse_alphas[i] <- model$cvsd[which(model$lambda==model$lambda.1se)]
       		cv_models[[i]] <- model
	}
        
	# select alpha with error within min-cv-error + sd                
	alpha_id <- which.min(cvcss_alphas)
	alpha_id <- which(cvcss_alphas > cvcss_alphas[alpha_id] & cvcss_alphas < (cvcss_alphas+cvse_alphas)[alpha_id])
       	alpha_id <- which.min(alphas[alpha_id])
	
	return( cv_models[[alpha_id]] ) 
}

doEnet <- function(data, ygroup, xgroup, alphas, cv_n, mccores){

	cv_coefficients <- mclapply(ygroup, predictY, data=data, xgroup=xgroup, alphas=alphas, cv_n=cv_n, mc.cores=mccores)
	names(cv_coefficients) <- ygroup

	return( cv_coefficients )
}

predictY <- function(cm, data, xgroup, alphas, cv_n){

	# current data
        y <- data[,cm]
        cur_xgroup <- xgroup[xgroup!=cm]
        X <- apply(as.matrix(data[,cur_xgroup]),2,as.numeric)

	# get the ids for cv
	ids <- makeCVids(length(y), cv_n=cv_n)

	# do a cv for getting stability of features
	cv_coeffs_mat <- matrix(NA, ncol=cv_n, nrow=length(cur_xgroup), dimnames=list(cur_xgroup, 1:10))
	for(i in 1:cv_n){
		# find best alpha
		model_fit_cv <- findBestParams(X[ids!=i,], y[ids!=i], alphas, cv_n)
		cv_coeffs_mat[,i] <- as.vector(coef(model_fit_cv))[-c(1)]
	}

	return( cv_coeffs_mat )
}

makeEnetTable <- function(cv_coefficients, sdfactor){

        # get the median coefficient over all cross-folds for each prediction 
        cv_coeffs_median <- lapply(cv_coefficients,function(x){apply(x,1,median)})

        # select the important features
	group <- names(cv_coeffs_median)
        important_features <- vector("list",length(cv_coeffs_median))
        for(i in 1:length(cv_coeffs_median)){

                # do the filtering
                imp <- 	cv_coeffs_median[[i]] <= mean(cv_coeffs_median[[i]])-(sdfactor*sd(cv_coeffs_median[[i]])) | 
			cv_coeffs_median[[i]] >= mean(cv_coeffs_median[[i]])+(sdfactor*sd(cv_coeffs_median[[i]]))

                # add if there is something to add and if any feature is predictive
                if( (!is.null(imp)) & any(round(cv_coeffs_median[[i]],5)!=0) ){
               		important_features[[i]] <- cbind( 	rep(group[i],sum(imp)), 
								names(cv_coeffs_median[[i]])[imp],
								round(cv_coeffs_median[[i]][imp],5),
								as.numeric(cv_coeffs_median[[i]][imp]>0)) 
		}
        }
	
	return( do.call(rbind, important_features) )
}

enetWrapper <- function(data, hms, cms, alphas, cv_n, sdfactor, mccores=1){

	enet_res_hms <- doEnet(data, hms, cms, alphas, cv_n, mccores)
	imp_res <- makeEnetTable(enet_res_hms, sdfactor)

	enet_res_cms <- doEnet(data, cms, cms, alphas, cv_n, mccores)
	imp_res <- rbind(imp_res, makeEnetTable(enet_res_cms, sdfactor))

	return( list(net=imp_res, coeffs=c(enet_res_hms,enet_res_cms)) )
}

combineAll <- function(enet, spcn, hms, cms){

	# prepare result tables for easy merging
	cm_cm_spcn <- spcn[which(spcn[,1]%in%cms & spcn[,2]%in%cms),]
	cm_cm_spcn <- rbind(cm_cm_spcn, cm_cm_spcn[,c(2,1,3,4)]) # the spcn is undirected
	hm_cm_spcn <- spcn[(spcn[,2]%in%hms & spcn[,1]%in%cms) | (spcn[,1]%in%hms & spcn[,2]%in%cms),]
        tmp <- hm_cm_spcn[,1]
        select <- tmp%in%cms
        tmp[select] <- hm_cm_spcn[select,2]
        hm_cm_spcn[select,2] <- hm_cm_spcn[select,1]
        hm_cm_spcn[,1] <- tmp

	# merge enet and spcn
        combined_all <- rbind(as.matrix(enet[,c(1,2)]),hm_cm_spcn[,c(1,2)],cm_cm_spcn[,c(1,2)])
        combined_all <- combined_all[!duplicated(combined_all),]
        combined_all <- cbind(combined_all,matrix(-1,nrow=nrow(combined_all),ncol=4))
        dimnames(combined_all) <- list(NULL,c("HM","CM","ENET_W","ENET_B","SPCN_W","SPCN_B"))
        for(i in 1:nrow(combined_all)){
		select <- apply(enet[,c(1,2)],1,function(x){all(x==combined_all[i,1:2])})
		if(any(select)){ combined_all[i,3:4] <- as.matrix(enet[select,3:4]) }
		select <- apply(hm_cm_spcn[,c(1,2)],1,function(x){all(x==combined_all[i,1:2])})
		if(any(select)){ combined_all[i,5:6] <- hm_cm_spcn[select,3:4] }
		select <- apply(cm_cm_spcn[,c(1,2)],1,function(x){all(x==combined_all[i,1:2])})
		if(any(select)){ combined_all[i,5:6] <- cm_cm_spcn[select,3:4] }		
	}

	# make the net symmetric
        symnet <- c() 
        tmp_net <- combined_all 
        for(i in 1:nrow(combined_all)){
		current_row <- combined_all[i,]
       
	 	# this is only neccessary for the cms	
		if(!any(current_row%in%hms)){
			# select the reverse interaction and check if both agree on the sign
			test <- tmp_net[,1]%in%current_row[2] & tmp_net[,2]%in%current_row[1] 
			test <- test & tmp_net[,4]%in%current_row[4] & tmp_net[,6]%in%current_row[6]

			# if there is any...
       		        if(sum(test)>0){
                		# add it to the symmetric net	
				symnet <- rbind(symnet, current_row[c(1,2,4,6)])
				
				# and take it out of the temporary matrix	
				tmp_net <- tmp_net[-c(which(test), which(tmp_net[,1]%in%current_row[1] & tmp_net[,2]%in%current_row[22])),]
                        }
                }else{
			symnet <- rbind(symnet, combined_all[i, c(1,2,4,6)])
		}
        }
	colnames(symnet) <- c("HM","CM","ENET_B","SPCN_B")
		
	return( symnet[(symnet[,"ENET_B"]!="-1" & symnet[,"SPCN_B"]!="-1" & symnet[,"ENET_B"]==symnet[,"SPCN_B"]),] )
}













