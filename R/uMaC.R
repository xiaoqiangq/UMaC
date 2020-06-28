# file:   based on ichorCNA.R
# author: Qiang Wei, Ph.D.
#               Vanderbilt Genetics Institute
# contact: <qiang.wei@vanderbilt.edu>
#
# need some R packages
# ichorCNA: https://github.com/broadinstitute/ichorCNA
# HMMcopy website: http://compbio.bccrc.ca/software/hmmcopy/ and https://www.bioconductor.org/packages/release/bioc/html/HMMcopy.html
# date:   June 24, 2019
# description: A hidden Markov model (HMM) to detect CNAs based on sequencing depth profiles in data with enriched ctDNA fragments

library(optparse);

option_list <- list(
	make_option(c("--BED"), type = "character", help = "Path to tumor BED file. Required."),
	make_option(c("--id"), type = "character", default="test", help = "Patient ID. Default: [%default]"),
	make_option(c("--normal"), type="character", default="0.5", help = "Initial normal contamination; can be more than one value if additional normal initializations are desired. Default: [%default]"),
	make_option(c("--ploidy"), type="character", default="c(2,3)", help = "Initial tumour ploidy; can be more than one value if additional ploidy initializations are desired. Default: [%default]"),
	make_option(c("--normalPanel"), type="character", default=NULL, help="Median corrected depth from panel of normals. Default: [%default]"),

	make_option(c("--altFracThreshold"), type="numeric", default=0.08, help="Minimum proportion of bins altered required to estimate tumor fraction; if below this threshold, then will be assigned zero tumor fraction. Default: [%default]"),
	make_option(c("--chrs"), type="character", default="c(1:22)", help = "Specify chromosomes to analyze. Default: [%default]"),
	#make_option(c("--chrs"), type="character", default="c(1:22,\"X\")", help = "Specify chromosomes to analyze. Default: [%default]"),
	make_option(c("--chrTrain"), type="character", default="c(1:22)", help = "Specify chromosomes to estimate params. Default: [%default]"),
	make_option(c("--coverage"), type="numeric", default=NULL, help = "PICARD sequencing coverage. Default: [%default]"),
	make_option(c("--estimateScPrevalence"), type="logical", default=TRUE, help = "Estimate subclonal prevalence. Default: [%default]"),  
	make_option(c("--lambda"), type="character", default="NULL", help="Initial Student's t precision; must contain 4 values (e.g. c(1500,1500,1500,1500)); if not provided then will automatically use based on variance of data. Default: [%default]"),
	make_option(c("--maxCN"), type="numeric", default=3, help = "Total clonal CN states. Default: [%default]"),
	make_option(c("--maxFracCNASubclone"), type="numeric", default=0.7, help="Exclude solutions with fraction of subclonal events greater than this value. Default: [%default]"),
	make_option(c("--maxFracGenomeSubclone"), type="numeric", default=0.5, help="Exclude solutions with subclonal genome fraction greater than this value. Default: [%default]"),
	make_option(c("--minSegmentBins"), type="numeric", default=50, help="Minimum number of bins for largest segment threshold required to estimate tumor fraction; if below this threshold, then will be assigned zero tumor fraction."),
	make_option(c("--plotYLim"), type="character", default="c(-2,2)", help = "ylim to use for chromosome plots. Default: [%default]"),
	make_option(c("--txnE"), type="numeric", default=0.99999, help = "Self-transition probability. Increase to decrease number of segments. Default: [%default]"),
	#make_option(c("--txnStrength"), type="numeric", default=1e7, help = "Transition pseudo-counts. Exponent should be the same as the number of decimal places of --txnE. Default: [%default]"),
	make_option(c("--txnStrength"), type="numeric", default=10000, help = "Transition pseudo-counts. Exponent should be the same as the number of decimal places of --txnE. Default: [%default]"),
	make_option(c("--outDir"), type="character", default=".", help = "Output Directory. Default: [%default]"),
	make_option(c("--includeHOMD"), type="logical", default=FALSE, help="If FALSE, then exclude HOMD state. Useful when using large bins (e.g. 1Mb). Default: [%default]"),
	make_option(c("--genomeBuild"), type="character", default="hg19", help="Geome build. Default: [%default]"),
	make_option(c("--genomeStyle"), type = "character", default = "NCBI", help = "NCBI or UCSC chromosome naming convention; use UCSC if desired output is to have \"chr\" string. [Default: %default]"),
	make_option(c("--libdir"), type = "character", default=NULL, help = "Script library path. Usually exclude this argument unless custom modifications have been made to the uMaC R package code and the user would like to source those R files. Default: [%default]")
)

parseobj <- OptionParser(option_list=option_list)
opt <- parse_args(parseobj)
print(opt)
options(scipen=0, stringsAsFactors=F)

##############
#library(mclust);
library(HMMcopy);
library(GenomicRanges)
library(GenomeInfoDb)

options(stringsAsFactors=FALSE)
options(bitmapType='cairo')

tumour_file <- opt$BED;
id <- opt$id
normal <- eval(parse(text = opt$normal))
ploidy <- eval(parse(text = opt$ploidy))
normal_panel <- opt$normalPanel

altFracThreshold <- opt$altFracThreshold
chrs <- eval(parse(text = opt$chrs))
chrTrain <- eval(parse(text=opt$chrTrain))
coverage <- opt$coverage
estimateScPrevalence <- opt$estimateScPrevalence
lambda <- eval(parse(text = opt$lambda))
maxCN <- opt$maxCN
maxFracCNASubclone <- opt$maxFracCNASubclone
maxFracGenomeSubclone <- opt$maxFracGenomeSubclone
minSegmentBins <- opt$minSegmentBins
plotYLim <- eval(parse(text=opt$plotYLim))
txnE <- opt$txnE
txnStrength <- opt$txnStrength

genomeBuild <- opt$genomeBuild
genomeStyle <- opt$genomeStyle
includeHOMD <- as.logical(opt$includeHOMD)

outDir <- opt$outDir
libdir <- opt$libdir

#######

gender <- NULL
#outImage <- paste0(outDir,"/", id,"/",id,".RData")

############
#

if (!is.null(libdir)){
	source(paste0(libdir,"/utils.R"))
	source(paste0(libdir,"/segmentation.R"))
	source(paste0(libdir,"/EM.R"))
	source(paste0(libdir,"/output.R"))
	source(paste0(libdir,"/plotting.R"))
} else {
	library(ichorCNA)
}

#####

#seqinfo <- getSeqInfo(genomeBuild, genomeStyle, chrTrain)
seqinfo <- getSeqInfo(genomeBuild, genomeStyle)
######################################################################################
##loading input and filter
######################################################################################

if (is.null(tumour_file) || tumour_file == "None" || tumour_file == "NULL") {
	stop("BED file is required");
} else {
	input <- read.table(tumour_file);
}

tumour_reads <- GenomicRanges::GRanges(ranges = IRanges(start = input$V2, end = input$V3),
				seqnames = factor(input$V1,levels=c(1:22)), reads = input$V6, normal = input$V7, rate = input$V8, noise=input$V4, noise.nor = input$V5)#, noise3 = input$V8)

tumour_reads$normal.nor <- tumour_reads$normal/max(tumour_reads$normal);

tumour_TFx <- sum(tumour_reads$reads)/(sum(tumour_reads$reads)+sum(tumour_reads$normal));

if(tumour_TFx<0.025){
	normal <- c(0.7,0.8,0.9,0.91);
	ploidy <- 2;
}else if(tumour_TFx<0.05){
	normal <- c(0.6,0.7,0.8,0.9,0.91);
	ploidy <- 2;
}else if (tumour_TFx<0.1){
	normal <- c(0.4,0.5,0.6,0.7,0.8);
}else{
	normal <- c(0.3,0.4,0.5,0.6,0.7);
}

tumour_counts <- list()
tumour_copy <- list()

## create output directories  ##
dir.create(paste0(outDir, "/", id, "/"), recursive = TRUE)
#########################
##filter and can use GMM cluster
##########################
sumNormal <- sum(tumour_reads$normal);
sumReads <- sum(tumour_reads$reads);
column <- length(tumour_reads$reads)

linearSlope <- sumNormal/sumReads;
linearY <- sumNormal/column;
linearX <- sumReads/column;

tumour_reads$valid <- TRUE
tumour_reads$ideal <- TRUE
tumour_reads$valid[tumour_reads$reads <=100  | tumour_reads$gc < 0] <- FALSE  #500
tumour_reads$ideal[tumour_reads$reads <=100  | tumour_reads$gc < 0] <- FALSE  #500


#tumour_reads$cor.normal <- tumour_reads$reads / (tumour_reads$normal/linearSlope);

routlier <- 0.01
range <- quantile(tumour_reads$reads[tumour_reads$valid], prob = c(routlier, 1 - routlier), na.rm = TRUE)
rangeNormal <- quantile(tumour_reads$normal[tumour_reads$valid], prob = c(routlier, 1 - routlier), na.rm = TRUE)
tumour_reads$ideal [!tumour_reads$valid | tumour_reads$reads <= range[1] |
tumour_reads$reads > range[2] | tumour_reads$normal <= rangeNormal[1] | tumour_reads$normal > rangeNormal[2] ] <- FALSE  ## | tumour_reads$gc < domain[1] | tumour_reads$gc > domain[2]

#plot(tumour_reads$normal[tumour_reads$ideal],tumour_reads$reads[tumour_reads$ideal])

###
samplesize=50000
set <- which(tumour_reads$ideal)
select <- sample(set, min(length(set), samplesize))
rough = loess(tumour_reads$reads[select] ~ tumour_reads$normal.nor[select], span = 0.03)

i <- seq(0, 1, by = 0.001)
final = loess(predict(rough, i) ~ i, span = 0.3)
tumour_reads$cor.normal <- tumour_reads$reads / predict(final, tumour_reads$normal.nor) ## gc divid into 1000, and predict i.

if(0){
	plot(tumour_reads$normal.nor,tumour_reads$reads,col ="red")

	for (i in seq(0.01, 1, length = 100)) {
		loessMod10 <- loess(tumour_reads$reads ~ tumour_reads$normal.nor, span=i)
		smoothed10 <- predict(loessMod10) 
		lines(tumour_reads$normal.nor,smoothed10, col =  gray(i), lwd=1.5, type = "pch")
		Sys.sleep(0.15)
	}
}


tumour_reads$cor.map <- tumour_reads$cor.normal

coutlier <- 0.01
range <- quantile(tumour_reads$cor.normal[which(tumour_reads$valid) ], prob = c(0, 1 - coutlier), na.rm = TRUE)
set <- which(tumour_reads$cor.normal < range[2])
select <- sample(set, min(length(set), samplesize))


final = approxfun(lowess(tumour_reads$noise[select]/tumour_reads$noise.nor[select], tumour_reads$cor.normal[select],f=0.5))
tumour_reads$cor.map <- tumour_reads$cor.normal / final(tumour_reads$noise/tumour_reads$noise.nor)

 #tumour_reads$cor.map <- tumour_reads$cor.normal

pdf(file = paste0(outDir,"/",id,"/",id,".reads.pdf"), width=8, height=8);
plot(tumour_reads$reads[which(tumour_reads$ideal)],tumour_reads$normal[which(tumour_reads$ideal)]);
dev.off();
#plot(tumour_reads$reads,tumour_reads$cor.normal)

tumour_reads$copy <- tumour_reads$cor.map
tumour_reads$copy[tumour_reads$copy <= 0] = NA
tumour_reads$copy <- log(tumour_reads$copy, 2)

###########################################################################
### OUTPUT FILE ###
### PUTTING TOGETHER THE COLUMNS IN THE OUTPUT ###
outMat <- as.data.frame(tumour_reads)
#outMat <- outMat[,c(1,2,3,12)]
outMat <- outMat[,c("seqnames","start","end","copy")]
colnames(outMat) <- c("chr","start","end","log2_TNratio_corrected")
outFile <- paste0(outDir,"/",id,"/",id,".correctedDepth.txt")
message(paste("Outputting to:", outFile))
write.table(outMat, file=outFile, row.names=F, col.names=T, quote=F, sep="\t")

########################################################################
##gender and use to identify in the furture
#########################################################################
gender <- "unknown" # chrX is not provided
chrYCov <- NA
chrXMedian <- NULL
gender = list(gender=gender, chrYCovRatio=chrYCov, chrXMedian=chrXMedian)

counts <- list(counts = tumour_reads, gender = gender)

################
#########
#tumour_copy[[id]] <- counts$counts #as(counts$counts, "GRanges")
gender <- counts$gender

normal_copy <- NULL

### DETERMINE GENDER ###
## if normal file not given, use chrY, else use chrX
#message("Determining gender...", appendLF = FALSE)
gender.mismatch <- FALSE
if (!is.null(normal_copy)){
	if (gender$gender != gender.normal$gender){ #use tumour # use normal if given
		# check if normal is same gender as tumour
		gender.mismatch <- TRUE
	}
}
#message("Gender ", gender$gender)

#########
#################
tumour_copy[[id]] <- counts$counts;

#
if(is.null(normal_panel) || normal_panel == "None" || normal_panel == "NULL"){
	} else {
	panel <- readRDS(normal_panel);
	tumour_copy[[id]]$copy <- tumour_copy[[id]]$copy - panel$Median;
}

valid <- tumour_copy[[id]]$valid

### RUN HMM ###
## store the results for different normal and ploidy solutions ##
ptmTotalSolutions <- proc.time() # start total timer
results <- list()
loglik <- as.data.frame(matrix(NA, nrow = length(normal) * length(ploidy), ncol = 7, 
				dimnames = list(c(), c("init", "n_est", "phi_est", "BIC", 
				"Frac_genome_subclonal", "Frac_CNA_subclonal", "loglik"))))


counter <- 1
compNames <- rep(NA, nrow(loglik))
mainName <- rep(NA, length(normal) * length(ploidy))
numSamples =1;


#libdir = "/scratch/cgg/weiq1/result/jinliang/ichorCNA/results/ichorCNA-master/R"
#source(paste0(libdir,"/utils.R"))
	#source(paste0(libdir,"/segmentation.R"))
	#source(paste0(libdir,"/EM.R"))
	#source(paste0(libdir,"/output.R"))
	#source(paste0(libdir,"/plotting.R"))


#### restart for purity and ploidy values ####
for (n in normal){
	for (p in ploidy){
		if (n == 0.95 & p != 2) {
			next
		}

		#logR <- as.data.frame(lapply(tumour_copy, "[[", "copy")) # NEED TO EXCLUDE CHR X #
		logR <- as.data.frame(lapply(tumour_copy, function(x) { x$copy }))
		param <- getDefaultParameters(logR[valid, , drop=F], maxCN = maxCN, includeHOMD = FALSE, 
				ct.sc=c(1,3), ploidy = floor(p), e = txnE, e.same = 50, strength=txnStrength)  ## 

##########
		param$phi_0 <- rep(p, numSamples)
		param$n_0 <- rep(n, numSamples)

		if (is.null(lambda)){
				logR.var <- 1 / ((apply(logR, 2, sd, na.rm = TRUE) / sqrt(length(param$ct))) ^ 2)
				param$lambda <- rep(logR.var, length(param$ct))
				param$lambda[param$ct %in% c(2)] <- logR.var 
				param$lambda[param$ct %in% c(1,3)] <- logR.var 
				param$lambda[param$ct >= 4] <- logR.var / 5
				param$lambda[param$ct == max(param$ct)] <- logR.var / 15
				param$lambda[param$ct.sc.status] <- logR.var / 10
		}

		param$alphaLambda <- rep(3, length(param$ct)) 

################ RUN HMM ####################
#############################################

		hmmResults.cor <- HMMsegment(tumour_copy, valid, dataType = "copy", 
					param = param, chrTrain = chrTrain, maxiter = 50,
					estimateNormal = as.logical('True'), estimatePloidy = as.logical('True'),
					estimateSubclone = estimateScPrevalence, verbose = TRUE)
		s = 1

		iter <- hmmResults.cor$results$iter
		id <- names(hmmResults.cor$cna)[s]

		correctedResults <- correctIntegerCN(cn = hmmResults.cor$cna[[s]],
		segs = hmmResults.cor$results$segs[[s]], 
		purity = 1 - hmmResults.cor$results$n[s, iter], ploidy = hmmResults.cor$results$phi[s, iter],
		cellPrev = 1 - hmmResults.cor$results$sp[s, iter], 
		maxCNtoCorrect.autosomes = maxCN, maxCNtoCorrect.X = maxCN, minPurityToCorrect = 0.03, 
		gender = gender$gender, chrs = chrs, correctHOMD = includeHOMD)

		hmmResults.cor$results$segs[[s]] <- correctedResults$segs
		hmmResults.cor$cna[[s]] <- correctedResults$cn
		## convert full diploid solution (of chrs to train) to have 1.0 normal or 0.0 purity
		## check if there is an altered segment that has at least a minimum # of bins
		segsS <- hmmResults.cor$results$segs[[s]]
		segsS <- segsS[segsS$chr %in% chrTrain, ]

		segAltInd <- which(segsS$event != "NEUT")

		maxBinLength = -Inf
		if (sum(segAltInd) > 0){
			maxInd <- which.max(segsS$end[segAltInd] - segsS$start[segAltInd] + 1)
			maxSegRD <- GRanges(seqnames=segsS$chr[segAltInd[maxInd]], 
								ranges=IRanges(start=segsS$start[segAltInd[maxInd]], end=segsS$end[segAltInd[maxInd]]))
			hits <- findOverlaps(query=maxSegRD, subject=tumour_copy[[s]][valid, ])
			maxBinLength <- length(subjectHits(hits))  ## need to ajust
		}


## check if there are proportion of total bins altered 
		# if segment size smaller than minSegmentBins, but altFrac > altFracThreshold, then still estimate TF
		cnaS <- hmmResults.cor$cna[[s]]
		altInd <- cnaS[cnaS$chr %in% chrTrain, "event"] == "NEUT"
		altFrac <- sum(!altInd, na.rm=TRUE) / length(altInd)   ## indel/normal,>0.05  need to ajust

		if ((maxBinLength <= minSegmentBins) & (altFrac <= altFracThreshold)) #altFracThreshold 0.05
		{
			hmmResults.cor$results$n[s, iter] <- 1.0
		}

		## plot solution ##
		outPlotFile <- paste0(outDir, "/", id, "/", id, "_genomeWide_", "n", n, "-p", p)
		mainName[counter] <- paste0(id, ", n: ", n, ", p: ", p, ", log likelihood: ", signif(hmmResults.cor$results$loglik[hmmResults.cor$results$iter], digits = 4))

		#plotGWSolution(hmmResults.cor, s=s, outPlotFile=outPlotFile, plotFileType="pdf", 
		         #plotYLim=plotYLim, estimateScPrevalence=estimateScPrevalence, main=mainName[counter])

		plotGWSolution(hmmResults.cor, s=s, outPlotFile=outPlotFile, plotFileType="pdf", 
					logR.column = "logR", call.column = "Corrected_Call",
					plotYLim=plotYLim, estimateScPrevalence=estimateScPrevalence, seqinfo=seqinfo, main=mainName[counter])

		iter <- hmmResults.cor$results$iter
		results[[counter]] <- hmmResults.cor
		loglik[counter, "loglik"] <- signif(hmmResults.cor$results$loglik[iter], digits = 4)
		subClonalBinCount <- unlist(lapply(hmmResults.cor$cna, function(x){ sum(x$subclone.status) }))
		fracGenomeSub <- subClonalBinCount / unlist(lapply(hmmResults.cor$cna, function(x){ nrow(x) }))
		fracAltSub <- subClonalBinCount / unlist(lapply(hmmResults.cor$cna, function(x){ sum(x$copy.number != 2) }))
		fracAltSub <- lapply(fracAltSub, function(x){if (is.na(x)){0}else{x}})
		loglik[counter, "Frac_genome_subclonal"] <- paste0(signif(fracGenomeSub, digits=2), collapse=",")
		loglik[counter, "Frac_CNA_subclonal"] <- paste0(signif(as.numeric(fracAltSub), digits=2), collapse=",")
		loglik[counter, "init"] <- paste0("n", n, "-p", p)
		loglik[counter, "n_est"] <- paste(signif(hmmResults.cor$results$n[, iter], digits = 2), collapse = ",")
		loglik[counter, "phi_est"] <- paste(signif(hmmResults.cor$results$phi[, iter], digits = 4), collapse = ",")

		counter <- counter + 1
	}
}


## get total time for all solutions ##
elapsedTimeSolutions <- proc.time() - ptmTotalSolutions
message("Total HMM Runtime: ", format(elapsedTimeSolutions[3] / 60, digits = 2), " min.")

### SAVE R IMAGE ###
#save.image(outImage)
#save(tumour_copy, results, loglik, file=paste0(outDir,"/",id,".RData"))

### SELECT SOLUTION WITH LARGEST LIKELIHOOD ###
if (estimateScPrevalence){ ## sort but excluding solutions with too large % subclonal 
	fracInd <- which(loglik[, "Frac_CNA_subclonal"] <= maxFracCNASubclone & 
							loglik[, "Frac_genome_subclonal"] <= maxFracGenomeSubclone)
	if (length(fracInd) > 0){ ## if there is a solution satisfying % subclonal
		ind <- fracInd[order(loglik[fracInd, "loglik"], decreasing=TRUE)]
	}else{ # otherwise just take largest likelihood
		ind <- order(as.numeric(loglik[, "loglik"]), decreasing=TRUE) 
	}
}else{#sort by likelihood only
	ind <- order(as.numeric(loglik[, "loglik"]), decreasing=TRUE) 
}

#new loop by order of solutions (ind)
outPlotFile <- paste0(outDir, "/", id, "/", id, "_genomeWide_all_sols")
for(i in 1:length(ind)) {
	hmmResults.cor <- results[[ind[i]]]
	turnDevOff <- FALSE
	turnDevOn <- FALSE
	if (i == 1){
		turnDevOn <- TRUE
	}
	if (i == length(ind)){
		turnDevOff <- TRUE
	}
	#plotGWSolution(hmmResults.cor, s=s, outPlotFile=outPlotFile, plotFileType="pdf", 
	#		plotYLim=plotYLim, estimateScPrevalence=estimateScPrevalence, 
	#		turnDevOn = turnDevOn, turnDevOff = turnDevOff, main=mainName[ind[i]])

	plotGWSolution(hmmResults.cor, s=s, outPlotFile=outPlotFile, plotFileType="pdf", 
				logR.column = "logR", call.column = "Corrected_Call",
				plotYLim=plotYLim, estimateScPrevalence=estimateScPrevalence, 
				seqinfo = seqinfo,
				turnDevOn = turnDevOn, turnDevOff = turnDevOff, main=mainName[ind[i]])
}


###############
##
hmmResults.cor <- results[[ind[1]]]
hmmResults.cor$results$loglik <- as.data.frame(loglik)
hmmResults.cor$results$gender <- gender$gender
hmmResults.cor$results$chrYCov <- gender$chrYCovRatio
hmmResults.cor$results$chrXMedian <- gender$chrXMedian
hmmResults.cor$results$coverage <- coverage

###Sample1.seg.txt ; Sample1.seg 
outputHMM(cna = hmmResults.cor$cna, segs = hmmResults.cor$results$segs, 
			results = hmmResults.cor$results, patientID = id, outDir=paste0(outDir,"/",id))
#############
##Sample1.params.txt
outFile <- paste0(outDir, "/", id, "/",id, ".params.txt")
outputParametersToFile(hmmResults.cor, file = outFile)

## plot solutions for all samples 
#plotSolutions(hmmResults.cor, tumour_copy, chrs, outDir, numSamples=numSamples,
#		plotFileType="pdf", plotYLim=plotYLim, 
#		estimateScPrevalence=estimateScPrevalence, maxCN=maxCN)

plotSolutions(hmmResults.cor, tumour_copy, chrs, outDir, numSamples=numSamples,
			logR.column = "logR", call.column = "Corrected_Call",
			plotFileType = "pdf", plotYLim=plotYLim, seqinfo = seqinfo,
			estimateScPrevalence=estimateScPrevalence, maxCN=maxCN)


