#sampling.R

# source("~/Dropbox/code/R/common_src/occFns.R")
# source("~/Dropbox/code/R/common_src/utils_marcot.R")
# source("~/Dropbox/code/R/common_src/taxonomicEv.R")

subsampleInterval_SQS <- function(thisInterval, thisOccs, occDates, quota=0.4) {
	if (quota <= 0.0 | quota > 1.0) stop("Value of quota for SQS should be between 0 and 1")
	intOccs <- thisOccs[occDates >= thisInterval["ageTop"] & occDates < thisInterval["ageBase"],]
	intCols <- unique(intOccs$collection_no)
	if (length(intCols)==0) { return(NA)
	} else if (length(intCols)==1) { return(intOccs$occurrence_no)
	} else {
		if (nrow(intOccs)>2) {
	
			noccs <- table(intOccs$taxon) 

			# ni <- sum(noccs==1)							# number of single occurrence taxa, by definition, this is the number of their occurrences
			p1 <- sum(noccs[sapply(names(noccs), function(x) length(unique(intOccs$occurrence.reference_no[intOccs$taxon==x]))) == 1])		# number of occurrences of single publication taxa (within the interval)
			# p1 <- sum(noccs[sapply(names(noccs), function(x) length(unique(thisOccs$occurrence.reference_no[thisOccs$taxon==x])))==1])	# number of occurrences of single publication taxa

			n1 <- max(noccs) 								# number of occurrences of the dominant taxon
			dominantTaxon <- names(noccs[noccs==n1])		# name(s) of the dominant taxon
			# if (length(dominantTaxon) > 1) { 
				# freq <- noccs/sum(noccs)
				# n1 <- 0
			# } else {
				freq <- noccs/sum(noccs[!names(noccs) %in% dominantTaxon])	# the most common taxon is ignored in frequency calculations
				freq[dominantTaxon] <- 0.0									# the most common taxon is ignored in frequency calculations
			# }
	
			maxCol <- intCols[which.max(sapply(intCols, function(x) sum(intOccs$collection_no==x)))]
			tmax <- sum(names(noccs[noccs==1]) %in% intOccs$taxon[intOccs$collection_no %in% maxCol])	#number of taxa ONLY found in collection(s) with maximum number of occurrences
	
			u <- (sum(noccs) - p1 - n1 + tmax) / (sum(noccs) - n1)
			quota <- quota/u
	
			thisCoverage=0.0
			subsample <- intCols[sample(length(intCols), length(intCols), replace=FALSE)]
			i <- 0
			while(thisCoverage < quota & i <= length(subsample)) {
				i <- i+1
				thisCoverage <- sum(freq[names(freq) %in% intOccs$taxon[intOccs$collection_no %in% subsample[1:i]]])
			}
			subsample <- subsample[1:i]
			list(coverage=u, occs=intOccs$occurrence_no[intOccs$collection_no %in% subsample])
		} else return(NA)
	}
}

subsampleInterval_R<-function(thisInterval, thisOccs, quota, quota.mismatch=FALSE) {
	intOccs <- thisOccs[thisOccs$ma_rand>thisInterval["ageTop"] & thisOccs$ma_rand<thisInterval["ageBase"],]
	if (nrow(intOccs) < quota) {
		if (quota.mismatch) return(NA) else return(intOccs$occurrence_no)
	} 
	sample(intOccs$occurrence_no, quota)
}

subsampleInterval_UW<-function(thisInterval, thisOccs, quota, quota.mismatch=FALSE) {
	intOccs<-thisOccs[thisOccs$ma_rand>thisInterval["ageTop"] & thisOccs$ma_rand<thisInterval["ageBase"],]
	intCols<-unique(intOccs$collection_no)
	if (length(intCols)<quota) {
		if (quota.mismatch) return(NA) else return(intOccs$occurrence_no)
	} 
	intOccs[intOccs$collection_no%in%sample(intCols, quota),]$occurrence_no
}

getWeightedFromOccList<-function(thisOccs, intList, method=c("OW", "O2W")) {
	occList<-lapply(intList, function(x, y) { thisOccs[thisOccs$ma_rand>x["ageTop"] & thisOccs$ma_rand<x["ageBase"],] }, thisOccs)
	if (method=="OW") { sapply(occList, function(x) { sum(table(x$collection_no)) } )
	} else sapply(occList, function(x) { sum(table(x$collection_no)^2) } )
}

subsampleInterval_Weighted<-function(thisInterval, thisOccs, method=c("OW", "O2W"), quota, quota.mismatch=FALSE) {
	intOccs<-thisOccs[thisOccs$ma_rand>thisInterval["ageTop"] & thisOccs$ma_rand<thisInterval["ageBase"],]
	intCols<-unique(intOccs$collection_no)
	subsample<-sample(intCols, length(intCols))
	i<-0
	noccs<-0
	while(noccs<quota & i<length(intCols)) {
		i<-i+1
		if (method=="OW") { noccs<-sum(intOccs$collection_no%in%subsample[1:i])
		} else if (method=="O2W") { noccs<-sum(intOccs$collection_no%in%subsample[1:i])^2 }
	}
	if (quota.mismatch & i==length(intCols)) return(NA)
	intOccs[intOccs$collection_no%in%subsample[1:i],]$occurrence_no
}

oneRep<-function(rep=1, occs, intList, method=c("R", "UW", "OW", "O2W", "SQ", "G"), do.G=FALSE, thisQuota=NULL, thisQuota.G=NULL, auto.quota=TRUE, age.determination=c("midpoint", "random"), doNotSample=NULL, do.parallel=FALSE) {
	method <- match.arg(method)
	age.determination <- match.arg(age.determination)
	if (method=="G") do.G <- FALSE
	if (method=="G" | do.G) source("~/Dropbox/Code/R/common_src/collectionGeography.R")

	thisOccs <- cbind(occs, ma_rand= getOccurrenceDatesFromOccs(occs, age.determination))					# redate collections every rep
	# thisOccs <- cbind(thisOccs, ma_mid=getCollectionDatesFromOccs(occs, age.determination="midpoint"))
	
	if (do.G) {
		if (auto.quota) {
			ACH <- sapply(intList, function(x) getACHFromOccsOneInterval(thisOccs[thisOccs$ma_rand > x["ageTop"] & thisOccs$ma_rand < x["ageBase"],]))
			thisQuota.G <- min(ACH[ACH>0][-doNotSample], na.rm=TRUE)
		}
		if (do.parallel) { thisSample <- mclapply(intList, subsampleInterval_G, thisOccs, method, quota=thisQuota.G, quota.mismatch=FALSE, mc.cores=detectCores()-2) 
		} else thisSample <- lapply(intList, subsampleInterval_G, thisOccs=thisOccs, quota=thisQuota.G, quota.mismatch=FALSE)
		thisOccs <- thisOccs[thisOccs$occurrence_no %in% unlist(thisSample),]
	}
	
	if (method=="SQ") { 
		if (do.parallel) { thisSample <- mclapply(intList, subsampleInterval_SQ, thisOccs, quota=thisQuota, mc.cores=detectCores()-2) 
		} else thisSample <- lapply(intList, subsampleInterval_SQ, thisOccs, quota=thisQuota)
	} else {
		if (auto.quota) {
			if (method=="R") {
				noccs <- sapply(intList, function(x, y) { sum(y$ma_rand>x["ageTop"] & y$ma_rand<x["ageBase"]) }, thisOccs)
				newQuota <- min(noccs[-doNotSample])
				obsRichness <- sapply(intList, function (x,y) { length(unique(y[y$ma_rand>x["ageTop"] & y$ma_rand<x["ageBase"],]$taxon)) }, thisOccs)
				if (newQuota < max(obsRichness[-doNotSample])) newQuota <- max(obsRichness[-doNotSample])
				print(paste("AutoQuota reset quota from", thisQuota, "to", newQuota))
				thisQuota <- newQuota
			} else if (method=="UW") {
				ncols <- sapply(intList, function(x, y) { length(unique(y[y$ma_rand>x["ageTop"] & y$ma_rand<x["ageBase"],]$collection_no)) }, thisOccs)
				newQuota <- min(ncols[-doNotSample])
				# noccs<-sapply(intList, function(x, y) { sum(y$ma_rand>x["ageTop"] & y$ma_rand<x["ageBase"]) }, thisOccs)
				# obsRichness<-sapply(intList, function (x,y) { length(unique(y[y$ma_rand>x["ageTop"] & y$ma_rand<x["ageBase"],]$taxon)) }, thisOccs)
				
				
				# min(ncols[which(noccs[-doNotSample]>max(obsRichness[-doNotSample]))])
				print(paste("AutoQuota reset quota from", thisQuota, "to", newQuota))
				thisQuota <- newQuota
				# thisQuota<-min(ncols[which(noccs[-doNotSample]>max(obsRichness[-doNotSample]))])
				# if (thisQuota<max(obsRichness)) thisQuota<-max(obsRichness)
			} else if (method=="OW" | method=="O2W") {
				weighted <- getWeightedFromOccList(thisOccs, intList, method)
				print(paste("AutoQuota reset quota from", thisQuota, "to", min(weighted)))
				thisQuota <- min(weighted)
			} else if (method=="G") {
				ACH <- sapply(intList, function(x) getACHFromOccsOneInterval(thisOccs[thisOccs$ma_rand > x["ageTop"] & thisOccs$ma_rand < x["ageBase"],]))
				thisQuota.G <- min(ACH[ACH>0][-doNotSample], na.rm=TRUE)
			}

		}
		if (method=="R") { 
			if (do.parallel) { thisSample <- mclapply(intList, subsampleInterval_R, thisOccs, quota=thisQuota, quota.mismatch=FALSE, mc.cores=detectCores()-2)  
			} else thisSample <- lapply(intList, subsampleInterval_R, thisOccs, quota=thisQuota, quota.mismatch=FALSE)
		} else if (method=="UW") { 
			if (do.parallel) { thisSample <- mclapply(intList, subsampleInterval_UW, thisOccs, quota=thisQuota, quota.mismatch=FALSE, mc.cores=detectCores()-2) 
			} else thisSample <- lapply(intList, subsampleInterval_UW, thisOccs, quota=thisQuota, quota.mismatch=FALSE)
		} else if (method=="OW" | method=="O2W") { 
			if (do.parallel) { thisSample <- mclapply(intList, subsampleInterval_Weighted, thisOccs, method, quota=thisQuota, quota.mismatch=FALSE, mc.cores=detectCores()-2) 
			} else thisSample <- lapply(intList, subsampleInterval_Weighted, thisOccs, method, quota=thisQuota, quota.mismatch=FALSE)
		} else if (method=="G") {
			if (do.parallel) { thisSample <- mclapply(intList, subsampleInterval_G, thisOccs, method, quota=thisQuota.G, quota.mismatch=FALSE, mc.cores=detectCores()-2) 
			} else thisSample <- lapply(intList, subsampleInterval_G, thisOccs=thisOccs, quota=thisQuota.G, quota.mismatch=FALSE)
		}
	}
	if (any(doNotSample>0)) for (i in seq_along(doNotSample)) thisSample[[doNotSample[i]]] <- subsampleInterval_R(intList[[doNotSample[i]]], thisOccs, 10e10, quota.mismatch=FALSE)
	thisSample
}

getSubsampledOccs <- function(occs, intList, method=c("R", "UW", "OW", "O2W", "SQ", "G"), do.G=FALSE, quota, quota.G=NULL, auto.quota=FALSE, reps=100,age.determination=c("midpoint","random"), doNotSample=NULL, restrict.col=NULL, restrict.class=NULL, do.parallel=FALSE) {
	method <- match.arg(method)
	age.determination <- match.arg(age.determination)

	# subsampledOccs <- replicate(reps, oneRep(rep=1, occs, intervals, method, quota, auto.quota, age.determination, doNotSample, do.parallel=FALSE), simplify=FALSE)
	if (do.parallel) { mclapply(1:reps, oneRep, occs=occs, intList=intList, method=method, do.G=do.G, thisQuota=quota, thisQuota.G=quota.G, auto.quota=auto.quota, age.determination=age.determination, doNotSample=doNotSample, do.parallel=FALSE, mc.cores=detectCores()-2) ## do parallel is set to FALSE to prevent parallel processes within parallel processes
	} else lapply(1:reps, oneRep, occs=occs, intList=intList, method=method, do.G=do.G, thisQuota=quota, thisQuota.G=quota.G, auto.quota=auto.quota, age.determination=age.determination, doNotSample=doNotSample, do.parallel=FALSE)
}

doResamplingAnalysis <- function(occs, method=c("R", "UW", "OW", "O2W", "SQ", "G"), do.G=FALSE, quota=NULL, quota.G=NULL, auto.quota=FALSE, reps=100, intervals, richness.method=c("SIB", "RT", "BC", "TT"), age.determination=c("midpoint","random"), doNotSample=NULL, logRichness=FALSE, restrict.col=NULL, restrict.class=NULL, do.parallel=FALSE) {
	method <- match.arg(method)
	age.determination <- match.arg(age.determination)
	richness.method <- match.arg(richness.method)
	if (is.null(doNotSample)) doNotSample <- -seq_len(nrow(intervals))
	if (do.parallel) require(parallel)
	
	intList <- listifyMatrixByRow(intervals)
	subsampledOccs <- getSubsampledOccs(occs=occs, intList=intList, method=method, quota=quota, quota.G=quota.G, auto.quota=auto.quota, reps=reps, age.determination=age.determination, doNotSample=doNotSample, restrict.col=restrict.col, restrict.class=restrict.class, do.parallel=do.parallel)
	# occs<-occs[match(unique(unlist(subsampledOccs)), occs$occurrence_no),]  # reduces the occ List down to only those in the complete subsample

	richnessMat <- getRichnessBoxFromOccBox(occBox=subsampledOccs, occs, method=richness.method, restrict.col, restrict.class, do.parallel)
	
	richnessMat <- t(apply(richnessMat, 1, quantile, probs=c(0.025, 0.50, 0.975), na.rm=TRUE))
	# richnessMat<-t(app ly(richnessMat, 1,  function(x) { m<-mean(x, na.rm=TRUE); se<-sd(x, na.rm=TRUE)/sqrt(sum(is.finite(x))); c(m-se, m, m+se) } ))
	# richnessMat<-t(apply(richnessMat, 1,  function(x) { m<-mean(x, na.rm=TRUE); se<-sd(x, na.rm=TRUE); c(m-se, m, m+se) } ))

	occs$ma_mid <- getOccurrenceDatesFromOccs(occs, age.determination="midpoint")
	if (do.parallel) { fvOccs <- mclapply(intList, function(x) { occs[rowMeans(occs[,c("ma_max", "ma_min")]) < x["ageBase"] & rowMeans(occs[,c("ma_max", "ma_min")]) > x["ageTop"],]$occurrence_no }, mc.cores=detectCores()-2) 
	} else fvOccs <- lapply(intList, function(x) { occs[rowMeans(occs[,c("ma_max", "ma_min")]) < x["ageBase"] & rowMeans(occs[,c("ma_max", "ma_min")]) > x["ageTop"],]$occurrence_no }) 
	fvRichness <- getSIBRichnessFromOneOccList(fvOccs, occs)
	fvRichnessRT <- getRangeThroughRichnessFromOneOccList(fvOccs, occs)
	
		par(mfrow=c(3,3))
	
		# plot(rowMeans(intervals), richnessMat[,3], type="n", xlim=c(max(intervals), min(intervals)), ylim=c(0, max(richnessMat)), xlab="Time (Ma)", ylab="Taxonomic Richness")
		# overlayCzTimescale()
		if (logRichness) thisY<-log(fvRichnessRT[-doNotSample]) else thisY<-fvRichnessRT[-doNotSample]
		plot(rowMeans(intervals)[-doNotSample], thisY, type="n", xlim=c(max(intervals), min(intervals)), ylim=c(0, max(thisY, na.rm=TRUE)), xlab="Time (Ma)", ylab="Taxonomic Richness", main="Taxonomic Richness", ylog=logRichness)
		overlayCzTimescale()
		if (!auto.quota & method=="R") abline(h=quota, lty=3, col="red")
	
		# if (logRichness) thisY<-log(richnessMat[-doNotSample,2]) else thisY<-richnessMat[-doNotSample,2]
		# text(x=rowMeans(intervals)[-doNotSample], y=thisY, labels=paste(rowMeans(intervals)[-doNotSample], "Ma", sep=""), pos=3, col="gray50", cex=0.5)
		thisY <- fvRichnessRT[-doNotSample]
		if (logRichness) thisY <- log(thisY)
		lines(rowMeans(intervals)[-doNotSample], thisY, lwd=2, lty=1, type="o", pch=21, col="black", bg="gray50")	# face-value range-through richness
		# thisY <- fvRichness[-doNotSample]
		# if (logRichness) thisY <- log(thisY)
		# lines(rowMeans(intervals)[-doNotSample], thisY, lwd=2, lty=2, type="o", pch=21, col="gray25", bg="gray75")	# face-value sampled-in-bin richness
		thisY <- getThreeTimerRichnessFromOneOccList(fvOccs, occs, doNotSample=doNotSample)[-doNotSample]
		if (logRichness) thisY <- log(thisY)
		lines(rowMeans(intervals)[-doNotSample], thisY, lwd=2, lty=3, type="o", pch=21, col="black", bg="gray50")	# rescaled sampled-in-bin richness

		if (length(subsampledOccs)>1) {
			thisY <- c(richnessMat[-doNotSample,1], richnessMat[-doNotSample,3][nrow(richnessMat[-doNotSample,]):1])
			if (logRichness) thisY <- log(thisY)
			polygon(x=c(rowMeans(intervals[-doNotSample,]), rowMeans(intervals[-doNotSample,])[nrow(intervals[-doNotSample,]):1]), y=thisY, col=alphaColor("darkgreen", 0.25), lwd=0.5, ylog=logRichness)
		}
		thisY <- richnessMat[-doNotSample,2]
		if (logRichness) thisY <- log(thisY)
		lines(rowMeans(intervals)[-doNotSample], thisY, lwd=3, col="darkgreen", type="o", pch=21, bg="green")		# resampled richness
	
	if (do.parallel) { mat <- mclapply(subsampledOccs, fourClassesFromOccs, occs, restrict.col, restrict.class, mc.cores=detectCores()-2) 
	} else mat <- lapply(subsampledOccs, fourClassesFromOccs, occs, restrict.col, restrict.class) 

	if (do.parallel) pc <- simplify2array(mclapply(mat, cmp_PCrates, mc.cores=detectCores()-2)) else pc <- sapply(mat, cmp_PCrates, simplify="array")
	pc[!(is.finite(pc))] <- NA
	pcMat <- apply(pc, c(1,2), quantile, probs=c(0.025, 0.5, 0.975), na.rm=TRUE)
	# pcMat<-apply(pc, c(1,2), function(x) { m<-mean(x, na.rm=TRUE); se<-sd(x, na.rm=TRUE)/sqrt(sum(is.finite(x))); c(m-se, m, m+se) } )
	# pcMat<-apply(pc, c(1,2), function(x) { m<-mean(x, na.rm=TRUE); se<-sd(x, na.rm=TRUE); c(m-se, m, m+se) } )
	colnames(pcMat) <- rownames(intervals)
	
	# quartz(width=4, height=9)
	# quartz(width=10, height=6.25)
	# par(mfrow=c(2,1), mar=c(2,4,1,2))
	
		plot(cbind(rowMeans(intervals)[-doNotSample], rowMeans(intervals)[-doNotSample]), pcMat[3,-doNotSample,], type="n", xlim=c(max(intervals), min(intervals)), ylim=c(0, max(pcMat[3,-doNotSample,], na.rm=TRUE)), xlab="Time (Ma)", ylab="Per-Capita Origination Rate", main="")
		overlayCzTimescale()
		# text(rowMeans(intervals)[-doNotSample], pcMat[2,-doNotSample,"orig"], paste(rowMeans(intervals)[-doNotSample], "Ma", sep=""), pos=3, col="gray50", cex=0.5)
		if (length(subsampledOccs)>1) polygon(x=c(rowMeans(intervals)[2:(nrow(intervals)-1)], rowMeans(intervals)[(nrow(intervals)-1):2]), y=c(pcMat[1,2:(dim(pcMat)[2]-1),"orig"], pcMat[3,(dim(pcMat)[2]-1):2,"orig"]), col=alphaColor("dodgerblue4", 0.25), border="dodgerblue4", lwd=0.5)
		lines(rowMeans(intervals)[-doNotSample], pcMat[2,-doNotSample,"orig"], lwd=2, col="dodgerblue1")
	
	# par(mar=c(4,4,1,1))
		plot(cbind(rowMeans(intervals)[-doNotSample], rowMeans(intervals)[-doNotSample]), pcMat[3,-doNotSample,], type="n", xlim=c(max(intervals), min(intervals)), ylim=c(0, max(pcMat[3,-doNotSample,], na.rm=TRUE)), xlab="Time (Ma)", ylab="Per-Capita Extinction Rate", , main="")
		overlayCzTimescale()
		# text(rowMeans(intervals)[-doNotSample], pcMat[2,-doNotSample,"ext"], paste(rowMeans(intervals)[-doNotSample], "Ma", sep=""), pos=3, col="gray50", cex=0.5)
		if (length(subsampledOccs)>1) polygon(x=c(rowMeans(intervals)[2:(nrow(intervals)-1)], rowMeans(intervals)[(nrow(intervals)-1):2]), y=c(pcMat[1,2:(dim(pcMat)[2]-1),"ext"], pcMat[3,(dim(pcMat)[2]-1):2,"ext"]), col=alphaColor("firebrick4", 0.25), border="firebrick4", lwd=0.5)
		lines(rowMeans(intervals)[-doNotSample], pcMat[2,-doNotSample,"ext"], lwd=2, col="firebrick1")
		
		# plot(pcMat[2,(1:(nrow(intervals)-1)),1], pcMat[2,(2:nrow(intervals)),2], xlab="Per-Capita Origination Rate", ylab="Per-Capita Extinction Rate", col="darkviolet")
		# abline(lm(pcMat[2,(2:nrow(intervals)),2] ~ pcMat[2,(1:(nrow(intervals)-1)),1], na.action="na.omit"), col="darkviolet")
	
		par(mar=c(5, 5, 4, 5)) 
		plot(rowMeans(intervals), sapply(fvOccs, length), type="n", xlim=c(max(intervals), min(intervals)), main="No. Collections and Occurrences", xlab="Time (Ma)", ylab="", yaxt="n")
		overlayCzTimescale()
		# text(rowMeans(intervals), sapply(fvOccs, length), paste(rowMeans(intervals), "Ma", sep=""), pos=3, col="darkorange", cex=0.5)
		lines(rowMeans(intervals), sapply(fvOccs, length), lwd=4, lty=1, col="orange")
		Axis(side=2, col="orange", col.ticks="orange")
		title(ylab="No. Occurrences", col.lab="orange")
		if (!auto.quota & method=="R") {
			# abline(h=quota, lty=3)
			arrows(x0=par()$usr[2], y0=quota, x1=par()$usr[1], y1=quota, length=0.025*par()$fin[1], lty=3, col="red")
		}
		# if (length(subsampledOccs)>1) polygon(x=c(rowMeans(intervals)[2:(nrow(intervals)-1)], rowMeans(intervals)[(nrow(intervals)-1):2]), y=c(pcMat[1,2:(dim(pcMat)[2]-1),"ext"], pcMat[3,(dim(pcMat)[2]-1):2,"ext"]), col=rgb(col2rgb("firebrick4")[1], col2rgb("firebrick4")[2], col2rgb("firebrick4")[3], 25, maxColorValue=256), border="firebrick4")
		par(new=TRUE)
		plot(rowMeans(intervals), sapply(fvOccs, function(x) length(unique(occs[occs$occurrence_no%in%x,]$collection_no))), type="l", lwd=4, lty=1, new=TRUE, xlim=c(max(intervals), min(intervals)), xlab="", yaxt="n", ylab="", col="sienna")
		# text(rowMeans(intervals), sapply(fvOccs, function(x) length(unique(occs[occs$occurrence_no%in%x,]$collection_no))), paste(rowMeans(intervals), "Ma", sep=""), pos=3, col="sienna4", cex=0.5)
		if (!auto.quota & method=="UW") {
			# abline(h=quota, lty=3)
			arrows(x0=par()$usr[1], y0=quota, x1=par()$usr[2], y1=quota, length=0.025*par()$fin[1], lty=3, col="red")
		}
		Axis(side=4, col="sienna", col.ticks="sienna")
		mtext(text="No. Collections", side=4, line=2, col="sienna", cex=par()$cex)
		# lines(rowMeans(intervals), sapply(fvOccs, length), lwd=4, lty=3)
		par(new=FALSE, mar=c(5,4,4,2), mgp=c(3,1,0))
		
	if (do.parallel) pt<-simplify2array(mclapply(mat, cmp_PTrates, mc.cores=detectCores()-2)) else pt<-sapply(mat, cmp_PTrates, simplify="array")
	pt[!(is.finite(pt))]<-NA
	ptMat<-apply(pt, c(1,2), quantile, probs=c(0.025, 0.5, 0.975), na.rm=TRUE)
	# ptMat<-apply(pt, c(1,2), function(x) { m<-mean(x, na.rm=TRUE); se<-sd(x, na.rm=TRUE)/sqrt(sum(is.finite(x))); c(m-se, m, m+se) } )
	# ptMat<-apply(pt, c(1,2), function(x) { m<-mean(x, na.rm=TRUE); se<-sd(x, na.rm=TRUE); c(m-se, m, m+se) } )
	colnames(ptMat)<-rownames(intervals)

		plot(cbind(rowMeans(intervals)[-doNotSample], rowMeans(intervals)[-doNotSample]), ptMat[3,-doNotSample,1:2], type="n", xlim=c(max(intervals), min(intervals)), ylim=c(min(ptMat[1,-doNotSample,], na.rm=TRUE), max(ptMat[3,-doNotSample,], na.rm=TRUE)), xlab="Time (Ma)", ylab="Per-Taxon Origination Rate", main="Per-Taxon Origination Rate")
		overlayCzTimescale()
		text(rowMeans(intervals)[-doNotSample], ptMat[2,-doNotSample,"orig"], paste(rowMeans(intervals)[-doNotSample], "Ma", sep=""), pos=3, col="gray50", cex=0.5)
		if (length(subsampledOccs)>1) polygon(x=c(rowMeans(intervals), rowMeans(intervals)[nrow(intervals):1]), y=c(ptMat[1,,"orig"], ptMat[3,dim(ptMat)[2]:1,"orig"]), col=alphaColor("dodgerblue4", 0.25), border="dodgerblue4", lwd=0.5)
		lines(rowMeans(intervals)[-doNotSample], ptMat[2,-doNotSample,"orig"], lwd=2, col="dodgerblue1")
	
		plot(cbind(rowMeans(intervals)[-doNotSample], rowMeans(intervals)[-doNotSample]), ptMat[3,-doNotSample,1:2], type="n", xlim=c(max(intervals), min(intervals)), ylim=c(min(ptMat[1,-doNotSample,], na.rm=TRUE), max(ptMat[3,-doNotSample,], na.rm=TRUE)), xlab="Time (Ma)", ylab="Per-Taxon Extinction Rate", main="Per-Taxon Extinction Rate")
		overlayCzTimescale()
		text(rowMeans(intervals)[-doNotSample], ptMat[2,-doNotSample,"ext"], paste(rowMeans(intervals)[-doNotSample], "Ma", sep=""), pos=3, col="gray50", cex=0.5)
		if (length(subsampledOccs)>1) polygon(x=c(rowMeans(intervals), rowMeans(intervals)[nrow(intervals):1]), y=c(ptMat[1,,"ext"], ptMat[3,dim(ptMat)[2]:1,"ext"]), col=alphaColor("firebrick4", 0.25), border="firebrick4", lwd=0.5)
		lines(rowMeans(intervals)[-doNotSample], ptMat[2,-doNotSample,"ext"], lwd=2, col="firebrick1")
		
	if (do.parallel) { mat <- mclapply(subsampledOccs, alroyCountsFromOccs, occs, restrict.col, restrict.class, mc.cores=detectCores()-2) 
	} else mat <- lapply(subsampledOccs, alroyCountsFromOccs, occs, restrict.col, restrict.class) 
	alrates <- sapply(mat, cmp_3Trates, simplify="array")
	alrates[!(is.finite(alrates))] <- NA
	arMat <- apply(alrates, c(1,2), quantile, probs=c(0.025, 0.5, 0.975), na.rm=TRUE)
	# arMat<-apply(alrates, c(1,2), function(x) { m<-mean(x, na.rm=TRUE); se<-sd(x, na.rm=TRUE)/sqrt(sum(is.finite(x))); c(m-se, m, m+se) } )
	# arMat<-apply(alrates, c(1,2), function(x) { m<-mean(x, na.rm=TRUE); se<-sd(x, na.rm=TRUE); c(m-se, m, m+se) } )
	colnames(arMat)<-rownames(intervals)
	unStandardizedAlroyPres <- cmp_3Trates(alroyCountsFromOccs(fvOccs, occs, restrict.col, restrict.class))[,"pres"]
	arMat[!is.finite(arMat)] <- unStandardizedAlroyPres[!is.finite(unStandardizedAlroyPres)] <- NA

		plot(rowMeans(intervals[-doNotSample,]), unStandardizedAlroyPres[-doNotSample], type="n", xlim=c(max(intervals[-doNotSample,]), min(intervals[-doNotSample,])), ylim=c(0,1), xlab="Time (Ma)", ylab="Preservation Probability", main="Preservation Probability")
		overlayCzTimescale()
		lines(rowMeans(intervals[-doNotSample,]), unStandardizedAlroyPres[-doNotSample], lwd=2, lty=3, col="mediumpurple4")
		text(rowMeans(intervals[-doNotSample,]), arMat[2,-doNotSample,"pres"], paste(rowMeans(intervals[-doNotSample,]), "Ma", sep=""), pos=3, col="gray50", cex=0.5)
		if (!any(doNotSample>0)) nas <- which(apply(arMat[c(1,3),,"pres"], 2, function(x) any(is.na(x)))) else nas <- sort(unique(c(doNotSample,which(apply(arMat[c(1,3),,"pres"], 2, function(x) any(is.na(x)))))))
		if (length(subsampledOccs)>1) polygon(x=c(rowMeans(intervals[-nas,]), rowMeans(intervals[-nas,])[nrow(intervals[-nas,]):1]), y=c(arMat[1,-nas,"pres"], arMat[3,-nas,"pres"][nrow(intervals[-nas,]):1]), col=alphaColor("mediumpurple4", 0.25), border="mediumpurple4", lwd=0.5)
		lines(rowMeans(intervals[-doNotSample,]), arMat[2,-doNotSample,"pres"], lwd=2, col="mediumpurple4")
		# lines(rowMeans(intervals[-doNotSample,]), getPresProbPerIntv(occs, intList)[-doNotSample], lwd=1, col="goldenrod4")
		
		plot(rep(rowMeans(intervals[-doNotSample,]), times=4), arMat[c(1,3),-doNotSample,2:3], type="n", xlim=c(max(intervals[-doNotSample,]), min(intervals[-doNotSample,])), xlab="Time (Ma)", ylab="Origination Rate", main="")
		overlayCzTimescale()
		# text(rowMeans(intervals[-doNotSample,]), arMat[2,-doNotSample,"orig"], paste(rowMeans(intervals[-doNotSample,]), "Ma", sep=""), pos=3, col="gray50", cex=0.5)
		if (!any(doNotSample>0)) nas <- which(apply(arMat[c(1,3),,"orig"], 2, function(x) any(is.na(x)))) else nas <- sort(unique(c(doNotSample,which(apply(arMat[c(1,3),,"orig"], 2, function(x) any(is.na(x)))))))
		if (length(subsampledOccs)>1) polygon(x=c(rowMeans(intervals[-nas,]), rowMeans(intervals[-nas,])[nrow(intervals[-nas,]):1]), y=c(arMat[1,-nas,"orig"], arMat[3,-nas,"orig"][nrow(intervals[-nas,]):1]), col=alphaColor("dodgerblue4", 0.25), border="dodgerblue4", lwd=0.5)
		lines(rowMeans(intervals[-doNotSample,]), arMat[2,-doNotSample,"orig"], lwd=2, col="dodgerblue4")
	
		plot(rep(rowMeans(intervals[-doNotSample,]), times=4), arMat[c(1,3),-doNotSample,2:3], type="n", xlim=c(max(intervals[-doNotSample,]), min(intervals[-doNotSample,])), xlab="Time (Ma)", ylab="Extinction Rate", main="")
		overlayCzTimescale()
		# text(rowMeans(intervals[-doNotSample,]), arMat[2,-doNotSample,"ext"], paste(rowMeans(intervals[-doNotSample,]), "Ma", sep=""), pos=3, col="gray50", cex=0.5)
		if (!any(doNotSample>0)) nas <- which(apply(arMat[c(1,3),,"ext"], 2, function(x) any(is.na(x)))) else nas <- sort(unique(c(doNotSample,which(apply(arMat[c(1,3),,"ext"], 2, function(x) any(is.na(x)))))))
		if (length(subsampledOccs)>1) polygon(x=c(rowMeans(intervals[-nas,]), rowMeans(intervals[-nas,])[nrow(intervals[-nas,]):1]), y=c(arMat[1,-nas,"ext"], arMat[3,-nas,"ext"][nrow(intervals[-nas,]):1]), col=alphaColor("firebrick4", 0.25), border="firebrick4", lwd=0.5)
		lines(rowMeans(intervals[-doNotSample,]), arMat[2,-doNotSample,"ext"], lwd=2, col="firebrick4")
}

