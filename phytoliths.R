source("~/Dropbox/code/R/common_src/strat.R")
source("~/Dropbox/code/R/common_src/utils_marcot.R")
phytoliths <- read.csv ("~/Dropbox/code/common_dat/Strömberg2006.csv")
ps <- c("wetland", "palm", "otherFI", "closed", "pooid", "pooid_non", "PACCAD", "otherGSSC", "nondiag")
# ps <- c("palm", "otherFI", "closed", "pooid", "pooid_non", "PACCAD", "otherGSSC")


doPhytolithRateAnalysis <- function(intervals, sig=0.01, reps=1) {
	source("~/Dropbox/code/R/amandaTeeth/src/amandaSrc.R")
	intervals <- intervals[intervals[,"ageTop"] <= max(phytoliths$ageMax) & intervals[,"ageBase"] > min(phytoliths$ageMin),]

	repIntPhy <- list()
	for (rep in seq_len(reps)) {
		phytoliths$ma_rand <- apply(phytoliths[,c("ageMax", "ageMin")], 1, function(x) runif(1, max=x[1], min=x[2]))
		intPhy <- list()
		for (int in seq_len(nrow(intervals))) {
			intPhy[[int]] <- which(is.finite(phytoliths$ma_rand) & (phytoliths$ma_rand >= intervals$ageTop[int]) & (phytoliths$ma_rand < intervals$ageBase[int]))
		}
		repIntPhy[[rep]] <- intPhy
	}
		
	thisCounts <- apply(sapply(repIntPhy, function(y) sapply(y, function(x) colSums(phytoliths[x,ps], na.rm=TRUE)), simplify="array"), c(1,2), median, na.rm=TRUE)
	# thisCounts <- rbind(closed=colSums(sapply(intPhy, function(x) colSums(phytoliths[x,ps], na.rm=TRUE))[1:4,]), open=colSums(sapply(intPhy, function(x) colSums(phytoliths[x,ps], na.rm=TRUE))[5:8,]))
	colnames(thisCounts) <- rownames(intervals)

	# plot(rowMeans(intervals), thisCounts[2,]/thisCounts[1,], xlim=c(40,0))

	zeros <- which(colSums(thisCounts)==0)
	if (length(zeros) > 0) {
		thisCounts <- thisCounts[,-zeros]
		intervals <- intervals[-zeros,]
	}

	par(mfrow=c(ncol(thisCounts), 1), mar=c(0.5,1,0.5,0))
	apply(thisCounts, 2, function(x) barplot(x, ylim=c(0,max(thisCounts))))
	
	# apply(thisCounts, 2, function(x) x/sum(x))
	# doHandleyTest(thisCounts, sig=sig, do.heuristic=TRUE, do.parallel=FALSE)
	optList <- doHandleyTest(thisCounts, sig=sig, do.heuristic=TRUE, do.parallel=FALSE)
	# intervals[optList[[length(optList)-1]]$optBreaks,]
	# intervals[sort(optList[[length(optList)-1]]$optBreaks),]
}

thisClust <- hclust(vegdist(phytoliths[apply(phytoliths[,ps], 1, function(x) all(is.finite(x))),ps], method="bray", na.rm=TRUE), method="ward")

plotTopesRateAnalysis <- function (optList_topes, intervals) {
	# plot(phytoliths$ma_mid, zachos2001$d18Oadj, xlim=c(max(intervals), min(intervals)), ylim=c(max(zachos2001$d18Oadj, na.rm=TRUE), min(zachos2001$d18Oadj, na.rm=TRUE)), type="n", xaxt="n", cex=0.5, ylab=expression(paste(delta^18,plain(O))))
	# overlayCzTimescale(do.subepochs=TRUE)
	
	# intColors <- rainbow(length(optList_topes[[length(optList_topes)-1]]$optBreaks))
	# for (i in length(optList_topes[[length(optList_topes)-1]]$optBreaks):2 ) { 
		# # rect(xleft=intervals[sort(optList_topes[[length(optList_topes)-1]]$optBreaks)[i],2], xright=intervals[sort(optList_topes[[length(optList_topes)-1]]$optBreaks)[i-1],2], ytop=10, ybottom=-10, col=alphaColor(intColors[i],0.25), border=intColors[i])
		# rect(xleft=intervals[sort(optList_topes[[length(optList_topes)-1]]$optBreaks)[i],2], xright=intervals[sort(optList_topes[[length(optList_topes)-1]]$optBreaks)[i-1],2], ytop=10, ybottom=-10, col=NA, border=alphaColor("dodgerblue1", 0.5), lwd=1)
	# }
	# rect(xleft=intervals[sort(optList_topes[[length(optList_topes)-1]]$optBreaks)[1],2], xright=min(intervals), ytop=10, ybottom=-10, col=NA, border=alphaColor("dodgerblue1", 0.5), lwd=1)
	# rect(xleft=max(intervals), xright=intervals[sort(optList_topes[[length(optList_topes)-1]]$optBreaks)[length(optList_topes[[length(optList_topes)-1]]$optBreaks)],2], ytop=10, ybottom=-10, col=NA, border=alphaColor("dodgerblue1", 0.5), lwd=1)
	
	# intTopes <- list()
	# for (int in seq_len(nrow(intervals))) {
		# intTopes[[int]] <- which(is.finite(phytoliths$ma_mid) & phytoliths$ma_mid >= intervals$ageTop[int] & phytoliths$ma_mid < intervals$ageBase[int] & is.finite(zachos2001$d18Oadj))
	# }

	# quants <- sapply(intTopes, function(x) quantile(zachos2001$d18Oadj[x], probs=c(0,0.25, 0.5, 0.75,1), na.rm=TRUE))
	# polygon(c(rowMeans(intervals), rev(rowMeans(intervals))), c(quants[1,], rev(quants[5,])), col=alphaColor("dodgerblue1", 0.25), border="dodgerblue1")
	# polygon(c(rowMeans(intervals), rev(rowMeans(intervals))), c(quants[2,], rev(quants[4,])), col=alphaColor("dodgerblue4", 0.25), border="dodgerblue4")
	
	# points(phytoliths$ma_mid, zachos2001$d18Oadj, col=alphaColor("gray0", 0.25), cex=0.35)
	
	# lines(rowMeans(intervals), quants[3,], col=alphaColor("white", 0.5), lwd=5,)
	# lines(rowMeans(intervals), quants[3,], col=alphaColor("dodgerblue4", 1.0), lwd=3)
	# points(rowMeans(intervals), quants[3,], col=alphaColor("dodgerblue1", 0.5), cex=0.5)
}

