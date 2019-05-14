source("~/Dropbox/code/R/common_src/strat.R")
source("~/Dropbox/code/R/common_src/utils_marcot.R")
zachos2001 <- read.csv ("~/Dropbox/code/common_dat/zachos2001.csv")

doTopesRateAnalysis <- function(intervals, thisBreaks=NULL, sig=0.01) {
	source("~/Dropbox/code/R/amandaTeeth/src/amandaSrc.R")

	# intervals_topes <- makeIntervals(max(zachos2001$Age), 0, 0.5)

	intTopes <- list()
	for (int in seq_len(nrow(intervals))) {
		intTopes[[int]] <- which(is.finite(zachos2001$Age) & zachos2001$Age >= intervals$ageTop[int] & zachos2001$Age < intervals$ageBase[int] & is.finite(zachos2001$d18Oadj))
	}

	if (is.null(thisBreaks)) thisBreaks <- hist(zachos2001$d18Oadj[unique(unlist(intTopes))], plot=FALSE)$breaks

	thisCounts <- sapply(intTopes, function(x) hist(zachos2001$d18Oadj[x], breaks=thisBreaks, plot=FALSE)$counts)
	doHandleyTest(thisCounts, n=nrow(zachos2001), sig=sig, do.heuristic=TRUE, do.parallel=FALSE)
	# intervals[sort(optList_topes[[length(optList_topes)-1]]$optBreaks),]
}

plotTopesRateAnalysis <- function (optList_topes, intervals, x.axis=TRUE) {
	intTopes <- which(is.finite(zachos2001$Age) & zachos2001$Age >= min(intervals) & zachos2001$Age < max(intervals) & is.finite(zachos2001$d18Oadj))
	# plot(zachos2001$Age, zachos2001$d18Oadj, xlim=c(max(intervals), min(intervals)), ylim=c(max(zachos2001$d18Oadj[intTopes], na.rm=TRUE), min(zachos2001$d18Oadj[intTopes], na.rm=TRUE)), type="n", xaxt="n", cex=0.5, xlab="Time (Ma)", ylab=expression(paste(delta^18,plain(O))))
	# if(x.axis) axis(side=1, at=seq(min(intervals), c(max(intervals)), by=5))
	plot(zachos2001$Age, zachos2001$d18Oadj, xlim=c(max(intervals), min(intervals)), ylim=c(max(zachos2001$d18Oadj[intTopes], na.rm=TRUE), min(zachos2001$d18Oadj[intTopes], na.rm=TRUE)), type="n", xaxt="n", cex=0.5, xlab="Time (Ma)", ylab=expression(paste(delta^18,plain(O))), col="gray75", fg="gray75", col.lab="gray75", col.axis="gray75", col.main="gray75", family="sans")
	plot(zachos2001$Age, zachos2001$d18Oadj, xlim=c(54,4), ylim=c(max(zachos2001$d18Oadj[intTopes], na.rm=TRUE), min(zachos2001$d18Oadj[intTopes], na.rm=TRUE)), type="n", xaxt="n", cex=0.5, xlab="Time (Ma)", ylab=expression(paste(delta^18,plain(O))), col="gray75", fg="gray75", col.lab="gray75", col.axis="gray75", col.main="gray75", family="sans")
	if(x.axis) axis(side=1, at=seq(min(intervals), c(max(intervals)), by=5), col="gray75", col.ticks="gray75")
	overlayCzTimescale(do.subepochs=TRUE)

	intColors <- rainbow(length(optList_topes[[length(optList_topes)-1]]$optBreaks))
	for (i in length(optList_topes[[length(optList_topes)-1]]$optBreaks):2 ) {
		# rect(xleft=intervals[sort(optList_topes[[length(optList_topes)-1]]$optBreaks)[i],2], xright=intervals[sort(optList_topes[[length(optList_topes)-1]]$optBreaks)[i-1],2], ytop=10, ybottom=-10, col=alphaColor(intColors[i],0.25), border=intColors[i])
		rect(xleft=intervals[sort(optList_topes[[length(optList_topes)-1]]$optBreaks)[i],2], xright=intervals[sort(optList_topes[[length(optList_topes)-1]]$optBreaks)[i-1],2], ytop=10, ybottom=-10, col=NA, border=alphaColor("dodgerblue1", 0.5), lwd=1)
	}
	rect(xleft=intervals[sort(optList_topes[[length(optList_topes)-1]]$optBreaks)[1],2], xright=min(intervals), ytop=10, ybottom=-10, col=NA, border=alphaColor("dodgerblue1", 0.5), lwd=1)
	rect(xleft=max(intervals), xright=intervals[sort(optList_topes[[length(optList_topes)-1]]$optBreaks)[length(optList_topes[[length(optList_topes)-1]]$optBreaks)],2], ytop=10, ybottom=-10, col=NA, border=alphaColor("dodgerblue1", 0.5), lwd=1)

	intTopes <- list()
	for (int in seq_len(nrow(intervals))) {
		intTopes[[int]] <- which(is.finite(zachos2001$Age) & zachos2001$Age >= intervals$ageTop[int] & zachos2001$Age < intervals$ageBase[int] & is.finite(zachos2001$d18Oadj))
	}

	quants <- sapply(intTopes, function(x) quantile(zachos2001$d18Oadj[x], probs=c(0,0.25, 0.5, 0.75,1), na.rm=TRUE))
	polygon(c(rowMeans(intervals), rev(rowMeans(intervals))), c(quants[1,], rev(quants[5,])), col=alphaColor("dodgerblue1", 0.25), border="dodgerblue1")
	# polygon(c(rowMeans(intervals), rev(rowMeans(intervals))), c(quants[2,], rev(quants[4,])), col=alphaColor("dodgerblue4", 0.25), border="dodgerblue4")

	points(zachos2001$Age, zachos2001$d18Oadj, col=alphaColor("dodgerblue4", 0.5), cex=1)

	lines(rowMeans(intervals), quants[3,], col=alphaColor("white", 0.5), lwd=5)
	lines(rowMeans(intervals), quants[3,], col=alphaColor("dodgerblue4", 1.0), lwd=3)
	points(rowMeans(intervals), quants[3,], col=alphaColor("dodgerblue1", 0.5), cex=0.5)
}

getAlroyStatisticsOneInterval <- function(thisInt, zachos2001) {
	intTopes <- which(is.finite(zachos2001$Age) & is.finite(zachos2001$d18Oadj) & zachos2001$Age >= thisInt["ageTop"] & zachos2001$Age < thisInt["ageBase"])
	thisLM <- lm(formula=zachos2001$d18Oadj[intTopes] ~ zachos2001$Age[intTopes])
	zachos2001[intTopes,]
	# plot(x=zachos2001$Age[intTopes], y=zachos2001$d18Oadj[intTopes])
	# abline(thisLM)

	# mean(zachos2001$d18Oadj[intTopes])
	c(n=length(intTopes), midpoint=(thisLM$coefficients[2] * mean(thisInt)) + thisLM$coefficients[1], stdev=sd(thisLM$residuals))
}

getAlroyStatistics <- function(intervals) {
	dat <- t(apply(intervals, 1, getAlroyStatisticsOneInterval,zachos2001))
  data.frame(n=dat[,1], midpoint=dat[,2], change=c(diff(dat[,2]), NA), volatility=c(abs(diff(dat[,2])), NA), stdev=dat[,3]) #det.stdev=c(diff(dat[,3]), NA)
}
