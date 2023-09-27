#richness/rates
require(compiler)

getFourClassesForOneTaxon<-function(x) {
	y<-array(FALSE, dim=c(base::length(x), 4), dimnames=list(NULL, c("Ft", "bt", "bL", "FL")))
	if (sum(x)==0) { return(y)
	} else if (sum(x)==1) { y[which(x),"FL"]<-TRUE
	} else {
		intv <- which(x)
		y[intv[1],"bL"] <- TRUE
		y[intv[-c(1, base::length(intv))], "bt"] <- TRUE
		y[intv[base::length(intv)],"Ft"] <- TRUE
	}
	y
}
cmp_GetFourClassesForOneTaxon <- cmpfun(getFourClassesForOneTaxon)

fourClassesFromOccs <- function(occList, this.occs, restrict.col=NULL, restrict.class=NULL) {
	# for (i in 1:length(occList)) if (!is.null(occList[[i]])) as.character(occs[occs$occurrence_no%in%occList[[i]],]$taxon) else NA 
	taxList <- lapply(occList, function(x) { if (!is.null(x)) as.character(occs$taxon[occs$occurrence_no %in% x]) else NA })	#gets the taxa in each (occurrence in each) interval
	taxVec <- unique(unlist(taxList))
	taxVec <- taxVec[taxVec !=""]
	if (!is.null(restrict.col) & !is.null(restrict.class)) taxVec <- taxVec[taxVec %in% occs$taxon[occs[,restrict.col] %in% restrict.class]]
	
	tax_out <- taxVec[taxVec %in% occs$taxon[!(occs[,restrict.col] %in% restrict.class)]]
	tax_in <- taxVec[taxVec %in% occs$taxon[occs[,restrict.col] %in% restrict.class]]	
	tax_in[tax_in %in% tax_out]
	# occs[occs$taxon=="Aderidae",]	
	
	oList <- sapply(taxVec, function(x) { sapply(taxList, function(y) { x%in%y }) })										#converts taxa in each interval into a boolean
	oList <- apply(oList, 2, function(x) { x[which(x)[1]:which(x)[length(which(x))]]<-TRUE;x }) 											# makes range throughs complete
	classMat <- array(dim=c(nrow(oList), 4, ncol(oList)), dimnames=list(NULL, c("Ft", "bt", "bL", "FL")))
	for (i in seq_len(ncol(oList))) classMat[,,i] <- cmp_GetFourClassesForOneTaxon(oList[,i])
	apply(classMat, c(1,2), sum)
}

getRichnessBoxFromOccBox <- function (occBox, occs, method=c("RT", "SIB", "BC", "TT"), restrict.col=NULL, restrict.class=NULL, doNotSample=NULL, do.parallel=FALSE) {
	if (do.parallel) require(parallel)
	method <- match.arg(method)
	if (method=="SIB") {
		if (do.parallel) { mat <- simplify2array(mclapply(occBox, getSIBRichnessFromOneOccList, occs, restrict.col, restrict.class, mc.cores=detectCores()-2))
		} else mat <- sapply(occBox, getSIBRichnessFromOneOccList, occs, restrict.col, restrict.class) 
	} else if (method=="RT") {
		if (do.parallel) { mat <- simplify2array(mclapply(occBox, getRangeThroughRichnessFromOneOccList, occs, restrict.col, restrict.class, mc.cores=detectCores()-2))
		} else mat <- sapply(occBox, getRangeThroughRichnessFromOneOccList, occs, restrict.col, restrict.class)
	} else if (method=="BC") {
		if (do.parallel) { mat <- simplify2array(mclapply(occBox, getRangeBoundaryCrosserRichnessFromOneOccList, occs, restrict.col, restrict.class, mc.cores=detectCores()-2))
		} else mat <- sapply(occBox, getRangeBoundaryCrosserRichnessFromOneOccList, occs, restrict.col, restrict.class) 
	} else if (method=="TT") {
		if (is.null(doNotSample)) doNotSample=-seq_along(occBox[[1]])
		if (do.parallel) { mat <- simplify2array(mclapply(occBox, getThreeTimerRichnessFromOneOccList, occs, restrict.col, restrict.class, doNotSample, mc.cores=detectCores()-2))
		} else mat <- sapply(occBox, getThreeTimerRichnessFromOneOccList, occs, restrict.col, restrict.class, doNotSample) 
	}
	mat
}

getSIBRichnessFromOneOccList <- function(occList, occs, restrict.col=NULL, restrict.class=NULL, do.parallel=FALSE) {
	if (do.parallel) require(parallel)
	if (is.null(restrict.col) | is.null (restrict.class)) { 
		if (do.parallel) { simplify2array(mclapply(occList, function (x) { length(unique(occs[occs$occurrence_no%in%x,]$taxon)) }, mc.cores=detectCores()-2))
		} else sapply(occList, function (x) { length(unique(occs[occs$occurrence_no%in%x,]$taxon)) })
	} else {
		if (do.parallel) { simplify2array(mclapply(occList, function (x) { length(unique(occs[occs$occurrence_no%in%x & occs[,restrict.col] %in% restrict.class,]$taxon)) }, mc.cores-detectCores()-2))
		} else sapply(occList, function (x) { length(unique(occs[occs$occurrence_no%in%x & occs[,restrict.col] %in% restrict.class,]$taxon)) })
	}
}

getRangeThroughRichnessFromOneOccList <- function(occList, occs, restrict.col=NULL, restrict.class=NULL) {
	rowSums(fourClassesFromOccs(occList, occs, restrict.col, restrict.class))
}

getRangeBoundaryCrosserRichnessFromOneOccList <- function(occList, occs, restrict.col=NULL, restrict.class=NULL) {
	apply(fourClassesFromOccs(occList, occs, restrict.col, restrict.class), 1, function (y) sum(y[c("bt", "bL")]))
}

alroyCounter<-function(intv, x) {
	y<-array(FALSE, dim=c(1,5))
	if (!x[intv[1]] & x[intv[2]] & !x[intv[3]]) { y[1]<-TRUE
	} else if (x[intv[1]] & x[intv[2]] & !x[intv[3]]) { y[2]<-TRUE
	} else if (!x[intv[1]] & x[intv[2]] & x[intv[3]]) { y[3]<-TRUE
	} else if (x[intv[1]] & !x[intv[2]] & x[intv[3]]) { y[4]<-TRUE
	} else if (all(x[intv])) y[5]<-TRUE
	y
}

alroyCountsForOneTaxon <- function(o) {
	counts <- array(FALSE, dim=c(length(o),5), dimnames=list(NULL, c("1T", "2Ti", "2Ti+1", "PT", "3T")))
	y <- suppressWarnings(cbind(3:length(o), 2:length(o), 1:length(o)))
	y <- y[-which(y[,1]<row(y)[,1]),] # cuts off last rows from the ordering thing (previous line)
	counts[2:(length(o)-1),] <- t(apply(y, 1, alroyCounter, o))
	counts
}

alroyCountsFromOccs <- function(occList, occs, restrict.col=NULL, restrict.class=NULL) {
	# for (i in 1:length(occList)) if (!is.null(occList[[i]])) as.character(occs[occs$occurrence_no%in%occList[[i]],]$taxon) else NA 
	taxList <- lapply(occList, function(x) { if (!is.null(x)) unique(as.character(occs[occs$occurrence_no%in%x,]$taxon)) else NA })	#gets the taxa in each (occurrence in each) interval
	taxVec <- unique(unlist(taxList))
	if (!is.null(restrict.col) & !is.null(restrict.class)) taxVec <- taxVec[taxVec %in% occs$taxon[occs[,restrict.col] %in% restrict.class]]
	oList<-sapply(taxVec, function(x) { sapply(taxList, function(y) { x%in%y }) })										#converts taxa in each interval into a boolean
	# oList<-apply(oList, 2, function(x) { x[which(x)[1]:which(x)[length(which(x))]]<-TRUE;x }) 											# makes range throughs complete
	classMat<-array(dim=c(nrow(oList), 5, ncol(oList)), dimnames=list(NULL, c("1T", "2Ti", "2Ti+1", "PT", "3T")))
	for (i in seq_along(taxVec)) classMat[,,i]<-alroyCountsForOneTaxon(oList[,i])
	apply(classMat, c(1,2), sum)
}

getThreeTimerRichnessFromOneOccList <- function(occList, occs, restrict.col=NULL, restrict.class=NULL, doNotSample=NULL) {
	# alroy's three-timer correction
	if (is.null(doNotSample)) doNotSample <- -seq_along(occList)
	arCounts <- alroyCountsFromOccs(occList, occs, restrict.col, restrict.class)
	pres <- arCounts[,"3T"]/(arCounts[,"3T"]+arCounts[,"PT"]) # sampling completeness statistic of Alroy et al. 2008
	grandEstimate <- sum(arCounts[-doNotSample, "3T"]) / (sum(arCounts[-doNotSample,"3T"]) + sum(arCounts[-doNotSample, "PT"]))
	pres[!(is.finite(pres))] <- 1.0
	rescaled <- getSIBRichnessFromOneOccList(occList, occs, restrict.col, restrict.class, do.parallel=FALSE)
	rescaled[-doNotSample] <- (rescaled[-doNotSample] / pres[-doNotSample]) * grandEstimate
	rescaled[!(is.finite(rescaled))] <- NA
	rescaled
}

perCapitaRatesFromFourClasses <- function(thisMat, div.rates=FALSE) {
	orig<-apply(thisMat, 1, function(x) { -log(x["bt"]/{x["bt"]+x["Ft"]}) })
	ext<-apply(thisMat, 1, function(x) { -log(x["bt"]/{x["bt"]+x["bL"]}) })
	if (div.rates) { cbind(orig, ext, div=orig-ext) 
	} else cbind(orig, ext)
}
cmp_PCrates <- cmpfun(perCapitaRatesFromFourClasses)

perTaxonRatesFromFourClasses<-function(thisMat) {
	orig<-apply(thisMat, 1, function(x) { {x["Ft"]+x["FL"]}/sum(x) })
	ext<-apply(thisMat, 1, function(x) { {x["bL"]+x["FL"]}/sum(x) })
	cbind(orig, ext)
}
cmp_PTrates <- cmpfun(perTaxonRatesFromFourClasses)

getThreeTimerRatesFromCounts <- function(arCounts) {
	#three timer rates Alroy (2008) PNAS
	rates <- array(NA, dim=c(nrow(arCounts), 3), dimnames=list(NULL, c("pres", "orig", "ext")))
	rates[c(1,nrow(rates)),"pres"] <- 1.0
	for (i in 2:(nrow(arCounts)-1)) rates[i,"pres"] <- arCounts[i,"3T"]/(arCounts[i,"3T"]+ arCounts[i,"PT"])
	for (i in 2:(nrow(arCounts)-1)) {
		rates[i,"ext"]  <- log(arCounts[i,"2Ti"]  /arCounts[i,"3T"]) + log(rates[i-1,"pres"]) 	# intervals are numbered from the youngest, so a smaller number (i-1) means the next interval
		rates[i,"orig"] <- log(arCounts[i,"2Ti+1"]/arCounts[i,"3T"]) + log(rates[i+1,"pres"])	# intervals are numbered from the youngest, so a larger number (i+1) means the previous interval
	}
	rates[c(1,nrow(rates)),"pres"] <- NA
	rates
}
cmp_3Trates <- cmpfun(getThreeTimerRatesFromCounts)

wobbleMachine<-function(index, richness) {
	abs(log((richness[index[2]]^2)/(richness[index[1]]*richness[index[3]])))
}

wobbleIndex<-function(richness) {
	y<-cbind(3:length(richness), 2:length(richness), 1:length(richness))
	y<-y[-which(y[,1]<row(y)[,1]),]
	median(apply(y, 1, wobbleMachine, richness))
}

proportionalRichnessOneInterval <- function(thisO, do.prop=TRUE) {
	# counts <- sapply(higherTax, function(x) sum(x %in% names(thisO)[thisO]) )
	# sapply(names(thisO),
	counts <- sapply(sort(unique(names(thisO))), function(x) sum(names(thisO[thisO]) %in% x))
	if (do.prop) counts/sum(counts) else counts
}

getRichnessMatrixFromOccs <- function(this.occs, intervals, occs.column.higher="family", occs.column.subtaxa="accepted_name", range.through=FALSE,  do.prop=FALSE, do.log=FALSE, random=FALSE, do.parallel=FALSE) {
		if (do.parallel) require(parallel)
		this.occs[,"ma_det"] <- getOccDatesFromColletions(this.occs, random=random)
		
		intList <- listifyMatrixByRow(intervals)
	
		if (do.parallel) { occList <- mclapply(intList, function(x) this.occs$occurrence_no[this.occs$ma_det <= x[2] & this.occs$ma_det > x[1]], mc.cores=detectCores()-2) 
		} else occList <- lapply(intList, function(x) this.occs$occurrence_no[this.occs$ma_det <= x[2] & this.occs$ma_det > x[1]])
		
		if (occs.column.subtaxa != "accepted_name") { subtaxa.list <- lapply(occList, function(x) { if (!is.null(x)) unique(as.character(this.occs[this.occs$occurrence_no %in% x, occs.column.subtaxa])) else NA })	# gets the taxa in each (occurrence in each) interval
		} else subtaxa.list <- lapply(occList, function(x) { if (!is.null(x)) unique(as.character(this.occs[(this.occs$occurrence_no %in% x) & grepl(pattern="species", x=this.occs$accepted_rank), occs.column.subtaxa])) else NA })	# gets the taxa in each (occurrence in each) interval
		
		oList <- sapply(unique(unlist(subtaxa.list)), function(x) { sapply(subtaxa.list, function(y) { x %in% y }) })								# converts taxa in each interval into a boolean
		if (range.through) oList <- apply(oList, 2, function(x) { x[which(x)[1]:which(x)[length(which(x))]] <- TRUE; x }) 					# makes range throughs complete
		
		tax.table <- unique(cbind(as.character(this.occs[,occs.column.higher]), as.character(this.occs[,occs.column.subtaxa])))
		tax.table[tax.table[,1] == "", 1] <- "indeterminate"
		prop <- t(apply(oList, 1, function(x) sapply(sort(unique(tax.table[,1])), function(this.taxon, this.names) sum(this.names==this.taxon), this.names=tax.table[tax.table[,2] %in% names(x)[x], 1])))
		rownames(prop) <- rownames(intervals)
	
		if (do.log) prop <- log(prop)
		if (do.prop) if (do.prop) prop <- prop/rowSums(prop) else counts
		
		prop
}

plotStackedRichness <- function(this.box, intervals, reorder.taxa = TRUE, do.log = FALSE, xlim = NULL, ylim = NULL, 
                                xaxp = NULL, yaxp = NULL, cex.axis = 1, cex.lab = 1, las = 0, plot.adj = NULL, 
                                col.axis = "black", col.lab = "black", poly.col = NULL,
                                xaxt = NULL, yaxt = NULL, 
                                ylab = "Richness (Number of Subtaxa)", xlab = "Time (Ma)",
                                add.legend = TRUE, prop.ylab = FALSE,
                                numbers.only = FALSE, overlay.labels = FALSE, overlay.color = TRUE, do.subepochs = FALSE, thisAlpha.text = 0.33, thisAlpha.intervals = 0.33, borderCol = "white", invertTime = FALSE, scale.cex = 0.75, scale.headers = 0.95, text.offset = 0.025) 
{
  # this.box <-this.box[,order(this.box[nrow(this.box)-1,], decreasing=TRUE)]
  # this.box <-this.box[,order(colMeans(this.box, na.rm=TRUE))]
  thisOrder <- apply(this.box, 2, function(x) max(which(x==max(x, na.rm=TRUE)))) #what does this line do? why?
  if (reorder.taxa) thisOrder <- sort(thisOrder, decreasing=TRUE)
  # if (do. this.box) { thisOrder <- sort(apply(this.box, 2, function(x) max(which(x==max(x, na.rm=TRUE)))), decreasing=TRUE) 
  # } else thisOrder <- sort(apply(this.box, 2, function(x) max(which(x==max(x, na.rm=TRUE)))), decreasing=TRUE)
  this.box <- this.box[,match(names(thisOrder), colnames(this.box))]
  this.box[!is.finite(this.box)] <- 0
  for (i in 2:ncol(this.box)) this.box[,i] <- this.box[,i] + this.box[,i-1]
  
  if(!is.null(poly.col) & is.character(poly.col))
  {
    this.colors <- poly.col
    this.fills <-alphaColor(poly.col, alpha = 0.5)
  } else {
    this.colors <-rainbow(ncol(this.box))
    this.fills <-rainbow(ncol(this.box), alpha=0.50)
  }
#   if(ylab)
#  {
    if (all(rowSums(this.box) == 1, na.rm=TRUE) | prop.ylab) 
    {  #changed prop to this.box in this line
      if(!is.null(ylim)) {y_lim <- ylim
      } else y_lim <- c(0,1)
#      y_lab<-"Proportion of Subtaxa"
    } else {
      if (do.log) y_lim <- c(min(this.box[this.box[] > 0], na.rm=TRUE), 1.05*max(this.box, na.rm=TRUE)) else y_lim<-c(0, max(this.box))
      if(!is.null(ylim)) y_lim <- ylim
#      y_lab<-"Richness (Number of Subtaxa)"
    }
#  } 
  
  if (is.null(xlim)) { plot(rowMeans(intervals), this.box[,1], xlim=c(max(rowMeans(intervals)), min(rowMeans(intervals))), 
                            ylim=y_lim, xaxp= xaxp, yaxp= yaxp, type="n", xaxt = xaxt, yaxt = yaxt, las = las,
                            usr=c(max(rowMeans(intervals), min(rowMeans(intervals)), 0, 1)), 
                            xaxs="i", yaxs="i", xlab = x_lab, ylab = ylab, 
                            col.axis = col.axis, col.lab = col.lab, cex.axis = cex.axis, cex.lab = cex.lab, adj = plot.adj)
  } else plot(rowMeans(intervals), this.box[,1], xlim= xlim, ylim= y_lim, xaxp= xaxp, yaxp= yaxp, las = las, type="n", 
              xaxt = xaxt, yaxt = yaxt, xlab = xlab, ylab = ylab, col.axis = col.axis, col.lab = col.lab, cex.axis = cex.axis, cex.lab = cex.lab, adj = plot.adj)
  
  overlayCzTimescale(do.subepochs=do.subepochs, color = overlay.color, thisAlpha.text = thisAlpha.text, thisAlpha.intervals = thisAlpha.intervals, borderCol = borderCol, invertTime = invertTime, scale.cex = scale.cex, scale.headers = scale.headers, text.offset = text.offset)
  
  polygon(c(rowMeans(intervals), rev(rowMeans(intervals))), c(this.box[,1], rep(0, nrow(intervals))), col=this.fills[1], border="gray33", xlab="Time (Ma)", ylab="Proportion of Taxa", lwd=0.5)
  for (i in 2:ncol(this.box)) polygon(c(rowMeans(intervals), rev(rowMeans(intervals))), c(this.box[,i], this.box[nrow(intervals):1,i-1]), col=this.fills[i], border="gray33", lwd=0.5)
  
  if (overlay.labels) 
  {
    if (numbers.only) 
    {
      translation.table <- cbind(names(thisOrder), seq_along(colnames(prop)))
      this.labels <- translation.table[,2]
    } else this.labels <- colnames(this.box)
    text(rowMeans(intervals)[thisOrder[1]], this.box[thisOrder[1],1]/2, labels=this.labels[1], adj=(0-(rowMeans(intervals)[thisOrder[1]]-par()$usr[1])/(par()$usr[1]-par()$usr[2])), col="white", cex=0.5, font=2)
    for (i in 2:ncol(this.box)) 
    {
      text(rowMeans(intervals)[thisOrder[i]], mean(c(this.box[thisOrder[i],i], this.box[thisOrder[i],i-1])), adj=(0-(rowMeans(intervals)[thisOrder[i]]-par()$usr[1])/(par()$usr[1]-par()$usr[2])), labels=this.labels[i], col="white", cex=0.5, font=2)
    }
  }
  
  if (add.legend) 
  {
    legend("topleft", legend=rev(colnames(this.box)), inset=0.055, fill=NA, border=NA, cex=0.5, box.lty=0, bg=rgb(0,0,0,0.5)) #this is the dropshadow on the legend
    legend("topleft", legend=rev(colnames(this.box)), inset=0.05, fill=this.fills[ncol(this.box):1], border=this.colors[ncol(this.box):1], cex=0.5, bg="white")
  }
  if (numbers.only) {translation.table}
}
