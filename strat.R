#strat.R

makeIntervals <- function(startDate, endDate, intervalLength=1) {
	# function begins intervals at start date, and includes an interval >= end date
	#this function works back in time, so if the start date is greater than the end date, it swaps them
	if (startDate>endDate) {
		holder<-startDate
		startDate<-endDate
		endDate<-holder
	}
	# (endDate+(intervalLength-(endDate%%intervalLength)))-intervalLength
	tops<-seq(from=startDate, to=endDate, by=intervalLength)
	mat<-data.frame(tops, tops+intervalLength)
	colnames(mat)<-c("ageTop", "ageBase")
	rownames(mat)<-paste(rowMeans(mat), "Ma")
	mat
}

getIntervals_Alroy2000<-function() {
	intv<-cbind(ageTop=c(0.0, 0.0117, 1.8, 4.9, 10.3, 13.6, 16.3, 20.6, 24.8, 26.3, 30.8, 33.3, 33.9, 38.0, 42.0, 46.2, 50.3, 55.4, 56.8, 60.2, 63.3), 
		ageBase=c(0.0117, 1.8, 4.9, 10.3, 13.6, 16.3, 20.6, 24.8, 26.3, 30.8, 33.3, 33.9, 38.0, 42.0, 46.2, 50.3, 55.4, 56.8, 60.2, 63.3, 65.5))
	rownames(intv)<-c("Holocene", "Irvingtonian", "Blancan", "Hemphillian", "Clarendonian", "Barstovian", "Hemingfordian", "Harrisonian", "Monroecreekian", "Geringian", "Whitneyan", "Orellan", "Chadrodnian", "Duchesnian", "Uintan", "Bridgerian", "Wasatchian", "Clarkforkian", "Tiffanian", "Torrejonian", "Puercan")
	intv
}

getIntervals_StdStages<-function() {
	intv<-cbind(ageTop <- c(0.0, 0.0117, 0.126, 0.781, 1.806, 2.588, 3.6, 5.332, 7.246, 11.608, 13.82, 15.97, 20.43, 23.03, 28.4, 33.9, 37.2, 40.4, 48.6, 55.8, 58.7, 61.1, 65.5, 70.6, 83.5, 85.8, 88.6, 93.6, 99.6, 112.0, 125.0, 130.0, 133.9, 140.2), 
		ageBase <- c(0.0117, 0.126, 0.781, 1.806, 2.588, 3.6, 5.332, 7.246, 11.608, 13.82, 15.97, 20.43, 23.03, 28.4, 33.9, 37.2, 40.4, 48.6, 55.8, 58.7, 61.1, 65.5, 70.6, 83.5, 85.8, 88.6, 93.6, 99.6, 112.0, 125.0, 130.0, 133.9, 140.2, 145.5))
	rownames(intv)<-c("Upper", "Ionian", "Calabrian", "Gelasian", "Piacenzian", "Zanclean", "Messinian", "Tortonian", "Serravallian", "Langhian", "Burdigalian", "Aquitanian", "Chattian", "Rupelian", "Priabonian", "Bartonian", "Lutetian", "Ypresian", "Thanetian", "Selandian", "Danian", "Maastrichtian", "Campanian", "Santonian", "Coniacian", "Turonian", "Cenomanian", "Albian", "Aptian", "Barremian", "Hauterivian", "Valanginian", "Berriasian")
	intv
}

############################################################################################################################################

getCollectionAgesFromOccs <- function(occs, random=FALSE, max.age.name="max_ma", min.age.name="min_ma") {
	cols <- unique(occs[ ,c("collection_no", max.age.name, min.age.name)])
	if (random) { cols <- data.frame(collection_no=cols$collection_no, collection_age=runif(n=nrow(cols), min=cols[ ,min.age.name], max=cols[ ,max.age.name]))
	} else cols <- data.frame(collection_no=cols$collection_no, collection_age=rowMeans(cols[ ,c(max.age.name, min.age.name)], na.rm=TRUE))
	cols
}

getOccDatesFromColletions <- function(occs, random=FALSE, max.age.name="max_ma", min.age.name="min_ma") {
	coll.dates <- getCollectionAgesFromOccs(occs, random=random)
	occ.dates <- coll.dates$collection_age[match(occs$collection_no,coll.dates$collection_no)]
	occ.dates
}

####
#### this is the function to use - dates collections, then assigns occurrences those dates, then dates taxa
####
getTaxonRangesFromOccs <- function(occs, rank="species", random=FALSE) {
	if (all(c("collection_no", "ma_max", "ma_min") %in% names(occs))) { 
		max.age.name <- "ma_max"
		min.age.name <- "ma_min"
	} else if (all(c("collection_no", "max_ma", "min_ma") %in% names(occs))) { 
		max.age.name <- "max_ma"
		min.age.name <- "min_ma"
	} else error()

	cols <- getCollectionAgesFromOccs(occs=occs, random=random, max.age.name=max.age.name, min.age.name= min.age.name)
	
	if(rank == "genus") 
	{
	  taxa <- sort(unique(occs$genus[occs$accepted_rank %in% c(rank, "species")]))
	  ranges <- array(NA, dim=c(length(taxa), 2), dimnames=list(taxa, c("FO", "LO")))
	  for (i in seq_len(nrow(ranges))) {
	    ranges[i,] <- rev(range(cols$collection_age[cols$collection_no %in% occs$collection_no[occs$genus==rownames(ranges)[i]]], finite=TRUE))
	    # thisDates <- cols[cols$collection_no %in% occs$collection_no[occs$accepted_name==rownames(ranges)[i]],"collection_age"]
	    # ranges[i,] <- c(max(thisDates, na.rm=TRUE), min(thisDates, na.rm=TRUE))
	  }
	} else if(rank == "species")
	  {
	  taxa <- sort(unique(occs$accepted_name[occs$accepted_rank==rank]))
	  ranges <- array(NA, dim=c(length(taxa), 2), dimnames=list(taxa, c("FO", "LO")))
	  for (i in seq_len(nrow(ranges))) {
	    ranges[i,] <- rev(range(cols$collection_age[cols$collection_no %in% occs$collection_no[occs$accepted_name==rownames(ranges)[i]]], finite=TRUE))
	    # thisDates <- cols[cols$collection_no %in% occs$collection_no[occs$accepted_name==rownames(ranges)[i]],"collection_age"]
	    # ranges[i,] <- c(max(thisDates, na.rm=TRUE), min(thisDates, na.rm=TRUE))
	  }
	} else {
	  print("Select species or genus for this.rank")
	  return()
	  }

	ranges
}

# getAges <- function(thisTaxon, occs) {
	# this <- occs[occs$taxon==thisTaxon,]
	# thisMat <- t(as.matrix(c(max(this$ageMax, na.rm=TRUE), max(this$ageMin, na.rm=TRUE), min(this$ageMax, na.rm=TRUE), min(this$ageMin, na.rm=TRUE))))
	# rownames(thisMat) <- thisTaxon
	# thisMat
# }

# getAges2 <- function(thisTaxon, occs) {
	# this <- occs[occs$taxon==thisTaxon,]
	# thisMat <- t(as.matrix(c(max(this$early_age, na.rm=TRUE), max(this$late_age, na.rm=TRUE), min(this$early_age, na.rm=TRUE), min(this$late_age, na.rm=TRUE))))
	# rownames(thisMat) <- thisTaxon
	# thisMat
# }

# makeRangesWithErrorFromOccurrences <- function(occs=NULL, filename=NULL) {
	# if (!is.null(filename)) occs <- read.csv(filename) else if (is.null(occs)) return(NA)
	# ranges <- t(sapply(unique(occs$taxon), getAges2, occs))
	# colnames(ranges) <- c("ageFO_max", "ageFO_min", "ageLO_max", "ageLO_min")
	# ranges
# }

# getDatesWithinError <- function(dates, random=TRUE) {
	# if (random) { array(cbind(FO=runif(nrow(dates), dates[,"ageFO_min"], dates[,"ageFO_max"]), LO=runif(nrow(dates), dates[,"ageLO_min"], dates[,"ageLO_max"])), dim=c(nrow(dates), 2), dimnames=list(rownames(dates), c("FO", "LO")))
	# } else cbind(FO=rowMeans(dates[,c("ageFO_min", "ageFO_max")]), LO=rowMeans(dates[,c("ageLO_min", "ageLO_max")]))
# }
	
getRandomRange <- function(species) {
	if (all(!is.finite(species))) {
		FO=NA
		LO=NA
	} else if (all(is.finite(species)) & species["ageFO_min"]==species["ageLO_min"] & species["ageFO_max"]==species["ageLO_max"]) {
		oneDate <- runif(n=1, min=species["ageFO_min"], max=species["ageLO_max"])
		FO <- oneDate
		LO <- oneDate
	} else {
		if (!is.finite(species["ageFO_min"]) & !is.finite(species["ageFO_max"])) { FO <- NA
		} else if (!is.finite(species["ageFO_min"])) { FO <- species["ageFO_max"]
		} else if (!is.finite(species["ageFO_max"])) { FO <- species["ageFO_min"]
		} else if (species["ageFO_min"]!=species["ageFO_max"]) { FO <- runif(n=1, min=species["ageFO_min"], max=species["ageFO_max"]) 
		} else FO <- species["ageFO_min"]
		
		if (!is.finite(species["ageLO_min"]) & !is.finite(species["ageLO_max"])) { LO <- NA
		} else if (!is.finite(species["ageLO_min"])) { LO <- species["ageLO_max"]
		} else if (!is.finite(species["ageLO_max"])) { LO <- species["ageLO_min"] 
		} else if (species["ageLO_min"]!=species["ageLO_max"]) { LO <- runif(n=1, min=species["ageLO_min"], max=species["ageLO_max"]) 
		} else LO <- species["ageLO_min"]
	}
	c(FO=FO, LO=LO)
}

getRandomRanges<-function(thisRanges) {
	thisRanges  <-  t(apply(thisRanges, 1, getRandomRange))
	dimnames(thisRanges) <- list(rownames(thisRanges), c("FO", "LO"))
	thisRanges
}

# getOneMidpointRange <- function(species) {
	# c(FO=mean(species["ageFO_min"], species["ageFO_max"]), LO=mean(species["ageLO_min"], species["ageLO_max"]))
# }

getMidpointRanges <- function(ranges) {
	# thisRanges <- t(apply(ranges[,c("ageFO_min", "ageFO_max", "ageLO_min", "ageLO_max")], 1, getOneMidpointRange))
	cbind(FO=rowMeans(ranges[,c("ageFO_min","ageFO_max")]), LO=rowMeans(ranges[,c("ageLO_min","ageLO_max")]))
	# dimnames(thisRanges) <- list(rownames(ranges), c("FO", "LO"))
	# thisRanges
}

getSingleTaxonRangeFromOccList <- function(taxon, thisOccs, indet.by.collection=FALSE) {
	if (all(c("ma_max","ma_min") %in% names(thisOccs))) {
		ages <- thisOccs[thisOccs$taxon==taxon, c("ma_max","ma_min")]
		if (indet.by.collection & grepl(pattern="sp.", taxon)) { return(array(data=cbind(ages[,1], ages[,1], ages[,2], ages[,2]), dim=c(nrow(ages),4), dimnames=list(paste(taxon, thisOccs[thisOccs$taxon==taxon, "collection_no"]), c("ageFO_max","ageFO_min","ageLO_max","ageLO_min"))))
		} else array(data=c(ageFO_max=max(ages, na.rm=TRUE), ageFO_min=max(ages$ma_min, na.rm=TRUE), ageLO_max=min(ages$ma_max, na.rm=TRUE), ageLO_min=min(ages, na.rm=TRUE)), dim=c(1,4), dimnames=list(taxon, c("ageFO_max","ageFO_min","ageLO_max","ageLO_min")))
	} else {
		ages <- thisOccs[thisOccs$taxon==taxon, c("max_ma","min_ma")]
		if (indet.by.collection & grepl(pattern="sp.", taxon)) { return(array(data=cbind(ages[,1], ages[,1], ages[,2], ages[,2]), dim=c(nrow(ages),4), dimnames=list(paste(taxon, thisOccs[thisOccs$taxon==taxon, "collection_no"]), c("ageFO_max","ageFO_min","ageLO_max","ageLO_min"))))
		} else array(data=c(ageFO_max=max(ages, na.rm=TRUE), ageFO_min=max(ages$min_ma, na.rm=TRUE), ageLO_max=min(ages$max_ma, na.rm=TRUE), ageLO_min=min(ages, na.rm=TRUE)), dim=c(1,4), dimnames=list(taxon, c("ageFO_max","ageFO_min","ageLO_max","ageLO_min")))
	# c(ageFO_max=max(ages, na.rm=TRUE), ageFO_min=max(ages$ma_min[ages$ma_max==max(ages, na.rm=TRUE)], na.rm=TRUE), ageLO_max=min(ages$ma_max[ages$ma_min==min(ages, na.rm=TRUE)], na.rm=TRUE), ageLO_min=min(ages, na.rm=TRUE), n.cols=length(unique(thisOccs$collection_no[thisOccs$taxon==taxon])))
	# c(FO=max(ages), LO=min(ages), dur=diff(range(ages)))
	}
}

list2array <- function(x) {
	m <- array(data=NA, dim=c(0,dim(x[[1]])[2]), dimnames=list(NULL, dimnames(x[[1]])[[2]]))
	for (i in seq_along(x)) m <- rbind(m, x[[i]])
	m
}

taxonRangesFromUndatedOccList <- function(thisOccs, indet.by.collection=FALSE) {
	# taxa <- unique(unlist(sapply(thisOccs, function(x) { if (!is.null(x)) as.character(occs[occs$occurrence_no%in%x,]$taxon) else NA })))
	ranges <- list2array(lapply(X=sort(unique(thisOccs$taxon)), FUN=getSingleTaxonRangeFromOccList, thisOccs=thisOccs, indet.by.collection=indet.by.collection))
	# rownames(ranges) <- sort(unique(thisOccs$taxon))
	ranges
}

taxonRangesFromDatedOccList <- function(thisOccs) {
	t(sapply(unique(thisOccs$taxon), function(x) c(FO=max(thisOccs$ma_rand[thisOccs$taxon==x], na.rm=TRUE), LO=min(thisOccs$ma_rand[thisOccs$taxon==x], na.rm=TRUE))))
}

############################################################################################################################################

getPresProbOneIntv <- function(thisInt, thisOccs, colSet=unique(thisOccs$collection_no)) {
	sampledBeforeIntv <- sort(unique(thisOccs$taxon[thisOccs$ma_mid > thisInt["ageBase"]]))
	sampledAfterIntv <- sort(unique(thisOccs$taxon[thisOccs$ma_mid <= thisInt["ageTop"]]))	# <= because if sampled before interval, and sampled at the last instant of this interval then it ranges through interval
	rangeThroughs <- sampledBeforeIntv[sampledBeforeIntv %in% sampledAfterIntv]
	sampledInIntv <- unique(thisOccs$taxon[thisOccs$collection_no %in% colSet & thisOccs$ma_mid > thisInt["ageTop"] & thisOccs$ma_mid <= thisInt["ageBase"]])
	sum(rangeThroughs %in% sampledInIntv)/length(rangeThroughs)
}

getPresProbPerIntv<-function(thisOccs, intList, colSet=unique(thisOccs$collection_no)) {
	sapply(intList, getPresProbOneIntv, thisOccs, colSet)
}

getPresProbOneTaxon<-function(thisTaxon, thisOccs, intList) {
	taxOccs<-thisOccs[thisOccs$taxon==thisTaxon,]
	taxIntv<-unique(sapply(taxOccs$ma_rand, function(x) which(sapply(intList, intervalContainsAge, date=x))))
	if (length(taxIntv)>2) (length(taxIntv)-2)/((max(taxIntv)-min(taxIntv))-1) else NA
}

getPresProbTaxa<-function(thisOccs, intList) {
	sapply(unique(thisOccs$taxon), getPresProbOneTaxon, thisOccs, intList)
}
