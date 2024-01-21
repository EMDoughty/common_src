server <- "https://paleobiodb.org/data1.2/"
# tax.out <- read.csv(URLencode(paste0("https://paleobiodb.org/data1.2/taxa/list.csv?taxon_name=", paste(gsub(pattern=" ", replacement="%20", x=tax.vec[119]), collapse=","), "&rel=current")))
# tax.out <- fromJSON(URLencode(paste0("https://paleobiodb.org/data1.2/taxa/list.json?taxon_name=", paste(gsub(pattern=" ", replacement="%20", x=tax.vec[43]), collapse=","), "&rel=current")))

getStuffWithoutWarnings <- function(this.URL) {
	# yy <- url(description=this.URL, method="libcurl")
	# on.exit(close(yy))

	on.exit(closeAllConnections())
	# url.con <- url(this.URL)
	dat <- try(read.csv(this.URL))
	if (class(dat)=="try-error") {
		closeAllConnections()
		return(NULL)
	}
	dat
}

getOneCurrentTaxon <- function(this.taxon, only.taxon_no = FALSE, keep.changed.rank=TRUE) {
	this.taxon <- as.character(this.taxon)
	if (!is.na(this.taxon) & this.taxon != "") {
		if (grepl(pattern=" sp[[:punct:]]", x=this.taxon)) this.taxon <- gsub(pattern=" sp[[:punct:]]", replacement="", x=this.taxon)
		if (grepl(pattern=" indet[[:punct:]]", x=this.taxon)) this.taxon <- gsub(pattern=" indet[[:punct:]]", replacement="", x=this.taxon)		# return (this.taxon)	# if this.taxon has punctuation (e.g., "Equus sp.") or "indet", then just return this.taxon
		this.URL <- URLencode(paste0(server, "taxa/single.csv?name=", this.taxon))
		this.names <- getStuffWithoutWarnings(this.URL)
		if (is.null(this.names$taxon_name)) {
			warning(paste0("The taxon \"", this.taxon, "\" cannot be found in the Paleobiology Database"), call.=FALSE)
			return(paste0("\"",this.taxon,"\""))
		} else {
			if (this.names$accepted_rank != this.names$taxon_rank & !keep.changed.rank) {
				warning(paste0("The rank taxon \"", this.taxon, "\" has been changed from ", this.names$taxon_rank," to ", this.names$accepted_rank, " (", this.names$accepted_name,")"), call.=FALSE)
				return(paste0("\"",this.taxon,"\""))
			} else if (only.taxon_no) this.names$accepted_no else as.character(this.names$accepted_name)
		}
	} else NA
}

getCurrentTaxa <- function(tax.vec, show.progress=TRUE, only.taxon_no=FALSE, keep.changed.rank=TRUE) { #, do.parallel=FALSE
	tax.vec.short <- unique(as.character(tax.vec))

	output <- vector()
	for (i in seq_along(tax.vec.short)) {
		if (show.progress) cat("\t**** Getting current taxon from PaleoBioDB\t", i, " of ", length(tax.vec.short), " (", round(i/length(tax.vec.short)*100, 1), "%)\r", sep="")
		output <- c(output, getOneCurrentTaxon(this.taxon=tax.vec.short[i], only.taxon_no=only.taxon_no, keep.changed.rank=keep.changed.rank))
	}
	if (show.progress) print("")
	name.mat <- cbind(input=tax.vec.short, output=output)
	if (!only.taxon_no) factor(name.mat[match(tax.vec, name.mat[,"input"]),"output"]) else as.numeric(name.mat[match(tax.vec, name.mat[,"input"]),"output"])
}


# getTaxonomyWithoutWarnings <- function(this.URL) {
	# yy <- url(description=this.URL, method="libcurl")
	# rawFile <- try(readLines(yy), silent=TRUE) # empty
	# close(yy)
	# if (class(rawFile)=="try-error") return(NULL) 
	# # if (any(sapply(rawFile, grepl, pattern="Warning:"))) {
		# # close(yy)
		# # warningMessages <- read.csv(textConnection(rawFile[grep("Warning:", rawFile)]), header=FALSE, stringsAsFactors=FALSE)
		# # sapply(warningMessages[,2], warning)
		# # dat <- read.csv(textConnection(rawFile[-(1:grep("Records:", rawFile))]))
	# # } else {
		# dat <- fromJSON(rawFile)$records
		# this.dat <- sapply (c("phl", "cll", "odl", "fml", "gnl", "nam"), function(x) if (x %in% colnames(dat)) dat[1,x] else NA)
	# # }			
	# this.dat
# }

getOneTaxonTaxonomy <- function(this.taxon) {
	this.taxon <- as.character(this.taxon)
	if (!is.na(this.taxon) & this.taxon != "") {
		if (grepl(pattern="[[:punct:]]", x=this.taxon)) return(array(data=NA, dim=c(1,6), dimnames=list(NULL, c("phylum", "class", "order", "family", "genus", "accepted_name"))))
		# if (grepl(pattern="[[:punct:]]", x=this.taxon)) return (array(data=NA, dim=c(1,6), dimnames=list(NULL, c("phl", "cll", "odl", "fml", "gnl", "nam"))))
		this.URL <- URLencode(paste0(server, "taxa/single.csv?taxon_name=", this.taxon, "&show=class"))
		this.names <- getStuffWithoutWarnings(this.URL)
		if (is.null(this.names$taxon_name)) {
			warning(paste0("The taxon \"", this.taxon, "\" cannot be found in the Paleobiology Database"), call.=FALSE)
			return(array(data=NA, dim=c(1,6), dimnames=list(NULL, c("phylum", "class", "order", "family", "genus", "accepted_name"))))
		} else as.matrix(this.names[,c("phylum", "class", "order", "family", "genus", "accepted_name")])
	} else return(array(data=NA, dim=c(1,6), dimnames=list(NULL, c("phylum", "class", "order", "family", "genus", "accepted_name"))))
	
}
# taxonomy <- t(sapply(gsub("_", " ", tax.vec), getOneTaxonTaxonomy))
# taxonomy <- t(sapply(gsub("_", " ", rownames(head(thisMat, 10))), getOneTaxonTaxonomy))
# for (this.taxon in gsub("_", " ", rownames(thisMat))) getOneTaxonTaxonomy(this.taxon)

getTaxonomyForTaxa <- function(tax.vec=NULL) {
	# require(jsonlite)
	# for (this.taxon in tax.vec) getOneTaxonTaxonomy(this.taxon)
	# m <- as.data.frame(t(sapply(head(tax.vec,32), getOneTaxonTaxonomy)), stringsAsFactors=TRUE)
	m <- as.data.frame(t(sapply(tax.vec, getOneTaxonTaxonomy)), stringsAsFactors=TRUE)
	colnames(m) <- c("phylum", "class", "order", "family", "genus", "species")
	m
}

getTaxonomyForOneBaseTaxon <- function(this.taxon) {
	this.URL <- URLencode(paste0(server, "taxa/list.csv?base_name=", this.taxon, "&taxon_status=all&show=class"))
	this.names <- getStuffWithoutWarnings(this.URL)
	this.names[,c("phylum", "class", "order", "family", "genus", "taxon_name")]
}

getFADateRange <- function(this.taxon) {
	occURL <- URLencode(paste0("http://testpaleodb.geology.wisc.edu/data1.2/occs/list.csv?limit=all&base_name=", this.taxon))
	occs <- getStuffWithoutWarnings(occURL)
	if (nrow(occs)>0) c(max(occs$early_age), max(occs$late_age)) else c(NA, NA)
}

getOccsFromCollection <- function(thisCol) {
	occURL <- URLencode(paste0("http://testpaleodb.geology.wisc.edu/data1.2/occs/list.csv?all_records&limit=all&coll_id=", thisCol))
	getStuffWithoutWarnings(occURL)
	# if (nrow(occs)>0) c(max(occs$early_age), max(occs$late_age)) else c(NA, NA)
}

appendTaxonNames <- function(occs, taxonomic.level=c("family", "subfamily", "tribe", "genus", "species"), keep.indet=FALSE, use.original=FALSE, thisSep=" ") {
	taxonomic.level <- match.arg(taxonomic.level)
	if (tolower(taxonomic.level)=="species") {
		if (!keep.indet) occs <- occs[occs$occurrence.species_name != "sp.",]
		occs$taxon <- factor(paste(occs$occurrence.genus_name, occs$occurrence.species_name, sep=thisSep))
	} else if (tolower(taxonomic.level)=="genus") { occs$taxon <- factor(occs$occurrence.genus_name) # for genus level 
	} else if (tolower(taxonomic.level)=="tribe") { occs$taxon <- factor(occs$Tribe) # for genus level 
	} else if (tolower(taxonomic.level)=="subfamily") { occs$taxon <- factor(occs$Subfamily) # for genus level 
	} else if (tolower(taxonomic.level)=="family") occs$taxon <- factor(occs$family_name) # for family level
	if (keep.indet) {
		levels(occs$taxon) <- c(levels(occs$taxon), paste("[taxon indet. ", seq_along(occs$taxon[occs$taxon==""]), "]", sep=""))
		occs$taxon[is.na(occs$taxon)] <- factor(paste("[taxon indet. ", seq_along(occs$taxon[occs$taxon==""]), "]", sep=""))
	} else occs <- occs[!is.na(occs$taxon),]
	occs
}

appendTaxonNames1.2 <- function(occs, taxonomic.level=c("family", "genus", "species"), keep.indet=FALSE, use.original=FALSE, thisSep=" ") {
	taxonomic.level <- match.arg(taxonomic.level)
	if (tolower(taxonomic.level)=="species") {
		if (!keep.indet) occs <- occs[grepl(pattern="species", x=occs$accepted_rank),]
		occs$taxon <- as.character(occs$accepted_name)
		occs$taxon[occs$accepted_rank == "subspecies"] <- sapply(strsplit(as.character(occs$accepted_name[occs$accepted_rank == "subspecies"]), split=" "), function(x) paste(x[1], x[2], sep=thisSep))
		occs$taxon[occs$accepted_rank == "genus"] <- paste(occs$accepted_name[occs$accepted_rank == "genus"], "sp.", sep=thisSep)
		occs$taxon[!(occs$accepted_rank %in% c("subspecies", "species", "genus"))] <- paste(occs$accepted_name[!(occs$accepted_rank %in% c("subspecies", "species", "genus"))], "indet.", sep=thisSep)
	} else if (tolower(taxonomic.level)=="genus") { occs$taxon <- factor(occs$genus) # for genus level 
	# } else if (tolower(taxonomic.level)=="tribe") { occs$taxon <- factor(occs$Tribe) # for genus level 
	# } else if (tolower(taxonomic.level)=="subfamily") { occs$taxon <- factor(occs$Subfamily) # for genus level 
	} else if (tolower(taxonomic.level)=="family") occs$taxon <- factor(occs$family) # for family level
	# if (keep.indet) {
		# levels(occs$taxon) <- c(levels(occs$taxon), paste("[taxon indet. ", seq_along(occs$taxon[occs$taxon==""]), "]", sep=""))
		# occs$taxon[is.na(occs$taxon)] <- factor(paste("[taxon indet. ", seq_along(occs$taxon[occs$taxon==""]), "]", sep=""))
	# } else occs <- occs[!is.na(occs$taxon),]
	occs$taxon <- factor(occs$taxon)
	occs
}

appendMissingOrders <- function(occs) {
	famList <- read.csv("~/Dropbox/Code/common_dat/taxonomy.csv", stringsAsFactors=TRUE, strip.white=TRUE)
	# missingOrders <- read.csv("~/Dropbox/Code/common_dat/PBDBMissingOrders.csv", stringsAsFactors=FALSE, strip.white=TRUE)
	occs$order_name[is.na(occs$order_name)] <- factor(famList$occurrences.order_name[match(occs$family_name[is.na(occs$order_name)], famList$occurrences.family_name)], levels=unique(c(levels(occs$order_name), unique(famList$occurrences.order_name))))
	occs$order_name[is.na(occs$order_name)] <- factor(famList$occurrences.order_name[match(occs$occurrence.genus_name[is.na(occs$order_name)], famList$occurrences.genus_name)], levels=unique(c(levels(occs$order_name), unique(famList$occurrences.order_name))))
	# occs[!is.na(occs$family_name) & 
	     # occs$family_name %in% missingOrders$Family[!is.na(missingOrders$Family)], ]$order_name <- 
	     # missingOrders[match(occs[!is.na(occs$family_name) & occs$family_name %in% missingOrders$Family, ]$family_name, missingOrders$Family), 1]
	# occs[!is.na(occs$occurrence.genus_name) & 
	     # occs$occurrence.genus_name %in% missingOrders$Genus[!is.na(missingOrders$Genus)], ]$order_name <- 
	     # missingOrders[match(occs[!is.na(occs$occurrence.genus_name) & occs$occurrence.genus_name %in% missingOrders$Genus, ]$occurrence.genus_name, missingOrders$Genus), 1]
	# # occs[occs$occurrence.genus_name!="" & 
	     # # occs$occurrence.genus_name %in% missingOrders$Genus, ]$order_name <- 
	     # # missingOrders[match(occs[occs$occurrence.genus_name!="" & occs$occurrence.genus_name%in%missingOrders$Genus,]$occurrence.genus_name, missingOrders$Genus), 1]
	occs
}

makeOneAcceptedRankFromNOW <- function(x) {
	if (grepl(pattern="indet.", x=x["SPECIES"]) | grepl(pattern="sp.", x=x["SPECIES"])) {
		if (grepl(pattern="indet.", x=x["genus"]) | grepl(pattern="gen.", x=x["genus"], fixed=TRUE)) {
			if (grepl(pattern="\\N", x=x["SUBFAMILY"])) {
				if (grepl(pattern="indet.", x=x["family"])) "order" else "family"
			} else "subfamily"
		} else "genus"
	} else "species"
}

makeOneAcceptedNameFromNOW <- function(x) {
	if (grepl(pattern="indet", x=x["SPECIES"]) | grepl(pattern="sp.", x=x["SPECIES"])) {
		if (grepl(pattern="indet", x=x["genus"]) | grepl(pattern="gen.", x=x["genus"])) {
			if (grepl(pattern="indet", x=x["family"])) x["order"] else x["family"]				
		} else x["genus"]
	} else paste(x["genus"], x["SPECIES"], sep=" ")
}

now2paleoDB <- function(db, maxCol=0, maxOcc=0, maxTax=0, fetch.pbdb.accepted_no=FALSE) {
	colnames(db)[colnames(db)=="LIDNUM"] <- "collection_no_NOW"
	colnames(db)[colnames(db)=="NAME"] <-  "collection_name"
	colnames(db)[colnames(db)=="LAT"] <- "lat"
	colnames(db)[colnames(db)=="LONG"] <- "lng"
	colnames(db)[colnames(db)=="MAX_AGE"] <- "max_ma"
	colnames(db)[colnames(db)=="MIN_AGE"] <- "min_ma"
	colnames(db)[colnames(db)=="COUNTRY"] <-  "cc"
	colnames(db)[colnames(db)=="STATE"] <- "state"
	colnames(db)[colnames(db)=="COUNTY"] <- "county"
	colnames(db)[colnames(db)=="LOC_SYNONYMS"] <- "collection_aka"
	colnames(db)[colnames(db)=="SIDNUM"] <- "accepted_no"
	colnames(db)[colnames(db)=="ORDER"] <- "order"
	colnames(db)[colnames(db)=="FAMILY"] <- "family"
	colnames(db)[colnames(db)=="GENUS"] <- "genus"
	  # need to check the following taxonomy for indeterminate assignments - not sure how cf. etc. are handled
	  # "UNIQUE", "TAXON_STATUS", "ID_STATUS", "ADD_INFO", "SOURCE_NAME", 
	  # "SVLENGTH", "BODYMASS", "SXDIMSZE", "SXDIMDIS", "TSHM",
	  # "TCRWNHT", "CROWNTYP", "DIET_1", "DIET_2", "DIET_3",
	  # "LOCOMO1", "LOCOMO2", "LOCOMO3",
	  # "SPCOMMENT", "SYNONYMS")

	db$occurrence_no <- seq(from=(maxOcc + 1), to=(maxOcc + nrow(db)), by=1)
	db$collection_no <- as.numeric(as.character(db$collection_no_NOW)) + maxCol # on col number is "We use Pinnipedia on order level for practical reasons." therefore the column is brought in as a factor
	db$accepted_rank <- apply(db[,c("order", "family", "SUBFAMILY", "genus", "SPECIES")], 1, makeOneAcceptedRankFromNOW)
	db$accepted_name <- apply(db[,c("order", "family", "genus", "SPECIES")], 1, makeOneAcceptedNameFromNOW)
	db$accepted_no <- db$accepted_no + maxTax			# need to check if these names are already in the PBDB
	if (fetch.pbdb.accepted_no) {
		pbdb.numbers <- getCurrentTaxa(db$accepted_name, only.taxon_no=TRUE)
		# cbind(db$accepted_no[is.finite(pbdb.numbers)], pbdb.numbers[is.finite(pbdb.numbers)])
		db$accepted_no[is.finite(pbdb.numbers)] <- pbdb.numbers[is.finite(pbdb.numbers)]
	}
	return(db)
}

# getCollectionDatesFromOccs <- function(occs, age.determination=c("midpoint", "random")) {
	# age.determination <- match.arg(age.determination)
	# cols <- unique(occs[, c("collection_no", "max_ma", "min_ma")])
	# if (age.determination=="random") thisCols <- cbind(cols[,1], apply(cols[,2:3], 1, function(x) { runif(1, max=x[1], min=x[2]) })) else thisCols <- cbind(cols[,1], rowMeans(cols[,2:3]))
	# array(thisCols[,2], dimnames=list(thisCols[,1]))
# }

# getOccurrenceDatesFromOccs <- function(occs, age.determination=c("midpoint", "random")) {
	# age.determination <- match.arg(age.determination)
	# thisCols <- getCollectionDatesFromOccs(occs=occs[, c("collection_no", "max_ma", "min_ma")], age.determination=age.determination)
	# thisCols[match(occs$collection_no, as.numeric(rownames(thisCols)))]
# }

dateOccsWithAEO <- function(occs) {
	aeo <- read.table("~/Dropbox/Code/java/data/PBDBDates/11Nov07_tcdm.collnoages.txt", header=TRUE, strip.white=TRUE)
	aeo <- aeo[aeo$ma_max >= 0,]
	for (i in seq_along(aeo$collection_no)) {
		if (aeo$collection_no[i] %in% occs$collection_no) {
			occs[occs$collection_no==aeo$collection_no[i], c("ma_max", "ma_min")] <- aeo[i, c("ma_max", "ma_min")]
		}
	}
	occs
}

dateColsWithAEO <- function(occs, age.determination=c("midpoint", "random")) {
	age.determination <- match.arg(age.determination)
	cols <- unique(occs[, c("collection_no", "ma_max", "ma_min")])
	aeo <- read.table("~/Dropbox/Code/java/data/PBDBDates/11Nov07_tcdm.collnoages.txt", header=TRUE, strip.white=TRUE)
	aeo <- aeo[aeo$ma_max >= 0,]
	cols[cols$collection_no %in% aeo$collection_no, 2:3] <- aeo[match(cols$collection_no[cols$collection_no %in% aeo$collection_no], aeo$collection_no), 2:3]
	if (age.determination=="random") thisCols <- cbind(cols[,1], apply(cols[,2:3], 1, function(x) { runif(1, max=x[1], min=x[2]) })) else thisCols <- cbind(cols[,1], rowMeans(cols[,2:3]))
	# thisCols[match(occs$collection_no, thisCols[,1]), 2]
	array(thisCols[,2], dimnames=list(thisCols[,1]))
}

getNOWLocalityCodesFromPBDBCollectionNo <- function(collection_no.vec, this.file="~/Dropbox/Code/R/dentalMeasurements/dat/pbdb_NOW_map.csv") {
	this.map <- read.csv(this.file)
	this.map[this.map==""] <- NA
	loc.vec <- this.map$NOW_loc[match(collection_no.vec, this.map$collection_no)]
	factor(loc.vec)
}

getNOWDatesMatFromLocCodes <- function(code.vec, this.file="~/Dropbox/Code/R/dentalMeasurements/dat/NOW_loc_dates.csv", dates.only=FALSE) {
  date.mat <- read.csv(this.file)
  date.mat$MAX_AGE[date.mat$MAX_AGE==64.1] <- 65.83	### makes any collection from Pu3 from combined Pu2/Pu3. as of 1/12/2024 the bound is set to 65.84 
  date.mat$MIN_AGE[date.mat$MIN_AGE==64.1] <- 64.81	### makes any collection from Pu2 from combined Pu2/Pu3. as of 1/12/2024 the bound is set to 64.75 
  missing.codes <- sort(unique(code.vec[!code.vec%in%date.mat$LOC_SYNONYMS]))
  if (length(missing.codes)>0) warning(paste("These collection codes were not found in the map", paste0(as.character(missing.codes), sep="\n")))
  this.dates <- date.mat[match(code.vec, date.mat$LOC_SYNONYMS), c("BFA_MAX", "BFA_MIN", "MAX_AGE", "MIN_AGE", "CHRON")]
  if (dates.only) return(this.dates[,c("MAX_AGE", "MIN_AGE")]) else return (this.dates)
}
