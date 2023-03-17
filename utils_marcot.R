# common source utilities - JDM
require(compiler)

listifyMatrixByRow <- function(m) {
	intList<-list()
	for (i in 1:nrow(m)) intList[[i]] <- data.matrix(m)[i,]
	intList
}

listifyMatrixByColumn <- function(m, do.parallel=FALSE) {
	if (do.parallel) {
		require (parallel)
		mclapply(as.data.frame(m), function(x) { x }, mc.cores=detectCores()-2)
	} else lapply(as.data.frame(m), function(x) x)
}

makeMatrixFromList <- function(this.list) {
  m <- this.list[[1]]
  for (i in seq(from=2, to=length(this.list))) m <- merge(m, this.list[[i]], all=TRUE, sort=FALSE)
  m
}

alphaColor <- function(color, alpha) {
	sapply(color, function(color) rgb(col2rgb(color)[1], col2rgb(color)[2], col2rgb(color)[3], alpha*255, maxColorValue=255))
}

getPosArray <- function(x, y) {
	posArray <- array (NA, dim=length(x))
	posArray[y < 0 & abs(y/x) >= 1] <- 1
	posArray[x < 0 & abs(y/x) < 1] <- 2
	posArray[y >= 0 & abs(y/x) >= 1] <- 3
	posArray[x >= 0 & abs(y/x) < 1] <- 4
	posArray
}

# getRateClassList <- function(nRateClasses, totalClasses, ranks=NULL) {
	# getRateClassListMachine(thisClass=1, index=1, nRateClasses, rateClass=vector(mode="double", length=totalClasses), classList=list(), ranks=ranks)
# }

# getRateClassListMachine <- function(thisClass, index, nRateClasses, rateClass, classList, ranks=NULL) {
	# for (i in index : (length(rateClass) - (nRateClasses-thisClass))) {
		
		# if (is.null(ranks)) { rateClass[i] <- thisClass 
		# } else rateClass[which(ranks==i)] <- thisClass
		
		# if (thisClass<nRateClasses) { classList <- getRateClassListMachine(thisClass=(thisClass+1), index=(i+1), nRateClasses=nRateClasses, rateClass=rateClass, classList=classList, ranks=ranks)
		# } else if (i==length(rateClass)) classList[[length(classList)+1]] <- rateClass
	# }	
	# classList
# }

# getIntervalBreaks <- function(nclasses, nbins) {
	# getBreaksMachine(nclasses=nclasses, nbins=nbins)	#array(dim=c(0, (nclasses-1)))
# }

# getAllBreaks <- function(nclasses, nbins=nbins, breakMat=data.frame(), thisBreak=1, index=1) {
	# if (nclasses > nbins) { simpleError("Number of classes cannot exceed the number of bins")
	# } else if (thisBreak < (nclasses-1)) { 
		# thisMat <- data.frame()
		# for (i in seq(from=index, to= ((nbins - 1) - ((nclasses - 1) - thisBreak)))) thisMat <- rbind(thisMat, cbind(i, getAllBreaks(nclasses=nclasses, nbins=nbins, breakMat=thisMat, thisBreak=thisBreak+1, index=i+1)))
		# return(thisMat)
	# }	
	# return(seq(from=index, to=(nbins - 1)))
# }

getAllBreaks_matrix <- function(nclasses, nbins=nbins, thisBreak=1, index=1) {
	if (nclasses > nbins) { simpleError("Number of classes cannot exceed the number of bins")
	} else if (thisBreak < (nclasses-1)) { 
		thisMat <- vector()
		for (i in index : ((nbins - 1) - ((nclasses - 1) - thisBreak))) thisMat <- rbind(thisMat, cbind(i, getAllBreaks_matrix(nclasses=nclasses, nbins=nbins, thisBreak=thisBreak+1, index=i+1)))
		return(thisMat)
	}	
	as.matrix(index : (nbins - 1), ncol=1)
}

cmp_getAllBreaks_matrix <- cmpfun(getAllBreaks_matrix)


# def choose_iter(elements, length):
    # for i in xrange(len(elements)):
        # if length == 1:
            # yield (elements[i],)
        # else:
            # for next in choose_iter(elements[i+1:len(elements)], length-1):
                # yield (elements[i],) + next


# incltxt <- '
  # int fibonacci(const int x) {
     # if (x == 0) return(0);
     # if (x == 1) return(1);
     # return (fibonacci(x - 1)) + fibonacci(x - 2);
  # }'

# ## now use the snippet above as well as one argument conversion
# ## in as well as out to provide Fibonacci numbers via C++
# fibRcpp <- cxxfunction(signature(xs="int"),
                       # plugin="Rcpp",
                       # incl=incltxt,
                       # body
                       # ='
   # int x = Rcpp::as<int>(xs);
   # return Rcpp::wrap( fibonacci(x) );
# ')

gmean <- function(x) {
	if (all(is.finite(x))) prod(x)^(1 / length(x)) else NA
}