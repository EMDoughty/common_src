#### occupancy_src
if("sp" %in% installed.packages()[,"Package"]) {require(sp)} else {install.packages("sp")}
if("geosphere" %in% installed.packages()[,"Package"]) {require(geosphere)} else {install.packages("geosphere")}
if("nloptr" %in% installed.packages()[,"Package"]) {require(nloptr)} else {install.packages("nloptr")}
if("pracma" %in% installed.packages()[,"Package"]) {require(pracma)} else {install.packages("pracma")}

######################################################################################################################################################
##### parameters for integration routines
######################################################################################################################################################
# this.tol <- max(50*.Machine$double.eps, 0.5e-28)	# .Machine$double.eps^0.25
this.tol <- .Machine$double.eps^0.25
max.intervals <- 1e2L

######################################################################################################################################################
##### function to place occs into grid cells labeled with unique identifying label
######################################################################################################################################################

getGridCellsFromOccs <- function(lat.increment = 1, lng.increment = 1) {
	cell.label <- 1
	grid.cells <- array(data = NA, dim = c(0, 5), dimnames = list(NULL, c("label", "lat.start", "lat.end", "lng.start", "lng.end")))
	for (this.lat in seq(from = ceiling(max(occs$lat)), to = floor(min(occs$lat)), by = -lat.increment)) {
		for (this.lng in seq(from = floor(min(occs$lng)), to = ceiling(max(occs$lng)), by = lng.increment)) {
			grid.cells <- rbind(grid.cells, c(label = cell.label, lat.start = this.lat, lat.end = this.lat + lat.increment, lng.start = this.lng, lng.end = this.lng + lng.increment))
			cell.label <- cell.label + 1
		}
	}
	grid.vec <- apply(X=occs[, c("lat", "lng")], MARGIN=1, FUN=function(x) grid.cells[grid.cells[, "lat.start"] <= x[1] & grid.cells[, "lat.end"] > x[1] & grid.cells[, "lng.start"] <= x[2] & grid.cells[, "lng.end"] > x[2], "label"])
	list(grid.mat=data.frame(grid.cells), grid.vec=grid.vec)
}

######################################################################################################################################################
##### functions to calculate MPWD among collections
######################################################################################################################################################

getMPWDColsOneInt <- function(this.intv, this.grid.mat) {
	this.cols <- sort(unique(occs$grid.cell[occs$occurrence_no %in% this.intv]))
	this.grid <- this.grid.mat[this.grid.mat$label %in% this.cols,]
	col.centroids <- cbind(rowMeans(this.grid[,c("lat.start", "lat.end")]), rowMeans(this.grid[,c("lng.start", "lng.end")]))
	m <- mean(spDists(x=col.centroids, longlat=TRUE))
	this.hull <- col.centroids[chull(x=col.centroids),]
	this.area <- areaPolygon(this.hull[,2:1])
	c(MPWD=m, area.hull=this.area)
}

getMPWDColsOneRep <- function(this.rep, this.grid.mat) {
	data.frame(t(sapply(this.rep, getMPWDColsOneInt, this.grid.mat=this.grid.mat)))
}

######################################################################################################################################################
##### functions to calculate occupancy (as k/n)
######################################################################################################################################################

getOccupancyOneInterval <- function(this.intv, site.type=match("PBDB", "NOW", "grid"), bigList=NULL) {
	if (settings$this.rank %in% names(occs)) tax.label <- settings$this.rank else tax.label <- "accepted_name"
	if (site.type=="NOW") { col.label <- "NOW_loc" 
	} else if (site.type=="PBDB") { col.label <- "collection_no"
	} else if (site.type=="grid") col.label <- "grid.cell"

	tax.vec <- sort(unique(occs[occs$occurrence_no %in% this.intv & occs$accepted_rank %in% rank.vec[1:which(rank.vec==settings$this.rank)], tax.label]))
	if (!is.null(bigList)) tax.vec <- tax.vec[tax.vec  %in% bigList[ , tax.label]]
	if (length(tax.vec) < 1) return(NA)

	nfinds.vec <- sapply(tax.vec, function(x) length(unique(occs[occs$occurrence_no %in% this.intv & occs[ ,tax.label]==x, col.label])))
	names(nfinds.vec) <- tax.vec
	ncols.total <- length(unique(occs[occs$occurrence_no %in% this.intv, col.label]))
	list(occupancy=nfinds.vec/ncols.total, n.finds=nfinds.vec, ncols.tot=ncols.total)
}

getOccupancyOneRep <- function(this.rep, site.type=match("PBDB", "NOW", "grid"), bigList=NULL) {
	lapply(this.rep, getOccupancyOneInterval, site.type=site.type, bigList=bigList)
}

######################################################################################################################################################
##### functions to calculate mean p from Foote 2016
######################################################################################################################################################

# p = occupancy probability
# f(p) = p varies
# k = number of sites occupied by a species
# n = total number of sites
# S = observe number of species
################################################################################

# A1 <- log((choose(n=n, k=k)*(p^k)*((1-p)^(n-k))))
# A1 <- dbinom(x = k, size = n, prob = p, log = TRUE)

################################################################################

# A2.log <- function(this.k=NULL, this.n=NULL, this.p=NULL) {
	# dbinom(x=this.k, size=this.n, prob=this.p, log = TRUE) - log(1 - ((1 - this.p)^this.n))
# }

################################################################################

# A3 <- sum(A2)

################################################################################

fpB <- function(this.p=NULL, this.par=NULL, this.k=NULL, this.n=NULL) {
	# must not be logarithm to be integrated, therefore exp() on sum of logs 
	# cat(this.par, this.k, "\n")
	exp(dlnorm(this.p, meanlog=this.par[1], sdlog=this.par[2], log=TRUE) + dbinom(x=this.k, size=this.n, prob=this.p, log=TRUE))
}

fpB.vectorized <- Vectorize(FUN=fpB, vectorize.args="this.p")

# x <- seq(from = 0, to = 1, length.out = 100001)
# curve(expr = fpB(x, this.par=this.par, this.k=this.k, this.n=this.n), from = 0, to=this.par[3], n=1e6)
# curve(expr = dlnorm(x, meanlog = this.par[1], sdlog = this.par[2], log = FALSE), from=0, to=this.par[3], n=1e6)

A6.log <- function(this.k=NULL, this.par=NULL, this.n=NULL) {
	# if(any(this.par==0)) return(NA_real_)
	# cat(this.par, this.k, "\n")

	A6.num <- try(expr=integral(fun=fpB.vectorized, xmin=0, xmax=this.par[3], no_intervals=max.intervals, reltol=this.tol, this.par=this.par, this.k=this.k, this.n=this.n), silent= TRUE)

	if (inherits(A6.num, "try-error")) {
		# cat("integral() failed on A6.num", this.par, "k=", this.k, "n=", this.n, "\n")
		A6.num <- try(expr=integrate(f=fpB.vectorized, lower=0, upper=this.par[3], this.par=this.par, this.k=this.k, this.n=this.n, subdivisions=max.intervals, rel.tol=this.tol), silent= TRUE)
		if (inherits(A6.num, "try-error")) {
			cat("integrate() failed on A6.num", this.par, "k=", this.k, "n=", this.n,"\n")
			warning(as.vector(A6.num))
			return(NA_real_)			### if neither integral worked, return NA)
		} else A6.num <- A6.num$value
	} 

	A6.denom <- try(expr=integral(fun=dlnorm, xmin=0, xmax=this.par[3], no_intervals=max.intervals, reltol= this.tol, meanlog=this.par[1], sdlog=this.par[2], log=FALSE), silent= TRUE)	# dlnorm must not be logarithm to be integrated

	if (inherits(A6.denom, "try-error")) {
		# cat("integral() failed on A6.denom", this.par, "k=", this.k, "n=", this.n, "\n")
		A6.denom <- try(expr=integrate(f=dlnorm, lower=0, upper=this.par[3], meanlog=this.par[1], sdlog=this.par[2], log=FALSE, subdivisions=max.intervals, rel.tol=this.tol), silent= TRUE)	# dlnorm must not be logarithm to be integrated
		if (inherits(A6.denom, "try-error")) {
			cat("integrate() failed on A6.denom", this.par, "k=", this.k, "n=", this.n,"\n")
			warning(as.vector(A6.denom))
			return(NA_real_)
		} else A6.denom <- A6.denom$value
	}
	
	log(A6.num) - log(A6.denom)
}

################################################################################

A7.log <- function(this.k=NULL, this.par=NULL, this.n=NULL) {
	A6.log(this.k=this.k, this.par=this.par, this.n=this.n) - log(1 - exp(A6.log(this.k=0, this.par=this.par, this.n=this.n)))
}

################################################################################

# A8 <- function(this.par=NULL, k.vec=NULL, this.n=NULL) {
	# if (length(this.par)==2) this.par <- c(this.par, 1)		#### if pmax is fixed (at 1.0), append pmax=1.0 to this.par, to integrate from 0 to pmax=1
	# lnL <- sum(sapply(k.vec, A7.log, this.par=this.par, this.n=this.n))
	# cat(this.par, lnL, "\n")
	# -lnL
# }

A9.log <- function(this.par=NULL, k.vec=NULL, this.n=NULL) {
	if (length(this.par)==2) this.par <- c(this.par, 1)		#### if pmax is fixed (at 1.0), append pmax=1.0 to this.par, to integrate from 0 to pmax=1
	this.table <- table(k.vec)
	p.vec <- vector("double", length=length(this.table))
	for (i in seq_along(this.table)) p.vec[i] <- A7.log(this.k=as.numeric(names(this.table))[i], this.par=this.par, this.n=this.n) * as.integer(this.table[i])
	lnL <- sum(p.vec)
	# cat("A9:", this.par, lnL, "\n")
	-lnL
}

################################################################################

fitTruncatedLogNormalOneIntv <- function(this.intv=NULL, fix.pmax=FALSE) {
	# meanlog <- mean(log(this.intv$occupancy))	# mean of log(k/n) or log of geometric mean of k/n
	# sdlog <- sd(log(this.intv$occupancy))		# sd of log(k/n)
	par.init <- c(mean(log(this.intv$occupancy)), sd(log(this.intv$occupancy)))
	# if (!fix.pmax) par.init <- c(par.init, max(this.intv$occupancy))
	this.lower <- c(-Inf, 1e-15)
	this.upper <- c(Inf, Inf)

	if (!fix.pmax) {
		par.init <- c(par.init, runif(n=1, min=max(this.intv$occupancy), max=1.0))
		this.lower <- c(this.lower, 1e-6)
		this.upper <- c(this.upper, 1)
	}
	# this.fit <- optim(par=par.init, fn=A9, k.vec=this.intv$n.finds, this.n=this.intv$ncols.tot, method="L-BFGS-B", lower=this.lower, upper=this.upper, control=list(trace=2))
	# this.fit <- nlminb(start= par.init, objective=A9, k.vec=this.intv$n.finds, this.n=this.intv$ncols.tot, lower=this.lower, upper= this.upper, control=list(eval.max=10000, iter.max=1000, trace=1))
	this.fit <- nloptr(x0=par.init, eval_f=A9.log, lb=this.lower, ub=this.upper, opts=list(algorithm = "NLOPT_LN_SBPLX", xtol_rel=this.tol, maxeval=10000), k.vec=this.intv$n.finds, this.n=this.intv$ncols.tot)
	this.fit
}

################################################################################
#### calculate mean value of probability distribution
################################################################################

pfp <- function(this.p=NULL, this.par=NULL) {
	this.p * dlnorm(this.p, meanlog = this.par[1], sdlog = this.par[2], log = FALSE)
}

A10 <- function(this.par=NULL) {
	if (length(this.par)==2) this.par <- c(this.par, 1)
	
	try(A10.num <- integrate(f=Vectorize(pfp, vectorize.args="this.p"), lower=0, upper=this.par[3], this.par=this.par), silent=TRUE)
	if (inherits(A10.num, "try-error")) {
		A10.num <- try(integral(fun = Vectorize(pfp, vectorize.args="this.p"), xmin = 0, xmax = this.par[3], this.par = this.par), silent=TRUE)
		if (inherits(A10.num , "try-error")) {
			warning(as.vector(A10.num))
			return(NA_real_)
		} 
	} else A10.num <- A10.num$value

	try(A10.denom <- integrate(f = dlnorm, lower = 0, upper = this.par[3], meanlog=this.par[1], sdlog=this.par[2], log=FALSE), silent=TRUE)
	if (inherits(A10.denom, "try-error")) {
		A10.denom <- try(integral(fun = dlnorm, xmin = 0, xmax = this.par[3], meanlog = this.par[1], sdlog = this.par[2], log = FALSE), silent=TRUE)
		if (inherits(A10.denom, "try-error")) {
			warning(as.vector(A10.denom))
			return(NA_real_)
		} 
	} else A10.denom <- A10.denom $value
	log(A10.num) - log(A10.denom)
}
	
################################################################################

getPMeanFromResultsOneRep <- function (rez.this.rep=NULL) {
	sapply(rez.this.rep, function(x) A10(this.par=x$solution))
}
