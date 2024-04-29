source("~/Dropbox/Code/R/common_src/occupancy_src.R")

######################################################################################################################################################
##### main
######################################################################################################################################################

# this.strat <- unique(occs[,c("collection_no", "NOW_loc", "stratgroup", "formation", "member")])
# # this.strat <- this.strat[order(this.strat$formation, this.strat$stratgroup, this.strat$member),]
# write.csv(this.strat, file="~/Desktop/pbdb_strat_all.csv")

# this.intv <- repIntOccs.master[[1]][[1]]
# this.rep <- repIntOccs.master[[1]]
# getOccupancyOneInterval(this.intv=this.intv, site.type="NOW", bigList=bigList)
# getOccupancyOneRep(this.rep=this.rep, site.type="NOW", bigList=bigList)

this.occs <- repIntOccs.master
# this.occs <- repIntOccs

this.grid.cells <- getGridCellsFromOccs(lat.increment=0.125, lng.increment=0.125) 
occs$grid.cell <- this.grid.cells$grid.vec

require(parallel)
# m.PBDB <- mclapply(this.occs, getOccupancyOneRep, site.type="PBDB", bigList=bigList, mc.cores=detectCores()-2)
# m.NOW <- mclapply(this.occs, getOccupancyOneRep, site.type="NOW", bigList=bigList, mc.cores=detectCores()-2)
system.time(m.grid <- mclapply(this.occs, getOccupancyOneRep, site.type="grid", bigList=bigList, mc.cores=detectCores()-2))

##### if bigList is NULL, then it will do for all taxa in this.occs (=all Mammalia)
# m.grid.mammals <- mclapply(this.occs, getOccupancyOneRep, site.type="grid", bigList=NULL, mc.cores=dete					ctCores()-2)

	# par(mfrow=c(ceiling(sqrt(nrow(intervals))), ceiling(sqrt(nrow(intervals)))), mar=c(0,0,0,0))
	# m.counts <- sapply(lapply(m[[1]], hist, breaks=seq(from=0, to=max(unlist(m)), by=0.01), xaxt="n", yaxt="n"), function(x) x$counts)
	# # m.counts/colSums(m.counts)

# par(mfrow=c(3,1))
##### switch to choose which occupancy metric gets plotted
# this.m <- m.PBDB
# this.m <- m.NOW
this.m <- m.grid
# this.m <- m.grid.mammals

######################################################################################################################################################
##### plot occupancy over time
######################################################################################################################################################
m.mean <- sapply(this.m, function(x) sapply(x, function(y) if(length(y)>1) mean(y$occupancy) else NA))
m.quart <- t(apply(m.mean, MARGIN=1, FUN=quantile, probs=c(0.025, 0.25, 0.50, 0.75, 0.975), na.rm=TRUE))
m.quart[is.na(m.quart)] <- 0
m.quart <- log(m.quart)

	plot(rowMeans(intervals), m.quart[,5], type="n", pch=21, col="dodgerblue4", bg="dodgerblue1", xlim=rev(range(intervals)), ylim<-c(-5, max(m.quart)), xlab="time (Ma)", ylab="mean interval occupancy")
	overlayCzTimescale(include.intervals=c("epoch"), plot.NALMA = TRUE, plot.NALMA_sub = TRUE)
	polygon(x=c(rowMeans(intervals), rev(rowMeans(intervals))), y=c(m.quart[,1], rev(m.quart[,5])), col=adjustcolor("dodgerblue4", alpha.f = 0.25), border=NA)
	lines(rowMeans(intervals), m.quart[,3], type="o", pch=21, col="dodgerblue4", bg="dodgerblue1")

######################################################################################################################################################
##### plot number of collections over time
######################################################################################################################################################
ncols.mean <- sapply(this.m, function(x) sapply(x, function(y) if(length(y)>1) y$ncols.tot else NA))
ncols.quart <- t(apply(ncols.mean, MARGIN=1, FUN=quantile, probs=c(0.025, 0.25, 0.50, 0.75, 0.975), na.rm=TRUE))
ncols.quart <- log(ncols.quart)
ncols.quart[is.na(ncols.quart)] <- 0

	plot(rowMeans(intervals), ncols.quart[,5], type="n", pch=21, col="firebrick4", bg="firebrick1", xlim=rev(range(intervals)), ylim=c(0,max(ncols.quart, na.rm=TRUE)), xlab="time (Ma)", ylab="no. collections")
	overlayCzTimescale(include.intervals=c("epoch"), plot.NALMA = TRUE, plot.NALMA_sub = TRUE)
	polygon(x=c(rowMeans(intervals), rev(rowMeans(intervals))), y=c(ncols.quart[,1], rev(ncols.quart[,5])), col=adjustcolor("firebrick4", alpha.f = 0.25), border=NA)
	lines(rowMeans(intervals), ncols.quart[,3], type="o", pch=21, col="firebrick4", bg="firebrick1")

######################################################################################################################################################
##### correlation between number of collections and mean occupancy
######################################################################################################################################################
cor.test(diff(rev(ncols.quart[,3])), diff(rev(m.quart[,3])), method="spearman")
plot(diff(rev(ncols.quart[,3])), diff(rev(m.quart[,3])))
abline(lm(diff(rev(m.quart[,3])) ~ diff(rev(ncols.quart[,3]))))

######################################################################################################################################################
##### plot mean number of taxa per collection
######################################################################################################################################################
ntax.mean <- sapply(this.m, function(x) sapply(x, function(y) if(length(y)>1) length(y$occupancy) else NA))
	plot(rowMeans(intervals), ntax.mean[,1]/ncols.mean[,1], type="o", pch=21, col="firebrick4", bg="firebrick1", xlim=rev(range(intervals)), ylim=c(0, max(ntax.mean[,1]/ncols.mean[,1], na.rm=TRUE)), xlab="time (Ma)", ylab="mean no. taxa per collection")

######################################################################################################################################################
##### plot MPWDist among collections and area of convex hull bounding collections
######################################################################################################################################################
geog <- getMPWDColsOneRep(this.occs[[1]], this.grid.mat = this.grid.cells$grid.mat)
	plot(rowMeans(intervals), geog$MPWD, xlim=c(66, 0), type="o", pch=21, col="orchid4", bg="orchid1", lwd=3)
	overlayCzTimescale(include.intervals=c("epoch"), plot.NALMA=TRUE, plot.NALMA_sub=TRUE)
	
	plot(rowMeans(intervals), log10(geog$area.hull/1e06), xlim=c(66, 0), type="o", pch=21, col="orchid4", bg="orchid1", lwd=3)
	overlayCzTimescale(include.intervals=c("epoch"), plot.NALMA=TRUE, plot.NALMA_sub=TRUE)
	
######################################################################################################################################################
##### plot p.mean from Foote (2016)
######################################################################################################################################################

# require (parallel)
# system.time(rez.this.rep <- mclapply(this.rep, fitTruncatedLogNormalOneIntv, fix.pmax=TRUE, mc.cores=detectCores()-2))
# which(sapply(rez.this.rep, class)== "try-error")

# system.time(rez.this.rep <- lapply(this.rep, fitTruncatedLogNormalOneIntv, fix.pmax=TRUE))

# rez.this.rep <- list()
# for (i in seq_along(this.rep)) {
	# rez.this.rep[[i]] <- fitTruncatedLogNormalOneIntv(this.rep[[i]])
# }

# mean.p <- getPMeanFromResultsOneRep(rez.this.rep)
# sapply(rez.this.rep, function(x) x$status)

	plot(rowMeans(intervals), mean.p, type="n", xlim=rev(range(intervals)), xlab="time (Ma)", ylab="mean interval occupancy")
	# overlayCzTimescale(include.intervals=c("epoch"), plot.NALMA = TRUE, plot.NALMA_sub = TRUE)
	
	lines(rowMeans(intervals), mean.p, type="o", pch=21, col="orchid4", bg="orchid1")
	text(rowMeans(intervals), mean.p, pos=4, cex=0.5)

system.time(rez <- lapply(this.m[1:25], FUN=function(x) lapply(X=x, FUN=fitTruncatedLogNormalOneIntv, fix.pmax=TRUE)))
# system.time(rez <- mclapply(this.m[1:10], FUN=function(x) lapply(X=x, FUN=fitTruncatedLogNormalOneIntv, fix.pmax=TRUE), mc.cores=detectCores()-2))
# system.time(rez <- lapply(this.m[1:10], FUN=function(x) mclapply(X=x, FUN=fitTruncatedLogNormalOneIntv, fix.pmax=TRUE, mc.cores=detectCores()-2)))
save(rez, file="~/Desktop/occupancy/occupancyResults.Rdata")

hist(sapply(rez, function(x) sapply(x, function(y) y$solution[3])))

################################################################################

# rez.this.rep <- rez[[1]]
p.mean <-  sapply(X=rez, FUN=getPMeanFromResultsOneRep)

p.quart <- t(apply(p.mean, MARGIN=1, FUN=quantile, probs=c(0.025, 0.25, 0.50, 0.75, 0.975), na.rm=TRUE))
p.quart[is.na(p.quart)] <- 0
# p.quart <- exp(p.quart)
# p.quart_Euungulata

	plot(rowMeans(intervals), p.quart[,5], type="n", pch=21, col="orchid4", bg="orchid1", xlim=rev(range(intervals)), ylim=range(p.quart), xlab="time (Ma)", ylab="mean interval occupancy")
	overlayCzTimescale(include.intervals=c("epoch"), plot.NALMA = TRUE, plot.NALMA_sub = TRUE)
	polygon(x=c(rowMeans(intervals), rev(rowMeans(intervals))), y=c(p.quart[,1], rev(p.quart[,5])), col=adjustcolor("orchid4", alpha.f = 0.25), border=NA)
	lines(rowMeans(intervals), p.quart[,3], type="o", pch=21, col="orchid4", bg="orchid1")


plot(log(sapply(this.m[[1]], function(x) length(x$occupancy))), p.quart[,3], xlab="log(no. taxa)", ylab="log Pr(occupancy)")	#x$ncols.tot , method="spearman"
abline(lm(p.quart[,3] ~ log(sapply(this.m[[1]], function(x) length(x$occupancy)))))
text(rowMeans(intervals), p.quart[,3], pos=4, cex=0.5)

cbind(ntaxa= rowMeans(sapply(this.m, function(y) sapply(y, function(x) length(x$occupancy)))), ncols= rowMeans(sapply(this.m, function(y) sapply(y, function(x) x$ncols.tot))))
matrix(rowMeans(sapply(this.occs, function(y) sapply(y, function(x) length(unique(occs$collection_no[occs$occurrence_no %in% x]))))), ncol=1)