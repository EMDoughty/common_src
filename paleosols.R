source("~/Dropbox/Code/R/common_src/strat.R")
source("~/Dropbox/Code/R/common_src/utils_marcot.R")
source("~/Dropbox/Code/R/common_src/CzTimescale.R")

intervals <- makeIntervals(1, 54.5, 1)

paleosol_OR <- read.csv("~/Dropbox/Code/common_dat/paleosol_OR.csv") #these common_dat files are not in any of my folders 11/1/2023
paleosol_MT <- read.csv("~/Dropbox/Code/common_dat/paleosol_MT.csv")
paleosol_NE <- read.csv("~/Dropbox/Code/common_dat/paleosol_NE.csv")
# paleosol <- rbind(paleosol, read.csv("~/Dropbox/code/common_dat/paleosol_MT.csv"))
# paleosol <- rbind(paleosol, read.csv("~/Dropbox/code/common_dat/paleosol_NE.csv"))

reps <- 1000
# randMAP <- sapply(seq_len(nrow(intervals)), function(x) sapply(randMAP, function(y, x) y[[x]], x=x))
# randMAP <- replicate(n=reps, expr=apply(intervals, 1, function(x) apply(paleosol[paleosol$Age..Ma. >= x[1] & paleosol$Age..Ma. < x[2], c("error.", "error..1")], 1, function(y) runif(n=1, min=y[2], max=y[1]))), simplify=FALSE)
randMAP <- replicate(n=reps, expr=apply(intervals, 1, function(x) quantile(apply(paleosol[paleosol$Age..Ma. >= x[1] & paleosol$Age..Ma. < x[2], c("error.", "error..1")], 1, function(y) runif(n=1, min=y[2], max=y[1])), probs=c(0.0, 0.025, 0.25, 0.50, 0.75, 0.975, 1.0), na.rm=TRUE)), simplify="array")

randMAP_OR <- replicate(n=reps, expr=apply(intervals, 1, function(x) quantile(apply(paleosol_OR[paleosol_OR$Age..Ma. >= x[1] & paleosol_OR$Age..Ma. < x[2], c("error.", "error..1")], 1, function(y) runif(n=1, min=y[2], max=y[1])), probs=c(0.0, 0.025, 0.25, 0.50, 0.75, 0.975, 1.0), na.rm=TRUE)), simplify="array")
randMAP_MT <- replicate(n=reps, expr=apply(intervals, 1, function(x) quantile(apply(paleosol_MT[paleosol_MT$Age..Ma. >= x[1] & paleosol_MT$Age..Ma. < x[2], c("error.", "error..1")], 1, function(y) runif(n=1, min=y[2], max=y[1])), probs=c(0.0, 0.025, 0.25, 0.50, 0.75, 0.975, 1.0), na.rm=TRUE)), simplify="array")
randMAP_NE <- replicate(n=reps, expr=apply(intervals, 1, function(x) quantile(apply(paleosol_NE[paleosol_NE$Age..Ma. >= x[1] & paleosol_NE$Age..Ma. < x[2], c("error.", "error..1")], 1, function(y) runif(n=1, min=y[2], max=y[1])), probs=c(0.0, 0.025, 0.25, 0.50, 0.75, 0.975, 1.0), na.rm=TRUE)), simplify="array")

# MAP <- sapply(randMAP[1:21], quantile, probs=c(0.0, 0.025, 0.25, 0.50, 0.75, 0.975, 1.0), na.rm=TRUE)
# MAP <- apply(randMAP[1:21, ], 1, quantile, probs=c(0.0, 0.025, 0.25, 0.50, 0.75, 0.975, 1.0), na.rm=TRUE)
MAP <- apply(randMAP, c(1,2), median, na.rm=TRUE)
MAP_OR <- apply(randMAP_OR, c(1,2), median, na.rm=TRUE)
MAP_MT <- apply(randMAP_MT, c(1,2), median, na.rm=TRUE)
MAP_NE <- apply(randMAP_NE, c(1,2), median, na.rm=TRUE)


index <- apply(MAP, 2, function(x) all(is.finite(x)))
# index <- seq_len(21)
plot(rowMeans(intervals[index,]), MAP[7,index], xlim=c(65,0), ylim=c(0, max(MAP, na.rm=TRUE)), type="n", xlab="Time (Ma)", ylab="MAP")
overlayCzTimescale(do.subepochs=TRUE)
polygon(c(rowMeans(intervals[index,]), rev(rowMeans(intervals[index,]))), c(MAP[1,index], rev(MAP[7,index])), border=alphaColor("burlywood4", 0.5), col=alphaColor("burlywood4", 0.2), lwd=0.5)
polygon(c(rowMeans(intervals[index,]), rev(rowMeans(intervals[index,]))), c(MAP[2,index], rev(MAP[6,index])), border=alphaColor("burlywood4", 0.5), col=alphaColor("burlywood4", 0.2), lwd=0.5)
polygon(c(rowMeans(intervals[index,]), rev(rowMeans(intervals[index,]))), c(MAP[3,index], rev(MAP[5,index])), border=alphaColor("burlywood4", 0.5), col=alphaColor("burlywood4", 0.2), lwd=0.5)
lines(rowMeans(intervals[index,]), MAP[4,index], type="o", lwd=2, col="burlywood4", pch=21, bg="burlywood1")
lines(rowMeans(intervals[index,]), MAP_OR[4,index], type="o", lwd=2, col="darkolivegreen4", pch=21, bg="darkolivegreen1")
lines(rowMeans(intervals[index,]), MAP_MT[4,index], type="o", lwd=2, col="firebrick4", pch=21, bg="firebrick1")
lines(rowMeans(intervals[index,]), MAP_NE[4,index], type="o", lwd=2, col="goldenrod4", pch=21, bg="goldenrod1")
legend("bottomleft", legend=c("Oregon", "Montana", "Nebraska"), bg="white", lty=1, pch=21, col=c("darkolivegreen4", "firebrick4", "goldenrod4"), pt.bg=c("darkolivegreen1", "firebrick1", "goldenrod1"))

