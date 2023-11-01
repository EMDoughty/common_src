# topes<-read.table("~/Dropbox/code/java/data/isotopes.txt")
# colnames(topes)<-c("ID", "date", "d18O", "d13C")
# topes[topes==-999]<-NA

# roundList<-unique(round(topes$d18O[is.finite(topes$d18O)], digits=1))
# fillList<-rainbow(n=length(roundList), start=0.6, end=1, alpha=0.5)
# borderList<-rainbow(n=length(roundList), start=0.6, end=1, alpha=1)

# quartz(width=7.5, height=3.5)
# par(mar=c(1,1,1,1))
# plot(topes$date, topes$d18O, ylim=c(max(topes$d18O, na.rm=TRUE),min(topes$d18O, na.rm=TRUE)), pch=21, col=borderList[match(round(topes$d18O, digits=1), roundList)], bg=fillList[match(round(topes$d18O, digits=1), roundList)])
# # plot(topes$date, topes$d18O, ylim=c(max(topes$d18O, na.rm=TRUE),min(topes$d18O, na.rm=TRUE)), pch=21, col="dodgerblue4", bg="dodgerblue1")
# lo<-loess(topes$d18O ~ topes$date, span=0.075)
# lo.predict<- predict(lo, data.frame(x=topes$date))
# lines(topes$date, lo.predict, lwd=3, col="cyan")


overlayIsotopes<-function(fill.alpha=0.25, border.alpha=1.0, plot.axis=TRUE, plot.loess=TRUE, loess.span=0.075) {
	holdPar<-par()
	topes<-read.table("~/Dropbox/Code/java/data/isotopes.txt") #not in any files currently 11/1/2023
	colnames(topes)<-c("ID", "date", "d18O", "d13C")
	topes[topes==-999]<-NA
	
	roundList<-sort(unique(round(topes$d18O[is.finite(topes$d18O)], digits=1)), decreasing=TRUE)
	# breaks<-pretty(topes$d18O, n=10)
	fillList<-rainbow(n=length(roundList), start=0.6, end=1, alpha=fill.alpha)
	borderList<-rainbow(n=length(roundList), start=0.6, end=1, alpha=border.alpha)
	
	par(new=TRUE)
	plot(topes$date, topes$d18O, usr=par()$usr, xlim=c((par()$usr[1]+((0.04/1.04)*par()$usr[2]))/(1-(0.04*(0.04/1.04))), (par()$usr[2]+((0.04/1.04)*par()$usr[1]))/(1-(0.04*(0.04/1.04)))), ylim=c(max(topes$d18O, na.rm=TRUE),min(topes$d18O, na.rm=TRUE)), pch=21, col=borderList[match(round(topes$d18O, digits=1), roundList)], bg=fillList[match(round(topes$d18O, digits=1), roundList)], xaxt="n", yaxt="n", xlab="", ylab="", main="")
	if (plot.axis) {
		Axis(side=4)
		mtext(text=expression(delta^18*O), side=4, line=2, cex=par()$cex)
	}
	if (plot.loess) {
		lo<-loess(topes$d18O ~ topes$date, span=loess.span)
		lo.predict<- predict(lo, data.frame(x=topes$date))
		lines(topes$date, lo.predict, lwd=3, col="cyan")
	}
	par(new=FALSE)
}