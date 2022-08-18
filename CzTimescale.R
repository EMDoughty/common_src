
overlayCzTimescale <- function(do.subepochs=FALSE, color=TRUE, thisAlpha.intervals=0.33, thisAlpha.text = 0.33, borderCol="white", invertTime=FALSE, scale.cex=0.75, scale.headers = 0.95, text.offset = 0.025) {
	textCol <- rgb(0,0,0,thisAlpha.text)
	# textShadowCol<-"gray50"
	old.cex<-par()$cex
	par(cex=old.cex * scale.cex)
	epochs=data.frame(name=c("Camb",
							"eO",
							"mO",
							"lO",
							"eS",
							"mS",
							"lS",
							"eD",
							"mD",
							"lD",
							"Miss",
							"Penn",
							"eP",
							"mP",
							"lP",
							"eT", 
							"mT",
							"lT",
							"eJ",
							"mJ",
							"lJ",
							"eK",
							"lK",
							"Paleo",
							"Eo",
							"Oligo",
							"Mio",
							"Plio",
							"Pl.",
							"Recent"),
					ageBase =c(542,
								486,
								472,
								461,
								444,
								428,
								423,
								416,
								398,
								385,
								359,
								318,
								299,
								271,
								260,

								251.0,
								245.0,
								235.0,
								201.6,
								176.0,
								161.0,
								145.5,
								99.6,
								65.5,
								55.8,
								33.9,
								23.03,
								5.33,
								2.58,
								0))
	if (color) { epochs<-data.frame(epochs, rgb =c(
		rgb(0.533, 0.671, 0.494, thisAlpha.intervals),
		rgb(0.0, 0.686, 0.565, thisAlpha.intervals),
		rgb(0.118, 0.737, 0.624, thisAlpha.intervals),
		rgb(0.459, 0.796, 0.690, thisAlpha.intervals),
		rgb(0.565, 0.835, 0.788, thisAlpha.intervals),
		rgb(0.675, 0.871, 0.831, thisAlpha.intervals),
		rgb(0.725, 0.894, 0.867, thisAlpha.intervals),
		rgb(0.906, 0.698, 0.471, thisAlpha.intervals),
		rgb(0.953, 0.788, 0.557, thisAlpha.intervals),
		rgb(0.953, 0.875, 0.710, thisAlpha.intervals),
		rgb(0.427, 0.624, 0.533, thisAlpha.intervals),
		rgb(0.584, 0.769, 0.780, thisAlpha.intervals),
		rgb(0.937, 0.463, 0.404, thisAlpha.intervals),
		rgb(0.988, 0.549, 0.478, thisAlpha.intervals),
		rgb(0.996, 0.702, 0.647, thisAlpha.intervals),

		rgb(0.643, 0.365, 0.627, thisAlpha.intervals),
		rgb(0.718, 0.510, 0.71, thisAlpha.intervals),
		rgb(0.745, 0.616, 0.776, thisAlpha.intervals),
		rgb(0.0, 0.718, 0.906, thisAlpha.intervals),
		rgb(0.392, 0.816, 0.918, thisAlpha.intervals),
		rgb(0.647, 0.882, 0.973, thisAlpha.intervals),
		rgb(0.5803922, 0.7960784, 0.4745098, thisAlpha.intervals),
		rgb(0.7803922, 0.8784314, 0.6156863, thisAlpha.intervals),
		rgb(0.9803922, 0.6980392, 0.4862745, thisAlpha.intervals),
		rgb(0.9843137, 0.7372549, 0.5294118, thisAlpha.intervals),
		rgb(0.9960784, 0.8588235, 0.6745098, thisAlpha.intervals),
		rgb(1, 0.945098, 0, thisAlpha.intervals),
		rgb(1, 0.9764706, 0.6823529, thisAlpha.intervals),
		rgb(1, 0.9411765, 0.7490196, thisAlpha.intervals),
		rgb(1, 0.9529412, 0.9333333, thisAlpha.intervals)), stringsAsFactors = FALSE) 
	} else { epochs<-data.frame(epochs, rgb =c(
		rgb(0.57, 0.57, 0.57, thisAlpha.intervals),
		rgb(0.51, 0.51, 0.51, thisAlpha.intervals),
		rgb(0.59, 0.59, 0.59, thisAlpha.intervals),
		rgb(0.68, 0.68, 0.68, thisAlpha.intervals),
		rgb(0.79, 0.79, 0.79, thisAlpha.intervals),
		rgb(0.79, 0.79, 0.79, thisAlpha.intervals),
		rgb(0.79, 0.79, 0.79, thisAlpha.intervals),
		rgb(0.69, 0.69, 0.69, thisAlpha.intervals),
		rgb(0.77, 0.77, 0.77, thisAlpha.intervals),
		rgb(0.85, 0.85, 0.85, thisAlpha.intervals),
		rgb(0.51, 0.51, 0.51, thisAlpha.intervals),
		rgb(0.68, 0.68, 0.68, thisAlpha.intervals),
		rgb(0.53, 0.53, 0.53, thisAlpha.intervals),
		rgb(0.61, 0.61, 0.61, thisAlpha.intervals),
		rgb(0.73, 0.73, 0.73, thisAlpha.intervals),

		rgb(0.39, 0.39, 0.39, thisAlpha.intervals),
		rgb(0.51, 0.51, 0.51, thisAlpha.intervals),
		rgb(0.59, 0.59, 0.59, thisAlpha.intervals),
		rgb(0.58, 0.58, 0.58, thisAlpha.intervals),
		rgb(0.71, 0.71, 0.71, thisAlpha.intervals),
		rgb(0.81, 0.81, 0.81, thisAlpha.intervals),
		rgb(0.69, 0.69, 0.69, thisAlpha.intervals),
		rgb(0.81, 0.81, 0.81, thisAlpha.intervals),
		rgb(0.71, 0.71, 0.71, thisAlpha.intervals),
		rgb(0.75, 0.75, 0.75, thisAlpha.intervals),
		rgb(0.85, 0.85, 0.85, thisAlpha.intervals),
		rgb(0.93, 0.93, 0.93, thisAlpha.intervals),
		rgb(0.96, 0.96, 0.96, thisAlpha.intervals),
		rgb(0.93, 0.93, 0.93, thisAlpha.intervals),
		rgb(1, 1, 1, thisAlpha.intervals)), stringsAsFactors = FALSE) }
	# } else { epochs <- data.frame(epochs, rgb = rep("000000", times=nrow(epochs)), stringsAsFactors = FALSE) }
	
	if (do.subepochs) subepochs=data.frame(modifier=c("e", "m", "l", "e", "m", "l", "e", "l", "e", "m", "l", "PP"), ageBase=c(65.5, 61.7, 58.7, 55.8, 48.6, 37.2, 33.9, 28.4, 23.03, 15.97, 11.61, 5.33), stringsAsFactors = FALSE)

	if (invertTime) {
		epochs$ageBase<- par()$usr[2]/1.137059-epochs$ageBase
		if (do.subepochs) subepochs$ageBase<- par()$usr[2]/1.137059-subepochs$ageBase
	}

	top=par()$usr[4]
	bottom=par()$usr[3]

#lays down colors
	for (i in 1:(nrow(epochs)-1)) {
		polygon(c(epochs[i,2], epochs[i,2], epochs[(i+1),2], epochs[(i+1),2]), c(bottom, top, top, bottom), col=epochs$rgb[i], border=borderCol, lty=0)
	}

#lays down Era boundaries
	abline(v=251,lwd=1.5, col=borderCol)
	abline(v=65.5,lwd=1.5, col=borderCol)

	if (do.subepochs) {
		top=(scale.headers*(top - bottom)) + bottom
		for (i in 1:(nrow(subepochs)-1)) {
			polygon(cbind(c(subepochs[i,2], subepochs[i,2], subepochs[(i+1),2], subepochs[(i+1),2]), c(bottom, top, top, bottom)), lwd=0.25, border=borderCol)
			text (x=mean(c(as.numeric(subepochs[i,2]), as.numeric(subepochs[(i+1),2]))), y=(((scale.headers + text.offset)*(top - bottom)) + bottom), label=subepochs[i,1], col=textCol)
		}
	}

#lays down borders and text
	top=par()$usr[4]
	for (i in 1:(nrow(epochs)-1)) {
		polygon(cbind(c(epochs[i,2], epochs[i,2], epochs[(i+1),2], epochs[(i+1),2]), c(bottom, top, top, bottom)), col=NA, border=borderCol, lwd=0.5)
		# text (mean(c(as.numeric(epochs[i,2]), as.numeric(epochs[i+1,2]))), (0.975*top), label=epochs[i,1], col=textShadowCol, adj=c(0.475,0.525), font=2) #cex=par()$cex+0.05
		# text (mean(c(as.numeric(epochs[i,2]), as.numeric(epochs[i+1,2]))), (0.975*top), label=epochs[i,1], col=textShadowCol, font=2) #cex=par()$cex+0.05
		text (x=mean(c(as.numeric(epochs[i,2]), as.numeric(epochs[i+1,2]))), y=(((scale.headers + text.offset)*(top - bottom)) + bottom), label=epochs[i,1], col=textCol, font=2)
	}
	abline(h=((scale.headers*(top - bottom)) + bottom), lwd=1.0, col=borderCol)

	par(cex=old.cex)
}

