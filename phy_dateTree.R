try(source("~/Dropbox/Code/R/common_src/strat.R"), silent=TRUE)
try(source("https://dl.dropbox.com/s/8jy9de5owxj72p7/strat.R"), silent=TRUE)

tipToNodeBranchID <- function(this.tip, this.tree) {	
	tipToNodeBranchID_machine(this.desc=this.tip, edge = this.tree$edge, root.node = length(this.tree$tip.label) + 1, branch.vec=vector())
}

tipToNodeBranchID_machine <- function (this.desc, edge, root.node, branch.vec) {
	branch.vec <- c(branch.vec, which(edge[,2]==this.desc))
	if (edge[edge[,2]==this.desc, 1]!=root.node) branch.vec <- tipToNodeBranchID_machine(this.desc= edge[edge[,2]==this.desc, 1], edge, root.node, branch.vec)
	branch.vec
}

getMaxShortBranchesRootToTip <- function(this.tree, min.bl=0.0) {
	branch.list <- lapply(X=seq_along(this.tree$tip.label), FUN=tipToNodeBranchID, this.tree=this.tree)
	max(sapply(branch.list, function(x) sum(this.tree$edge.length[x]<=min.bl)))
}

makeBLFromNodeDates <- function(this.tree) {
	this.tree$edge.length <- this.tree$node.date[this.tree$edge[ , 1]] - this.tree$node.date[this.tree$edge[ , 2]]
	# this.tree$edge.length[this.tree$edge.length < min.bl] <- min.bl									min.bl shouldn't be added here, but in the node dates themselves...
	this.tree
}

# lengthToRootMachine <- function(edge, this.tree, thisLength) {
	# thisLength <- thisLength + this.tree$edge.length[edge]
	# if (this.tree$edge[edge,1]!=length(this.tree$tip.label)+1) {
		# thisLength <- lengthToRootMachine(which(this.tree$edge[ ,2]==this.tree$edge[edge,1]), this.tree, thisLength)
	# }
	# return(thisLength)
# }

lengthToRoot <- function(this.node, this.tree) {
	ltr <- 0
	while (this.node != (length(this.tree$tip.label)+1)) {
		ltr = ltr + this.tree$edge.length[this.tree$edge[ ,2]==this.node]
		this.node <- this.tree$edge[this.tree$edge[ ,2]==this.node, 1]
	}
	ltr
}

makeNodeDatesFromBL <- function(this.tree, root.age=NULL, youngestTipAge=0.0) {
	n.tips <- length(this.tree$tip.label)
	this.tree$node.date <- array(NA, dim=n.tips+this.tree$Nnode)
	ltr <- sapply(seq_along(this.tree$node.date), lengthToRoot, this.tree)
	if (is.null(root.age)) {
		root.age <- youngestTipAge + max(ltr)
	}
	this.tree$node.date <- root.age - ltr
	this.tree$tsr <- max(this.tree$node.date) - this.tree$node.date
	this.tree
}


setNodeDatesMachine <- function(this.tree, nodeID, min.bl=5e-5) {
	desc <- this.tree$edge[this.tree$edge[,1] == nodeID, 2]				#get all desc of this node
	noDates <- desc[which(!is.finite(this.tree$node.date[desc]))] 	#find desc without dates
	if (length(noDates) > 0) {
		for (i in noDates) {
			if (i > length(this.tree$tip.label)) this.tree <- setNodeDatesMachine(this.tree=this.tree, nodeID=i, min.bl) # recursively check those nodes, then get back to me
		}
	}

	if (any(is.finite(this.tree$node.date[desc]))) {
		this.tree$node.date[which(!is.finite(this.tree$node.date[desc]))] <- max(this.tree$node.date[desc], na.rm=TRUE)	# if any descendants are STILL NA, then assign them the max desc age -> STILL ALLOWS FOR NA DATES IF SISTERS ALL DESC OF A DESC ARE MISSING
		this.tree$node.date[nodeID] <- max(this.tree$node.date[desc], na.rm=TRUE) + min.bl
	}
	this.tree$tsr <- max(this.tree$node.date) - this.tree$node.date
	this.tree
}

updateNodeDatesUsingTipNodeDates<-function(this.tree, min.bl=5e-5) {
	#requires that this.tree$node.date for tips are correct and will not be changed
	this.tree$node.date <- c(this.tree$node.date[this.tree$tip.label], array(NA, dim=this.tree$Nnode))
	setNodeDatesMachine(this.tree, nodeID=length(this.tree$tip.label)+1, min.bl=min.bl)
}

dropTipKeepDates <- function(this.tree, tip=NULL, min.bl=0.5) {
	if (length(tip) < 1) return(this.tree)
	this.tree <- ape::drop.tip(phy=this.tree, tip=tip, rooted=FALSE)
	this.tree$node.date <- this.tree$node.date[this.tree$tip.label]
	makeNodeDatesFromBL(this.tree, youngestTipAge=min(this.tree$node.date))
}

# scaleNodeDates<-function(cladeDuration=1, this.tree) {
	# this.tree$node.date<-this.tree$node.date*(cladeDuration/max(sapply(1:length(this.tree$tip.label), lengthToRoot, this.tree)))
	# return(this.tree)
# }

# scaleBL<-function(cladeDuration=1, this.tree) {
	# this.tree$edge.length<-this.tree$edge.length*(cladeDuration/max(sapply(1:length(this.tree$tip.label), lengthToRoot, this.tree)))
	# return(this.tree)
# }

### function ensures that anc nodes are at least [min.bl] older than their oldest descendant
applyMinBL <- function(this.tree, terminals.only=FALSE, min.bl=0.0) {
	for (thisNode in seq(from=length(this.tree$node.date), to=length(this.tree$tip.label)+1)) { # going from youngest internal node to the root node
		this.desc <- this.tree$edge[this.tree$edge[,1] == thisNode, 2]
		if (terminals.only) {
			if (any(this.desc <= length(this.tree$tip.label))) {
				this.max <- max(this.tree$node.date[this.desc[this.desc <= length(this.tree$tip.label)]], na.rm=TRUE)
			} else this.max = 0
		} else this.max <- max(this.tree$node.date[this.desc], na.rm=TRUE)
		if (this.max > (this.tree$node.date[thisNode] - min.bl)) {
			this.tree$node.date[thisNode] <- this.max + min.bl
		}
	}
	this.tree <- makeBLFromNodeDates(this.tree)
	if (any(this.tree$edge.length < 0)) this.tree <- applyMinBL(this.tree, min.bl=0.0)
	this.tree$tsr <- max(this.tree$node.date) - this.tree$node.date
	this.tree
}

smoothNodeAges <- function(this.tree, method=c("uniform", "exponential"), root.max=max(this.tree$node.date), min.bl=0.0) {
	# randomly pushes internal nodes back in time, beginning with root and working toward tips 
	method <- match.arg(method)
	ntaxa <- length(this.tree$tip.label)
	if (root.max > max(this.tree$node.date)) {
		thisNode <- ntaxa + 1	# this is ID of the root node
		if (method=="uniform") this.tree$node.date[thisNode] <- runif(n=1, min=this.tree$node.date[thisNode], max=root.max) #adjust root node age
		# if (method=="uniform") this.tree$node.date[thisNode] <- runif(n=1, min=this.tree$node.date[thisNode], max=(this.tree$node.date[thisNode] + rootAdj)) #adjust root node age
		# if (method=="exponential") this.tree$node.date[thisNode]<-rexp(1), min=this.tree$node.date[thisNode], max=(this.tree$node.date[thisNode]+rootAdj)) #adjust root node age
	} else if (root.max < max(this.tree$node.date)) {
		warning(paste0("In smoothNodeAges: oldest node is already older (", max(this.tree$node.date), ") than root.max (",root.max,") ***"))
	}
	for (thisNode in seq(from=(ntaxa + 2), to=length(this.tree$node.date))) { # root node (ntaxa+1) is already fixed, other nodes can be smoothed
		ancNode <- this.tree$edge[this.tree$edge[,2]==thisNode, 1]
		if ((this.tree$node.date[ancNode] - this.tree$node.date[thisNode]) > min.bl) this.tree$node.date[thisNode] <- runif(1, min=this.tree$node.date[thisNode], max=(this.tree$node.date[ancNode] - min.bl))
	}
	this.tree$tsr <- max(this.tree$node.date) - this.tree$node.date
	makeBLFromNodeDates(this.tree)
}

backfillMachine <- function(this.tree, nodeID, min.bl=0.0, min.date=0.0) {
	if (is.finite(this.tree$node.date[nodeID])) return(this.tree) 
	ancNode <- this.tree$edge[this.tree$edge[,2]==nodeID,1]
	if (!is.finite(this.tree$node.date[ancNode])) {
		if (ancNode == length(this.tree$tip.label) + 1) { print("backfillMachine error:	root node date is missing")
		} else this.tree <- backfillMachine(this.tree, ancNode, min.bl, min.date)
	} 
	
	if ((this.tree$node.date[ancNode] - min.bl) < min.date) this.tree$node.date[nodeID] <- min.date else this.tree$node.date[nodeID] <- (this.tree$node.date[ancNode] - min.bl)
	this.tree
}

backfillMissingNodeDates <- function(this.tree, min.bl=0.0, min.date=0.0) {
	for (i in sort(which(!is.finite(this.tree$node.date)), decreasing=FALSE)) this.tree <- backfillMachine(this.tree=this.tree, nodeID=i, min.bl=min.bl, min.date=min.date)
	this.tree$tsr <- max(this.tree$node.date) - this.tree$node.date
	this.tree
}

dateTreeWithRanges <- function(this.tree, strat.ranges, within.error=FALSE, random.ages=TRUE, min.date=0.0, min.bl=0.0, smooth.node.ages=FALSE, root.max=max(strat.ranges), keep.missing=FALSE, missing.date=NULL, include.ranges=FALSE) {
	if (keep.missing) {
		strat.ranges <- strat.ranges[match(this.tree$tip.label, rownames(strat.ranges)),]		#reorder taxa in strat.ranges by the order in the this.tree
		rownames(strat.ranges) <- this.tree$tip.label	# the match above erases unmatched species from the rownames - this restablishes them
		if (!is.null(missing.date)) strat.ranges[!is.finite(strat.ranges)] <- missing.date
	} else if (any(!(this.tree$tip.label %in% rownames(strat.ranges)))) {
		this.tree <- drop.tip(this.tree, this.tree$tip.label[!(this.tree$tip.label %in% rownames(strat.ranges))]) # cut taxa that are not in ranges
		# strat.ranges <- strat.ranges[match(this.tree$tip.label, rownames(strat.ranges)),]		#reorder taxa in strat.ranges by the order in the this.tree
	}

	strat.ranges <- strat.ranges[rownames(strat.ranges) %in% this.tree$tip.label,] 		#cull taxa in strat.ranges list that aren't in the this.tree
	strat.ranges <- strat.ranges[match(this.tree$tip.label, rownames(strat.ranges)),]		#reorder taxa in strat.ranges by the order in the this.tree
	
	# set node dates of tips
	if (within.error) {
		if (random.ages) { this.tree$node.date <- getRandomRanges(strat.ranges)[,"FO"] 
		} else this.tree$node.date <- rowMeans(strat.ranges[,1:2]) #if without error, use the midpoint of FA dates
	} else this.tree$node.date <- strat.ranges[,"FO"]	# if error on the FO and LO is not provided, a column of FOs should be
		
	# make space for internal nodes
	this.tree$node.date <- c(this.tree$node.date, array(data=NA, dim=this.tree$Nnode)) #add NA dates for internal nodes
	this.tree <- setNodeDatesMachine(this.tree=this.tree, nodeID=length(this.tree$tip.label)+1, min.bl=0)
	this.tree$tsr <- max(this.tree$node.date) - this.tree$node.date
	this.tree <- makeBLFromNodeDates(this.tree)

	# apply min.bl to terminal branches ONLY
 	this.tree <- applyMinBL(this.tree, terminals.only=TRUE, min.bl=min.bl) # apply min.bl to terminal branches ONLY

 	# this determines whether to reduce min.bl to stay under max.root
 	this.min <- min.bl
 	if (any(this.tree$edge.length <= min.bl)) {
		this.min <- (root.max - max(this.tree$node.date)) / getMaxShortBranchesRootToTip(this.tree, min.bl=min.bl)
		if (this.min < min.bl) {
			message <- paste0("in dateTreeWithRanges: To keep the root younger than root.max (",root.max,"), minimum length of internal branches set to ",round(this.min,3)," (shorter length than min.bl [",round(min.bl,3),"])")
			warning(message)
		}
		this.min <- min(this.min, min.bl)
		this.tree <- applyMinBL(this.tree, min.bl=this.min)
	}
	
	# if (any(!is.finite(this.tree$node.date))) this.tree <- backfillMissingNodeDates(this.tree, min.date = min.bl)
	if (any(!is.finite(this.tree$node.date))) {
		this.tree <- backfillMissingNodeDates(this.tree, min.date = this.min)
		this.tree$tsr <- max(this.tree$node.date) - this.tree$node.date
	}

	# if (smooth.node.ages) this.tree <- smoothNodeAges(this.tree, method="uniform", root.max=root.max, min.bl=min.bl)
	if (smooth.node.ages) {
		this.tree <- smoothNodeAges(this.tree, method="uniform", root.max=root.max, min.bl=this.min)
		this.tree$tsr <- max(this.tree$node.date) - this.tree$node.date
	}

	this.tree <- makeBLFromNodeDates(this.tree) # called within smoothNodeAges, but has to be here in case smooth.node.ages is FALSE
	if (include.ranges) {
		if (within.error) { this.tree$strat.range <- cbind(this.tree$node.date[seq_along(this.tree$tip.label)], apply(strat.ranges[,c("ageLO_max", "ageLO_min")], 1, function(x) runif(1, min=x[2], max=x[1]) ))
		} else this.tree$strat.range <- cbind(this.tree$node.date[seq_along(this.tree$tip.label)], strat.ranges[,"LO"])
		this.tree$strat.range[this.tree$strat.range[,1]<this.tree$strat.range[,2],] <- this.tree$strat.range[this.tree$strat.range[,1]<this.tree$strat.range[,2],c(1,1)]
	}
	this.tree
}

getProportionsOfBranchesWithinIntervals<-function(this.tree, intervals) {
	oneBranchOneInterval<-function(thisInterval, thisBranch) {
		if (thisBranch[1] > thisInterval["ageBase"] & thisBranch[2] < thisInterval["ageTop"]) {														return ((thisInterval["ageBase"] - thisInterval["ageTop"]) / (thisBranch[1]-thisBranch[2])) # bt
		} else if (thisBranch[1] <= thisInterval["ageBase"] & thisBranch[2] >= thisInterval["ageTop"]) {											return (1.0)  #FL
		} else if (thisBranch[1] <= thisInterval["ageBase"] & thisBranch[1] > thisInterval["ageTop"] & thisBranch[2] < thisInterval["ageTop"]) {	return ((thisBranch[1] - thisInterval["ageTop"]) / (thisBranch[1] - thisBranch[2]))	#Ft
		} else if (thisBranch[1] > thisInterval["ageBase"] & thisBranch[2] < thisInterval["ageBase"] & thisBranch[2] >= thisInterval["ageTop"]) {	return ((thisInterval["ageBase"]-thisBranch[2]) / (thisBranch[1] - thisBranch[2]))	#bL
		} else return(0.0)
	}
	
	# ancAge=this.tree$node.date[this.tree$edge[i,1]]
	# descAge=this.tree$node.date[this.tree$edge[i,2]]
	intvList<-listifyMatrixByRow(intervals)
	t(apply(this.tree$edge, 1, function(x) sapply(intvList, oneBranchOneInterval, c(this.tree$node.date[x[1]], this.tree$node.date[x[2]]))))
}

getLMAForSingleInterval<-function(thisInt, thisBranch) {
	if (thisBranch["ageDesc"]>=thisInt["ageBase"] | thisBranch["ageAnc"]<=thisInt[ "ageTop"]) { return(0)
	} else if (thisBranch["ageAnc"]>=thisInt["ageBase"] & thisBranch["ageDesc"]<=thisInt["ageTop"]) { return(thisInt["ageBase"] - thisInt["ageTop"])		#bt
	} else if (thisBranch["ageAnc"]<=thisInt["ageBase"] & thisBranch["ageDesc"]>=thisInt["ageTop"]) { return(thisBranch["ageAnc"]-thisBranch["ageDesc"] )	#FL
	} else if (thisBranch["ageAnc"]<=thisInt["ageBase"] & thisBranch["ageDesc"]<=thisInt["ageTop"]) { return(thisBranch["ageAnc"] - thisInt["ageTop"])		#Ft
	} else if (thisBranch["ageAnc"]>=thisInt["ageBase"] & thisBranch["ageDesc"]>=thisInt["ageTop"]) { return(thisInt["ageBase"] - thisBranch["ageDesc"]) }	#bL
}

getLMAForSingleBranch<-function(thisBranch, intervals) {
	apply(intervals, 1, getLMAForSingleInterval, thisBranch)
}

getBranchAges<-function(thisBranch, this.tree, dates, use.ranges=FALSE) { 
	if (use.ranges) {
		if (thisBranch[2]<=length(this.tree$tip.label)) return(c(ageAnc=this.tree$node.date[thisBranch[1]], ageDesc=unname(this.tree$strat.range[thisBranch[2],2])) )
	} #else 
	c(ageAnc=this.tree$node.date[thisBranch[1]], ageDesc=this.tree$node.date[thisBranch[2]])
}

