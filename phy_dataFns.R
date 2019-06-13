require(compiler)
require(stringr)

convertOneChrStateVecToUL <- function(state.vec) {
	if (grepl(pattern="/", x= state.vec) | grepl(pattern="&", x= state.vec)) {
		require(stringr)
		state.vec <- unlist(str_split(state.vec, pattern="[[:punct:]]"))
	}
	stateListToUL(as.numeric(state.vec))
}

convertChrVecToUnsignedLong <- function(chr.vec, missing.symbols=c("?","-")) {
	chr.vec[toupper(chr.vec)=="A"]<-"10"
	chr.vec[toupper(chr.vec)=="B"]<-"11"
	chr.vec[toupper(chr.vec)=="C"]<-"12"
	chr.vec[toupper(chr.vec)=="D"]<-"13"
	chr.vec[toupper(chr.vec)=="E"]<-"14"
	chr.vec[toupper(chr.vec)=="F"]<-"15"
	chr.vec[toupper(chr.vec)=="G"]<-"16"
	chr.vec[toupper(chr.vec)=="H"]<-"17"
	chr.vec[toupper(chr.vec)=="J"]<-"18"
	chr.vec[toupper(chr.vec)=="K"]<-"19"
	chr.vec[toupper(chr.vec)=="L"]<-"20"
	
	chr.vec[chr.vec %in% missing.symbols] <- NA
	
	chr.vec <- sapply(chr.vec, convertOneChrStateVecToUL, USE.NAMES=FALSE)
	# chr.vec <- 2^(as.numeric(as.character(factor(as.vector(chr.vec),)))) ## adds 1, so state 0 becomes state 1 
	if (any(is.finite(chr.vec)) & !all(is.finite(chr.vec))) chr.vec[!is.finite(chr.vec)] <- stateListToUL(0:log2(max(chr.vec, na.rm=TRUE))) # a list of states from 0 (not zero-indexed) to the maximum observed state for the chr - replaces missing data in matrix
	chr.vec
}

convertAPEDataToUL<-function(dat) {
	thisDat <- apply(t(as.data.frame(dat, check.names=TRUE, check.rows=TRUE)), 2, convertChrVecToUnsignedLong)
	rownames(thisDat)<-names(dat)
	return(thisDat)
}

ULToStateList <- function(a) {
# returns the true states from a bitwise unsigned long integeter state set
	if (!is.finite(a) | a==0) return(NA)
	states <- vector(mode="numeric")
	while (a>0) {
		states<-c(states, trunc(log2(a)))
		a<-a-(2^trunc(log2(a)))
	}
	return(sort(states))
}
ULToStateList_vec <- Vectorize(ULToStateList)

ULToStateListChr <- function(a) {
# returns the true states from a bitwise unsigned long integeter state set - as a single chr
	if (!is.finite(a) | a==0) return(NA)
	states <- vector(mode="numeric")
	while (a > 0) {
		states <- c(states, trunc(log2(a)))
		a <- a-(2^trunc(log2(a)))
	}
	paste(sort(states), collapse="/")
}
ULToStateListChr_vec <- Vectorize(ULToStateListChr)

stateListToUL <- function(state.vec) {
	sum(2^state.vec)
}

getNStatesULVec <- function(chr) {
	length(sort(unique(unlist(ULToStateList_vec(chr)))))
}


isULMissing <- function(chr.vec) {
	missing.value = -1
	max.state <- max(unlist(sapply(chr.vec, ULToStateList_vec)))
	if (max(chr.vec) == 2^(max.state + 1) - 1) missing.value <- 2^(max.state + 1) - 1
	if (missing.value > 0) { chr.vec==missing.value
	} else rep(FALSE, length(chr.vec))
}

# recodeChr <- function(chr) {
	# states <- sort(unique(unlist(ULToStateList_vec(chr))))
	# if ((max(states)) > (length(states) - 1)) {
		
	# }
# }

# ULToStateList<-function(a) {
# # returns the true states from a bitwise unsigned long integeter state set
	# state.vec<-vector()
	# while (a>0) {
		# state.vec<-c(state.vec, trunc(log2(a)))
		# a<-a-(2^state.vec[length(state.vec)])
	# }
	# return(sort(state.vec))
# }

# convertULListToBit<-function(chr.vec) {
	# require(bit)
	# lapply(sapply(chr.vec, ULToStateList), as.bit.which, length=trunc(log2(max(chr.vec, na.rm=TRUE))))
# }

# convertBitListToUL<-function(bitList) {
	# require(bit)
	# sapply(sapply(bitList, as.which), stateListToUL)
# }

bitMinUL<-function(a) {
# returns the minimum (UL) state from a bitwise unsigned long integeter state set
	# state<-NA
	while (a>0) {
		state<-trunc(log2(a))
		a<-a-(2^state)
	}
	2^state
}

bitMaxUL<-function(a) {
# returns the maximum (UL) state from a bitwise unsigned long integeter state set
	2^(trunc(log2(a)))
}

bitMaxTrue<-function(a) {
# returns the maximum (UL) state from a bitwise unsigned long integeter state set
	trunc(log2(a))
}
cbmaxt<-cmpfun(bitMaxTrue)

bitMinTrue<-function(a) {
# returns the minimum (true) state from a bitwise unsigned long integeter state set
	# state<-NA
	while (a>0) {
		state<-trunc(log2(a))
		a<-a-(2^state)
	}
	state
}
cbmint<-cmpfun(bitMinTrue)

simulateOneCharacter<-function(tree, delta) {
	tree<-reorder(tree, order="c")
	dat<-matrix(0, nrow=nrow(tree$edge)+1, ncol=1)
	dat[tree$edge[1,1],1]<-1 #sets the root node's state as O (1 is the unsigned long representation of state O)
	rownames(dat)<-c(tree$tip.label, paste("node", (length(tree$tip.label)+1:tree$Nnode)))
	chReal<-vector(mode="numeric", length=nrow(tree$edge))
	for (i in 1:nrow(tree$edge)) {
		chReal[i]<-rpois(1, delta*tree$edge.length[i])
		if ((chReal[i]%%2)>0) {
			if (dat[tree$edge[i,1],1]==1) dat[tree$edge[i,2],1]<-2
			else if (dat[tree$edge[i,1],1]==2) dat[tree$edge[i,2],1]<-1
		} else dat[tree$edge[i,2],1]<-dat[tree$edge[i,1],1]
	}
	chObs<-sapply(1:nrow(tree$edge), function(x, tree, dat) { if (dat[tree$edge[x,][1],1]!=dat[tree$edge[x,][2],1]) return(1) else return(0) }, tree, dat)
	return(list(dat=dat, chReal=chReal, chObs=chObs))
}

simulateMatrix<-function(nchrs, tree, delta) {
	datMat<-matrix(0, nrow=nrow(tree$edge)+1, ncol=0)
	obsMat<-matrix(0, nrow=nrow(tree$edge), ncol=0)
	realMat<-matrix(0, nrow=nrow(tree$edge), ncol=0)
	for (i in 1:nchrs) {
		thisChr<-simulateOneCharacter(tree, delta)
		datMat<-cbind(datMat, thisChr$dat)
		obsMat<-cbind(obsMat, thisChr$chObs)
		realMat<-cbind(realMat, thisChr$chReal)
	}
	return (list(dat=datMat, chObs=obsMat, chReal=realMat))
}