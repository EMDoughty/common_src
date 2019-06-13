getDescNodeList<-function(thisNode, tree, includeThis=TRUE) {
	nodeList<-getDescNodeListMachine(thisNode, tree)
	if (includeThis) nodeList<-c(thisNode, nodeList)
	return (nodeList)
}

getDescNodeListMachine<-function(thisNode, tree) {
	thisDesc<-tree$edge[tree$edge[,1]%in%thisNode,2]
	if (any(thisDesc>length(tree$tip.label))) thisDesc<-c(thisDesc, getDescNodeListMachine(thisDesc[thisDesc>length(tree$tip.label)], tree))
	return(thisDesc)
}

is.polytomy<-function(node, tree) {
	if (is.vector(node)) thisAnc<-node[1] else thisAnc<-node
	length(grep(thisAnc,tree$edge[,1]))>2
}

getAllPolyNodes<-function(tree) {
	ancList<-unique(tree$edge[,1])
	return(ancList[which(sapply(ancList, is.polytomy, tree))])
}

getAllPolyEdges<-function(tree) {
	return(which(apply(tree$edge, 1, is.polytomy, tree)))
}

permutationMachine<-function(thisVec, thisMat, n) {
	x<-seq(1:n)
	x<-x[!x%in%thisVec]
	for (i in 1:length(x)) {
		thisVec[n-length(x)+1]=x[i]
		if (length(x)>1) { thisMat<-permutationMachine(thisVec, thisMat, n) 
		} else thisMat<-rbind(thisMat, thisVec)
	}
	return(thisMat)
}

perm<-function(n, k) {
	combMat<-combn(n,k)
	thisMat<-matrix(nrow=0, ncol=k)
	thisVec<-vector(mode="numeric", length=k)
	thisMat<-permutationMachine(thisVec, thisMat, k)

	permMat<-matrix(nrow=0, ncol=k)
	for (i in 1:ncol(combMat)) {
		for (j in 1:nrow(thisMat)) {
			permMat<-rbind(permMat, combMat[,i][thisMat[j,]])
		}
	}
	return(permMat)
}

permFull<-function(n) {
	thisMat<-matrix(nrow=0, ncol=n)
	thisVec<-vector(mode="numeric", length=n)
	permutationMachine(thisVec, thisMat, n)
}

# masterResList<-vector(mode="numeric", length=0)
# polyList<-getAllPolyNodes(tree)
# for (i in 1:length(polyList)) {
	# nodeList<-which(tree$edge[,1]==polyList[i])
	# perms<-permFull(length(nodeList))
	# if (length(masterResList>0)) {
		# dummy<-vector(mode="numeric", length=0)
		# for (j in 1:nrow(masterResList)) {
			# for (k in 1:nrow(perms)) {
				# dummy<-rbind(dummy, c(masterResList[j,], nodeList[perms[k,]]))
			# }
		# }
		# masterResList<-dummy
	# } else masterResList<-matrix(nodeList[perms], ncol=length(nodeList))
# }

getMaxPairs <- function(tree) {
	spPairs <- array(dim=c(0,2))
	while(length(tree$tip.label) > 1) {
		nodesWithSisterTips <- sort(unique(sapply(seq(from=(length(tree$tip.label)+1), to=max(tree$edge[,1])), FUN=function(x) if(all(tree$edge[tree$edge[,1]==x,2] <= length(tree$tip.label))) x else NA)))
		thisPairs <- t(sapply(nodesWithSisterTips, function(y) tree$tip.label[tree$edge[tree$edge[,1]==y,2]]))
		spPairs <- rbind(spPairs, thisPairs)
		if (length(tree$tip.label) - length(thisPairs) > 1) tree <- drop.tip(tree, as.vector(thisPairs)) else break
		thisPairs <- NULL
	}
	spPairs
}
