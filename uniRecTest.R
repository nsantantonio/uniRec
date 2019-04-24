# Killington skiing trip !!!

recombine <- function(chrom, trackRecomb = FALSE, force1recomb = FALSE){
	cM <- nrow(chrom)
	recomb <- c(as.logical(rpois(cM - 1, (cM -1) / 1e4)), FALSE)
	if(force1recomb) while(!any(recomb)) recomb <- c(as.logical(rpois(cM - 1, (cM -1) / 1e4)), FALSE)
	rsite <- which(recomb)
	sisters <- chrom	
	for(i in rsite) sisters[(i + 1):nrow(sisters), ] <- sisters[(i + 1):nrow(sisters), c(2, 1)]
	if(trackRecomb) sisters <- list(sisters = sisters, rsites = rsite)
	sisters
}

meiosis <- function(geno){
	chromatids <- lapply(geno, recombine)
	gamete <- lapply(chromatids, function(x) x[, sample(1:2,1)])
	gamete
}

pollinate <- function(gamete1, gamete2 = NULL){
	seed <- list()
	if(is.null(gamete2)) gamete2 <- lapply(gamete1, function(x) rep(FALSE, length(x))) 
	for (i in 1:length(gamete1)) seed[[i]] <- cbind(gamete1[[i]], gamete2[[i]])
	seed
}

cross <- function(geno1, geno2 = NULL, DH = FALSE, BC = FALSE){
	gamete1 <- meiosis(geno1)
	if(DH) {
		gamete2 <- gamete1
	} else if(!BC) {
		gamete2 <- if(is.null(geno2)) meiosis(geno1) else meiosis(geno2)
	} else {
		gamete2 <- NULL
	} 
	pollinate(gamete1, gamete2)
}

randomIntermate <- function(genolist, newpopsize = length(genolist)){
	newpop <- list()
	for (k in 1:newpopsize){
		flowers <- sample(1:length(genolist),2)
		newpop[[k]] <- cross(geno1 = genolist[[flowers[1]]], geno2 = genolist[[flowers[2]]])
	}
	return(newpop)
}


seqInd <- function(ind) lapply(ind, rowSums) # need wrapper function

getSeqMatrix <- function(alCounts, chr = NULL){
	if(is.null(chr)) chr <- 1:length(alCounts[[1]])
	gmat <- list()
	for(i in chr){
		gmat[[i]] <- do.call(rbind, lapply(alCounts, "[[", i))
	}
	gmat
}

genotype <- function(markers, seqMat) {
	Ml <- list()
	for(i in 1:length(markers)) Ml[[i]] <- seqMat[[i]][, markers[[i]]]
	Ml
}

getCO <- function(X){ # need wrapper function
	COl <- list() 
	for(i in 2:ncol(X)) COl[[i - 1]] <- X[, i - 1] != X[, i]
	do.call(cbind, COl)
}

haldane <- function(r) -0.5 * log(1 - 2*r)

recSel <- function(R, N, method = "uniRec", tR = TRUE){
	if(is.list(R)){
		if(tR) R <- lapply(R, t)
		R <- do.call(rbind, R)
	} else {
		if(tR) R <- t(R)
	}
	cM <- haldane(rowSums(R) / ncol(R)) 
	if(is.logical(R)) class(R) <- "numeric"
	cat(nrow(R), "markers scored on", ncol(R), "projeny\n")
	if(method == "maxRec"){
		cdotj <- colSums(R)
		sel <- order(-cdotj)[1:N]
	} else if (method == "uniRec"){
		m <- (matrix(1, 1, ncol(R)) %x% cM)
		d <- R / m
		ddotj <- colSums(d)
		sel <- which.max(ddotj)
		n <- 1
		while(n < N){
			S <- d[,sel, drop = FALSE]
			didot <- rowSums(S)
			p <- max(didot)
			R[, sel] <- 0
			tj <- (p - didot) %*% R
			sel <- c(sel, which.max(tj))
			n <- length(sel) 
		}
	} else {
		cat("Please specify 'uniRec' or 'maxRec' as the argument to method!\n")
	}
	sel
}


sampleMarkerPos <- function(nMarkers, cM = 100, nChrom = NULL, method = "uniform"){
	if(length(nMarkers) > 1 & is.null(nChrom)) nChrom <- length(nMarkers)
	if(length(cM) != nChrom & length(cM) == 1) cM <- rep(cM, nChrom)
	if(length(nMarkers) == 1) rem <- nMarkers %% nChrom else rem = 0
	if(length(nMarkers) == 1 & nChrom > 1) nMarkers <- rep(floor(nMarkers / nChrom), nChrom)
	if (rem > 0) nMarkers[[1]] <- nMarkers[[1]] + rem
	mPos <- list()
	if(method == "uniform"){
		for(i in 1:nChrom) mPos[[i]] <- floor(seq(1, cM[[i]], length.out = nMarkers[[i]])) # random sampling
	} else {
		for(i in 1:nChrom) mPos[[i]] <- sample(1:cM[[i]], nMarkers[[i]]) # random sampling
	}
	mPos
}

makeF1 <- function(markerPos, nChrom = NULL, cM = NULL){
	if(!is.list(markerPos)) markerPos <- list(markerPos)
	if(is.null(cM)) cM <- lapply(markerPos, max)
	if(is.null(nChrom)) nChrom <- length(markerPos)
	if(length(cM) != nChrom & length(cM) == 1) cM <- rep(cM, nChrom)
	if (length(markerPos) != nChrom & length(markerPos) == 1) markerPos <- rep(markerPos, nChrom)
	f1 <- list()
	for (i in 1:nChrom) f1[[i]] <- cbind(rep(TRUE, cM[i]), rep(FALSE, cM[i]))
	f1
}

makePop <- function(f1, popSize = 100, type = "BC", gen = 6){
	dh <- if(type == "DH") TRUE else FALSE
	bc <- if(type == "BC") TRUE else FALSE
	pop <- lapply(1:popSize, function(x) cross(f1, DH = dh, BC = bc))
	if(type %in% c("RIL", "outcross")){
		f <- 2
		while(f <= gen){
			pop <- if(type == "RIL") lapply(pop, cross) else randomIntermate(pop)
			f <- f + 1
		}	
	}
	pop
}




markers <- sampleMarkerPos(nMarkers = c(11, 11, 11, 6, 6, 6), cM = 100)
f1 <- makeF1(markers)
BCpop <- makePop(f1, 100)
DHpop <- makePop(f1, 200, type = "DH")
RILpop <- makePop(f1, 1000, type = "RIL")

BCgeno <- lapply(BCpop, seqInd)
BCmat <- getSeqMatrix(BCgeno)

M <- genotype(markers, BCmat)
CO <- lapply(M, getCO)

selected <- recSel(CO, N = 20)
length(selected)


RILgeno <- lapply(RILpop, seqInd) # need wrapper function
RILmat <- getSeqMatrix(RILgeno)

M <- genotype(markers, RILmat)
CO <- lapply(M, getCO) # need wrapper function

selected <- recSel(CO, N = 20)
length(selected)


addNA <- function(M, propMiss = 0.3){
	np <- prod(dim(M))
	M[sample(1:np, np * propMiss)] <- NA 
	M
}

M <- genotype(markers, RILmat)
M <- lapply(M, addNA)
CO <- lapply(M, getCO) # need wrapper function

selected <- recSel(CO, N = 20)
length(selected)



# chromosomes <- rep(101, 6)
# nchrom <- length(chromosomes)
# markers <- c(rep(list(0:10 * 10 + 1), 3), rep(list(0:5 * 20 + 1), 3))


# Make Parents, where alleles are tracked through Parent 1 (designated as TRUE, verses FALSE for Parent 2)
# cM <- 100
# nChrom <- 6
# markerPos <- markers[[1]]
# nMarkers <- 63




# seel(RILpop)
# nBC <- 100
# BC <- list()
# for(i in 1:nBC) BC[[i]] <- cross(f1, p2)




# cM <- lapply(CO, colMeans)
# cM <- lapply(CO, colSums)
# lapply(cM, mean)
# R <- CO[[1]]
# R <- CO
# dim(R)

# CO <- lapply(CO, t)

# selected <- lapply(CO, recSel)
# lapply(selected, length)




	# for (i in 1:nchrom) p1[[i]] <- matrix(rep(TRUE, 2* chromosomes[i]), ncol = 2)
	# p2 <-lapply(p1, function(x) !x)
	# f1 <-lapply(p1, function(x) {x[, 2] <- !x[, 2]; x})
	# return(list(p1 = p1, p2 = p2, f1 = f1))
