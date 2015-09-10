hclustPH <- function(Z) {
	

	n.vertices <- dim(Z)[1]
	v.list <- 1:n.vertices
	
	# compute the delaunay triangulation
	triang <- delaunayn(Z, options="Qt Pp")
	
	# need vertices small to large
	triang <- apply(triang,1,sort)
	
	# extract some guesses
	edges <- matrix(as.vector(apply(triang,2,function(x) combn(x,2))),ncol=2,byrow=TRUE)
	
	edges <- unique(edges)

	
	n.edges <- dim(edges)[1]
	
	
	simplices <- vector("list",length=1+n.vertices)
	

	
	# compute the edge lengths
	p1 <- Z[edges[,1],]
	p2 <- Z[edges[,2],]
	ds <- sqrt(rowSums((p1 - p2)^2))
	
	

	
	edge.list <- lapply(1:(dim(edges)[1]), function(i) return(edges[i,]))

	
	# reorder the list according to birth times
	ds.order <- order(ds)
	edge.list.ordered <- edge.list[ds.order]
	
	ds <- ds[ds.order] 
	
	
	simplices[[1]] <- vector("integer",length=0)
	simplices[v.list+1] <- 1
	simplices <- c(simplices, edge.list.ordered)
	
	
	birth <- c(rep(0,n.vertices+1), ds)

	reduced.list <- simplices
	
	# function that returns the lowest 1 (highest number) in a simplex
	top <- function(vec) {
		if (is.na(vec[1])) -1 else max(vec)
	}

	
	# initialise
	top.simplexes <- sapply(reduced.list,top)
	
	
	# add two vectors together mod 2
	add.simplices <- function(x,y) {
		c(x,y)[!c(match(x,y,0L),match(y,x,0L))]
	}
	
	dups <- duplicated(top.simplexes,incomparables=-1)
	indexvec <- 1:length(top.simplexes)
	
	while (any(dups)) {
		
		adds <- match(top.simplexes[dups], top.simplexes)
		reps <- indexvec[dups]
		
		# solve the addition
		replacement <- mapply(add.simplices,reduced.list[reps],reduced.list[adds],SIMPLIFY=FALSE)
				
		# put the solutions back in the correct place 
		reduced.list[reps] <- replacement
		
		top.simplexes[reps] <- sapply(replacement, top)
		dups <- duplicated(top.simplexes,incomparables=-1)
	}
	
	class(reduced.list) <- "reduced.matrix"
	
	n <- dim(Z)[1]
	
	hcass <- makeHclustMerge(reduction=reduced.list, n=n)
	
	hc <- list()
	hc$merge <- cbind(hcass$iia[1L:(n - 1)], hcass$iib[1L:(n - 1)])
	
	
	# compute the heights
	# find the size of each simplex
	simplex.size <- sapply(reduced.list[-2],length)
	# orig.simplex.size <- sapply(simplices,length)
	
	# death <- intersect(which(simplex.size >= 1),which(orig.simplex.size == 2))
	death <- which(simplex.size >= 1)
	hc$height <- birth[death]
	
	
	
	# hc$height <- findHeights(reduced.list,simplices)
	hc$order <- hcass$order
	class(hc) <- "hclust"
	
	return(hc)

}

makeHclustMerge <- function(reduction,n) {
	
	Z <- matrix(unlist(reduction[-2]), nrow=2)
	Y <- apply(Z,2,sort.int)
	
	left <- as.integer(Y[1,])
	right <- as.integer(Y[2,])
	
	recurseSearch <- function(x, i, l, r) {

		nexti <- which(r[1:(i-1)]==x)
		# nexti <- match(x, r[1:(i-1)])
		nextx <- l[nexti]
		# test <- r[1:(i-1)]==x
		# nextx <- (l[1:(i-1)])[test]
		if (length(nexti) == 0) return(x) else return(recurseSearch(nextx, i, l, r))
		# if (is.na(nexti)) return(x) else return(recurseSearch(nextx, i, l, r)) 
	}
	# recurseSearch <- cmpfun(recurseSearch)
	
	
	dups <- duplicated(as.vector(Y, mode="integer"))[seq(1,2*n-2,by=2)]
	ia <- left
	ia[dups] <- rev(mapply(recurseSearch, x=rev(left[dups]), i=rev((1:(n-1))[dups]), MoreArgs=list(l=left, r=right)))

	hcass <- .Fortran(stats:::C_hcass2, n = n, ia = as.integer(c(ia,0)), ib = as.integer(c(right,0)), order = integer(n), iia = integer(n), iib = integer(n))
	
	return(hcass)

}

