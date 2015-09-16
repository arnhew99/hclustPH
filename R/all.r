hclustPH <- function(Z) {
	
	n.vertices <- dim(Z)[1]
	
	data.d <- dim(Z)[2]
	if (data.d >= 9) stop(cat("Data dimension must be 8 or lower\n"))
	
	# Compute the delaunay triangulation, and extract a matrix with all the edges
	# [NB (a,b) == (b,a)]
	# Construct a list containing the edges, which represents the simplicial complex. 
	

	# delaunayn comes from the geometry package - uses the QHull algorithm and is
	# well-behaved for sets of points in 8 dimensions or fewer...
	triang <- delaunayn(Z, options="Qt Pp")
	triang <- apply(triang,1,sort.int)
	
	edges <- matrix(as.vector(apply(triang,2,function(x) combn(x,2))),ncol=2,byrow=TRUE)
	edges <- unique(edges)
	n.edges <- dim(edges)[1]
	
	
	
	# Compute the edge lengths. These are referred to as the edges "birth times" in the 
	# PH literature.
	p1 <- Z[edges[,1],]
	p2 <- Z[edges[,2],]
	ds <- sqrt(rowSums((p1 - p2)^2))
	
	# reorder the list according to birth times
	ds.order <- order(ds)
	birth <- ds[ds.order] 
	edges <- edges[ds.order,]
	
	simplices <- lapply(1:n.edges, function(i) return(as.integer(edges[i,])))
	
	
	# The following reduction is essentially a sparse implementation of Gaussian
	# elimination over the field Z_2.  If any pair of edges share the same 
	# vertex with highest index, then we eliminate that vertex from the edge with
	# later birth time by modulo 2 addition of the two sets of vertices. This 
	# repeats until no pair of edges meet the condition described.	
	# The order of the edges matters, which is why we sorted by birth time.

	### REDUCTION ###
	top <- function(vec) {
		if (is.na(vec[1])) -1 else max(vec)
	}

	
	# initialise
	top.simplexes <- unlist(lapply(simplices,top))
	
	
	# # add two vectors together mod 2
	# add.simplices <- function(x,y) {
		# c(x,y)[!c(match(x,y,0L),match(y,x,0L))]
	# }
	
	
	
	dups <- duplicated(top.simplexes,incomparables=-1)
	indexvec <- 1:length(top.simplexes)
	
	while (any(dups)) {
		
		adds <- match(top.simplexes[dups], top.simplexes)
		reps <- indexvec[dups]
		
		# solve the addition
		replacement <- mapply(addSimplices,simplices[reps],simplices[adds],SIMPLIFY=FALSE)
				
		# put the solutions back in the correct place 
		simplices[reps] <- replacement
		
		top.simplexes[reps] <- unlist(lapply(replacement, top))
		dups <- duplicated(top.simplexes,incomparables=-1)
	}
	
	class(simplices) <- "reduced.matrix"
	### END OF REDUCTION ###

	
	# the reduced form needs forcing into the right format so that we can plot it
	# using the generic plot.hclust() method
	# hcass <- makeHclustMerge(reduction=simplices, n=n.vertices)
	
	Y <- matrix(unlist(simplices), nrow=2)
	Y <- apply(Y,2,sort.int)
	
	left <- as.integer(Y[1,])
	right <- as.integer(Y[2,])
	
	
	ia <- hclustIntermediate(left, right)
	hcass <- .Fortran("hcass2", n = n.vertices, ia = as.integer(c(ia,0)), ib = as.integer(c(right,0)), order = integer(n.vertices), iia = integer(n.vertices), iib = integer(n.vertices))
	
	
	# compute the hierachical clustering "heights", which are the birth times of any 
	# edge left in the reduced form
	simplex.size <- sapply(simplices,length)
	death <- which(simplex.size >= 1)
	
	
	
	# output a list which is compatible with the generic hclust method
	hc <- list()
	hc$merge <- cbind(hcass$iia[1L:(n.vertices - 1)], hcass$iib[1L:(n.vertices - 1)])
	hc$height <- birth[death]
	hc$order <- hcass$order
	hc$method <- "single"
	hc$labels <- NULL  # we could handle this at some point
	hc$dist.method <- NULL
	hc$call <- NULL
	class(hc) <- "hclust"
	
	return(hc)

}

makeHclustMerge <- function(reduction,n) {

	# in the stats:::hclust Fortran routines, an intermediate form of the
	# clustering result is formed by the "HCLUST" function, which is then
	# converted into the output form by "HCASS2", which is what plot.hclust
	# expects. This function turns our reduced form into the intermediate 
	# form that would come from HCLUST, and then uses the hclustPH copy of 
	# the HCASS2 Fortran code to convert this into the output form.
	
	Z <- matrix(unlist(reduction), nrow=2)
	Y <- apply(Z,2,sort.int)
	
	left <- as.integer(Y[1,])
	right <- as.integer(Y[2,])
	
	recurseSearch <- function(x, i, l, r) {

		nexti <- which(r[1:(i-1)]==x)
		nextx <- l[nexti]
		if (length(nexti) == 0) return(x) else return(recurseSearch(nextx, i, l, r))
	}
	
	
	dups <- duplicated(as.vector(Y, mode="integer"))[seq(1,2*n-2,by=2)]
	ia <- left
	ia[dups] <- rev(mapply(recurseSearch, x=rev(left[dups]), i=rev((1:(n-1))[dups]), MoreArgs=list(l=left, r=right)))
	
	# hcass2 is a copy of the Fortran routine from the base package "stats"
	# it outputs an ordering of the leaves in the plot so that tree branches
	# don't overlap, which is the "order" output, and a different description
	# of the sequence of merges which is used by the plot.hclust() method
	hcass <- .Fortran("hcass2", n = n, ia = as.integer(c(ia,0)), ib = as.integer(c(right,0)), order = integer(n), iia = integer(n), iib = integer(n))
	
	return(hcass)

}

