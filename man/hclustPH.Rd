\name{hclustPH}
\alias{hclustPH}
\title{
Single linkage clustering using persistent homology
}
\description{
Hiearchical cluster analysis using the single linkage agglomeration method, computed using persistent homology. This method avoids computing the distance matrix directly. 
}
\usage{
hclustPH(X)
}
\arguments{
  \item{X}{Data matrix with observations as rows and column dimension of 8 or less.}
}
\details{
}
\value{
An object of class \bold{hclust} which describes the tree produced by the clustering process. The object is a list with components: 
\item{merge}{an \eqn{n-1} by 2 matrix. Row \eqn{i} of \code{merge} describes the merging of clusters at step \eqn{i} of the clustering. If an element \eqn{j} in the row is negative, then observation \eqn{-j} was merged at this stage. If \eqn{j} is positive then the merge was with the cluster formed at the (earlier) stage \eqn{j} of the algorithm. Thus negative entries in merge indicate agglomerations of singletons, and positive entries indicate agglomerations of non-singletons.}
\item{height}{a set of \eqn{n-1} real values (non-decreasing for ultrametric trees). The clustering height: that is, the value of the criterion associated with the single linkage method for the particular agglomeration.}
\item{order}{a vector giving the permutation of the original observations suitable for plotting, in the sense that a cluster plot using this ordering and matrix \code{merge} will not have crossings of the branches.}
\item{labels}{Currently NULL}
\item{call}{Currently NULL}
\item{method}{Currently NULL}
}
\references{
}
\author{
Matthew Arnold
}
\examples{
X <- matrix(runif(1000), ncol=2)
clu <- hclustPH(X)
}
