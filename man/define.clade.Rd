\encoding{UTF-8}
\name{define.clade}
\alias{define.clade}

\title{Define clade}

\description{
Function to define groups (clades) in a phylogenetic tree. 
}

\usage{
define.clade(tree, threshold, time, method = c("threshold", "time"))
}

\arguments{
\item{tree}{Phylogenetic tree.}
\item{threshold}{A threshold value to form the groups.}
\item{time}{A cutting height (age) to form the groups.}
\item{method}{Method to define the clades, "threshold" or "time".}
}

\details{
In the method threshold the total length of phylogenetic tree is used as cutting factor. If threshold is near to zero the cutting is near the root, if threshold near to one cutting is near the tips.

The phylogenetic tree must contain the node labels for the function work. Use the \code{\link{makeNodeLabel}} for defining node labels in a flexible way.

}
\value{
\item{clades}{Tips and their clades.}
\item{height}{The cutting height.}
}

\author{
Vanderlei Júlio Debastiani <vanderleidebastiani@yahoo.com.br>
}

\seealso{ \code{\link{makeNodeLabel}}
}

\examples{
require(ape)	
tree<-makeNodeLabel(rcoal(20))
clades<-define.clade(tree, threshold = 0.8, method = "threshold")
clades
plot.phylo(tree, show.node.label = TRUE)
abline(v = clades$height)
}

\keyword{PCPS}