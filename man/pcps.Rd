\encoding{UTF-8}
\name{pcps}
\alias{pcps}
\alias{print.pcps}
\alias{summary.pcps}
\alias{scores.pcps}
\alias{plot.pcps}

\title{Principal Coordinates of Phylogenetic Structure}

\description{

Function to generate Principal Coordinates of Phylogenetic Structure (PCPS).

}

\usage{
pcps(comm, dist.spp, method = "bray", squareroot = TRUE)
\method{print}{pcps}(x, ...)
\method{summary}{pcps}(object, ...)
\method{scores}{pcps}(x, choices = c(1, 2), ...)
\method{plot}{pcps}(x, display = c("text", "points"), groups, showlabel = TRUE, ...)
}

\arguments{
\item{comm}{Community data, with species as columns and sampling units as rows. This matrix can contain either presence/absence or abundance data.}
\item{dist.spp}{Matrix containing phylogenetic distances between species.}
\item{method}{Dissimilarity index, as accepted by \code{\link{vegdist}} (Default dist="bray").}
\item{squareroot}{Logical argument (TRUE or FALSE) to specify if use square root of dissimilarity index (Default squareroot = TRUE).}
\item{object,x}{An object of class pcps.}
\item{choices}{ Axes for re-scaling. Choices must have length equal to two (Default choices = c(1,2)).}
\item{display}{Display text or points for the sampling units.}
\item{groups}{Factor giving the groups (Clades) for each species.}
\item{showlabel}{Label the groups by their names in the centroid of the object.}
\item{...}{Other parameters for the respective functions.}
}

\details{

The function obtains a matrix containing phylogeny-weighted species composition (\code{\link{matrix.p}}) and is submitted to principal coordinates analysis (PCoA). This method generates the principal coordinates of phylogenetic structure (PCPS) (Duarte, 2011).

The function scores.pcps re-scales the correlation values for \code{\link{biplot}} graphics. The function plot.pcps draws a simple biplot and represent clades as 'spider" graphs (see \code{\link{ordispider}}).
}

\value{
\item{P}{Phylogeny-weighted species composition matrix.}
\item{values}{The eigenvalues, relative eigenvalues and cumulative relative eigenvalues.}
\item{vectors}{The principal coordinates of phylogenetic structure (PCPS).}
\item{correlations}{Correlations between a PCPS axis and phylogenetically weighted species abundances or frequencies.}
}

\note{
\strong{IMPORTANT}: The sequence species show up in the community data matrix MUST be the same as they show up in the phylogenetic distance matrix. See \code{\link{organize.syncsa}}.
}

\references{
	
Duarte, L.S. (2011). Phylogenetic habitat filtering influences forest nucleation in grasslands. Oikos, 120, 208–215.

}

\author{
Vanderlei Júlio Debastiani <vanderleidebastiani@yahoo.com.br>
}

\seealso{ \code{\link{matrix.p}}, \code{\link{wcmdscale}} 
}

\examples{
data(flona)
res<-pcps(flona$community, flona$phylo)
res
summary(res)
scores(res)
plot(res, display = "text", groups = c(rep(1, 30), rep(2, 29)))
}

\keyword{PCPS}