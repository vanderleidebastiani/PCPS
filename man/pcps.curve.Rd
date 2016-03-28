\encoding{UTF-8}
\name{pcps.curve}
\alias{pcps.curve}
\alias{print.pcpscurve}
\alias{summary.pcpscurve}
\alias{plot.pcpscurve}

\title{Curve of phylogenetic signal at metacommunity level}

\description{
The function estimate the phylogenetic signal at metacommunity level and draws a representation curve.
}

\usage{
pcps.curve(comm, dist.spp, trait, method = "bray", squareroot = TRUE,
        null.model.ts = FALSE, null.model.bm = FALSE, tree, runs = 99,
        progressbar = FALSE)
\method{print}{pcpscurve}(x, ...)
\method{summary}{pcpscurve}(object, probs = c(0.025, 0.975), ...)
\method{plot}{pcpscurve}(x, type = "b", draw.model = c("none","ts","bm"),
        probs = c(0.025, 0.975), col = "black", model.col = "black", ...)
}

\arguments{

\item{comm}{Community data, with species as columns and sampling units as rows. This matrix can contain either presence/absence or abundance data.}
\item{dist.spp}{Matrix containing phylogenetic distances between species.}
\item{trait}{Matrix data of species described by traits, with traits as columns and species as rows.}
\item{method}{Dissimilarity index, as accepted by \code{\link{vegdist}} (Default dist = "bray").}
\item{squareroot}{Logical argument (TRUE or FALSE) to specify if use square root of dissimilarity index (Default squareroot = TRUE).}
\item{null.model.ts}{Logical argument (TRUE or FALSE) to specify if use null model that shuffles terminal tips across the phylogenetic tree to generate null curves. See details (Default null.model.ts = FALSE).}
\item{null.model.bm}{Logical argument (TRUE or FALSE) to specify if use null model that simulate trait evolving under Brownian motion to generate null curves. See details (Default null.model.bm = FALSE).}
\item{tree}{Phylogenetic tree, as phylo object.}
\item{runs}{Number of randomizations.}
\item{progressbar}{Logical argument (TRUE or FALSE) to specify if display a progress bar on the R console (Default progressbar = FALSE).}
\item{object, x}{An object of class pcpscurve.}
\item{probs}{Numeric vector of probabilities used by \code{\link{quantile}}. (Default probs = c(0.025, 0.975)).}
\item{type}{Type of the plot to be drawn (Default type = "b").}
\item{draw.model}{Type of null model to draw; none (none), taxa shuffle (ts), browian motion model (bm).}
\item{col}{Plot color.}
\item{model.col}{Color of lines of null models.}
\item{...}{Further graphical parameters for points.}
}

\details{

The PCPS are used, in a sequential manner, as predictors in a linear regression to model the trait averages across the metacommunity. The curve is drawn as the percentage of cumulative eigenvalues in the abscissa and as the determination coefficient of regressions in the ordinate.

Two null models are available. The first one (ts), the null curves are generated shuffling terminal tips across the phylogenetic tree, generates a set of random PCPS and recalculates the curves. The second (bm), the null curves are generated with simulate traits evolving under Brownian motion model. 
}

\value{
\item{curve.obs}{The cumulative PCPS eigenvalues and the coefficient of determination.}
\item{curve.null.ts}{The cumulative PCPS eigenvalues and the coefficient of determination for each randomization using the taxa shuffle null model.}
\item{curve.null.bm}{The cumulative PCPS eigenvalues and the coefficient of determination for each randomization using the Brownian motion null model.}
}

\references{
	
Duarte, L.S. (2011). Phylogenetic habitat filtering influences forest nucleation in grasslands. Oikos, 120, 208–215.

}

\author{
Vanderlei Júlio Debastiani <vanderleidebastiani@yahoo.com.br>
}

\seealso{ \code{\link{matrix.p}}, \code{\link{pcps}}
}

\examples{
data(flona)
res_curve<-pcps.curve(flona$community, flona$phylo, flona$trait[,1], method = "bray",
        squareroot = TRUE, null.model.ts = TRUE, runs = 9, progressbar = FALSE)
res_curve
summary(res_curve)
plot(res_curve, type = "b", draw.model = "ts", col = "red")
}

\keyword{PCPS}