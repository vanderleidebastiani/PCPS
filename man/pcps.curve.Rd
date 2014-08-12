\encoding{UTF-8}
\name{pcps.curve}
\alias{pcps.curve}
\alias{print.pcpscurve}
\alias{plot.pcpscurve}

\title{Curve of phylogenetic signal at metacommunity level}

\description{
The function estimate the phylogenetic signal at metacommunity level and draws a representation curve.
}

\usage{
pcps.curve(comm, dist.spp, trait, method = "bray", squareroot = TRUE,
        null.model = TRUE, runs = 99, progressbar = FALSE)
\method{print}{pcpscurve}(x, ...)        
\method{plot}{pcpscurve}(x, type = "b", errorbars = c("none", "sd", "se", "quantile"),
        probs = c(0.025, 0.975), col = "black", errbar.col = "black",
        errbar.pch = 16, ...)
}

\arguments{

\item{comm}{Community data, with species as columns and sampling units as rows. This matrix can contain either presence/absence or abundance data.}
\item{dist.spp}{Matrix containing phylogenetic distances between species.}
\item{trait}{Matrix data of species described by traits, with traits as columns and species as rows.}
\item{method}{Dissimilarity index, as accepted by \code{\link{vegdist}} (Default dist="bray").}
\item{squareroot}{Logical argument (TRUE or FALSE) to specify if use square root of dissimilarity index (Default squareroot = TRUE).}
\item{null.model}{Logical argument (TRUE or FALSE) to specify if use null model to generate null curves (Default null.model = TRUE)}
\item{runs}{Number of randomizations.}
\item{progressbar}{Logical argument (TRUE or FALSE) to specify if display a progress bar on the R console (Default progressbar = FALSE).}
\item{x}{An object of class pcpscurve.}
\item{type}{Type of the plot to be drawn (Default type="b").}
\item{errorbars}{Plot type of error bars; none (none), standard deviation (sd), standard error (se)  and quantiles corresponding to the given probabilities (quantile).}
\item{probs}{Numeric vector of probabilities used by \code{\link{quantile}}. (Default probs=c(0.025,0.975)).}
\item{col}{Plot color.}
\item{errbar.col}{Color of error bars.}
\item{errbar.pch}{Symbol to use of error bars.}
\item{...}{Further graphical parameters for points.}
}

\details{

The PCPS are used, in a sequential manner, as predictors in a linear regression to model the trait averages across the metacommunity. The curve is drawn as the percentage of cumulative eigenvalues in the abscissa and as the determination coefficient of regressions in the ordinate.

The null curves are generated via a null model that shuffles terminal tips across the phylogenetic tree, generates a set of random PCPS and recalculates the curves.

}

\value{
\item{curve.obs}{The cumulative PCPS eigenvalues and the coefficient of determination.}
\item{curve.null}{The cumulative PCPS eigenvalues and the coefficient of determination for each randomization.}
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
        squareroot = TRUE, null.model = TRUE, runs = 9, progressbar = FALSE)
res_curve
plot(res_curve, type = "b", errorbars = "se", col = "red",
        errbar.col = "black", errbar.pch = 21,pch = 19)
}

\keyword{PCPS}