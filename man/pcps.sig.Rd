\encoding{UTF-8}
\name{pcps.sig}
\alias{pcps.sig}

\title{Association between PCPS axes and environmental predictors}

\description{
Analyze the association between PCPS axes and environmental predictors via null model.
}

\usage{
pcps.sig(comm, dist.spp, envir, method = "bray", 
        squareroot = TRUE, formula, runs = 999, AsFactors = NULL)
}

\arguments{
\item{comm}{Community data, with species as columns and sampling units as rows. This matrix can contain either presence/absence or abundance data.}
\item{dist.spp}{Matrix containing phylogenetic distances between species.}
\item{envir}{Environmental variables for each community, with variables as columns and sampling units as rows.}
\item{method}{Dissimilarity index, as accepted by \code{\link{vegdist}} (Default dist="bray").}
\item{squareroot}{Logical argument (TRUE or FALSE) to specify if use square root of dissimilarity index (Default squareroot = TRUE).}
\item{formula}{An object of class \code{\link{formula}} quotation marks. See Details.}
\item{runs}{Number of permutations for assessing significance.}
\item{AsFactors}{Encode an environmental variable as factor. See Details.}
}

\details{

The function generates principal coordinates of phylogenetic structure (\code{\link{pcps}}) and run a generalized linear model  (\code{\link{glm}}) with a gaussian error distribution (see \code{\link{family}}). The significance is obtained via a null model that shuffles terminal tips across the phylogenetic tree and generates a set of null PCPS. Each null PCPS is submitted to a procrustean adjustment (see \code{\link{procrustes}}), and the fitted values between observed PCPS and null PCPS is obtained. The adjusted null PCPS is used to rerun the model, generating a null F value. The observed F value is compared with the null F values to generate a probability value of the original F value being generated merely by chance.

The item formula is an expression of the form pcps.1 ~ model. The response term must be the pcps name, for example pcps.1, pcps.2, pcps.12.

The item AsFactors changes a environmental variable for the class \code{\link{factor}}. The sequence is the same that in the environmental data matrix. Use \code{\link{c}} to combine more that one variable.
}

\value{
\item{model}{The model, an object of class glm.}
\item{Envir_class}{The class of each variable in environmental data.}
\item{formula}{The formula used.}
\item{f.obs}{Observed F value in GLM.}
\item{p}{The p value.}
}

\references{
	
Duarte, L.S. (2011). Phylogenetic habitat filtering influences forest nucleation in grasslands. Oikos, 120, 208–215.

}

\author{
Vanderlei Júlio Debastiani <vanderleidebastiani@yahoo.com.br>
}

\seealso{ \code{\link{matrix.p}}, \code{\link{pcps}}, \code{\link{procrustes}}, \code{\link{glm}} 
}

\examples{
data(flona)
pcps.sig(flona$community, flona$phylo, flona$environment,
        formula = "pcps.1~alt", runs = 100)
}

\keyword{PCPS}