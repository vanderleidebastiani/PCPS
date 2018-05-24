#' @title Principal Coordinates of Phylogenetic Structure
#' 
#' @description Function to generate Principal Coordinates of Phylogenetic Structure (PCPS).
#' 
#' @details The function obtains a matrix containing phylogeny-weighted species composition 
#' (\code{\link{matrix.p}}) and is submitted to principal coordinates analysis (PCoA). 
#' This method generates the principal coordinates of phylogenetic structure 
#' (PCPS) (Duarte, 2011).
#'
#' The function summary or the function scores.pcps re-scales the correlation values for obtain the scores for \code{\link{biplot}} 
#' graphics. The function plot draws a simple biplot and represent clades as 
#' "spider" graphs (see \code{\link{ordispider}}).
#' 
#' @encoding UTF-8
#' @import SYNCSA
#' @importFrom stats cor
#' @importFrom vegan ordispider wcmdscale ordilabel vegdist
#' @importFrom graphics plot points text
#' @aliases pcps print.pcps summary.pcps print.summarypcps plot.pcps scores.pcps
#' @param comm Community data, with species as columns and sampling units as rows. 
#' This matrix can contain either presence/absence or abundance data.
#' @param phylodist Matrix containing phylogenetic distances between species.
#' @param method Dissimilarity index, as accepted by \code{\link{vegdist}} (Default dist="bray").
#' @param squareroot Logical argument (TRUE or FALSE) to specify if use square root of 
#' dissimilarity index (Default squareroot = TRUE).
#' @param correlations Logical argument (TRUE or FALSE) to specify if are calculed the correlations
#' between each PCPS and each species in matrix P (Default correlations = TRUE).
#' @param object An object of class pcps.
#' @param x An object of class pcps.
#' @param choices Axes for re-scaling. Choices must have length equal to two (Default choices = c(1, 2)).
#' @param display Display text or points for the sampling units, partial match to "text" or "points" (Default display = "text").
#' @param groups Factor giving the groups (Clades) for each species  (Default groups = NULL).
#' @param showlabel Label the groups by their names in the centroid of the object.
#' @param ... Other parameters for the respective functions.
#' @return \item{P}{Phylogeny-weighted species composition matrix.} \item{values}{The eigenvalues, 
#' relative eigenvalues and cumulative relative eigenvalues.} \item{vectors}{The principal coordinates
#' of phylogenetic structure (PCPS).} \item{correlations}{Correlations between a PCPS axis and 
#' phylogenetically weighted species abundances or frequencies.} \item{scores}{Scores for biplot graphics.}
#' @note \strong{IMPORTANT}: The sequence species show up in the community data matrix MUST be the 
#' same as they show up in the phylogenetic distance matrix. See \code{\link{organize.syncsa}}.
#' @author Vanderlei Julio Debastiani <vanderleidebastiani@@yahoo.com.br>
#' @seealso \code{\link{matrix.p}}, \code{\link{wcmdscale}}, \code{\link{ordispider}}, \code{\link{ordilabel}}
#' @references Duarte, L.S. (2011). Phylogenetic habitat filtering influences forest nucleation 
#' in grasslands. Oikos, 120, 208:215.
#' @keywords PCPS
#' @examples
#'
#' data(ADRS)
#' res<-pcps(ADRS$community, ADRS$phylo)
#' res
#' summary(res)
#' summary(res, choices = c(1, 2))$scores
#' plot(res, display = "text", groups = c(rep("Clade-A", 2), rep("Clade-B", 4)))
#'
#' @export
pcps <- function(comm, phylodist, method = "bray", squareroot = TRUE, correlations = TRUE){
	P <- SYNCSA::matrix.p(comm, phylodist, notification = FALSE)$matrix.P
	ord.P <- wcmdscale.org(P, method = method, squareroot = squareroot, eig = TRUE, correlations = correlations)
	res <- ord.P
	res$call <- match.call()
	res$P <- P
	colnames(res$vectors) <- paste("pcps.", seq_len(ncol(res$vectors)), sep = "")
	row.names(res$values) <- colnames(res$vectors)
	if(correlations){
	  colnames(res$correlations) <- paste("pcps.", seq_len(ncol(res$vectors)), sep = "")
		rownames(res$correlations) <- rownames(res$correlations, do.NULL = FALSE, prefix = "spp.")
	}
	class(res) <- "pcps"
	return(res)
}