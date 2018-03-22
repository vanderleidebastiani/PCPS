#' @title Curve of phylogenetic signal at metacommunity level
#' 
#' @description The function estimate the phylogenetic signal at metacommunity level and draws
#' a representation curve.
#' 
#' @details The PCPS are used, in a sequential manner, as predictors in a linear regression
#' to model the trait averages across the metacommunity. The curve is drawn as the
#' percentage of cumulative eigenvalues in the abscissa and as the determination 
#' coefficient of regressions in the ordinate.
#' 
#' Two null models are available. The first one (ts), the null curves are generated
#' shuffling terminal tips across the phylogenetic tree, generates a set of random PCPS
#' and recalculates the curves. The second (bm), the null curves are generated with 
#' simulate traits evolving under Brownian motion model. 
#'
#' @encoding UTF-8
#' @importFrom ape rTraitCont
#' @importFrom vegan vegdist
#' @importFrom stats quantile
#' @importFrom RcppArmadillo fastLm
#' @importFrom graphics plot points segments
#' @importFrom parallel makeCluster clusterExport clusterApply stopCluster
#' @aliases pcps.curve print.pcpscurve summary.pcpscurve plot.pcpscurve
#' @param comm Community data, with species as columns and sampling units as rows. This 
#' matrix can contain either presence/absence or abundance data.
#' @param phylodist Matrix containing phylogenetic distances between species.
#' @param trait Matrix data of species described by traits, with traits as columns and species as rows.
#' @param method Dissimilarity index, as accepted by \code{\link{vegdist}} (Default dist = "bray").
#' @param squareroot Logical argument (TRUE or FALSE) to specify if use square root of dissimilarity
#' index (Default squareroot = TRUE).
#' @param ranks Logical argument (TRUE or FALSE) to specify if ordinal variables are 
#' convert to ranks (Default ranks = TRUE).
#' @param null.model.ts Logical argument (TRUE or FALSE) to specify if use null model that shuffles
#' terminal tips across the phylogenetic tree to generate null curves. See details (Default null.model.ts = FALSE).
#' @param null.model.bm Logical argument (TRUE or FALSE) to specify if use null model that simulate 
#' trait evolving under Brownian motion to generate null curves. See details (Default null.model.bm = FALSE).
#' @param tree Phylogenetic tree, as phylo object.
#' @param runs Number of randomizations.
#' @param progressbar Logical argument (TRUE or FALSE) to specify if display a progress bar 
#' on the R console (Default progressbar = FALSE).
#' @param parallel Number of parallel processes.  Tip: use detectCores() (Default parallel = NULL).
#' @param newClusters Logical argument (TRUE or FALSE) to specify if make new parallel 
#' processes or use predefined socket cluster. Only if parallel is different of NULL (Default newClusters = TRUE).
#' @param CL A predefined socket cluster done with parallel package.
#' @param values The eigenvalues, relative eigenvalues and cumulative relative eigenvalues returned by \code{\link{pcps}}. 
#' @param vectors The principal coordinates of phylogenetic structure returned by \code{\link{pcps}}.
#' @param mt Matrix containing trait average at community level for one trait.
#' @param object An object of class pcpscurve.
#' @param x An object of class pcpscurve.
#' @param probs Numeric vector of probabilities used by \code{\link{quantile}}. (Default probs = c(0.025, 0.975)).
#' @param type Type of the plot to be drawn (Default type = "b").
#' @param draw.model Type of null model to draw; none (none), taxa shuffle (ts), browian motion model (bm).
#' @param col Plot color.
#' @param model.col Color of lines of null models.
#' @param ... Further graphical parameters for points.
#' @return \item{curve.obs}{The cumulative PCPS eigenvalues and the coefficient of determination.}
#' \item{curve.null.ts}{The cumulative PCPS eigenvalues and the coefficient of determination for 
#' each randomization using the taxa shuffle null model.} \item{curve.null.bm}{The cumulative PCPS 
#' eigenvalues and the coefficient of determination for each randomization using the Brownian motion null model.}
#' @author Vanderlei Julio Debastiani <vanderleidebastiani@@yahoo.com.br>
#' @seealso \code{\link{matrix.p}}, \code{\link{pcps}}
#' @references Duarte, L.S. (2011). Phylogenetic habitat filtering influences forest nucleation
#' in grasslands. Oikos, 120, 208:215.
#' @keywords PCPS
#' @examples
#' 
#' data(flona)
#' res<-pcps.curve(flona$community, flona$phylo, flona$trait[,1,drop = FALSE], 
#'        null.model.ts = TRUE, runs = 9)
#' res
#' summary(res)
#' plot(res, draw.model = "ts", type = "b", col = "red")
#'
#' @export
pcps.curve<-function(comm, phylodist, trait, method = "bray", squareroot = TRUE, ranks = TRUE, null.model.ts = FALSE, null.model.bm = FALSE, tree, runs = 99, progressbar = FALSE, parallel = NULL, newClusters = TRUE, CL =  NULL){
	RES <- list(call= match.call())
	if(ncol(trait)!=1){
		stop("\n Only one trait is allowed\n")
	}
	MT <- SYNCSA::matrix.t(comm, trait, scale = FALSE, ranks = ranks, notification = FALSE)$matrix.T
	res.pcps <- pcps(comm, phylodist, method = method, squareroot = squareroot, correlations = FALSE)
	res.values <- res.pcps$values
	res.vectors <- res.pcps$vectors
	curve.obs <- pcpc.curve.calc(res.values, res.vectors, MT)
	rownames(curve.obs) <- rownames(res.values)
	RES$curve.obs <- curve.obs
	if(progressbar){
		if(null.model.ts & null.model.bm){
			BarRuns <- runs*2
		}else{
			BarRuns <- runs
		}
	}
	if(!is.null(CL)){
		parallel <- length(CL)
	}
	ptest.ts <- function(samp, comm, phylodist, method, squareroot, mt){
		pcps.null <- PCPS::pcps(comm, phylodist[samp, samp], method = method, squareroot = squareroot, correlations = FALSE)
		res <- PCPS::pcpc.curve.calc(pcps.null$values, pcps.null$vectors, mt)
		return(res)
	}
	ptest.bm <- function(samp, tree, comm, values, vectors, ranks){
		trait.null <- cbind(ape::rTraitCont(phy = tree, model = "BM"))
		match.names <- match(colnames(comm), rownames(trait.null))
		MT.null <- SYNCSA::matrix.t(comm, trait.null[match.names,,drop = FALSE], scale = FALSE, ranks = ranks, notification = FALSE)$matrix.T
        res <- PCPS::pcpc.curve.calc(values, vectors, MT.null)
		return(res)
	}
	if((null.model.ts | null.model.bm) & !is.null(parallel) & newClusters){
		CL <- parallel::makeCluster(parallel, type = "PSOCK")
	}	
	if(null.model.ts){
		seqpermutation <- SYNCSA::permut.vector(ncol(phylodist), nset = runs)
		seqpermutation <- lapply(seq_len(nrow(seqpermutation)), function(i) seqpermutation[i,])
		if(is.null(parallel)){
	    	res.curve.null.ts<-vector("list",runs)
		    for (i in 1:runs) {
	    	    res.curve.null.ts[[i]] <- ptest.ts(samp = seqpermutation[[i]], comm = comm, phylodist = phylodist, method = method, squareroot = squareroot, mt = MT)   
	    	    if(progressbar){
					    SYNCSA::ProgressBAR(i,BarRuns,style=3)
				}
    		}
		} else {
			res.curve.null.ts <- parallel::clusterApply(CL, seqpermutation, ptest.ts, comm = comm, phylodist = phylodist, method = method, squareroot = squareroot, mt = MT)		
		}	
		RES$curve.null.ts <- res.curve.null.ts
	}
	if(null.model.bm){
		seqpermutation <- vector("list",runs)
		if(is.null(parallel)){
			res.curve.null.bm <- vector("list",runs)
		    for (i in 1:runs) {
	    	    res.curve.null.bm[[i]] <- ptest.bm(NULL, tree, comm, res.values, res.vectors, ranks = ranks)
	    	    if(progressbar){
					SYNCSA::ProgressBAR(i+runs,BarRuns,style=3)
				}
    		}
		} else {
			res.curve.null.bm <- parallel::clusterApply(CL, seqpermutation, ptest.bm, tree = tree, comm = comm, values = res.values, vectors = res.vectors, ranks = ranks)		
		}
		RES$curve.null.bm <- res.curve.null.bm
	}
	if((null.model.ts | null.model.bm) & !is.null(parallel) & newClusters){
		parallel::stopCluster(CL)
	}
	class(RES) <- "pcpscurve"
	return(RES)	
}