#' @title Association between phylogeny-weighted species composition and environmental predictors
#' 
#' @description Analyses to relate an environmental gradient to the phylogenetic assembly of species 
#' across a metacommunity by means of phylogenetic fuzzy weighting.
#' 
#' @details Each metacommunity is submitted to phylogenetic fuzzy weighting, generating a matrix
#' that describing the phylogeny-weighted species composition of the communities
#' (\code{\link{matrix.p}}). The function matrix.p.sig test directly the association 
#' this matrix with the environmental predictors. The pairwise dissimilarities are 
#' submitted to Mantel test (\code{\link{mantel}}) or ADONIS test (\code{\link{adonis}})
#' to evaluate the influence of an environmental gradient on species dispersion across 
#' the communities. The function pcps.sig generates principal coordinates of phylogenetic
#' structure (\code{\link{pcps}}) and use a single axis for run a generalized linear 
#' model (GLM, \code{\link{glm}}) or use set of axis for run a distance-based redundancy
#' analysis (db-RDA, \code{\link{rda}}).
#' 
#' The significance is obtained via two null models, one that shuffles sites across the
#' environmental gradient and another that shuffles terminal tips (taxa) across the phylogenetic
#' tree. The first null model (site shuffle) shuffles the site position across the environmental
#' gradient and rerun the same model, generating a null F value (or r value in Mantel test). The
#' second null model (taxa shuffle), shuffles terminal tips across the phylogenetic tree and 
#' generates a null matrix containing phylogeny-weighted species composition and rerun the same
#' model, generating another null F value. In the pcps.sig function are generate set of null PCPS
#' and each null PCPS (or set of PCPS in RDA) is submitted to a procrustean adjustment 
#' (see \code{\link{procrustes}}), and the fitted values between observed PCPS and null PCPS is 
#' obtained. The adjusted null PCPS is used to rerun the model, generating another null F value. 
#' The observed F value (or r value) is compared independently with both null sets of F values 
#' (or r value) to generate a probability value of the original F value being generated merely by
#' chance according to each null model.
#' 
#' The item formula is an expression of the form pcps.1 ~ model. The response term must be the 
#' pcps name, for example pcps.1, pcps.2, pcps.12.
#' 
#' The item AsFactors changes a environmental variable for the class \code{\link{factor}}. The 
#' sequence is the same that in the environmental data matrix. Use \code{\link{c}} to combine 
#' more that one variable.
#' 
#' @encoding UTF-8
#' @import SYNCSA
#' @importFrom vegan procrustes rda adonis mantel vegdist
#' @importFrom parallel makeCluster clusterApply stopCluster parRapply clusterExport
#' @importFrom stats glm summary.lm as.formula fitted gaussian
#' @aliases pcps.sig matrix.p.sig
#' @param comm Community data, with species as columns and sampling units as rows. This matrix 
#' can contain either presence/absence or abundance data.
#' @param phylodist Matrix containing phylogenetic distances between species.
#' @param envir Environmental variables for each community, with variables as columns and 
#' sampling units as rows.
#' @param analysis Type of analysis. For the function pcps.sig \code{\link{glm}} or 
#' \code{\link{rda}}, for matrix.p.sig function \code{\link{adonis}} or \code{\link{mantel}}.
#' See Details.
#' @param method Dissimilarity index, as accepted by \code{\link{vegdist}} (Default dist = "bray").
#' @param squareroot Logical argument (TRUE or FALSE) to specify if use square root of 
#' dissimilarity index (Default squareroot = TRUE).
#' @param formula An object of class \code{\link{formula}} quotation marks used in GLM analysis. 
#' See Details.
#' @param family A description of the error distribution to be used in used in GLM analysis. 
#' See \code{\link{family}} (Dafault family = gaussian).
#' @param AsFactors Encode an environmental variable as factor used in GLM analysis. See Details.
#' @param pcps.choices PCPS used in RDA analysis (Default pcps.choices = c(1, 2, 3, 4)).
#' @param runs Number of permutations for assessing significance.
#' @param method.envir Resemblance index between communities based on environmental variables, 
#' as accepted by vegdist used in Mantel analysis (Default method.envir = "euclidean")
#' @param parallel Number of parallel processes.  Tip: use detectCores() (Default parallel = NULL).
#' @param newClusters Logical argument (TRUE or FALSE) to specify if make new parallel 
#' processes or use predefined socket cluster. Only if parallel is different of NULL (Default newClusters = TRUE).
#' @param CL A predefined socket cluster done with parallel package.
#' @return \item{model}{The model, an object of class glm, rda, adonis or mantel.}
#' \item{envir_class}{The class of each variable in environmental data in glm.}
#' \item{formula}{The formula used in glm.} \item{statistic.obs}{Observed F value or r value.}
#' \item{p.site.shuffle}{The p value for the site shuffle null model.}
#' \item{p.taxa.shuffle}{The p value for the taxa shuffle null model.}
#' @author Vanderlei Julio Debastiani <vanderleidebastiani@@yahoo.com.br>
#' @seealso \code{\link{matrix.p}}, \code{\link{pcps}}, \code{\link{procrustes}}, 
#' \code{\link{glm}}, \code{\link{rda}}, \code{\link{adonis}}, \code{\link{mantel}} 
#' @references Duarte, L.S. (2011). Phylogenetic habitat filtering influences forest 
#' nucleation in grasslands. Oikos, 120, 208:215.
#' @keywords PCPS
#' @examples
#' 
#' data(flona)
#' pcps.sig(flona$community, flona$phylo, flona$environment, analysis = "glm",
#'         formula = "pcps.1~alt", runs = 99)
#' matrix.p.sig(flona$community,flona$phylo,flona$environment[, 2, drop = FALSE],
#'         analysis = "adonis", runs = 99)
#' 
#' 
#' @export
pcps.sig<-function(comm, phylodist, envir, analysis = c("glm", "rda"), method = "bray", squareroot = TRUE, formula, family = stats::gaussian, AsFactors = NULL, pcps.choices = c(1, 2, 3, 4), runs = 999, parallel = NULL, newClusters = TRUE, CL =  NULL){
  RES <- list(call = match.call())
  F.rda <- function (x){
    Chi.z <- x$CCA$tot.chi
    q <- x$CCA$qrank
    Chi.xz <- x$CA$tot.chi
    r <- nrow(x$CA$Xbar) - x$CCA$QR$rank - 1
    F.0 <- (Chi.z/q)/(Chi.xz/r)
    F.0 <- round(F.0, 12)
    return(F.0)
  }
  Analysis <- c("glm", "rda")
  analysis <- pmatch(analysis, Analysis)
  if (length(analysis) != 1 | (is.na(analysis[1]))) {
    stop("\n Invalid analysis. Only one argument is accepted in analysis \n")
  }
  if(!is.null(CL)){
    parallel <- length(CL)
  }
  if(!is.null(parallel) & newClusters){
    CL <- parallel::makeCluster(parallel, type = "PSOCK")
  }
  if(analysis == 1){
    envir <- as.data.frame(envir)
    if(!is.null(AsFactors)){
      for(i in AsFactors){
        envir[,i] <- as.factor(envir[,i])
      }
    }
    envir.class <- SYNCSA::var.type(envir)
    RES$envir.class <- envir.class
  }
  pcps.obs <- pcps(comm, phylodist, method = method, squareroot = squareroot, correlations = FALSE)
  if(analysis == 1){
    mod.obs <- stats::glm(formula, data = data.frame(envir, pcps.obs$vectors), family = family)
    f.obs <- stats::summary.lm(mod.obs)$fstatistic[1]
    y.name <- as.character(stats::as.formula(formula)[[2]])
    vectors.obs <- pcps.obs$vectors[, y.name, drop = FALSE]
    RES$formula <- formula
  }
  if(analysis == 2){
    envir <- as.matrix(envir)
    vectors.obs <- pcps.obs$vectors[, pcps.choices, drop = FALSE]
    mod.obs <- vegan::rda(vectors.obs~envir)
    f.obs <- F.rda(mod.obs)
  }
  RES$model <- mod.obs
  RES$statistic.obs <- f.obs
  seqpermutation <- SYNCSA::permut.vector(ncol(phylodist), nset = runs)
  seqpermutation2 <- SYNCSA::permut.vector(nrow(vectors.obs), nset = runs)
  ptest <- function(i, seqperm1, seqperm2, comm, phylodist, envir, method, squareroot, analysis, vectors.obs, formula, y.name, family, pcps.choices){
    F.rda <- function (x){
      Chi.z <- x$CCA$tot.chi
      q <- x$CCA$qrank
      Chi.xz <- x$CA$tot.chi
      r <- nrow(x$CA$Xbar) - x$CCA$QR$rank - 1
      F.0 <- (Chi.z/q)/(Chi.xz/r)
      F.0 <- round(F.0, 12)
      return(F.0)
    }
    samp <- seqperm1[[i]]
    samp2 <- seqperm2[[i]]
    pcps.null <- PCPS::pcps(comm, phylodist[samp, samp], method = method, squareroot = squareroot, correlations = FALSE)
    if(analysis == 1){
      vector.null.taxa <- stats::fitted(vegan::procrustes(vectors.obs, pcps.null$vectors[, y.name, drop = FALSE], symmetric = TRUE, choices = 1))
      colnames(vector.null.taxa) <- y.name
      mod.null.taxa <- stats::glm(formula, data = data.frame(vector.null.taxa, envir), family = family)
      mod.null.site <- stats::glm(formula, data = data.frame(vectors.obs[samp2, , drop = FALSE], envir), family = family)
      F.null.taxa <- stats::summary.lm(mod.null.taxa)$fstatistic[1]
      F.null.site <- stats::summary.lm(mod.null.site)$fstatistic[1]
    }
    if(analysis == 2){
      if(length(pcps.choices)==1){
        vectors.null.taxa <- stats::fitted(vegan::procrustes(vectors.obs, pcps.null$vectors[, pcps.choices, drop = FALSE], symmetric = TRUE, choices = 1))
      } else {
        vectors.null.taxa <- stats::fitted(vegan::procrustes(vectors.obs, pcps.null$vectors[, pcps.choices, drop = FALSE], symmetric = TRUE))
      }
      mod.null.taxa <- vegan::rda(vectors.null.taxa~envir)
      mod.null.site <- vegan::rda(vectors.obs[samp2, , drop = FALSE]~envir)
      F.null.taxa <- F.rda(mod.null.taxa)
      F.null.site <- F.rda(mod.null.site)
    }
    return(cbind(F.null.taxa, F.null.site))
  }
  seqpermutation <- lapply(seq_len(nrow(seqpermutation)), function(i) seqpermutation[i,])
  seqpermutation2 <- lapply(seq_len(nrow(seqpermutation2)), function(i) seqpermutation2[i,])
  seqpermutation0 <- as.list(seq_len(runs))
  if(is.null(parallel)){
    res.F.null<-vector("list",runs)
    for (i in 1:runs) {
      if(analysis == 1){
        res.F.null[[i]] <- ptest(i = seqpermutation0[[i]], seqperm1 = seqpermutation, seqperm2 = seqpermutation2, comm = comm, phylodist = phylodist, envir = envir, method = method, squareroot = squareroot, analysis = 1, vectors.obs = vectors.obs, formula = formula, y.name = y.name, family = family)   
      }
      if(analysis == 2){
        res.F.null[[i]] <- ptest(i = seqpermutation0[[i]], seqperm1 = seqpermutation, seqperm2 = seqpermutation2,  comm = comm, phylodist = phylodist, envir = envir, method = method, squareroot = squareroot, analysis = 2, vectors.obs = vectors.obs, pcps.choices = pcps.choices)   
      }
    }
  } else {
    if(analysis == 1){
      res.F.null<-parallel::clusterApply(CL, seqpermutation0, ptest, seqperm1 = seqpermutation, seqperm2 = seqpermutation2, comm = comm, phylodist = phylodist, envir = envir, method = method, squareroot = squareroot, analysis = 1, vectors.obs = vectors.obs, formula = formula, y.name = y.name, family = family)
    }
    if(analysis == 2){
      res.F.null <- parallel::clusterApply(CL, seqpermutation0, ptest, seqperm1 = seqpermutation, seqperm2 = seqpermutation2, comm = comm, phylodist = phylodist, envir = envir, method = method, squareroot = squareroot, analysis = 2, vectors.obs = vectors.obs, pcps.choices = pcps.choices) 
    }
  } 
  F.null.taxa <- sapply(seq_len(runs), function(i) res.F.null[[i]][1,1])
  F.null.site <- sapply(seq_len(runs), function(i) res.F.null[[i]][1,2]) 
  if(!is.null(parallel) & newClusters){
    parallel::stopCluster(CL)
  }
  p.taxa.shuffle <- (sum(ifelse(F.null.taxa>=f.obs, 1, 0))+1)/(runs+1)
  p.site.shuffle <- (sum(ifelse(F.null.site>=f.obs, 1, 0))+1)/(runs+1)
  RES$p.taxa.shuffle <- p.taxa.shuffle
  RES$p.site.shuffle <- p.site.shuffle
  return(RES)
}