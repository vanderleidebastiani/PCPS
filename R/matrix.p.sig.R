#' @rdname pcps.sig
#' @encoding UTF-8
#' @export
matrix.p.sig<-function(comm, phylodist, envir, analysis = c("adonis", "mantel"), method = "bray", squareroot = TRUE, method.envir = "euclidean", runs = 999, parallel = NULL, newClusters = TRUE, CL =  NULL){
  RES <- list(call = match.call())
  Analysis <- c("adonis", "mantel")
  analysis <- pmatch(analysis, Analysis)
  if (length(analysis) != 1 | (is.na(analysis[1]))) {
    stop("\n Invalid analysis. Only one argument is accepted in analysis \n")
  }
  p.matrix <- SYNCSA::matrix.p(comm, phylodist, notification = FALSE)$matrix.P
  p.dist <- vegan::vegdist(p.matrix, method = method)
  if (squareroot) {
    p.dist <- sqrt(p.dist)
  }
  if(!is.null(CL)){
    parallel <- length(CL)
    if(analysis == 1){
      mod.obs <- vegan::adonis(p.dist~envir,permutations=runs, parallel = CL)
      statistic.obs <- mod.obs$aov.tab$F.Model[1]
      p.site.shuffle <- mod.obs$aov.tab$"Pr(>F)"[1]
    }
    if(analysis == 2){
      env.dist <- vegan::vegdist(envir, method = method.envir)
      mod.obs <- vegan::mantel(p.dist, env.dist, permutations = runs, parallel = CL)
      statistic.obs <- mod.obs$statistic
      p.site.shuffle <- mod.obs$signif
    }
  } else {
    if(is.null(parallel)){
      if(analysis == 1){
        mod.obs <- vegan::adonis(p.dist~envir, permutations = runs)
        statistic.obs <- mod.obs$aov.tab$F.Model[1]
        p.site.shuffle <- mod.obs$aov.tab$"Pr(>F)"[1]
      }
      if(analysis == 2){
        env.dist <- vegan::vegdist(envir, method = method.envir)
        mod.obs <- vegan::mantel(p.dist, env.dist, permutations = runs)
        statistic.obs <- mod.obs$statistic
        p.site.shuffle <- mod.obs$signif
      }
    } else {
      if(analysis == 1){
        mod.obs <- vegan::adonis(p.dist~envir,permutations=runs, parallel = parallel)
        statistic.obs <- mod.obs$aov.tab$F.Model[1]
        p.site.shuffle <- mod.obs$aov.tab$"Pr(>F)"[1]
      }
      if(analysis == 2){
        env.dist <- vegan::vegdist(envir, method = method.envir)
        mod.obs <- vegan::mantel(p.dist, env.dist, permutations = runs, parallel = parallel)
        statistic.obs <- mod.obs$statistic
        p.site.shuffle <- mod.obs$signif
      }
    }
  }
  RES$model <- mod.obs
  RES$statistic.obs <- statistic.obs
  RES$p.site.shuffle <- p.site.shuffle
  seqpermutation <- SYNCSA::permut.vector(ncol(phylodist), nset = runs)
  ptest <- function(samp, comm, phylodist, envir, method, squareroot, analysis){
    MP.null <- SYNCSA::matrix.p(comm, phylodist[samp, samp], notification = FALSE)$matrix.P
    dist.MP.null <- vegan::vegdist(MP.null, method = method)
    if(squareroot){
      dist.MP.null <- sqrt(dist.MP.null)
    }
    if(analysis == 1){
      mod.null <- vegan::adonis(dist.MP.null~envir, permutations = 0)
      res <- mod.null$aov.tab$F.Model[1]
    } 
    if(analysis == 2){
      mod.null <- vegan::mantel(dist.MP.null, envir, permutations = 0)
      res <- mod.null$statistic
    }
    return(res)
  }
  if(!is.null(parallel) & newClusters){
    CL <- parallel::makeCluster(parallel, type = "PSOCK")
  }
  if(is.null(parallel)){
    res.null <- matrix(NA,runs,1)
    for (i in 1:runs) {
      if(analysis == 1){
        res.null[i,1] <- ptest(samp = seqpermutation[i,], comm = comm, phylodist = phylodist, envir = envir, method = method, squareroot = squareroot, analysis = analysis)   
      }
      if(analysis == 2){
        res.null[i,1] <- ptest(samp = seqpermutation[i,], comm = comm, phylodist = phylodist, envir = env.dist, method = method, squareroot = squareroot, analysis = analysis)
      }
    }
  } else {
    if(analysis == 1){
      res.null <- parallel::parRapply(CL, seqpermutation, ptest, comm = comm , phylodist = phylodist, envir = envir, method = method, squareroot = squareroot, analysis = analysis)
    }
    if(analysis == 2){
      res.null <- parallel::parRapply(CL, seqpermutation, ptest, comm = comm , phylodist = phylodist, envir = env.dist, method = method, squareroot = squareroot, analysis = analysis)
    }
  }
  if(!is.null(parallel) & newClusters){
    parallel::stopCluster(CL)
  }
  p.taxa.shuffle <- (sum(ifelse(res.null>=statistic.obs, 1, 0))+1)/(runs+1)
  RES$p.taxa.shuffle <- p.taxa.shuffle
  return(RES)
}