#' @rdname pcps.sig
#' @include pcps.sig.R
#' @include matrix.p.sig.R
#' @encoding UTF-8
#' @export
FUN.ADONIS2.global <- function(x, envir, method.p, formula, sqrt.p = TRUE, return.model = FALSE){
  p.dist <- vegan::vegdist(x, method = method.p)
  if(sqrt.p){
    p.dist <- sqrt(p.dist)
  }
  assign("p.dist", p.dist, envir = globalenv())
  mod.obs <- vegan::adonis2(formula, data = data.frame(envir), permutations = 0, by = NULL, parallel = NULL)
  rm(p.dist, envir = globalenv())
  statistic.obs <- mod.obs$F[1]
  if(return.model){
    res <- list()
    res$mod.obs <- mod.obs
    res$statistic.obs <- statistic.obs
  } else{
    res <- statistic.obs
  }
  return(res)
}