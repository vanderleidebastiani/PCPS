#' @rdname pcoa.sig
#' @encoding UTF-8
#' @export
print.summarypcoasig<-function(x, ...){
  cat("Call:\n")
  cat(deparse(x$call), "\n\n")
  cat("Values:\n")
  print(as.matrix(x$values), ...)
  cat("\n Vectors:\n")
  print(as.matrix(x$vectors), ...)
  cat("\n Correlations:\n")
  print(as.matrix(x$correlations), ...)
  cat("\n Probabilities:\n")
  print(as.matrix(x$probabilities), ...)
  invisible(x)
}