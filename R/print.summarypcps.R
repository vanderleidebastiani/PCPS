#' @rdname pcps
#' @encoding UTF-8
#' @export
print.summarypcps<-function(x, ...){
  cat("Call:\n")
  cat(deparse(x$call), "\n\n")
  cat("PCPS values:\n")
  print(as.matrix(x$values), ...)
  cat("\n PCPS vectors:\n")
  print(as.matrix(x$vectors), ...)
  cat("\n Correlations:\n")
  print(as.matrix(x$correlations), ...)
  invisible(x)
}