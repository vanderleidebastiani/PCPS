#' @rdname pcps
#' @encoding UTF-8
#' @export
print.summarypcps<-function(x, ...){
  cat("Call:\n")
  cat(deparse(x$call), "\n\n")
  cat("\n$values:\n")
  print(as.matrix(x$values), ...)
  cat("\n$vectors:\n")
  print(as.matrix(x$vectors), ...)
  cat("\n$correlations:\n")
  print(as.matrix(x$correlations), ...)
  invisible(x)
}