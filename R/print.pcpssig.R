#' @rdname pcps.sig
#' @encoding UTF-8
#' @export
print.pcpssig <- function(x, ...){
  cat("$call:\n")
  cat(deparse(x$call), "\n\n")
  cat("$model:\n")
  print(x$model, ...)
  cat("\n$obs.statistic:\n")
  print(x$obs.statistic, ...)
  cat("\n$p.site.shuffle:\n")
  print(x$p.site.shuffle, ...)
  cat("\n$p.taxa.shuffle:\n")
  print(x$p.taxa.shuffle, ...)
  invisible(x)
}
