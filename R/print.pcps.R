print.pcps<-function(x , ...){
	cat("Call:\n")
	cat(deparse(x$call), "\n\n")
	cat("PCPS values:\n")
	print(as.matrix(x$values))
	invisible(x)
}