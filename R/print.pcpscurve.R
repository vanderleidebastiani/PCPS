print.pcpscurve<-function(x, ...){
	cat("Call:\n")
	cat(deparse(x$call), "\n\n")
	cat("PCPS curve observed:\n")
	res<-as.matrix(x$curve.obs)
	print(res)
	if(!is.null(x$curve.null)){
		N<-length(x$curve.null)
		X<-matrix(NA,N,dim(x$curve.obs)[1])
		Y<-X
		for(i in 1:N){
			X[i,]<-x$curve.null[[i]][,1]
			Y[i,]<-x$curve.null[[i]][,2]
		}
		mean_X<-apply(X,2,mean,na.rm=T)
		mean_Y<-apply(Y,2,mean,na.rm=T)
		resN<-cbind(mean_X,mean_Y)
		rownames(resN)=rownames(res)
		colnames(resN)=c("Cumulative_PCPS_eigenvalues","Coefficient_of_determination")
		cat("\n")
		cat("Mean PCPS curve null:\n")
		print(as.matrix(resN))
	}
	invisible(x)
}