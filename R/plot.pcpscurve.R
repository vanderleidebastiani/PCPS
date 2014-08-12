plot.pcpscurve<-function(x,type="b",errorbars=c("none","sd","se","quantile"),probs=c(0.025,0.975),col="black",errbar.col="black",errbar.pch=16,...){
	if (length(errorbars) > 1) {
		stop("\n Only one argument is accepted in errorbars \n")
	}
	plot(-1,-1,xlim=c(0,1),ylim=c(0,1),xlab="Cumulative PCPS eigenvalues (%)", ylab="Coefficient of determination (R2)")
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
		if (errorbars=="sd") {
			error_X<-apply(X,2,sd,na.rm=T)
			error_Y<-apply(Y,2,sd,na.rm=T)
			plotCI(mean_X,mean_Y,error_Y,error_Y,err="y",add=TRUE,col= errbar.col,pch=errbar.pch)
			plotCI(mean_X,mean_Y,error_X,error_X,err="x",add=TRUE,col= errbar.col,pch=errbar.pch)
		}
		if(errorbars=="se"){
			error_X<-apply(X,2,sd,na.rm=T)
			error_Y<-apply(Y,2,sd,na.rm=T)
			error_X<-error_X/sqrt(N)
			error_Y<-error_Y/sqrt(N)	
			plotCI(mean_X,mean_Y,error_Y,error_Y,err="y",add=TRUE,col= errbar.col,pch=errbar.pch)
			plotCI(mean_X,mean_Y,error_X,error_X,err="x",add=TRUE,col= errbar.col,pch=errbar.pch)
		}
		if (errorbars=="quantile") {
			if(length(probs)!=2){
				stop("\n Only two values are accepted in probs \n")
			}
			error_X<-apply(X,2,quantile,na.rm=T,probs=probs)
			error_Y<-apply(Y,2,quantile,na.rm=T,probs=probs)
			plotCI(mean_X,mean_Y,li=error_Y[1,],ui=error_Y[2,],err="y",add=TRUE,col= errbar.col,pch=errbar.pch)
			plotCI(mean_X,mean_Y,li=error_X[1,],ui=error_X[2,],err="x",add=TRUE,col= errbar.col,pch=errbar.pch)
		}
		if(errorbars=="none"){
			points(mean_X,mean_Y,col= errbar.col,pch=errbar.pch)
		}
	}
	segments(0,0,1,1,lty=2)
	points(x$curve.obs[,1],x$curve.obs[,2],type=type,col=col,...)
}