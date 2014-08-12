pcps.curve<-function(comm, dist.spp,trait,method = "bray", squareroot = TRUE,null.model=TRUE,runs=99,progressbar=FALSE){
	dis<-dist.spp
	m_t_obs<-matrix.t(comm,trait,scale=FALSE,notification=FALSE)$matrix.T
	ord<-pcps(comm,dis, method = method, squareroot = squareroot)
	values<-ord$values
	vectors<-ord$vectors
	calc.pcpc.curve<-function(values,vectors,matrixT){
		use<-1:(dim(vectors)[2])
		x<-vectors[,use]
		y<-matrixT[,1]
		fac<-length(use)
		xnam <- paste("x[,", 1:fac,"]", sep="")
		res.y<-matrix(NA,nrow=fac,ncol=1)
		for (j in 1:fac){
			res.y[j,1]<-as.numeric(summary(lm(as.formula(paste("y ~ ", paste(xnam[1:j], collapse= "+")))))$r.squared)
		}	
		colnames(res.y)="Coefficient_of_determination"
		res.x<-as.matrix(values[1:fac,3])
		colnames(res.x)="Cumulative_PCPS_eigenvalues"
		result<-cbind(res.x,res.y)
	return(result)
	}
	curve_obs<-calc.pcpc.curve(values,vectors,m_t_obs)
	if(null.model){
		res_curve_null<-vector("list",runs)
		for(k in 1:runs){
			dist_null<-taxaShuffle(dis)
			match.names <- match(colnames(comm), colnames(dist_null))
			m_p_null<-matrix.p(comm,as.matrix(dist_null[match.names, match.names]))$matrix.P
			dist_p_null <- vegdist(m_p_null, method = method)
		    if (squareroot == TRUE) {
    		    dist_p_null <- sqrt(dist_p_null)
    		}
			ord_null<-pcoa(dist_p_null)
			values_null<-ord_null$values[,c(1,2,4)]
			vectors_null<-ord_null$vectors
			res_curve_null[[k]]<-calc.pcpc.curve(values_null,vectors_null,m_t_obs)
			if(progressbar){
				ProgressBAR(k,runs,style=3)
			}
		}
	}
	ReTuRn<-list(curve_obs=curve_obs)
	if(null.model){
		ReTuRn<-list(curve.obs=curve_obs,curve.null=res_curve_null)	
	}
	class(ReTuRn) <- "pcpscurve"
	return(ReTuRn)	
}