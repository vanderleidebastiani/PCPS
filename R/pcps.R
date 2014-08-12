pcps<-function(comm,dist.spp,method="bray",squareroot=TRUE){
	dis<-dist.spp
	P<-matrix.p(comm,dis,notification=FALSE)$matrix.P
	P.dist<-vegdist(P,method=method)
	if(squareroot==TRUE){
		P.dist<-sqrt(P.dist)
	}
	ordi.P<-wcmdscale(P.dist,eig=TRUE)
	vectors<-ordi.P$points
	colnames(vectors)<-NULL
	colnames(vectors)<-colnames(vectors,do.NULL=FALSE,prefix="pcps.")
	values<-ordi.P$eig[which((ordi.P$eig>=0)==TRUE)]
	if(length(unique(ordi.P$eig<0))>1){
		warning("Warning: Negative eigenvalues are present in the decomposition result, but only positive eigenvalues were considered",call.=FALSE)
	}
	relative<-values/sum(values)
	cumulative<-as.vector(rep(NA,length(values)))
	for (i in 1:length(values)){
		cumulative[i]<-sum((values/sum(values))[1:i])
	}
	Values<-cbind(values,relative,cumulative)
	colnames(Values)=c("Eigenvalues","Relative_eig","Cumul_eig")
	rownames(Values)=1:length(values)
	n.col<-dim(P)[2]
	n.axis<-dim(vectors)[2]
	correlations<-matrix(NA,nrow=n.col,ncol=n.axis)
	for (i in 1:n.col){
		for (j in 1:n.axis){
			correlations[i,j]<-cor(P[,i],vectors[,j])
		}
	}		
	colnames(correlations)<-colnames(vectors)
	rownames(correlations)<-colnames(P,do.NULL=FALSE,prefix="Uni.")
	Res<-list(call= match.call(), P=P, values=Values, vectors=vectors, correlations=correlations)
	class(Res) <- "pcps"
	return(Res)
}