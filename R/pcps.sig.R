pcps.sig<-function(comm, dist.spp, envir, method = "bray", squareroot = TRUE, formula,runs=999,AsFactors=NULL){
	dis<-dist.spp
	envir<-as.data.frame(envir)
	if(!is.null(AsFactors)){
		for(i in AsFactors){
			envir[,i]<-as.factor(envir[,i])
		}
	}
	envir_class<-matrix(NA,dim(envir)[2],1)
	rownames(envir_class)=colnames(envir)
	colnames(envir_class)=c("Class")
	for(j in 1:dim(envir)[2]){
		envir_class[j,1]<-class(envir[,j])	
	}
	ord<-pcps(comm,dis, method = method, squareroot = squareroot)
	vectors<-ord$vectors
	colnames(vectors)=NULL
	colnames(vectors)=colnames(vectors,do.NULL=FALSE,prefix="pcps.")
	data_obs<-as.data.frame(cbind(vectors,envir))
	mod_obs<-glm(formula,data=data_obs)
	f_obs<-summary.lm(mod_obs)$fstatistic[1]
	y_name<-substr(formula,1,gregexpr("~",formula)[[1]][1]-1)
	res_F_null<-matrix(NA,runs,1)
		for(k in 1:runs){
			dist_null<-taxaShuffle(dis)
			match.names <- match(colnames(comm), colnames(dist_null))
			m_p_null<-matrix.p(comm,as.matrix(dist_null[match.names, match.names]))$matrix.P
			dist_p_null <- vegdist(m_p_null, method = method)
		    if (squareroot == TRUE) {
    		    dist_p_null <- sqrt(dist_p_null)
    		}
			ord_null<-pcoa(dist_p_null)
			vectors_null<-ord_null$vectors
			colnames(vectors_null)=NULL
			colnames(vectors_null)=colnames(vectors_null,do.NULL=FALSE,prefix="pcps.")
			res_pro<-procrustes(vectors[,y_name],vectors_null[,y_name],symmetric = TRUE, choices=1)
			vector_null<-fitted(res_pro)
			colnames(vector_null)=y_name
			data_null<-as.data.frame(cbind(vector_null,envir))
			mod_null<-glm(formula,data=data_null)
			res_F_null[k,]<-summary.lm(mod_null)$fstatistic[1]
		}
	p<-(sum(ifelse(res_F_null>=f_obs,1,0))+1)/(runs+1)
return(list(model=mod_obs, Envir_class=envir_class,formula=formula, f.obs=f_obs,p=p))
}