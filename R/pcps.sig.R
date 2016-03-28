pcps.sig<-function(comm, dist.spp, envir, analysis = c("glm", "rda"), method = "bray", squareroot = TRUE, formula, family = gaussian, AsFactors = NULL, pcps.choices=c(1,2,3,4), runs = 999){
	F.rda<-function (x){
		Chi.z <- x$CCA$tot.chi
		q <- x$CCA$qrank
		Chi.xz <- x$CA$tot.chi
		r <- nrow(x$CA$Xbar) - x$CCA$QR$rank - 1
		F.0 <- (Chi.z/q)/(Chi.xz/r)
		F.0 <- round(F.0, 12)
	return(F.0)
	}
	Analysis <- c("glm", "rda")
    analysis <- pmatch(analysis, Analysis)
    if (length(analysis) > 1) {
        stop("\n Only one argument is accepted in analysis \n")
    }
    if (is.na(analysis)) {
        stop("\n Invalid analysis \n")
    }
    if(analysis == 1){
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
	}
	pcps_obs <-pcps(comm, dist.spp, method = method, squareroot = squareroot)
	vectors<-pcps_obs$vectors
	if(analysis == 1){
		data_obs<-as.data.frame(cbind(vectors,envir))
		mod_obs<-glm(formula,data=data_obs,family=family)
		f_obs<-summary.lm(mod_obs)$fstatistic[1]
		y_name<-substr(formula,1,gregexpr("~",formula)[[1]][1]-1)
	}
	if(analysis == 2){
		vectors_obs<-pcps_obs$vectors[,pcps.choices,drop=FALSE]
		mod_obs<-rda(vectors_obs~envir)
		f_obs<-F.rda(mod_obs)
	}
	F_null_site<-matrix(NA,runs,1)
	F_null_taxa<-matrix(NA,runs,1)
	for(k in 1:runs){
		dist_null<-taxaShuffle(dist.spp)
		match.names <- match(colnames(comm), colnames(dist_null))
		pcps_null<-pcps(comm,as.matrix(dist_null[match.names, match.names]),method=method, squareroot= squareroot)
		if(analysis == 1){
			vector_null_taxa<-fitted(procrustes(vectors[,y_name], pcps_null$vectors[,y_name],symmetric = TRUE, choices=1))
			colnames(vector_null_taxa)=y_name
			data_null_taxa<-as.data.frame(cbind(vector_null_taxa,envir))
			mod_null_taxa<-glm(formula,data=data_null_taxa,family=family)
			F_null_taxa[k,]<-summary.lm(mod_null_taxa)$fstatistic[1]
			vectors_null_site<-cbind(vectors[sample(1:dim(vectors)[1]),y_name,drop=FALSE])
			data_null_site<-as.data.frame(cbind(vectors_null_site,envir))
			mod_null_site<-glm(formula,data=data_null_site,family=family)
			F_null_site[k,]<-summary.lm(mod_null_site)$fstatistic[1]
		}
		if(analysis == 2){
			vectors_null_taxa<-pcps_null$vectors[,pcps.choices,drop=FALSE]
			if(length(pcps.choices)==1){
				vector_null_taxa<-fitted(procrustes(vectors_obs,vectors_null_taxa,symmetric = TRUE,choices=1))
			}else{
				vector_null_taxa<-fitted(procrustes(vectors_obs,vectors_null_taxa,symmetric = TRUE))
			}
			mod_null_taxa<-rda(vectors_null_taxa~envir)
			F_null_taxa[k,]<-F.rda(mod_null_taxa)
			vectors_null_site<-cbind(vectors_obs[sample(1:dim(vectors_obs)[1]),,drop=FALSE])
			mod_null_site<-rda(vectors_null_site~envir)
			F_null_site[k,]<-F.rda(mod_null_site)
		}
	}
	p_taxa<-(sum(ifelse(F_null_taxa>=f_obs,1,0))+1)/(runs+1)
	p_site<-(sum(ifelse(F_null_site>=f_obs,1,0))+1)/(runs+1)
	if(analysis == 1){	
		res<-list(model=mod_obs, Envir_class=envir_class, formula=formula, statistic.obs=f_obs, p.site.shuffle = p_site, p.taxa.shuffle = p_taxa)
	}
	if(analysis == 2){
		res<-list(model=mod_obs, statistic.obs=f_obs, p.site.shuffle = p_site, p.taxa.shuffle = p_taxa)
	}
return(res)
}