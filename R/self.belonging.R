self.belonging<-function (dis,standardize=TRUE){
	diag.matrix<-diag(belonging(dis,standardize=standardize))
	return(diag.matrix)
}