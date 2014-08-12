define.clade<-function(tree,threshold,time,method=c("threshold","time")){
	if(is.null(tree$node.label)){
		stop("\n Node labels not found. Use the function makeNodeLabel\n")
	}
	if(tree$Nnode!=length(tree$node.label)){
		stop("\n Some nodes not found among all nodes in tree\n")
	}
	tree1<-phylo4(tree)
	NoDes<-node.depth.edgelength(tree) 
	names(NoDes)=c(tree$tip.label,tree$node.label)
	NoDes<-NoDes[-match(tree$tip.label,names(NoDes))]
	N1<-length(NoDes)
	if(method=="threshold"){
		NoDes2<-NoDes[NoDes/node.depth.edgelength(tree)[1]>=(threshold)]
		height<-node.depth.edgelength(tree)[1]*(threshold)
	}
	else{
		NoDes2<-NoDes[NoDes>=time]
		height<-time
	}
	clades<-vector(length=length(tree$tip.label))	
	names(clades)=tree$tip.label
	clades[]=tree$tip.label	
	N2<-length(NoDes2)
	if(N1==N2){
		N2=N2-1
	}
	n=1
	if(!N2==0){
		for (n in 1:N2){
			RM=n
			Descendants<-descendants(tree1,names(sort(NoDes,decreasing=T))[RM])
			NoDe<-names(sort(NoDes,decreasing=T))[RM]
			clades[Descendants]=NoDe
		}
	}
return(list(clades=clades,height=height))
}