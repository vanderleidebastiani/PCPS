plot.pcps<-function(x,display=c("text","points"),groups,showlabel=TRUE,...){
	sco<-scores(x,...)
	plot(sco$scores.sites,type="n",ylim=c(min(sco$scores.sites[,2],sco$scores.species[,2],na.rm=TRUE)-0.05, max(sco$scores.sites[,2],sco$scores.species[,2],na.rm=TRUE)+0.05),xlim=c(min(sco$scores.sites[,1],sco$scores.species[,1],na.rm=TRUE)-0.05,max(sco$scores.sites[,1],sco$sco[,1],na.rm=TRUE)+0.05),...)
	if(display=="text"){
		text(sco$scores.sites,labels=rownames(sco$scores.sites),...) 
	}
	if(display=="points"){
		points(sco$scores.sites,...)
	}
	ordispider(sco$scores.species,groups=groups,label=showlabel,...) 
}