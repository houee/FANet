MatBioPrep = function(matrixBio, filtersize=TRUE, nbgenes.min=50, nbgenes.max=200, optimal.nbclasses=TRUE,  nbclasses.min=5, nbclasses.max=20, nbclasses=NULL){
    if(filtersize==TRUE){
        matGO.temp=matrixBio
        sum_element_colonne=apply(matGO.temp,2,FUN=sum)
        i=max(nbgenes.min, min(sum_element_colonne))
        j=min(nbgenes.max, max(sum_element_colonne))
        ind=which((sum_element_colonne>=i)&(sum_element_colonne<=j))
        matGO.selec=matGO.temp[,ind]
    }
    if(filtersize==FALSE){
        matGO.selec=matrixBio
    }
    tree= hclust(dist(t(matGO.selec)))
    if (optimal.nbclasses==TRUE){
        bestcut=best.cutree(tree, min=nbclasses.min,max=nbclasses.max)
        cut=cutree(hclust(dist(t(matGO.selec))),k=bestcut)
    }
    if (optimal.nbclasses==FALSE){
        cut=cutree(hclust(dist(t(matGO.selec))),k=nbclasses)
        bestcut="Function Bestcut has not been called"
    }
    biosynt=matrix(0,dim(matGO.selec)[1],length(levels(as.factor(cut))))
    for (i in 1:length(levels(as.factor(cut)))){
        if (ncol(as.data.frame(matGO.selec[,cut==i]))==1) {
            biosynt[,i]=matGO.selec[,cut==i]}
        if (ncol(as.data.frame(matGO.selec[,cut==i]))!=1){
            biosynt[,i]=as.numeric(apply(matGO.selec[,cut==i],1,sum)!=0)
        }
    }
    zeros=biosynt
    rownames(zeros)=rownames(matrixBio)
    GOClassMembership =as.data.frame(cor(matGO.selec, zeros,  use ="p"));
    indice.max=function(x){
        ind=which(x==max(x))
        if(length(ind)>1){
            ind=ind[1]}
        return(ind)
    }
    maxvalue=apply(GOClassMembership, 2, indice.max)
    GOclustcloser=rownames(GOClassMembership)[maxvalue]
    return(list(zeros=zeros, tree = tree, bestcut = bestcut, MaxClassMembership=GOclustcloser))
}
