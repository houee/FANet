BioContribution = function(dta, Sbio, zeros, beta=6, cutHeight = 0.998){
    S = cor(t(dta))
    ADJ = abs(S)^beta
    dissTOM=TOMdist(ADJ)
    hierTOM = hclust(as.dist(dissTOM),method="average");
    classlabels = dynamicTreeCut::cutreeDynamic(hierTOM, distM= dissTOM , cutHeight = cutHeight,deepSplit=2, pamRespectsDendro = FALSE)
    colorDynamicHybridTOM = labels2colors(classlabels)
    
    MEList = moduleEigengenes(t(dta),colors= colorDynamicHybridTOM)
    MEs = MEList$eigengene
    ADJ.bio = abs(Sbio)^beta
    dissTOM.bio=TOMdist(ADJ.bio)
    hierTOM.bio = hclust(as.dist(dissTOM.bio),method="average");
    classlabels.bio = dynamicTreeCut::cutreeDynamic(hierTOM.bio, distM= dissTOM.bio , cutHeight =cutHeight,deepSplit=2, pamRespectsDendro = FALSE)
    colorDynamicHybridTOM.bio = labels2colors(classlabels.bio)
    
    MEList.bio = moduleEigengenes(t(dta),colors= colorDynamicHybridTOM.bio)
    MEs.bio = MEList.bio$eigengene
    SGO=cor(t(zeros))
    ADJ.GO = abs(SGO)^beta
    dissTOM.GO=TOMdist(ADJ.GO)
    hierTOM.GO = hclust(as.dist(dissTOM.GO),method="average");
    classlabels.GO = dynamicTreeCut::cutreeDynamic(hierTOM.GO, distM= dissTOM.GO , cutHeight = cutHeight,deepSplit=2, pamRespectsDendro = FALSE)
    colorDynamicHybridTOM.GO = labels2colors(classlabels.GO)
    
    MEList.GO = moduleEigengenes(t(dta),colors= colorDynamicHybridTOM.GO)
    MEs.GO = MEList.GO$eigengene
    
    dev.new()
    sizeGrWindow(11,6 );
    
    res.GO.bio = CompaClust(colorDynamicHybridTOM1=colorDynamicHybridTOM.GO, colorDynamicHybridTOM2=colorDynamicHybridTOM.bio, MEs1=MEs.GO, MEs2=MEs.bio, names1="ONLY BioInfo", names2="Prior BioInfo", intitle1="using BioInfo only", intitle2="using prior BioInfo", plot=T)
    
    res.GO.without = CompaClust(colorDynamicHybridTOM1=colorDynamicHybridTOM.GO, colorDynamicHybridTOM2=colorDynamicHybridTOM, MEs1=MEs.GO, MEs2=MEs, names1="ONLY BioInfo", names2="Without BioInfo", intitle1="using BioInfo only", intitle2="without using BioInfo", plot=F)
    
    dev.new()
    sizeGrWindow(11,6 );
    
    res.bio.without = CompaClust(colorDynamicHybridTOM1=colorDynamicHybridTOM, colorDynamicHybridTOM2=colorDynamicHybridTOM.bio, MEs1=MEs, MEs2=MEs.bio, names1="Without BioInfo", names2="Prior BioInfo", intitle1="without using BioInfo", intitle2="using prior BioInfo")
    BioContrib.Score = res.GO.bio$RandIndex/res.bio.without$RandIndex
    
    return(list(BioContrib.Score=BioContrib.Score, Clust.OnlyBioinfo.PriorBioInfo=res.GO.bio, Clust.WithoutBioinfo.PriorBioInfo=res.bio.without))
}
