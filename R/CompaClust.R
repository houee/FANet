CompaClust = function(colorDynamicHybridTOM1,colorDynamicHybridTOM2, MEs1, MEs2, names1, names2, intitle1, intitle2, plot=T){
    Modules1 = substring(names( MEs1), 3)
    Modules2 = substring(names(MEs2), 3)
    nMods1 =length(Modules1)
    nMods2 =length(Modules2)
    pTable =matrix(0,nrow= nMods1, ncol= nMods2);
    CountTbl =matrix(0,nrow= nMods1, ncol= nMods2);
    for(fmod in 1:nMods1){
        for(cmod in 1:nMods2){
            Members1 = (colorDynamicHybridTOM1 == Modules1[fmod]);
            Members2 = (colorDynamicHybridTOM2 == Modules2[cmod]);
            pTable[fmod, cmod] = -log10(fisher.test(Members1, Members2, alternative ="greater")$p.value);
            CountTbl[fmod, cmod] =sum(colorDynamicHybridTOM1 == Modules1[fmod]&colorDynamicHybridTOM2 ==Modules2[cmod])
        }
    }
    RI=randIndex(CountTbl ,adjust=F)
    pTable[is.infinite(pTable)] = 1.3*max(pTable[is.finite(pTable)]);
    pTable[pTable>50 ] = 50 ;
    ModTotals1 =apply(CountTbl, 1,sum)
    ModTotals2 =apply(CountTbl, 2,sum)
    if (plot==T){
        sizeGrWindow(10,7 );
        par(mfrow=c(1,1));
        par(cex = 1.0);
        par(mar=c(8, 10.4, 2.7, 1)+0.3);
        labeledHeatmap(Matrix = pTable,
        xLabels =paste(" ", Modules2),
        yLabels =paste(" ", Modules1),
        colorLabels = TRUE, xSymbols =paste(names2, " ", Modules2,": ",ModTotals2, sep=""), ySymbols =paste(names1," ", Modules1,": ", ModTotals1, sep=""),
        textMatrix = CountTbl, colors= greenWhiteRed(100)[50:100], main =paste("Clustering obtained", intitle2, "or", intitle1," RandIndex=", round(RI,2)),cex.text= 1.2, cex.lab = 0.8, setStdMargins = FALSE)
    }
    return(list(contingency=CountTbl,pvalueTable=pTable,RandIndex=RI))
}
