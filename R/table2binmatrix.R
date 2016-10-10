table2binmatrix = function(tableGeneGO){
    gene.symbol=tableGeneGO[,1]
    GOterm=tableGeneGO[,2]
    transf_pre=function(x){
        logic=match(GOterm,x,nomatch=0)
        datatemp=cbind.data.frame(gene.symbol,logic)
        res=table(datatemp)[,"1"]
        return(res)}
    data=cbind.data.frame(as.character(gene.symbol),as.character(GOterm))
    GOtermunique=unique(GOterm)
    datares=apply(as.data.frame(GOtermunique),1,transf_pre)
    colnames(datares)=GOtermunique
    return(datares)
}
