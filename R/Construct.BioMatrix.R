Construct.BioMatrix = function(IDhgnc, GO.classe="biological_process"){
    mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
    hgnc <- as.character(IDhgnc[,2])
    results <- getBM(attributes = c("hgnc_symbol", "go_id","name_1006","definition_1006","namespace_1003"), filters = "hgnc_symbol", values = hgnc, mart = mart)
    data.GO=results[,]
    indiceGO=which(results[,"namespace_1003"]==GO.classe)
    gene.symbol=results[indiceGO,"hgnc_symbol"]
    GOterm=results[indiceGO,"go_id"]
    tableGeneGO=cbind.data.frame(gene.symbol,GOterm)
    matrixGO=table2binmatrix(tableGeneGO)
    return(matrixGO)
}
