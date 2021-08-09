library(Seurat)


MyModuleScore <- function(obj, gene.list, save=F, filename=NULL){
n1 = ncol(obj@meta.data) + 1
obj <- AddModuleScore(object = obj, features = gene.list, name = names(gene.list))
n2 = ncol(obj@meta.data)
meta = obj@meta.data[,n1:n2]
colnames(meta) <- gsub("[0-9]+$", "", colnames(meta))
if (save) {
saveRDS(meta, filename)
}
for(x in paste0(names(gene.list), seq(1:length(gene.list)))) {
obj[[x]] <- NULL
}
obj <- AddMetaData(obj, meta)
return(obj)
}


