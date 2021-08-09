# Replicate the analysis pipeline presented in the original MerFISH paper
# https://science.sciencemag.org/content/sci/suppl/2018/10/31/science.aau5324.DC1/aau5324-Moffitt-SM.pdf


library(Seurat)
library(dplyr)
source("MerFISH_var_genes.R")


# Input parameters:
# 1. Seurat object
# 2. Array of variable genes (from "MerFISH_var_genes.R")
# 3. Boolean - whether to perform normalization
# 4. List of clusters to remove from the analysis
# 5. Clustering resolution
ProcessMerFISH <- function(seur, VG=var_genes, norm=T, LQ=NULL, nPC=NULL, clusterRes=0.5) {

# Normalize without log-transformation then compute z-scores
if (norm) {
seur <- NormalizeData(seur, normalization.method='RC', scale.factor=1)
}

# remove LQ clusters
if (!is.null(LQ)) {
seur <- subset(seur, idents=LQ, invert=T)
}

# Compute Z score
zscore=t(apply(seur@assays$RNA@data, 1, function(X){
return((X-mean(X))/sd(X))
}))
seur@assays$RNA@scale.data = zscore

# Specify variable genes to use
VariableFeatures(seur) <- VG  # from "MerFISH_var_genes.R"

# PCA
seur <- RunPCA(seur, npcs=10)

if (is.null(nPC)) {
# Determine number of PCs to use
seur <- JackStraw(seur, prop.freq=.5)
seur <- ScoreJackStraw(seur, dims = 1:10)
nPC=max(seur@reductions$pca@jackstraw$overall.p.values[,"PC"][seur@reductions$pca@jackstraw$overall.p.values[,'Score']<0.05])
print(paste0("using ", nPC, " PCs (JackStraw pval < 0.05)"))
}

# downstream
seur <- FindNeighbors(seur, dims=1:nPC, k.param=10)  # default k is 20, in the MerFISH paper, they used 10-12
seur <- FindClusters(seur, resolution=0.5)
seur <- RunTSNE(seur, dims=1:nPC, check_duplicates=F)
saveRDS(seur, paste0("seur_processed_custom_varGenes_",nPC,"pcs.rds"))

}
