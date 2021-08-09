library(Seurat)
library(stringr)
library(plyr)
library(dplyr)
library(pracma)
library(Hmisc)
options(future.globals.maxSize = 80 * 1024^3)


# Process Seurat object using SCTransform
# Input parameters:
# 1. Seurat object
# 2. array of low quality (LQ) clusters to be excluded
ProcessSCT <- function(seur, LQ=NULL, clusterRes=0.5, regressVars=NULL) {

# remove LQ clusters
if (!is.null(LQ)) {
seur <- subset(seur, idents=LQ, invert=T)
}

if (!is.null(regressVars)) {
regressVars = c("percent.mt", "Replicate")
}

# Process using SCTransform
seur <- seur %>%
        SCTransform(vars.to.regress=regressVars) %>%
        RunPCA(verbose=F) %>%
        RunTSNE(assay='SCT', reduction='pca', dims=1:30) %>%
        FindNeighbors(dims = 1:30) %>%
        FindClusters(resolution=clusterRes)

saveRDS(seur, "seur_sct_clustered.rds")

return(seur)
}


# Used to process Fezf2 Control and KO scRNA dataset
# This function additionally removes mitochondrial and ribosomal genes from the variable genes
# to address batch effect
ProcessSCT2 <- function(seur, LQ=NULL, clusterRes=0.5, regressVars=NULL) {

# remove LQ clusters
if (!is.null(LQ)) {
seur <- subset(seur, idents=LQ, invert=T)
}

if (!is.null(regressVars)) {
regressVars = c("percent.mt", "Replicate")
}

# Process using SCTransform
seur <- SCTransform(seur, vars.to.regress=regressVars)
vg = VariableFeatures(seur)
VariableFeatures(seur) = setdiff(vg, grep("^mt-|^Rp[sl]", vg, value=T))

seur <- seur %>%
        RunPCA(verbose=F) %>%
        RunTSNE(assay='SCT', reduction='pca', dims=1:30) %>%
        FindNeighbors(dims = 1:30) %>%
        FindClusters(resolution=clusterRes)

saveRDS(seur, "seur_sct2_clustered.rds")

return(seur)
}


