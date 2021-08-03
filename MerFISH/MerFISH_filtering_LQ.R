# Replicate the analysis pipeline presented in the original MerFISH paper
# https://science.sciencemag.org/content/sci/suppl/2018/10/31/science.aau5324.DC1/aau5324-Moffitt-SM.pdf


library(Seurat)
library(dplyr)
source("cell_type_specific_genes.R")
source("~/kwanho/src/seurat_tools.R")

seur <- readRDS("seur_processed_custom_varGenes_8pcs.rds")

print("Remove LQ cells!")
seur <- subset(seur, idents=9, invert=T)
print(table(Idents(seur)))

print("Compute Z scores")
zscore=t(apply(seur@assays$RNA@data, 1, function(X){
return((X-mean(X))/sd(X))
}))
seur@assays$RNA@scale.data = zscore

# Use cell type-specific genes from the gene probe list as the variable genes
# (all genes except L-R targets, and Ctgf which is in the Excel Jeff shared, but not in the codebook)
print("Selected variable genes:")
VariableFeatures(seur) <- var_genes

print("Run PCA")
seur <- RunPCA(seur, npcs=10)

print("Select nPCs")
seur <- JackStraw(seur, prop.freq=.5)
seur <- ScoreJackStraw(seur, dims = 1:10)

pdf("figures/pc_noLQ_choose_nPCs.pdf")
JackStrawPlot(seur, dims = 1:10)
ElbowPlot(seur)
dev.off()

pdf("figures/pc_noLQ_dim_loadings.pdf", height=14)
VizDimLoadings(seur, dims=1:10, reduction='pca')
dev.off()

pdf("figures/pc_noLQ_dimplots.pdf")
Idents(seur) <- 'sample'
for (i in 1:9) {
print(DimPlot(seur, reduction="pca", dims=c(i,i+1)))
}
dev.off()

npcs=max(seur@reductions$pca@jackstraw$overall.p.values[,"PC"][seur@reductions$pca@jackstraw$overall.p.values[,'Score']<0.05])
print(paste0("using ", npcs, " PCs (JackStraw pval < 0.05)"))


outf = paste0("seur_noLQ_processed_custom_varGenes_", npcs, "pcs.rds")

seur <- FindNeighbors(seur, dims=1:npcs, k.param=10)  # default k is 20, in the MerFISH paper, they used 10-12	
seur <- FindClusters(seur)
seur <- RunTSNE(seur, dims=1:npcs, check_duplicates=F)
saveRDS(seur, outf)

pdf(paste0("figures/tsne_noLQ_custom_", npcs, "pcs.pdf"))
DimPlot(seur, reduction='tsne', label=T)
DimPlot(seur, reduction='tsne', group.by='sample')
DimPlot(seur, reduction='tsne', group.by='run')
DimPlot(seur, reduction='tsne', group.by='slice')
plot_feature(seur, reduction='tsne', features=c("nCount_RNA", "nFeature_RNA"))
dev.off()

pdf(paste0("figures/heatmap_noLQ_custom_", npcs, "pcs.pdf"), height=15, width=30)
sseur <- subset(seur, downsample=500)
print(DoHeatmap(sseur, features=var_genes))
dev.off()

pdf("figures/heatmap_noLQ_major_celltype.pdf")
DoHeatmap(sseur, features=rev(c(
"Slc6a20a",
"Aldh1l1",
"Itpr2",
"Sox10",
"Gabbr2",
"Erbb4",
"Neurod2",
"Slc17a7",
"Fcrls",
"Tmem119"
)))
dev.off()

