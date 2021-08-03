# Replicate the analysis pipeline presented in the original MerFISH paper
# https://science.sciencemag.org/content/sci/suppl/2018/10/31/science.aau5324.DC1/aau5324-Moffitt-SM.pdf


library(Seurat)
library(dplyr)
library(ggplot2)
source("cell_type_specific_genes.R")
source("~/kwanho/src/seurat_tools.R")

print("Load Seurat object!")
seur <- readRDS("seur_init.rds")

dir.create('figures')

pdf("figures/violin_qc_raw_nCount_nGene_per_sample.pdf", height=10, width=5)
p1 = VlnPlot(seur, features=c('nCount_RNA'), pt.size=0) + NoLegend() + stat_summary(fun=median, geom='crossbar')
p2 = VlnPlot(seur, features=c('nFeature_RNA'), pt.size=0) + NoLegend() + stat_summary(fun=median, geom='crossbar')
p3 = VlnPlot(seur, features=c('area'), pt.size=0) + NoLegend() + stat_summary(fun=median, geom='crossbar')
p1+p2+p3
dev.off()

# Normalize without log-transformation then compute z-scores
seur <- NormalizeData(seur, normalization.method='RC', scale.factor=1)

print("compute Z scores")
zscore=t(apply(seur@assays$RNA@data, 1, function(X){
return((X-mean(X))/sd(X))
}))
seur@assays$RNA@scale.data = zscore

print("Run PCA")
seur <- RunPCA(seur, npcs=10)

print("Select nPCs")
seur <- JackStraw(seur, prop.freq=.5)
seur <- ScoreJackStraw(seur, dims = 1:10)

npcs=max(seur@reductions$pca@jackstraw$overall.p.values[,"PC"][seur@reductions$pca@jackstraw$overall.p.values[,'Score']<0.05])
print(paste0("using ", npcs, " PCs (JackStraw pval < 0.05)"))

outf = paste0(paste0("seur_processed_custom_varGenes_",npcs,"pcs.rds"))

seur <- FindNeighbors(seur, dims=1:npcs, k.param=10)  # default k is 20, in the MerFISH paper, they used 10-12	
seur <- FindClusters(seur, resolution=0.5)
seur <- RunTSNE(seur, dims=1:npcs, check_duplicates=F)
saveRDS(seur, outf)

pdf(paste0("tsne_final_", npcs, "pcs.pdf"))
DimPlot(seur, reduction='tsne', label=T)
DimPlot(seur, reduction='tsne', group.by='sample')
DimPlot(seur, reduction='tsne', group.by='run')
DimPlot(seur, reduction='tsne', group.by='slice')
plot_feature(seur, reduction='tsne', features=c("nCount_RNA", "nFeature_RNA"))
dev.off()

sseur <- subset(seur, idents=names(cols.mg)[12:19])
p1=DimPlot(sseur, cols=cols.n, split.by='sample', ncol=3, reduction='spatial', pt.size=.5)
sseur <- subset(seur, idents=names(cols.mg)[1:7])
p2=DimPlot(sseur, cols=cols.mg, split.by='sample', ncol=3, reduction='spatial', pt.size=.5)
sseur <- subset(seur, idents=names(cols.mg)[8:11])
p3=DimPlot(sseur, cols=cols.other, split.by='sample', ncol=3, reduction='spatial', pt.size=.5)

p1[[1]]$layers[[1]]$aes_params$alpha=0.7
p2[[1]]$layers[[1]]$aes_params$alpha=0.7
p3[[1]]$layers[[1]]$aes_params$alpha=0.7

pdf("spatial_final_neuronal.pdf", height=10, width=10)
p1
dev.off()
pdf("spatial_final_microglia.pdf", height=10, width=10)
p2
dev.off()
pdf("spatial_final_other.pdf", height=10, width=10)
p3
dev.off()

# plot all cells
p4=DimPlot(seur, cols=cols.n, split.by='sample', ncol=3, reduction='spatial', pt.size=.4)
p5=DimPlot(seur, cols=cols.mg, split.by='sample', ncol=3, reduction='spatial', pt.size=.4)
p6=DimPlot(seur, cols=cols.other, split.by='sample', ncol=3, reduction='spatial', pt.size=.4)

p4[[1]]$layers[[1]]$aes_params$alpha=0.7
p5[[1]]$layers[[1]]$aes_params$alpha=0.7
p6[[1]]$layers[[1]]$aes_params$alpha=0.7

pdf("spatial_all_cells_final_neuronal.pdf", height=10, width=10)
p4
dev.off()
pdf("spatial_all_cells_final_microglia.pdf", height=10, width=10)
p5
dev.off()
pdf("spatial_all_cells_final_other.pdf", height=10, width=10)
p6
dev.off()

