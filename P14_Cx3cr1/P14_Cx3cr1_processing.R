library(Seurat)
library(stringr)
library(dplyr)
library(pracma)
library(Hmisc)
library(plyr)
source("~/kwanho/src/seurat_tools.R")
options(future.globals.maxSize = 80 * 1024^3)

seur <- readRDS("../seur_sct_clustered.rds")
seur$prev_anno <- readRDS("../old_metadata_annotation.rds")

# 4, 12: Low Quality
# 7: BAMs
# 11: oligo
seur <- subset(seur, idents=c(4,7,11,12), invert=T)

# Process using SCTransform
seur <- seur %>%
        SCTransform(vars.to.regress=c("percent.mt", "Replicate")) %>%
        RunPCA(verbose=F) %>%
        RunUMAP(assay='SCT', reduction='pca', dims=1:30) %>%
        RunTSNE(assay='SCT', reduction='pca', dims=1:30) %>%
        FindNeighbors(dims = 1:30) %>%
        FindClusters(resolution=0.5)

saveRDS(seur, "seur_subset_sct_clustered.rds")

dir.create("figures")
pdf("figures/umap_subset_sct_0.5res.pdf")
DimPlot(seur, reduction='umap', pt.size=.1, label=T) + NoAxes()
DimPlot(seur, reduction='umap', pt.size=.1, group.by='Replicate') + NoAxes()
DimPlot(seur, reduction='umap', pt.size=.1, group.by='Layer') + NoAxes()
DimPlot(seur, reduction='umap', pt.size=.1, group.by='prev_anno', label=T) + NoAxes()
plot_feature(seur, reduction='umap', features=c("nCount_RNA", "nFeature_RNA","percent.mt", "percent.rb"), title="QC on UMAP")
DimPlot(seur, reduction='umap', pt.size=.1, group.by='Phase') + NoAxes()
dev.off()

pdf(paste0("figures/umap_subset_sct_split_0.5res.pdf"), height=8, width=16)
DimPlot(seur, reduction='umap', pt.size=.1, split.by="Replicate", ncol=2) + NoAxes()
DimPlot(seur, reduction='umap', pt.size=.1, split.by="Replicate", ncol=2, group.by='Layer') + NoAxes()
DimPlot(seur, reduction='umap', pt.size=.1, split.by="Layer", ncol=3) + NoAxes()
DimPlot(seur, reduction='umap', pt.size=.1, split.by="dataset", ncol=3) + NoAxes()
dev.off()

pdf("figures/tsne_subset_sct_0.5res.pdf")
DimPlot(seur, reduction='tsne', pt.size=.1, label=T) + NoAxes()
DimPlot(seur, reduction='tsne', pt.size=.1, group.by='Replicate') + NoAxes()
DimPlot(seur, reduction='tsne', pt.size=.1, group.by='Layer') + NoAxes()
DimPlot(seur, reduction='tsne', pt.size=.1, group.by='prev_anno', label=T) + NoAxes()
plot_feature(seur, reduction='tsne', features=c("nCount_RNA", "nFeature_RNA","percent.mt", "percent.rb"), title="QC on UMAP")
DimPlot(seur, reduction='tsne', pt.size=.1, group.by='Phase') + NoAxes()
dev.off()

pdf(paste0("figures/tsne_subset_sct_split_0.5res.pdf"), height=8, width=16)
DimPlot(seur, reduction='tsne', pt.size=.1, split.by="Replicate", ncol=2) + NoAxes()
DimPlot(seur, reduction='tsne', pt.size=.1, split.by="Replicate", ncol=2, group.by='Layer') + NoAxes()
DimPlot(seur, reduction='tsne', pt.size=.1, split.by="Layer", ncol=3) + NoAxes()
DimPlot(seur, reduction='tsne', pt.size=.1, split.by="dataset", ncol=3) + NoAxes()
dev.off()

pdf(paste0("figures/vlnplot_subset_cluster_specific_qc_0.5res.pdf"), height=20, width=10)
VlnPlot(seur, features=c('nCount_RNA', 'nFeature_RNA', 'percent.mt', 'percent.rb'), pt.size=0, ncol=1) + NoLegend()
dev.off()

genes <- c('Jun','Junb','Fos','Egr1','Hspa1a','Dusp1')
pdf("figures/vlnplot_subset_activation.pdf", height=8, width=16)
VlnPlot(seur, features=genes, split.by='Replicate', ncol=3, assay='RNA', slot='data', split.plot=T, pt.size=0)
dev.off()

# Find marker genes
library(xlsx)
markers = FindAllMarkers(seur, test.use='MAST', only.pos=T, latent.vars=c('percent.mt', 'nCount_RNA'), assay='RNA')
write.table(markers, file="markerSubsetListsMAST.txt", sep='\t', row.names=F, quote=F)
outfname = "markers_subset_MAST.xlsx"
for (i in 0:max(as.numeric(levels(Idents(seur))))) {
  cur_markers = subset(markers, cluster==i)
  if (nrow(cur_markers)>0) {
  if (i==0) {write.xlsx(cur_markers[order(-cur_markers$avg_logFC),], file=outfname, sheetName=paste0("cluster",i,"_markers"))}
  else {write.xlsx(cur_markers[order(-cur_markers$avg_logFC),], file=outfname, sheetName=paste0("cluster",i,"_markers"), append=T)}
  }
}

mg_subtype <- c('C1qa', 'Csf1r', 'Tmem119', 'Sall1', 'P2ry12', 'Ifi27l2a', 'Ifitm3', 'Cd63', 'Cd9', 'Apoe', 'Lyz2', 'Ccr1', 'Mcm6', 'Mki67', 'Gm26870', 'Gm10800')
plot_feature2(seur, reduction='tsne', features=mg_subtype, size=4, filename="figures/featureplot_subtype_genes.pdf")

print("DONE")



