library(Seurat)
library(stringr)
library(dplyr)
library(Hmisc)
source("~/kwanho/src/seurat_tools.R")
options(future.globals.maxSize = 80 * 1024^3)

seur <- readRDS("../seur_all_merged_initial_merged_logNorm.rds")

dir.create("figures")

# Filtering cells
seur <- subset(seur, percent.mt<10 & nFeature_RNA>500)

pdf("figures/qc_filtered_all_merged.pdf", width=10, height=5)
print(VlnPlot(seur, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), pt.size=0, ncol = 4, group.by='dataset'))
plot1 <- FeatureScatter(seur, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by='dataset')
plot2 <- FeatureScatter(seur, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by='dataset')
print(plot1 + plot2)
dev.off()

# Process using SCTransform
seur <- SCTransform(seur, vars.to.regress=c("percent.mt", "Replicate"))
vg = VariableFeatures(seur)
print(length(vg))
VariableFeatures(seur) = setdiff(vg, grep("^mt-|^Rp[sl]", vg, value=T))
print(length(VariableFeatures(seur)))

seur <- seur %>%
        RunPCA() %>%
        RunTSNE(assay='SCT', reduction='pca', dims=1:30) %>%
        FindNeighbors(dims = 1:30) %>%
        FindClusters(resolution=0.5)

print("Saving!")
saveRDS(seur, "seur_all_merged_sct_clustered.rds")

pdf("figures/tsne_all_merged_sct.pdf")
DimPlot(seur, reduction='tsne', pt.size=1, label=T) + NoAxes()
DimPlot(seur, reduction='tsne', pt.size=1, group.by='Layer') + NoAxes()
DimPlot(seur, reduction='tsne', pt.size=1, group.by='Replicate') + NoAxes()
plot_feature(seur, reduction='tsne', features=c("nCount_RNA", "nFeature_RNA","percent.mt", "percent.rb"), title="QC")
DimPlot(seur, reduction='tsne', pt.size=1, group.by='Phase') + NoAxes()
dev.off()

pdf(paste0("figures/tsne_all_merged_sct_split.pdf"), height=8, width=16)
DimPlot(seur, reduction='tsne', pt.size=1, group.by='Layer',  split.by="Replicate", ncol=3) + NoAxes()
DimPlot(seur, reduction='tsne', pt.size=1, group.by='Replicate',  split.by="Treat", ncol=3) + NoAxes()
DimPlot(seur, reduction='tsne', pt.size=1, split.by="Replicate", ncol=2) + NoAxes()
DimPlot(seur, reduction='tsne', pt.size=1, split.by="dataset", ncol=3) + NoAxes()
dev.off()

pdf(paste0("figures/vlnplot_all_merged_cluster_specific_qc.pdf"), height=20, width=10)
VlnPlot(seur, features=c('nCount_RNA', 'nFeature_RNA', 'percent.mt', 'percent.rb'), pt.size=0, ncol=1) + NoLegend()
dev.off()

colors = readRDS("../colors.rds")
seur$prev_anno <- readRDS("../metadata_final_Fezf2WTKO_scRNA_SCT_CellType.rds")
pdf("figures/tsne_prev_final_anno.pdf")
DimPlot(seur, reduction='tsne', group.by='prev_anno', cols=colors)
dev.off()

# Find marker genes
library(xlsx)
markers = FindAllMarkers(seur, test.use='MAST', only.pos=T, latent.vars=c('percent.mt', 'nCount_RNA'))
saveRDS(markers, "markers_MAST_all_merged.txt")
outfname = "markers_MAST_all_merged.xlsx"
for (i in 0:max(as.numeric(levels(Idents(seur))))) {
  cur_markers = subset(markers, cluster==i)
  if (nrow(cur_markers)>0) {
  if (i==0) {write.xlsx(cur_markers[order(-cur_markers$avg_logFC),], file=outfname, sheetName=paste0("cluster",i,"_markers"))}
  else {write.xlsx(cur_markers[order(-cur_markers$avg_logFC),], file=outfname, sheetName=paste0("cluster",i,"_markers"), append=T)}
  }
}

print("DONE")



