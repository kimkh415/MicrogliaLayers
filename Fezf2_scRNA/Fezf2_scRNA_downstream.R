library(Seurat)
library(stringr)
library(dplyr)
library(Hmisc)
source("~/kwanho/src/seurat_tools.R")
options(future.globals.maxSize = 80 * 1024^3)

seur <- readRDS("seur_all_merged_sct_clustered.rds")
print(table(Idents(seur)))

seur$sample = paste(seur$Layer, seur$Treat, seur$Replicate, sep='_')
seur$sample = factor(seur$sample, levels=c('L1-4_WT_rep2','L5_WT_rep2','L6_WT_rep2','L1-4_WT_rep3','L5_WT_rep3','L6_WT_rep3','L1-4_KO_rep1','L5_KO_rep1','L6_KO_rep1','L1-4_KO_rep3','L5_KO_rep3','L6_KO_rep3'))
seur$replicate = as.numeric(seur$Replicate == "rep3") #0=Rep1, 1=Rep2

pdf("figures/tsne_split_sct_cluster_by_sample.pdf", height=12, width=10)
DimPlot(seur, reduction='tsne',split.by='sample', ncol=3, pt.size=1.5)
dev.off()

# 4: 1256 activated
# 5,7: 1137, 529 LQ
# 11: 308 neuron
# 12: 191 dying cells
# 14: 61 bam
seur <- subset(seur, idents=c(4,5,7,11,12,14), invert=T)
print(table(Idents(seur)))

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
saveRDS(seur, "seur_final_Fezf2_WTKO_scRNA_sct_clustered.rds")

pdf("figures/tsne_split_final_by_sample.pdf", height=12, width=10)
DimPlot(seur, reduction='tsne',split.by='sample', ncol=3)
dev.off()

pdf("figures/tsne_final_Fezf2_WTKO_scRNA_sct.pdf")
DimPlot(seur, reduction='tsne', pt.size=1, label=T) + NoAxes()
DimPlot(seur, reduction='tsne', pt.size=1, group.by='Layer') + NoAxes()
DimPlot(seur, reduction='tsne', pt.size=1, group.by='Replicate') + NoAxes()
plot_feature(seur, reduction='tsne', features=c("nCount_RNA", "nFeature_RNA","percent.mt", "percent.rb"), title="QC")
DimPlot(seur, reduction='tsne', pt.size=1, group.by='Phase') + NoAxes()
dev.off()

pdf(paste0("figures/tsne_split_final_Fezf2_WTKO_scRNA_sct.pdf"), height=8, width=16)
DimPlot(seur, reduction='tsne', pt.size=1, group.by='Layer',  split.by="Replicate", ncol=3) + NoAxes()
DimPlot(seur, reduction='tsne', pt.size=1, group.by='Replicate',  split.by="Treat", ncol=3) + NoAxes()
DimPlot(seur, reduction='tsne', pt.size=1, split.by="Replicate", ncol=2) + NoAxes()
DimPlot(seur, reduction='tsne', pt.size=1, split.by="dataset", ncol=3) + NoAxes()
dev.off()

vp1 = VlnPlot(seur, features='nCount_RNA', pt.size=0, ncol=1) + NoLegend() + stat_summary(fun=median, geom='crossbar')
vp2 = VlnPlot(seur, features='nFeature_RNA', pt.size=0, ncol=1) + NoLegend() + stat_summary(fun=median, geom='crossbar')
vp3 = VlnPlot(seur, features='percent.mt', pt.size=0, ncol=1) + NoLegend() + stat_summary(fun=median, geom='crossbar')
pdf("figures/violin_final_Fezf2_WTKO_scRNA_cluster_specific_qc.pdf", height=15, width=10)
vp1+vp2+vp3
dev.off()

vp1 = VlnPlot(seur, group.by='sample', features='nCount_RNA', pt.size=0, ncol=1) + NoLegend() + stat_summary(fun=median, geom='crossbar')
levels(seur$sample)
vp2 = VlnPlot(seur, group.by='sample', features='nFeature_RNA', pt.size=0, ncol=1) + NoLegend() + stat_summary(fun=median, geom='crossbar')
vp3 = VlnPlot(seur, group.by='sample', features='percent.mt', pt.size=0, ncol=1) + NoLegend() + stat_summary(fun=median, geom='crossbar')
pdf("figures/violin_final_Fezf2_WTKO_scRNA_sample_qc.pdf", height=15, width=10)
vp1+vp2+vp3
dev.off()

print("exit")
q()

colors = readRDS("../colors.rds")
prev_anno <- readRDS("../metadata_Fezf2WTKO_scRNA_harmony_CellType.rds")
seur$prev_anno <- prev_anno
pdf("figures/tsne_final_prev_anno.pdf")
DimPlot(seur, reduction='tsne', group.by='prev_anno', cols=colors)
dev.off()

# Find marker genes
library(xlsx)
markers = FindAllMarkers(seur, test.use='MAST', only.pos=T, latent.vars=c('percent.mt', 'nCount_RNA', 'replicate'))
markers <- markers %>% filter(p_val_adj<=0.05)
saveRDS(markers, "markers_MAST_final_Fezf2_WTKO_scRNA.txt")
outfname = "markers_MAST_final_Fezf2_WTKO_scRNA.xlsx"
for (i in 0:max(as.numeric(levels(Idents(seur))))) {
  cur_markers = subset(markers, cluster==i)
  if (nrow(cur_markers)>0) {
  if (!file.exists(outfname)) {write.xlsx(cur_markers[order(-cur_markers$avg_logFC),], file=outfname, sheetName=paste0("cluster",i))}
  else {write.xlsx(cur_markers[order(-cur_markers$avg_logFC),], file=outfname, sheetName=paste0("cluster",i), append=T)}
  }
}


print("DONE")



