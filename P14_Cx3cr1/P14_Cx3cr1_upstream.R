library(Seurat)
library(stringr)
library(dplyr)
library(pracma)
library(Hmisc)
library(plyr)
source("~/kwanho/src/seurat_tools.R")
options(future.globals.maxSize = 80 * 1024^3)

setwd("/stanley/levin_dr_storage/kwanho/jeff_microglia/data/P14_Cx3cr1_layers/reprocess/counts")

samples = list.files()
print(samples)
dirlist_10x = paste0(samples, "/filtered_feature_bc_matrix/")
print(dirlist_10x)
stopifnot(all(sapply(dirlist_10x, FUN=dir.exists)))

# Create mtx
mtx.list = lapply(dirlist_10x, FUN=Read10X)
for (i in 1:length(mtx.list)) {
colnames(mtx.list[[i]]) = paste0(samples[i], '_', colnames(mtx.list[[i]]))
}

# Create Seurat objects
print("creating seurate objects")
seur.list = lapply(mtx.list, FUN=CreateSeuratObject, min.cells=3, min.features=200)
names(seur.list) <- samples
print(seur.list)

for (i in 1:length(seur.list)) {
seur.list[[i]] <- seur.list[[i]] %>%
	PercentageFeatureSet(pattern="^mt-", col.name="percent.mt") %>%
	PercentageFeatureSet(pattern="^Rp[sl]", col.name="percent.rb")%>%
	NormalizeData() %>%
	FindVariableFeatures() %>%
	ScaleData(vars.to.regress=c('percent.mt', 'percent.rb'))
}

# Add cell cycle score
s.genes <- upFirst(tolower(cc.genes$s.genes))
g2m.genes <- upFirst(tolower(cc.genes$g2m.genes))
for (i in 1:length(seur.list)) {
seur.list[[i]] <- CellCycleScoring(seur.list[[i]], s.features = s.genes, g2m.features = g2m.genes, set.ident=F)
}

setwd("../analysis/")
saveRDS(seur.list, "initial_seur_list.rds")

# initial_qc
dir.create("figures")

for (i in 1:length(seur.list)) {
obj = seur.list[[i]]
obj.name = names(seur.list)[i]
pdf(paste0("figures/", obj.name, "_initial_qc.pdf"), width=10, height=5)
print(VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), pt.size=0, ncol = 4))
plot1 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
print(plot1 + plot2)
dev.off()
}

pdf("figures/initial_hist_qc.pdf", height=12, width=8)
par(mfrow=c(6,4))
for (i in 1:length(seur.list)) {
obj = seur.list[[i]]
obj.name = names(seur.list)[i]
hist(obj$nCount_RNA, breaks=linspace(0, max(obj$nCount_RNA), 20), main=paste0(obj.name, " nUMI"))
hist(obj$nFeature_RNA, breaks=linspace(0, max(obj$nFeature_RNA), 20), main="nGene")
hist(obj$percent.mt, breaks=linspace(0, max(obj$percent.mt), 20), main="percent mito")
hist(obj$percent.rb, breaks=linspace(0, max(obj$percent.rb), 20), main="percent ribo")
}
dev.off()


# Filtering cells
for (i in 1:length(seur.list)) {
seur.list[[i]] <- subset(seur.list[[i]], percent.mt<10 & nFeature_RNA>500)  # 4 cells with >1e5 nCount
}

for (i in 1:length(seur.list)) {
obj = seur.list[[i]]
obj.name = names(seur.list)[i]
pdf(paste0("figures/", obj.name, "_filtered_qc.pdf"), width=10, height=5)
print(VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), pt.size=0, ncol = 4))
plot1 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
print(plot1 + plot2)
dev.off()
}

pdf("figures/filtered_hist_qc.pdf", height=12, width=8)
par(mfrow=c(6,4))
for (i in 1:length(seur.list)) {
obj = seur.list[[i]]
obj.name = names(seur.list)[i]
hist(obj$nCount_RNA, breaks=linspace(0, max(obj$nCount_RNA), 20), main=paste0(obj.name, " nUMI"))
hist(obj$nFeature_RNA, breaks=linspace(0, max(obj$nFeature_RNA), 20), main="nGene")
hist(obj$percent.mt, breaks=linspace(0, max(obj$percent.mt), 20), main="percent mito")
hist(obj$percent.rb, breaks=linspace(0, max(obj$percent.rb), 20), main="percent ribo")
}
dev.off()

# Merge datasets
seur = merge(seur.list[[1]], y=seur.list[2:length(seur.list)])

# Add metadata
sample_names = str_extract(colnames(seur), pattern="P14_[a-zA-Z]+_Rep[1-2]")
seur$dataset = sample_names
layer = as.factor(str_split(sample_names, "_", simplify=T)[,2])
layer = revalue(layer, c("Lower"="L6", "Middle"="L5", "Upper"="L1-4"))
layer = ordered(layer, levels=c("L1-4","L5","L6"))
seur$Layer = layer
replicate = as.factor(str_split(sample_names, "_", simplify=T)[,3])
seur$Replicate = replicate
seur$age = "P14"

# Noramlize RNA assay
seur <- NormalizeData(seur)

saveRDS(seur, "seur_merged_logNorm.rds")

# Process using SCTransform
seur <- seur %>%
        SCTransform(vars.to.regress=c("percent.mt", "Replicate")) %>%
	RunPCA(verbose=F) %>%
        RunUMAP(assay='SCT', reduction='pca', dims=1:30) %>%
	RunTSNE(assay='SCT', reduction='pca', dims=1:30) %>%
        FindNeighbors(dims = 1:30) %>%
        FindClusters(resolution=0.5)

saveRDS(seur, "seur_sct_clustered.rds")

pdf("figures/umap_sct_0.5res.pdf")
DimPlot(seur, reduction='umap', pt.size=.1, label=T) + NoAxes()
DimPlot(seur, reduction='umap', pt.size=.1, group.by='Replicate') + NoAxes()
DimPlot(seur, reduction='umap', pt.size=.1, group.by='Layer') + NoAxes()
plot_feature(seur, features=c("nCount_RNA", "nFeature_RNA","percent.mt", "percent.rb"), title="QC on UMAP")
DimPlot(seur, reduction='umap', pt.size=.1, group.by='Phase') + NoAxes()
dev.off()

pdf(paste0("figures/umap_sct_split_0.5res.pdf"), height=8, width=16)
DimPlot(seur, reduction='umap', pt.size=.1, split.by="Replicate", ncol=2) + NoAxes()
DimPlot(seur, reduction='umap', pt.size=.1, split.by="Replicate", ncol=2, group.by='Layer') + NoAxes()
DimPlot(seur, reduction='umap', pt.size=.1, split.by="Layer", ncol=3) + NoAxes()
DimPlot(seur, reduction='umap', pt.size=.1, split.by="dataset", ncol=3) + NoAxes()
dev.off()

pdf("figures/tsne_sct_0.5res.pdf")
DimPlot(seur, reduction='tsne', pt.size=.1, label=T) + NoAxes()
DimPlot(seur, reduction='tsne', pt.size=.1, group.by='Replicate') + NoAxes()
DimPlot(seur, reduction='tsne', pt.size=.1, group.by='Layer') + NoAxes()
plot_feature(seur, reduction='tsne', features=c("nCount_RNA", "nFeature_RNA","percent.mt", "percent.rb"), title="QC on UMAP")
DimPlot(seur, reduction='tsne', pt.size=.1, group.by='Phase') + NoAxes()
dev.off()

pdf(paste0("figures/tsne_sct_split_0.5res.pdf"), height=8, width=16)
DimPlot(seur, reduction='tsne', pt.size=.1, split.by="Replicate", ncol=2) + NoAxes()
DimPlot(seur, reduction='tsne', pt.size=.1, split.by="Replicate", ncol=2, group.by='Layer') + NoAxes()
DimPlot(seur, reduction='tsne', pt.size=.1, split.by="Layer", ncol=3) + NoAxes()
DimPlot(seur, reduction='tsne', pt.size=.1, split.by="dataset", ncol=3) + NoAxes()
dev.off()

pdf(paste0("figures/vlnplot_cluster_specific_qc_0.5res.pdf"), height=20, width=10)
VlnPlot(seur, features=c('nCount_RNA', 'nFeature_RNA', 'percent.mt', 'percent.rb'), pt.size=0, ncol=1) + NoLegend()
dev.off()

# Map previous annotations
seur$prev_anno <- readRDS("old_metadata_annotation.rds")
pdf("figures/umap_sct_prev_anno.pdf")
DimPlot(seur, reduction='umap', pt.size=.1, label=T, group.by='prev_anno') + NoAxes()
DimPlot(seur, reduction='tsne', pt.size=.1, label=T, group.by='prev_anno') + NoAxes()
dev.off()

# Find marker genes
library(xlsx)
markers = FindAllMarkers(seur, test.use='MAST', only.pos=T, latent.vars=c('percent.mt', 'nCount_RNA'), assay='RNA')
write.table(markers, file="markerListsMAST.txt", sep='\t', row.names=F, quote=F)
outfname = "markers_MAST.xlsx"
for (i in 0:max(as.numeric(levels(Idents(seur))))) {
  cur_markers = subset(markers, cluster==i)
  if (nrow(cur_markers)>0) {
  if (i==0) {write.xlsx(cur_markers[order(-cur_markers$avg_logFC),], file=outfname, sheetName=paste0("cluster",i,"_markers"))}
  else {write.xlsx(cur_markers[order(-cur_markers$avg_logFC),], file=outfname, sheetName=paste0("cluster",i,"_markers"), append=T)}
  }
}

## Update annotation
#lq1 = WhichCells(seur, idents='3')
#lq2 = WhichCells(seur, idents='5')
#lq3 = WhichCells(seur, idents='12')
#bam = WhichCells(seur, idents='9')
#
#seur$prev_anno2 = as.character(seur$prev_anno)
#seur$prev_anno2[lq1] <- "LQ1"
#seur$prev_anno2[lq2] <- "LQ2"
#seur$prev_anno2[lq3] <- "LQ3"
#seur$prev_anno2[bam] <- "BAM"
#unknown = names(which(is.na(seur$prev_anno2)))
#seur$prev_anno2[unknown] <- "Unknown"
#
#pdf("figures/umap_sct_prev_anno2.pdf")
#DimPlot(seur, reduction='umap', pt.size=.1, group.by='prev_anno2') + NoAxes()
#dev.off()
#
#saveRDS(seur$prev_anno2, "metadata_prev_anno2.rds")

print("DONE")



