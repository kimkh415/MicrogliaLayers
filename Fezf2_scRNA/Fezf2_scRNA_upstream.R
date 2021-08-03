library(Seurat)
library(stringr)
library(dplyr)
library(Hmisc)
source("~/kwanho/src/seurat_tools.R")
options(future.globals.maxSize = 80 * 1024^3)

# rep1 is only KO
# use only WT from rep2
# rep3 has both

rep1_dir = "/stanley/levin_dr_storage/kwanho/jeff_microglia/data/fezf2_wt_ko/scRNA/scRNA_rep1/counts/"
rep2_dir = "/stanley/levin_dr_storage/kwanho/jeff_microglia/data/fezf2_wt_ko/scRNA/scRNA_rep2/counts/"
rep3_dir = "/stanley/levin_dr_storage/kwanho/jeff_microglia/data/fezf2_wt_ko/scRNA/scRNA_rep3/counts/"

rep1_samples = dir(rep1_dir)[file.info(paste0(rep1_dir, dir(rep1_dir)))$isdir]
rep2_samples = dir(rep2_dir)[file.info(paste0(rep2_dir, dir(rep2_dir)))$isdir]
rep3_samples = dir(rep3_dir)[file.info(paste0(rep3_dir, dir(rep3_dir)))$isdir]
rep2_samples = rep2_samples[grep("WT", rep2_samples)]

samples = c(rep1_samples, rep2_samples, rep3_samples)
print(samples)

rep1_dirlist_10x = paste0(rep1_dir, rep1_samples, "/filtered_feature_bc_matrix/")
rep2_dirlist_10x = paste0(rep2_dir, rep2_samples, "/outs/filtered_feature_bc_matrix/")
rep3_dirlist_10x = paste0(rep3_dir, rep3_samples, "/raw_feature_bc_matrix/")
dirlist_10x = c(rep1_dirlist_10x, rep2_dirlist_10x, rep3_dirlist_10x)
print(dirlist_10x)
stopifnot(all(sapply(dirlist_10x, FUN=dir.exists)))

# for rep3, takes cells from reanalysis
cells_dir = "/stanley/levin_dr_storage/kwanho/jeff_microglia/data/fezf2_wt_ko/scRNA/scRNA_rep3/rerun_cell_calling/"

# Create mtx
print("Loading counts matrices")
mtx.list = lapply(dirlist_10x, FUN=Read10X)
for (i in 1:length(mtx.list)) {
print(samples[i])
if (grepl("rep3", samples[i])) {
print("rep3: loading cell barcodes from reanalysis result")
cells_file = paste0(cells_dir, samples[i], "/outs/analysis/umap/2_components/projection.csv")
cells.to.use = data.table::fread(cells_file)
cells.to.use = cells.to.use$Barcode
mtx.list[[i]] = mtx.list[[i]][,cells.to.use]
}
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
	NormalizeData()
}

# Add cell cycle score
s.genes <- upFirst(tolower(cc.genes$s.genes))
g2m.genes <- upFirst(tolower(cc.genes$g2m.genes))
for (i in 1:length(seur.list)) {
seur.list[[i]] <- CellCycleScoring(seur.list[[i]], s.features = s.genes, g2m.features = g2m.genes, set.ident=F)
}

saveRDS(seur.list, "initial_seur_list_all_merged.rds")

# Merge datasets
seur = merge(seur.list[[1]], y=seur.list[2:length(seur.list)])

# Add metadata
sample_names = str_extract(colnames(seur), pattern="P14_Fezf2_[WTKO]+_L.+_rep[1-3]")
seur$dataset = sample_names
layer = as.factor(str_split(sample_names, "_", simplify=T)[,4])
#layer = revalue(layer, c("Lower"="L6", "Middle"="L5", "Upper"="L1-4"))
layer = ordered(layer, levels=c("L1-4","L5","L6"))
seur$Layer = layer
replicate = as.factor(str_split(sample_names, "_", simplify=T)[,5])
seur$Replicate = replicate
treat = as.factor(str_split(sample_names, "_", simplify=T)[,3])
seur$Treat = treat  # as.factor(paste0(treat, "_", replicate))
seur$age = "P14"

print(table(seur$Treat, seur$Replicate))
saveRDS(seur, "seur_all_merged_initial_merged_logNorm.rds")

# initial_qc
dir.create("figures")

pdf("figures/qc_initial_all_merged.pdf", width=10, height=5)
print(VlnPlot(seur, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), pt.size=0, ncol = 4, group.by='dataset'))
plot1 <- FeatureScatter(seur, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by='dataset')
plot2 <- FeatureScatter(seur, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by='dataset')
print(plot1 + plot2)
dev.off()

# Filtering cells
seur <- subset(seur, percent.mt<10 & nFeature_RNA>500)

pdf("figures/qc_filtered_all_merged.pdf", width=10, height=5)
print(VlnPlot(seur, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), pt.size=0, ncol = 4, group.by='dataset'))
plot1 <- FeatureScatter(seur, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by='dataset')
plot2 <- FeatureScatter(seur, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by='dataset')
print(plot1 + plot2)
dev.off()

# Process using SCTransform
seur <- seur %>%
        SCTransform(vars.to.regress=c("percent.mt")) %>%
        RunPCA(verbose=F) %>%
        RunUMAP(assay='SCT', reduction='pca', dims=1:30) %>%
        RunTSNE(assay='SCT', reduction='pca', dims=1:30) %>%
        FindNeighbors(dims = 1:30) %>%
        FindClusters()

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
DimPlot(seur, reduction='tsne', pt.size=1, split.by="Replicate", ncol=2) + NoAxes()
DimPlot(seur, reduction='tsne', pt.size=1, split.by="dataset", ncol=3) + NoAxes()
dev.off()

pdf(paste0("figures/vlnplot_all_merged_cluster_specific_qc.pdf"), height=20, width=10)
VlnPlot(seur, features=c('nCount_RNA', 'nFeature_RNA', 'percent.mt', 'percent.rb'), pt.size=0, ncol=1) + NoLegend()
dev.off()

# Find marker genes
library(xlsx)
markers = FindAllMarkers(seur, test.use='MAST', only.pos=T, latent.vars=c('percent.mt', 'nCount_RNA'), assay='RNA')
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



