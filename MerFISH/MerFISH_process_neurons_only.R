library(Seurat)
source("cell_type_specific_genes.R")
source("~/kwanho/src/seurat_tools.R")

seur <- readRDS("seur_noLQ_processed_custom_varGenes_9pcs.rds")

plot_feature2(seur, features=c('Tmem119','Fcrls','Slc17a7','Neurod2','Sox10','Itpr2','Erbb4','Gabbr2','Aldh1l1','Slc6a20a'), title="Major cell type markers", size=4, filename="figures/featureplot_major_cellTypeMarkers.pdf", reduction='tsne')

a = c('PNs','PNs','PNs','PNs','INs','Oligo', 'Unknown','Astrocyte','Mg','PNs','Mg','Endothelia','Astrocyte','Oligo','PNs')
names(a) = levels(seur)

seur <- RenameIdents(seur, a)
Idents(seur) <- factor(Idents(seur), levels=c("Mg","PNs","INs","Oligo","Astrocyte","Endothelia","Unknown"))

my.cols = c('#53868B','#696969','#71C671','#FFA500','#CD7054','#8B7500', "gray")
names(my.cols)=c('PNs','INs','Mg','Astrocyte','Oligo','Endothelia', 'Unknown')

pdf("figures/violin_qc_majorCellType.pdf", height=10, width=5)
p1=VlnPlot(seur, features=c("nFeature_RNA"), cols=my.cols, pt.size=0) + stat_summary(fun=median, geom='crossbar') + NoLegend()
p2=VlnPlot(seur, features=c("nCount_RNA"), cols=my.cols, pt.size=0) + stat_summary(fun=median, geom='crossbar') + NoLegend()
p3=VlnPlot(seur, features=c("area"), cols=my.cols, pt.size=0) + stat_summary(fun=median, geom='crossbar') + NoLegend()
p1+p2+p3
dev.off()

pdf("figures/tsne_major_CellTypes.pdf")
DimPlot(seur, cols=my.cols)
dev.off()

sseur <- subset(seur, downsample=500)

pdf("figures/heatmap_major_cellTypes.pdf", height=6, width=8)
DoHeatmap(sseur, group.colors=my.cols[levels(sseur)], features=c('Tmem119','Fcrls','Slc17a7','Neurod2','Erbb4','Gabbr2','Sox10','Itpr2','Aldh1l1','Slc6a20a'))
dev.off()

pdf(paste0("figures/heatmap_major_CellTypes_var_genes.pdf"), height=8, width=15)
print(DoHeatmap(sseur, group.colors=my.cols, features=var_genes))
dev.off()

saveRDS(Idents(seur), "metadata_major_CellType.rds")

cycling_genes = c('Mcm4','Ube2c','Mki67')

for (ct in c("Neuronal")) {
print(paste0("Processing ", ct))
if (ct == "Neuronal") {
obj <- subset(seur, idents=c("PNs", "INs"))
} else {
obj <- subset(seur, idents=c("Mg", "Oligo", "Astrocyte","Endothelia","Unknown"))
}
print(table(Idents(obj)))
cur_vg = get(paste0(ct, "_genes"))
cur_vg = c(cur_vg, cycling_genes)
print(cur_vg)
# Set new variable genes
VariableFeatures(obj) = cur_vg
# Run PCA
obj <- RunPCA(obj, npcs=10)
obj <- JackStraw(obj, prop.freq=.9)
obj <- ScoreJackStraw(obj, dims = 1:10)
pdf(paste0("figures/pc_",ct,"_choose_nPCs.pdf"))
print(JackStrawPlot(obj, dims = 1:10))
print(ElbowPlot(obj))
dev.off()
pdf(paste0("figures/pc_",ct,"_dim_loadings.pdf"), height=14)
print(VizDimLoadings(obj, dims=1:10, reduction='pca'))
dev.off()
pdf(paste0("figures/pc_",ct,"_dimplots.pdf"))
for (i in 1:9) {
print(DimPlot(obj, reduction="pca", dims=c(i,i+1)))
}
dev.off()

#npcs=max(obj@reductions$pca@jackstraw$overall.p.values[,"PC"][obj@reductions$pca@jackstraw$overall.p.values[,'Score']<0.05])
#print(npcs)
#if (npcs < 2) {npcs=2}
if (ct == "Neuronal") {
npcs=7
} else {
npcs=
}
print(paste0("using the top ", npcs, "PCs"))

obj <- FindNeighbors(obj, dims=1:npcs, k.param=10)
obj <- FindClusters(obj, resolution=0.4)
obj <- RunTSNE(obj, dims=1:npcs, check_duplicates=F)
print("Saving!")
saveRDS(obj, paste0("seur_", ct, "_processed.rds"))

pdf(paste0("figures/tsne_", ct, ".pdf"))
print(DimPlot(obj, reduction='tsne', label=T))
print(DimPlot(obj, reduction='tsne', group.by='sample'))
print(DimPlot(obj, reduction='tsne', group.by='run'))
print(DimPlot(obj, reduction='tsne', group.by='slice'))
plot_feature(obj, reduction='tsne', features=c("nCount_RNA", "nFeature_RNA"))
dev.off()

pdf(paste0("figures/heatmap_", ct, ".pdf"), height=8, width=15)
sobj <- subset(obj, downsample=300)
print(DoHeatmap(sobj, features=cur_vg))
dev.off()
}

obj <- readRDS("seur_Neuronal_processed.rds")

a = c('L2-3 CPN','CThPN','L4 Stellate','INs','L4 Stellate','DL CPN','L5 CPN/CStrPN','L5 PT','L5 NP','L6b Subplate')
names(a) <- levels(obj)



# plot spatial
pdf("figures/spatial.pdf", height=10, width=10)
DimPlot(seur, reduction='spatial', split.by='sample',cols=cols.n, ncol=3)
DimPlot(seur, reduction='spatial', split.by='sample',cols=cols.other, ncol=3)
dev.off()

# Add mg subtype anno
# remove duplicates
a = readRDS("../../mg_labels/mg_labels.rds")
aa = unlist(a)
names(aa)[grep("Homeostatic1", names(aa))] <- "Homeostatic1"
names(aa)[grep("Homeostatic2", names(aa))] <- "Homeostatic2"
names(aa)[grep("Innate Immune", names(aa))] <- "Innate Immune"
names(aa)[grep("Inflammatory", names(aa))] <- "Inflammatory"
names(aa)[grep("Apoe", names(aa))] <- "Apoe+"
names(aa)[grep("Ccr1", names(aa))] <- "Ccr1+"
names(aa)[grep("Dividing", names(aa))] <- "Dividing"

mg <- subset(seur, idents='Mg')
mg$mg_subtype <- "Mg"

b = aa[which(!is.na(match(aa, mg$location)))]
mg@meta.data[match(b, mg$location), "mg_subtype"] <- names(b)

