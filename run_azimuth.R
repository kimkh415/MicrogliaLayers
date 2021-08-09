# This was run in R v4.0.3 and Seurat v4

library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(patchwork)


# input arguments:
# 1. path to downsampled Seurat dataset to be used as the reference
#	ref data must have cell types defined in 'subclass_label' column of the metadata
# 2. path to the query Seurat dataset
args = commandArgs(trailinOnly=T)

# Reference data
ref <- readRDS(args[1])
# query data
query <- readRDS(args[2])


print("Find achors!")
anchors <- FindTransferAnchors(
  reference = ref,
  query = query,
  normalization.method = "SCT",
  reference.reduction = "spca",
  dims = 1:50
)

print("Transfer data!")
query <- TransferData(
  anchorset = anchors, 
  reference = ref,
  query = query,
  refdata = list(
    celltype = "subclass_label")
)

print("Integrate embeddings!")
query <- IntegrateEmbeddings(
  anchorset = anchors,
  reference = ref,
  query = query, 
  new.reduction.name = "ref.spca"
)


print("Saving!")
saveRDS(query$predicted.celltype, "metadata_azimuth_preds.rds")
saveRDS(query$predicted.celltype.score, "metadata_azimuth_scores.rds")

print("Plotting!")
# Manually annotated cell types for the query dataset is stored in the 'CellType' column of the metadata

cols = readRDS("SeuratAnalysis/mg_colors.rds")
cols = cols[levels(ref$subclass_label)]

p3 = DimPlot(query, reduction = "tsne", group.by = "predicted.celltype", label.size = 3, repel = TRUE, cols=cols) + ggtitle("Azimuth prediction")
p4 = FeaturePlot(query, reduction = "tsne", features = "predicted.celltype.score") + ggtitle("Prediction score")
p5 = DimPlot(query, reduction = "tsne", group.by = "CellType", label.size = 3, cols=cols, repel = TRUE) + ggtitle("Manual annotation")
pdf("tsne_azimuth_pred.pdf", height=5, width=10)
print(p3 + p4)
print(p3 + p5)
dev.off()

tab = table(query$CellType, query$predicted.celltype)
prop = tab/rowSums(tab)
prop = prop[which(!rowSums(tab)==0),]
prop = apply(prop,2,rev)
prop = prop[,names(cols)]

pdf("barplot_pred_props.pdf", height=8, width=12)
par(mar=c(5,8,5,11), xpd=T, mgp=c(3,0.5,0))
barplot(t(prop), main="Azimuth prediction", horiz=T, xlab="Proportion", col=cols, las=1)
legend("right", inset=c(-0.225,0), legend=colnames(prop), fill=cols)
dev.off()

pdf("violin_manual_CellType_pred_score.pdf")
VlnPlot(query, group.by='CellType', features='predicted.celltype.score', pt.size=0, cols=cols) + NoLegend() + stat_summary(fun=median, geom='crossbar') + ggtitle("Prediction score for each manually labeled cell type")
dev.off()

query$predicted.celltype = factor(query$predicted.celltype, levels=names(cols))
pdf("violin_pred_CellType_pred_score.pdf")
VlnPlot(query, group.by='predicted.celltype', features='predicted.celltype.score', pt.size=0, cols=cols) + NoLegend() + stat_summary(fun=median, geom='crossbar') + ggtitle("Prediction score for each predicted cell type")
dev.off()


