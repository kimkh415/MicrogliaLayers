.libPaths(.libPaths()[2])

library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(patchwork)


ref <- readRDS("/stanley/levin_dr_storage/kwanho/jeff_microglia/paper/seur_P14_azimuth_ref_ds.rds")
query.all <- readRDS("../seur_final_Fezf2_WTKO_SCT_annotated.rds")
print(table(query.all$Treat))

for (gt in c('KO')) {
print(gt)
dir.create(gt)
setwd(gt)

query = subset(query.all, Treat==gt)
print(table(query$Treat))

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
saveRDS(query$predicted.celltype, paste0("metadata_",gt,"_azimuth_preds.rds"))
saveRDS(query$predicted.celltype.score, paste0("metadata_",gt,"_azimuth_scores.rds"))

print("Plotting!")

cols = readRDS("../../colors.rds")
cols = cols[levels(ref$subclass_label)]

p3 = DimPlot(query, reduction = "tsne", group.by = "predicted.celltype", label.size = 3, repel = TRUE, cols=cols) + ggtitle("Azimuth prediction")
p4 = FeaturePlot(query, reduction = "tsne", features = "predicted.celltype.score") + ggtitle("Prediction score")
p5 = DimPlot(query, reduction = "tsne", group.by = "CellType", label.size = 3, cols=cols, repel = TRUE) + ggtitle("Manual annotation")
pdf(paste0("tsne_",gt,"_azimuth_pred.pdf"), height=5, width=10)
print(p3 + p4)
print(p3 + p5)
dev.off()

tab = table(query$CellType, query$predicted.celltype)
prop = tab/rowSums(tab)
#prop = prop[which(!rownames(prop) %in% c('L5 NP','L5 PT','L6 CThPN 1','L6 CThPN 2','L6b Subplate')),]
prop = prop[which(!rowSums(tab)==0),]

prop = apply(prop,2,rev)
prop = prop[,names(cols)]
pdf(paste0("barplot_",gt,"_pred_props.pdf"), height=8, width=12)
par(mar=c(5,8,5,11), xpd=T, mgp=c(3,0.5,0))
barplot(t(prop), main="Azimuth prediction", horiz=T, xlab="Proportion", col=cols, las=1)
legend("right", inset=c(-0.225,0), legend=colnames(prop), fill=cols)
dev.off()

pdf(paste0("violin_",gt,"_manual_CellType_pred_score.pdf"))
VlnPlot(query, group.by='CellType', features='predicted.celltype.score', pt.size=0, cols=cols) + NoLegend() + stat_summary(fun=median, geom='crossbar') + ggtitle("Prediction score for each manually labeled cell type")
dev.off()

query$predicted.celltype = factor(query$predicted.celltype, levels=names(cols))
pdf(paste0("violin_",gt,"_pred_CellType_pred_score.pdf"))
VlnPlot(query, group.by='predicted.celltype', features='predicted.celltype.score', pt.size=0, cols=cols) + NoLegend() + stat_summary(fun=median, geom='crossbar') + ggtitle("Prediction score for each predicted cell type")
dev.off()

setwd("../")

}

