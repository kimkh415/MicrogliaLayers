library(Seurat)
library(ggplot2)




pdf("heatmap.pdf", height=15, width=30)
sseur <- subset(seur, downsample=500)
print(DoHeatmap(sseur, features=var_genes))
dev.off()

pdf("heatmap_major_celltype_markers.pdf")
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

