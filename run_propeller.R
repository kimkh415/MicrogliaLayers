# Draw cell type composition bar chart
# Compute cell type composition difference between/among groups
# Information about speckle package can be found here: https://github.com/Oshlack/speckle

library(Seurat)
library(speckle)
library(limma)
library(ggplot2)


# Input:
# 1. path to a Seurat object
# 2. path to an RDS file containing an array of colors
args = commandArgs(trailingOnly=T)

seur <- readRDS(args[1])
cols <- readRDS(args[2])

# Clusters are defined in the active idents slot of the Seurat object
# Layers of origination are the group
seur$group = seur$Layer
clusters=Idents(seur)
sample = factor(seur$sample)

prop.list <- getTransformedProps(cluster=Idents(seur), sample=seur$sample)
Proportions <- as.vector(t(prop.list$Proportions))
Samples <- rep(colnames(prop.list$Proportions), nrow(prop.list$Proportions))
Clusters <- rep(rownames(prop.list$Proportions), each = ncol(prop.list$Proportions))

plotdf <- data.frame(Samples = Samples, Clusters = Clusters, Proportions = Proportions)
plotdf$Clusters <- factor(plotdf$Clusters, levels=levels(seur))
plotdf$Samples <- factor(plotdf$Samples, levels=levels(sample))
mybar = ggplot(plotdf, aes(x = Samples, y = Proportions, fill = Clusters)) +
        geom_bar(stat = "identity") +  #, position=position_dodge()) +
        theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), axis.title = element_text(size = 14), legend.text = element_text(size = 12), legend.title = element_text(size = 14))

pdf("CellTypeProps.pdf")
mybar + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_fill_manual(values=cols[levels(seur)])
dev.off()

# compute sig for composition difference across all three layers for each cell type
# output one pvalue for each cell type
res = propeller(seur)
print(res)
write.table(res, file=paste0("speckle_CellType_composition_pval.tsv"), sep='\t', quote=F, row.names=T, col.names=NA)

