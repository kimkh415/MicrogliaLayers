library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(viridis)



# Draw and save Seurat's feature plot
# Most parameters are defined in the Seurat's FeaturePlot
# Specifying 'title' puts a title to the plot
MyPlotFeature <- function(seur, features, reduction='umap', title="", cols='none', pt.size=.1, min.cutoff=NA, max.cutoff=NA, filename="test.pdf", size=4, slot='data')
{
if(cols=='none'){cols <- rev(viridis(50))}
p <- FeaturePlot(seur, features=features, reduction=reduction, slot=slot, pt.size=pt.size, min.cutoff=min.cutoff, max.cutoff=max.cutoff, combine=F)
for (i in 1:length(p))
{
p[[i]] <- p[[i]] + NoAxes() + scale_colour_gradientn(colours=cols)
}
nc=ceiling(sqrt(length(features)))
print(paste0("n_columns=", nc))
plot.width = nc*size
plot.height = ceiling(length(features)/nc)*size
print(paste0("plot_dimensions=", plot.width, " x ", plot.height))
pdf(filename, width=plot.width, height=plot.height)
print(ggpubr::annotate_figure(
p = cowplot::plot_grid(plotlist = p, ncol = nc),
top = ggpubr::text_grob(label=title, face='bold', size=20)))
dev.off()
}


# Violin plots of Sam Marsh module scores and Homeostatic Mg signatures
# Takes a Seurat object as an input
# Assums that the module scores are already in the metadata
ViolinPlotModule <- function(seur) {
p1 = VlnPlot(seur, features=c('MG_ID_score'), pt.size=0) + NoLegend() + stat_summary(fun=median, geom='crossbar')
p2 = VlnPlot(seur, features=c('CNS_ACTIV_score'), pt.size=0) + NoLegend() + stat_summary(fun=median, geom='crossbar')
p3 = VlnPlot(seur, features=c('MG_ACTIV_score'), pt.size=0) + NoLegend() + stat_summary(fun=median, geom='crossbar')
pdf("Sam_Marsh_module_score_per_cluster.pdf", height=10, width=5)
p1+p2+p3
dev.off()

p1 = VlnPlot(seur, features=c('Hom1Final'), pt.size=0, idents=c('Homeostatic1','Homeostatic2'), cols=cols) + NoLegend() + stat_summary(fun=median, geom='crossbar')
p2 = VlnPlot(seur, features=c('Hom2Final'), pt.size=0, idents=c('Homeostatic1','Homeostatic2'), cols=cols) + NoLegend() + stat_summary(fun=median, geom='crossbar')
p3 = VlnPlot(sseur, features='Hom1Final', group.by='Layer', pt.size=0) + stat_summary(fun=median, geom='crossbar') + NoLegend()
p4 = VlnPlot(sseur, features='Hom2Final', group.by='Layer', pt.size=0) + stat_summary(fun=median, geom='crossbar') + NoLegend()
pdf("Hom_module_score.pdf", height=5, width=10)
p1|p2
p3|p4
dev.off()

}

