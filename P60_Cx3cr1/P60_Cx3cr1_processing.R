library(Seurat)
library(stringr)
library(dplyr)
library(pracma)
library(Hmisc)
library(plyr)
source("~/kwanho/src/seurat_tools.R")
options(future.globals.maxSize = 80 * 1024^3)

seur <- readRDS("../analysis/seur_sct_clustered.rds")
seur@reductions$umap <- NULL

seur <- RunTSNE(seur, assay='SCT', reduction='pca', dims=1:30)

pa = readRDS("prev_anno.rds")
seur$prev_anno = NA
seur$prev_anno[names(pa)] = pa

dir.create("figures")

pdf("figures/tsne_prev_anno.pdf")
DimPlot(seur, reduction='tsne', group.by='prev_anno', pt.size=.1, label=T) + NoAxes()
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

# activation score
source("~/microglia/paper/gene_modules_from_Sam_Marsh.R")  # loads variable "modules"
seur <- MyModuleScore(seur, gene.list=modules, save=T, filename="metadata_preQC_P60_Cx3cr1_SCT_Sam_Marsh_module_score.rds")
p1 = VlnPlot(seur, features=c('MG_ID_score'), pt.size=0) + NoLegend() + stat_summary(fun=median, geom='crossbar')
p2 = VlnPlot(seur, features=c('CNS_ACTIV_score'), pt.size=0) + NoLegend() + stat_summary(fun=median, geom='crossbar')
p3 = VlnPlot(seur, features=c('MG_ACTIV_score'), pt.size=0) + NoLegend() + stat_summary(fun=median, geom='crossbar')
pdf("figures/violin_preQC_P60_Cx3cr1_sam_marsh_module_score_per_cluster0.5.pdf", height=10, width=5)
p1+p2+p3
dev.off()

plot_feature2(seur, reduction='tsne', features=c('Lyve1', 'F13a1', 'Flt3', 'Ccr2'), size=3, filename="figures/featureplot_bam.pdf")
plot_feature2(seur, reduction='tsne', features=c('Snap25', 'Sox11', 'Cpa3', 'Mcm3', 'Mcm6', 'Lig1', 'Hells', 'Cald1', 'Vtn'), size=3, filename="figures/featureplot_other.pdf")

# filter LQ clusters
# 3, 5, 12: LQ  1329, 409, 84
# 8: Active
# 9: BAM  205
# 10: Neuronal
seur <- subset(seur, idents=c(3,5,8,9,10,12), invert=T)
print(table(Idents(seur)))

seur <- seur %>%
	SCTransform(vars.to.regress=c("percent.mt", "Replicate")) %>%
	RunPCA() %>%
        RunTSNE(assay='SCT', reduction='pca', dims=1:30) %>%
        FindNeighbors(dims = 1:30) %>%
        FindClusters(resolution=0.5)

saveRDS(seur, "seur_subset_clustered.rds")

pdf("figures/tsne_subset_prev_anno.pdf")
DimPlot(seur, reduction='tsne', group.by='prev_anno', pt.size=.1, label=T) + NoAxes()
dev.off()

pdf("figures/tsne_subset_sct_0.5res.pdf")
DimPlot(seur, reduction='tsne', pt.size=.1, label=T) + NoAxes()
DimPlot(seur, reduction='tsne', pt.size=.1, group.by='Replicate') + NoAxes()
DimPlot(seur, reduction='tsne', pt.size=.1, group.by='Layer') + NoAxes()
plot_feature(seur, reduction='tsne', features=c("nCount_RNA", "nFeature_RNA","percent.mt", "percent.rb"), title="QC on tSNE")
DimPlot(seur, reduction='tsne', pt.size=.1, group.by='Phase') + NoAxes()
dev.off()

pdf(paste0("figures/tsne_subset_sct_split_0.5res.pdf"), height=8, width=16)
DimPlot(seur, reduction='tsne', pt.size=.1, split.by="Replicate", ncol=2) + NoAxes()
DimPlot(seur, reduction='tsne', pt.size=.1, split.by="Replicate", ncol=2, group.by='Layer') + NoAxes()
DimPlot(seur, reduction='tsne', pt.size=.1, split.by="Layer", ncol=3) + NoAxes()
DimPlot(seur, reduction='tsne', pt.size=.1, split.by="dataset", ncol=3) + NoAxes()
dev.off()

pdf(paste0("figures/vlnplot_subset_cluster_specific_qc_0.5res.pdf"), height=20, width=10)
VlnPlot(seur, features=c('nCount_RNA', 'nFeature_RNA', 'percent.mt', 'percent.rb'), pt.size=0, ncol=1) + NoLegend()
dev.off()

library(xlsx)
markers = FindAllMarkers(seur, test.use='MAST', only.pos=T, latent.vars=c('percent.mt', 'nCount_RNA', 'replicate'))
markers <- markers %>% filter(p_val_adj<=0.05)
saveRDS(markers, "markers_MAST_P60_subset.rds")
outfname = "markers_MAST_P60_subset.xlsx"
for (i in 0:max(as.numeric(levels(Idents(seur))))) {
  cur_markers = subset(markers, cluster==i)
  if (nrow(cur_markers)>0) {
  if (!file.exists(outfname)) {write.xlsx(cur_markers[order(-cur_markers$avg_logFC),], file=outfname, sheetName=paste0("cluster",i))}
  else {write.xlsx(cur_markers[order(-cur_markers$avg_logFC),], file=outfname, sheetName=paste0("cluster",i), append=T)}
  }
}

# Annotate dataset
anno = c('Homeostatic1','Homeostatic2','Homeostatic2','Apoe+','Innate Immune','Inflammatory','Ccr1+','Homeostatic3')
names(anno) <- levels(seur)
seur <- RenameIdents(seur, anno)
Idents(seur) <- factor(Idents(seur), levels=c('Homeostatic1','Homeostatic2','Homeostatic3','Apoe+','Ccr1+','Innate Immune','Inflammatory'))
seur$CellType <- Idents(seur)

cols = readRDS("~/microglia/paper/new_colors.rds")

pdf("figures/tsne_final_annotated.pdf")
DimPlot(seur, reduction='tsne',cols=cols) + NoAxes()
dev.off()

pdf(paste0("figures/tsne_split_final_annotated.pdf"), height=8, width=16)
DimPlot(seur, reduction='tsne', pt.size=.1, split.by="Replicate", cols=cols, ncol=2) + NoAxes()
DimPlot(seur, reduction='tsne', pt.size=.1, split.by="Replicate", ncol=2, group.by='Layer') + NoAxes()
DimPlot(seur, reduction='tsne', pt.size=.1, split.by="Layer", ncol=3) + NoAxes()
DimPlot(seur, reduction='tsne', pt.size=.1, split.by="dataset", ncol=3) + NoAxes()
dev.off()

p1 = VlnPlot(seur, features='nCount_RNA', pt.size=0, cols=cols) + NoLegend() + stat_summary(fun=median, geom='crossbar')
p2 = VlnPlot(seur, features='nFeature_RNA', pt.size=0, cols=cols) + NoLegend() + stat_summary(fun=median, geom='crossbar')
p3 = VlnPlot(seur, features='percent.mt', pt.size=0, cols=cols) + NoLegend() + stat_summary(fun=median, geom='crossbar')
pdf(paste0("figures/vlnplot_final_celltype_qc.pdf"), height=10, width=8)
p1+p2+p3
dev.off()

seur$sample = paste(seur$age, seur$Replicate, seur$Layer, sep='_')
p1 = VlnPlot(seur, features='nCount_RNA', pt.size=0, group.by='sample') + NoLegend() + stat_summary(fun=median, geom='crossbar')
p2 = VlnPlot(seur, features='nFeature_RNA', pt.size=0, group.by='sample') + NoLegend() + stat_summary(fun=median, geom='crossbar')
p3 = VlnPlot(seur, features='percent.mt', pt.size=0, group.by='sample') + NoLegend() + stat_summary(fun=median, geom='crossbar')
pdf(paste0("figures/vlnplot_final_P60_sample_qc.pdf"), height=10, width=8)
p1+p2+p3
dev.off()

saveRDS(seur, "seur_final_P60_annotated_072221.rds")

seur$replicate = as.numeric(seur$Replicate != "Rep1")

library(xlsx)
markers = FindAllMarkers(seur, test.use='MAST', only.pos=T, latent.vars=c('percent.mt', 'nCount_RNA', 'replicate'))
markers <- markers %>% filter(p_val_adj<=0.05)
saveRDS(markers, "markers_MAST_P60_final.rds")
outfname = "markers_MAST_P60_final.xlsx"
for (i in levels(seur)) {
  cur_markers = subset(markers, cluster==i)
  if (nrow(cur_markers)>0) {
  if (!file.exists(outfname)) {write.xlsx(cur_markers[order(-cur_markers$avg_logFC),], file=outfname, sheetName=i)}
  else {write.xlsx(cur_markers[order(-cur_markers$avg_logFC),], file=outfname, sheetName=i, append=T)}
  }
}

hom.markers <- FindMarkers(seur, ident.1='Homeostatic1', ident.2='Homeostatic2', test.use='MAST', latent.vars=c('percent.mt','nCount_RNA','replicate'), only.pos=F)
hm = hom.markers %>% filter(p_val_adj<=0.05)
saveRDS(hm, "markers_MAST_P60_final_hom1_v_hom2.rds")

m1 = hm %>% subset(avg_logFC>0)
m2 = hm %>% subset(avg_logFC<0)
m2$avg_logFC = m2$avg_logFC*-1

outfname = "markers_MAST_P60_final_hom1_v_hom2.xlsx"
write.xlsx(m1[order(-m1$avg_logFC),], file=outfname, sheetName="Homeostatic1_markers")
write.xlsx(m2[order(-m2$avg_logFC),], file=outfname, sheetName="Homeostatic2_markers", append=T)

# hom1 v hom2 v apoe v ccr1
hom1.markers <- FindMarkers(seur, ident.1="Homeostatic1", ident.2="Homeostatic2", test.use='MAST', latent.vars=c('percent.mt','nCount_RNA','replicate'), only.pos=T)
hom2.markers <- FindMarkers(seur, ident.1="Homeostatic2", ident.2="Homeostatic1", test.use='MAST', latent.vars=c('percent.mt','nCount_RNA','replicate'), only.pos=T)
apoe.markers <- FindMarkers(seur, ident.1="Apoe+", ident.2=c("Homeostatic1","Homeostatic2"), test.use='MAST', latent.vars=c('percent.mt','nCount_RNA','replicate'), only.pos=T)
ccr1.markers <- FindMarkers(seur, ident.1="Ccr1+", ident.2=c("Homeostatic1","Homeostatic2"), test.use='MAST', latent.vars=c('percent.mt','nCount_RNA','replicate'), only.pos=T)
hom1.markers <- hom1.markers %>% filter(p_val_adj<=0.05)
hom2.markers <- hom2.markers %>% filter(p_val_adj<=0.05)
apoe.markers <- apoe.markers %>% filter(p_val_adj<=0.05)
ccr1.markers <- ccr1.markers %>% filter(p_val_adj<=0.05)

filename = "markers_MAST_P60_Cx3cr1_homs_apoe_ccr1.xlsx"
write.xlsx(hom1.markers[order(-hom1.markers$avg_logFC),], file=filename, sheetName="Homeostatic1")
write.xlsx(hom2.markers[order(-hom2.markers$avg_logFC),], file=filename, sheetName="Homeostatic2", append=T)
write.xlsx(apoe.markers[order(-apoe.markers$avg_logFC),], file=filename, sheetName="Apoe+", append=T)
write.xlsx(ccr1.markers[order(-ccr1.markers$avg_logFC),], file=filename, sheetName="Ccr1+", append=T)




# Cluster composition
library(speckle)
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

pdf("figures/speckle_P60_CellTypeProps.pdf")
mybar + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_fill_manual(values=cols[levels(seur)])
dev.off()

# compute sig for composition difference across all three layers for each cell type
# output one pvalue for each cell type
res = propeller(seur)
print(res)
write.table(res, file=paste0("speckle_P60_CellType_composition_pval.tsv"), sep='\t', quote=F, row.names=T, col.names=NA)

# state score
sig.module = readRDS("~/microglia/analysis_071321/Hom_Signature/nebula/mg_hom_state_signature_genes_final.rds")
seur <- MyModuleScore(seur, gene.list=sig.module, save=T, filename="metadata_P60_mg_state_module_score.rds")
p1 = VlnPlot(seur, features=c('Hom1Final'), pt.size=0, idents=c('Homeostatic1','Homeostatic2'), cols=cols) + NoLegend() + stat_summary(fun=median, geom='crossbar')
p2 = VlnPlot(seur, features=c('Hom2Final'), pt.size=0, idents=c('Homeostatic1','Homeostatic2'), cols=cols) + NoLegend() + stat_summary(fun=median, geom='crossbar')
p3 = VlnPlot(sseur, features='Hom1Final', group.by='Layer', pt.size=0) + stat_summary(fun=median, geom='crossbar') + NoLegend()
p4 = VlnPlot(sseur, features='Hom2Final', group.by='Layer', pt.size=0) + stat_summary(fun=median, geom='crossbar') + NoLegend()
pdf("figures/violin_P60_final_Hom_module_score.pdf", height=5, width=10)
p1|p2
p3|p4
dev.off()


##########################
# Used r v4.0.3 for the below analysis
library(reshape2)
age='P60'
# Composition analysis - Dirichlet mixed effect model
counts <- table(seur$sample, seur$CellType)
percent <- counts/rowSums(counts)
percent <- percent + 1e-6

my.data = as.data.frame.matrix(percent)
my.data$layer = str_split(rownames(my.data), "_", simplify=T)[,3]
my.data$Y = DR_data(percent)
saveRDS(my.data, paste0("data_", age, "_DirichReg.rds"))

fit1 = DirichReg(Y ~ layer, my.data, model='common')
u = summary(fit1)
pvals = u$coef.mat[grep('Intercept', rownames(u$coef.mat), invert=T), 4]
v = names(pvals)
pvals = matrix(pvals, ncol=length(u$varnames))
rownames(pvals) = gsub('layer', '', v[1:nrow(pvals)])
colnames(pvals) = u$varnames
fit1$pvals = pvals
saveRDS(fit1, paste0("fit_", age, "_DirichReg.rds"))

tab <- melt(fit1$pvals)
tab = tab[order(tab$Var1),]
colnames(tab) <- c("Layer", "CellType", "pval")
tab$padj <- p.adjust(tab$pval, method='BH')
tab = tab %>% select(-pval)
tab=reshape(tab, idvar='Layer', timevar='CellType', direction='wide')
write.table(tab, file=paste0(age, "_cluster_composition_pval_DirichReg.tsv"), sep='\t', quote=F, row.names=F, col.names=T)

###################
# module score significance
seur <- AddMetaData(seur, readRDS("metadata_P60_mg_state_module_score.rds"))
sseur <- subset(seur, idents=c("Homeostatic1", "Homeostatic2"))
df = sseur@meta.data[, c('Hom1Final','Hom2Final','CellType','Replicate','sample','Layer')]
df$CellType <- as.factor(as.character(df$CellType))
df$sample <- as.factor(df$sample)

library(lme4)
ret=c()
for (i in 1:2) {
colnames(df)[i] = "col"
res=lmer(col ~ CellType + (1|Replicate) + (1|sample),REML=F,data=df)
colnames(df)[i] = paste0('Hom',i,"Final")
res2 = summary(res)
#print(res2)
saveRDS(res2, paste0("P60_hom", i, "_lmerTest_summary.rds"))
ret = c(ret, res2$coefficients[2,5])
}
ret.padj = p.adjust(ret, method='BH')
names(ret.padj) <- c('Hom1Final', 'Hom2Final')
print(ret.padj)
write.table(ret.padj,file="P60_hom_signature_lmerTest_pvals.tsv",sep="\t",quote=F,row.names=T,col.names=F)

# module score sig by layer
ret=c()
for (i in 1:2) {
colnames(df)[i] = "col"
res=lmer(col ~ Layer + (1|Replicate),REML=F,data=df)
res2 = lmer(col ~ (1|Replicate),REML=F,data=df)
colnames(df)[i] = paste0('Hom',i,"Final")
anov=anova(res,res2)
ret=c(ret,anov$"Pr(>Chisq)"[2])
}
names(ret)=c('Hom1Final', 'Hom2Final')
ret.padj = p.adjust(ret, method="BH")
print(ret.padj)
write.table(ret.padj,file="P60_hom_signature_by_layer_lmm_anova_pvals.tsv",sep="\t",quote=F,row.names=T)


print("DONE")



