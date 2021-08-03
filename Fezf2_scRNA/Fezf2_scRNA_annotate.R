library(Seurat)
library(dplyr)
library(ggplot2)
source("~/kwanho/src/seurat_tools.R")
options(future.globals.maxSize = 80 * 1024^3)

seur <- readRDS("seur_final_Fezf2_WTKO_scRNA_sct_clustered.rds")

anno = c('Homeostatic1','Homeostatic2','Apoe+','Ccr1+','Innate Immune','Inflammatory','Homeostatic3','Dividing','Dividing')
names(anno) <- levels(seur)
seur <- RenameIdents(seur, anno)

cols = readRDS("../colors.rds")
cols = cols[-8]  # remove activated
Idents(seur) <- factor(Idents(seur), levels=names(cols))

saveRDS(Idents(seur), "metadata_final_Fezf2WTKO_scRNA_SCT_CellType.rds")
saveRDS(seur, "seur_final_Fezf2_WTKO_SCT_annotated.rds")

pdf("figures/tsne_split_final_colored_by_sample.pdf", height=12, width=10)
DimPlot(seur, reduction='tsne',split.by='sample', ncol=3, cols=cols)
dev.off()

seur$treat_rep <- paste(seur$Treat, seur$Replicate, sep='_')
seur$treat_layer <- paste(seur$Treat, seur$Layer, sep='_')
pdf("figures/tsne_final_colored_Fezf2_WTKO_scRNA_sct.pdf")
DimPlot(seur, reduction='tsne', pt.size=1, cols=cols) + NoAxes()
DimPlot(seur, reduction='tsne', pt.size=1, group.by='treat_layer') + NoAxes()
DimPlot(seur, reduction='tsne', pt.size=1, group.by='treat_rep') + NoAxes()
plot_feature(seur, reduction='tsne', features=c("nCount_RNA", "nFeature_RNA","percent.mt", "percent.rb"), title="QC")
DimPlot(seur, reduction='tsne', pt.size=1, group.by='Phase') + NoAxes()
dev.off()

vp1 = VlnPlot(seur, features='nCount_RNA', pt.size=0, ncol=1, cols=cols) + NoLegend() + stat_summary(fun=median, geom='crossbar')
vp2 = VlnPlot(seur, features='nFeature_RNA', pt.size=0, ncol=1, cols=cols) + NoLegend() + stat_summary(fun=median, geom='crossbar')
vp3 = VlnPlot(seur, features='percent.mt', pt.size=0, ncol=1, cols=cols) + NoLegend() + stat_summary(fun=median, geom='crossbar')
pdf("figures/violin_final_Fezf2_WTKO_scRNA_cluster_specific_qc.pdf", height=15, width=10)
vp1+vp2+vp3
dev.off()

# Cluster composition
library(speckle)
seur$group = paste(seur$Treat, seur$Layer, sep='_')
clusters=Idents(seur)
sample = seur$sample
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

pdf("figures/speckle_final_CellTypeProps.pdf")
mybar + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_fill_manual(values=cols[levels(seur)])
dev.off()

# compute sig for composition difference across all three layers for each cell type
# output one pvalue for each cell type
for (gt in c('WT','KO')) {
print(gt)
sseur <- subset(seur, Treat==gt)
res = propeller(sseur)
print(res)
write.table(res, file=paste0("speckle_",gt,"_CellType_composition_pval.tsv"), sep='\t', quote=F, row.names=T, col.names=NA)
}

for (lay in levels(seur$Layer)) {
print(lay)
sseur <- subset(seur, Layer==lay)
print(table(sseur$group))
print(table(sseur$sample))
res = propeller(sseur)
print(res)
write.table(res, file=paste0("speckle_",lay,"_CellType_composition_pval.tsv"), sep='\t', quote=F, row.names=T, col.names=NA)
}


library(xlsx)
markers = FindAllMarkers(seur, test.use='MAST', only.pos=T, latent.vars=c('percent.mt', 'nCount_RNA', 'replicate'))
markers <- markers %>% filter(p_val_adj<=0.05)
saveRDS(markers, "markers_MAST_final_CellType_Fezf2_WTKO_scRNA.rds")
outfname = "markers_MAST_final_CellType_Fezf2_WTKO_scRNA.xlsx"
for (i in levels(seur)) {
  cur_markers = subset(markers, cluster==i)
  if (nrow(cur_markers)>0) {
  if (!file.exists(outfname)) {write.xlsx(cur_markers[order(-cur_markers$avg_logFC),], file=outfname, sheetName=paste0("cluster",i))}
  else {write.xlsx(cur_markers[order(-cur_markers$avg_logFC),], file=outfname, sheetName=paste0("cluster",i), append=T)}
  }
}

#"Homeostatic1","Homeostatic2","Apoe+","Ccr1+"
hom1.markers <- FindMarkers(seur, ident.1="Homeostatic1", ident.2="Homeostatic2", test.use='MAST', latent.vars=c('percent.mt','nCount_RNA','replicate'), only.pos=T, logfc.threshold=0)
hom2.markers <- FindMarkers(seur, ident.1="Homeostatic2", ident.2="Homeostatic1", test.use='MAST', latent.vars=c('percent.mt','nCount_RNA','replicate'), only.pos=T, logfc.threshold=0)
apoe.markers <- FindMarkers(seur, ident.1="Apoe+", ident.2=c("Homeostatic1","Homeostatic2"), test.use='MAST', latent.vars=c('percent.mt','nCount_RNA','replicate'), only.pos=T, logfc.threshold=0)
ccr1.markers <- FindMarkers(seur, ident.1="Ccr1+", ident.2=c("Homeostatic1","Homeostatic2"), test.use='MAST', latent.vars=c('percent.mt','nCount_RNA','replicate'), only.pos=T, logfc.threshold=0)
hom1.markers <- hom1.markers %>% filter(p_val_adj<=0.05)
hom2.markers <- hom2.markers %>% filter(p_val_adj<=0.05)
apoe.markers <- apoe.markers %>% filter(p_val_adj<=0.05)
ccr1.markers <- ccr1.markers %>% filter(p_val_adj<=0.05)
filename = "markers_MAST_final_hom1_v_hom2_logFC0.xlsx"
write.xlsx(hom1.markers[order(-hom1.markers$avg_logFC),], file=filename, sheetName="Homeostatic1")
write.xlsx(hom2.markers[order(-hom2.markers$avg_logFC),], file=filename, sheetName="Homeostatic2", append=T)

sseur <- subset(seur, idents="Homeostatic1")
Idents(sseur) <- 'Treat'
hom1.markers <- FindMarkers(sseur, ident.1="WT", ident.2="KO", test.use='MAST', latent.vars=c('percent.mt','nCount_RNA','replicate'), logfc.threshold=0)
wt1.markers <- hom1.markers %>% filter(avg_logFC>0) %>% filter(p_val_adj<=0.05)
ko1.markers <- hom1.markers %>% filter(avg_logFC<0) %>% filter(p_val_adj<=0.05)
ko1.markers$avg_logFC = ko1.markers$avg_logFC * -1
sseur <- subset(seur, idents="Homeostatic2")
Idents(sseur) <- 'Treat'
hom2.markers <- FindMarkers(sseur, ident.1="WT", ident.2="KO", test.use='MAST', latent.vars=c('percent.mt','nCount_RNA','replicate'), logfc.threshold=0)
wt2.markers <- hom2.markers %>% filter(avg_logFC>0) %>% filter(p_val_adj<=0.05)
ko2.markers <- hom2.markers %>% filter(avg_logFC<0) %>% filter(p_val_adj<=0.05)
ko2.markers$avg_logFC = ko2.markers$avg_logFC * -1
filename = "markers_MAST_final_diff_by_genotype_logFC0.xlsx"
write.xlsx(wt1.markers[order(-wt1.markers$avg_logFC),], file=filename, sheetName="WT_HOM1")
write.xlsx(ko1.markers[order(-ko1.markers$avg_logFC),], file=filename, sheetName="KO_HOM1", append=T)
write.xlsx(wt2.markers[order(-wt2.markers$avg_logFC),], file=filename, sheetName="WT_HOM2", append=T)
write.xlsx(ko2.markers[order(-ko2.markers$avg_logFC),], file=filename, sheetName="KO_HOM2", append=T)



for (gt in levels(seur$Treat)) {
print(gt)
sseur <- subset(seur, Treat==gt)
print(table(sseur$Treat))
hom.markers <- FindMarkers(sseur, ident.1="Homeostatic1", ident.2="Homeostatic2", test.use='MAST', latent.vars=c('percent.mt','nCount_RNA','replicate'), only.pos=F)
hm = hom.markers %>% filter(p_val_adj<=0.05)
print(dim(hm))
saveRDS(hm, paste0("markers_MAST_",gt,"_only_Hom1_v_Hom2.rds"))
m1 = hm %>% subset(avg_logFC>0)
m2 = hm %>% subset(avg_logFC<0)
m2$avg_logFC = m2$avg_logFC*-1

outfname = paste0("markers_MAST_",gt,"_only_Hom1_v_Hom2.xlsx")
write.xlsx(m1[order(-m1$avg_logFC),], file=outfname, sheetName="Homeostatic1")
write.xlsx(m2[order(-m2$avg_logFC),], file=outfname, sheetName="Homeostatic2", append=T)

}

# Mg state module score
newm = readRDS("../../Hom_Signature/gene_module_Homs_new_MAST.rds")
oldm = readRDS("../../Hom_Signature/gene_module_Homs_old_MAST.rds")
modules = c(newm, oldm)
names(modules) = c('Hom1new','Hom2new','Hom1old','Hom2old')
seur <- MyModuleScore(seur, gene.list=modules, save=T, filename="metadata_final_old_and_new_Hom_module_score.rds")
p1 = VlnPlot(seur, features=c('Hom1new'), pt.size=0, idents=c('Homeostatic1','Homeostatic2'), split.by='Treat') + NoLegend() + stat_summary(fun=median, geom='crossbar')
p2 = VlnPlot(seur, features=c('Hom2new'), pt.size=0, idents=c('Homeostatic1','Homeostatic2'), split.by='Treat') + NoLegend() + stat_summary(fun=median, geom='crossbar')
p3 = VlnPlot(seur, features=c('Hom1old'), pt.size=0, idents=c('Homeostatic1','Homeostatic2'), split.by='Treat') + NoLegend() + stat_summary(fun=median, geom='crossbar')
p4 = VlnPlot(seur, features=c('Hom2old'), pt.size=0, idents=c('Homeostatic1','Homeostatic2'), split.by='Treat') + stat_summary(fun=median, geom='crossbar')
library(patchwork)
pdf("figures/violin_final_Hom_modules.pdf", height=10, width=10)
p1+p2|p3+p4
dev.off()


# activation score
source("~/microglia/paper/gene_modules_from_Sam_Marsh.R")  # loads variable "modules"
seur <- MyModuleScore(seur, gene.list=modules, save=T, filename="metadata_final_Fezf2WTKO_scRNA_SCT_Sam_Marsh_module_score.rds")
p1 = VlnPlot(seur, features=c('MG_ID_score'), pt.size=0) + NoLegend() + stat_summary(fun=median, geom='crossbar')
p2 = VlnPlot(seur, features=c('CNS_ACTIV_score'), pt.size=0) + NoLegend() + stat_summary(fun=median, geom='crossbar')
p3 = VlnPlot(seur, features=c('MG_ACTIV_score'), pt.size=0) + NoLegend() + stat_summary(fun=median, geom='crossbar')
pdf("figures/violin_final_sam_marsh_module_score_per_CellType.pdf", height=10, width=5)
p1+p2+p3
dev.off()



# Composition analysis - Dirichlet mixed effect model
library(reshape2)
for (gt in c("WT", "KO")) {
sseur <- subset(seur, Treat==gt)
sseur$sample <- as.character(sseur$sample)

counts <- table(sseur$sample, sseur$CellType)
percent <- counts/rowSums(counts)
percent <- percent + 1e-6

my.data = as.data.frame.matrix(percent)
my.data$layer = str_split(rownames(my.data), "_", simplify=T)[,1]
my.data$Y = DR_data(percent)
saveRDS(my.data, paste0("data_", gt, "_DirichReg.rds"))

fit1 = DirichReg(Y ~ layer, my.data, model='common')
u = summary(fit1)
pvals = u$coef.mat[grep('Intercept', rownames(u$coef.mat), invert=T), 4]
v = names(pvals)
pvals = matrix(pvals, ncol=length(u$varnames))
rownames(pvals) = gsub('layer', '', v[1:nrow(pvals)])
colnames(pvals) = u$varnames
fit1$pvals = pvals
saveRDS(fit1, paste0("fit_", gt, "_DirichReg.rds"))

tab <- melt(fit1$pvals)
tab = tab[order(tab$Var1),]
colnames(tab) <- c("Layer", "CellType", "pval")
tab$padj <- p.adjust(tab$pval, method='BH')
tab = tab %>% select(-pval)
tab=reshape(tab, idvar='Layer', timevar='CellType', direction='wide')
write.table(tab, file=paste0("Fezf2_scRNA_", gt, "_cluster_composition_pval_DirichReg.tsv"), sep='\t', quote=F, row.names=F, col.names=T)
}


library(lme4)
library(lmerTest)
# Mertk difference pval. Within each Mg state by genotype
sseur <- subset(seur, idents=c("Homeostatic1","Homeostatic2"))
df = sseur@meta.data[,c('sample','Replicate','CellType','Treat')]
df$CellType <- factor(as.character(df$CellType))

df$Mertk <- sseur@assays$SCT@data["Mertk",]

ret = c()
for (mgstate in levels(df$CellType)) {
print(mgstate)
sdf <- df %>% filter(CellType == mgstate)
print(table(sdf$CellType))
fit1 = lmer(Mertk ~ Treat + (1|Replicate), data = sdf)
fit2 = lmer(Mertk ~ (1|Replicate), data = sdf)
anov = anova(fit1, fit2)
print(anov)
ret=c(ret, anov$"Pr(>Chisq)"[2])
}
names(ret) = levels(df$CellType)
ret.padj = p.adjust(ret, method="BH")
print(ret.padj)
write.table(ret.padj,file="Mertk_exp_diff_by_genotype_within_each_Mgstate_lmm_anova_pvals.tsv",sep="\t",quote=F,row.names=T)

# Mertk diff between Hom1 and Hom2 (irrespective of genotype)
fit1 = lmer(Mertk ~ CellType + (1|Replicate), data=df)
summary(fit1)




print("DONE")



