library(Seurat)
library(PLIER)
library(dplyr)

ref <- readRDS("/stanley/levin_dr_storage/kwanho/jeff_microglia/paper/seur_P14_final_030221.rds")

# import genes to use as the prior
new_homs = readRDS("/stanley/levin_dr_storage/kwanho/jeff_microglia/analysis_071321/Hom_Signature/nebula/mg_hom_state_signature_genes_final.rds")
genes.hom1 <- new_homs$Hom1Final
genes.hom2 <- new_homs$Hom2Final

gl.comb = c(genes.hom1, genes.hom2)
prior = matrix(0L, nrow=length(gl.comb), ncol=2)
rownames(prior) <- gl.comb
colnames(prior) <- c("Homeostatic1","Homeostatic2")
prior[genes.hom1, "Homeostatic1"] = 1
prior[genes.hom2, "Homeostatic2"] = 1
saveRDS(prior, "prior_new_hom_genes.rds")

# Prepare gene expression matrix
gex = as.matrix(ref@assays$RNA@data)

cm.genes = commonRows(gex, prior)

prior = prior[cm.genes,]
gex = gex[cm.genes,]
write.table(prior, "prior.tsv", sep='\t', quote=F, row.names=T, col.names=NA)
write.table(gex, "gex.tsv", sep='\t', quote=F, row.names=T, col.names=NA)

# Precompute
gex.svd = svd(gex)  # makes re-running PLIER with different parameters (e.g. k) faster
saveRDS(gex.svd, "svd_of_gex_matrix.rds")
Chat = computeChat(prior)  # a ridge inverse of prior matrix
saveRDS(Chat, "inverse_of_prior_matrix.rds")


prior = as.matrix(read.table("prior.tsv",sep='\t',header=T, row.names=1))
gex = as.matrix(read.table("gex.tsv", sep='\t', header=T, row.names=1))
print(str(prior))
print(str(gex))
print("all rownames match?")
print(all(rownames(prior)==rownames(gex)))

zscore=t(apply(gex, 1, function(X){
return((X-mean(X))/sd(X))
}))
#zscore[is.na(zscore)]=0

gex.svd = readRDS("svd_of_gex_matrix.rds")
Chat = readRDS("inverse_of_prior_matrix.rds")

# Run PLIER
res = PLIER(zscore, prior, scale=F, svdres=gex.svd, Chat=Chat)
saveRDS(res, "plier_result.rds")

library(patchwork)
library(stringr)

# Find the max LV from selected LVs (ones with auc>0.5 and FDR < 0.05)
pval.cutoff = max(res$summary[res$summary[, 5] <= 0.05, 4])  # fifth col is FDR
print(paste0("p value cutoff = ", pval.cutoff, " (this gives FDR < 0.05)"))
selectLV=which(apply(res$Uauc*(res$Up<=pval.cutoff),2,max)>0.7)
print(selectLV)

pdf("PLIER_REF_selected_LVs.pdf")
plotU(res, indexCol=selectLV)
dev.off()

# Visualize gene loading for each significant LV
library(ggplot2)
gene.loading = res$Z
colnames(gene.loading) = paste0("LV", 1:ncol(gene.loading))
gene.loading = gene.loading[,names(selectLV)]
pl=list()
for (i in paste0("LV", rownames(cell.loading[selectLV,]))) {
print(i)
x=sort(gene.loading[, str_extract(i, "LV[0-9]+")], decreasing=T)  # [1:50]
x=sort(x)
y=factor(names(x), levels=names(x))
#print(qplot(x, y, xlab="Gene loading", main=i) + theme(axis.title.y=element_blank()))
df = as.data.frame(x)
df$y = y
pl[[i]] = ggplot(df, aes(x, y)) + geom_point() + theme(axis.title.y=element_blank()) + ggtitle(i)
}

pdf("PLIER_REF_gene_loading_sig_LVs.pdf", height=10, width=10)
wrap_plots(pl, ncol=4)
dev.off()

cols = readRDS("../colors.rds")

# Plot cell loading
cell.loading = res$B
cell.loading = cell.loading[selectLV,]
rownames(cell.loading) <- paste0("LV", str_split(rownames(cell.loading), ",", simplify=T)[,1], "_loading")
lv_scores = as.data.frame(t(cell.loading))
rownames(lv_scores) <- gsub("\\.", "-", rownames(lv_scores))
ref <- AddMetaData(ref, lv_scores)
source("~/kwanho/src/seurat_tools.R")
plot_feature2(ref, reduction='tsne', features=c('LV1_loading','LV4_loading','LV11_loading'), size=4, title="PLIER cell loading", filename="cell_loading_P14Cx3cr1_data.pdf")

p1=VlnPlot(ref, features='LV1_loading', idents=c('Homeostatic1', 'Homeostatic2'), cols=cols, pt.size=0) + stat_summary(fun=median, geom='crossbar') + NoLegend()
#p2=VlnPlot(ref, features='LV2_loading', idents=c('Homeostatic1', 'Homeostatic2'), cols=cols, pt.size=0) + stat_summary(fun=median, geom='crossbar') + NoLegend()
p3=VlnPlot(ref, features='LV4_loading', idents=c('Homeostatic1', 'Homeostatic2'), cols=cols, pt.size=0) + stat_summary(fun=median, geom='crossbar') + NoLegend()
p4=VlnPlot(ref, features='LV11_loading', idents=c('Homeostatic1', 'Homeostatic2'), cols=cols, pt.size=0) + stat_summary(fun=median, geom='crossbar') + NoLegend()
pdf("violin_refdata_PLIER_cell_loading.pdf", height=8, width=8)
(p1|plot_spacer())/(p3|p4)
dev.off()

# Select LVs for hom1 and hom2
umat = res$U[, selectLV]
pmat = res$Up[, selectLV]
umat[pmat>pval.cutoff] = 0
a = rownames(umat)[apply(umat, 2, which.max)]
hom1 = names(selectLV)[which(a=="Homeostatic1")]
hom1 = c("LV1")
hom2 = names(selectLV)[which(a=="Homeostatic2")]

# group LVs - collect genes with top loading scores
lvmodule = list()

for (lv in hom1) {
val = sort(gene.loading[,lv], decreasing=T)
cur_genes=names(val[val>0])
lvmodule[[paste0(lv, "_module")]] = cur_genes
#saveRDS(cur_genes, paste0("Hom1_REF_", lv, "_pos_genes_PLIER.rds"))
}

for (lv in hom2) {
val = sort(gene.loading[,lv], decreasing=T)
cur_genes=names(val[val>0])
lvmodule[[paste0(lv, "_module")]] = cur_genes
#saveRDS(cur_genes, paste0("Hom2_REF_", lv, "_pos_genes_PLIER.rds"))
}

seur <- readRDS("../seur_final_Fezf2_WTKO_SCT_annotated.rds")

ref <- MyModuleScore(ref, gene.list=lvmodule, save=T, "metadata_P14Cx3cr1_PLIER_module_score.rds")
seur <- MyModuleScore(seur, gene.list=lvmodule, save=T, "metadata_Fezf2WTKO_PLIER_module_score.rds")

# add our curated mg state signature module score
sigm = readRDS("../../../Hom_Signature/nebula/mg_hom_state_signature_genes_final.rds")
ref <- MyModuleScore(ref, gene.list=sigm, save=T, "metadata_P14Cx3cr1_our_final_hom_state_score.rds")
seur <- MyModuleScore(seur, gene.list=sigm, save=T, "metadata_Fezf2WTKO_our_final_hom_state_score.rds")

Idents(seur) <- 'Treat'
wt <- subset(seur, idents="WT")
ko <- subset(seur, idents="KO")
Idents(wt) <- 'CellType'
Idents(ko) <- 'CellType'
swt <- subset(wt, idents=c('Homeostatic1', 'Homeostatic2'))
sko <- subset(ko, idents=c('Homeostatic1', 'Homeostatic2'))
sref <- subset(ref, idents=c('Homeostatic1', 'Homeostatic2'))

# Plot individual modules
rplist = list()
rplist2 = list()
wplist = list()
wplist2 = list()
kplist = list()
kplist2 = list()
mins = c(0.28, -0.25, 0.2, -0.1, -0.16)
maxs = c(1.52,0.81, 1.03, 0.87, 0.8)
i=0
for (lv in c(names(sigm), names(lvmodule))) {
i=i+1
rplist[[lv]] = VlnPlot(ref, features=lv, idents=c('Homeostatic1', 'Homeostatic2'), cols=cols, pt.size=0) + stat_summary(fun=median, geom='crossbar') + NoLegend()
rplist2[[lv]] = VlnPlot(sref, features=lv, group.by='Layer', pt.size=0) + stat_summary(fun=median, geom='crossbar') + NoLegend()
wplist[[lv]] = VlnPlot(wt, features=lv, idents=c('Homeostatic1', 'Homeostatic2'), cols=cols, pt.size=0) + stat_summary(fun=median, geom='crossbar') + NoLegend() + scale_y_continuous(limits=c(mins[i],maxs[i]))
wplist2[[lv]] = VlnPlot(swt, features=lv, group.by='Layer', pt.size=0) + stat_summary(fun=median, geom='crossbar') + NoLegend() + scale_y_continuous(limits=c(mins[i],maxs[i]))
kplist[[lv]] = VlnPlot(ko, features=lv, idents=c('Homeostatic1', 'Homeostatic2'), cols=cols, pt.size=0) + stat_summary(fun=median, geom='crossbar') + NoLegend() + scale_y_continuous(limits=c(mins[i],maxs[i]))
kplist2[[lv]] = VlnPlot(sko, features=lv, group.by='Layer', pt.size=0) + stat_summary(fun=median, geom='crossbar') + NoLegend() + scale_y_continuous(limits=c(mins[i],maxs[i]))
}

pdf("violin_P14Cx3cr1_PLIER_modules.pdf", height=7, width=12)
wrap_plots(c(rplist, rplist2), ncol=5)
dev.off()

pdf("violin_Fezf2WT_PLIER_modules_matching_y.pdf", height=7, width=12)
wrap_plots(c(wplist, wplist2), ncol=5)
dev.off()

pdf("violin_Fezf2KO_PLIER_modules_matching_y.pdf", height=7, width=12)
wrap_plots(c(kplist, kplist2), ncol=5)
dev.off()


# compute module sig
# this is done on R v4.0.3
.libPaths(.libPaths()[2])
library(Seurat)
library(dplyr)
library(tidyr)
library(lme4)
library(lmerTest)

ref <- readRDS("/stanley/levin_dr_storage/kwanho/jeff_microglia/paper/seur_P14_final_030221.rds")
seur <- readRDS("../seur_final_Fezf2_WTKO_SCT_annotated.rds")

ref.meta.plier = readRDS("metadata_P14Cx3cr1_PLIER_module_score.rds")
ref.meta.hom = readRDS("metadata_P14Cx3cr1_our_final_hom_state_score.rds")
seur.meta.plier = readRDS("metadata_Fezf2WTKO_PLIER_module_score.rds")
seur.meta.hom = readRDS("metadata_Fezf2WTKO_our_final_hom_state_score.rds")
ref.meta = cbind(ref.meta.hom, ref.meta.plier)
seur.meta = cbind(seur.meta.hom, seur.meta.plier)

ref <- AddMetaData(ref, ref.meta)
seur <- AddMetaData(seur, seur.meta)

sref <- subset(ref, idents=c("Homeostatic1", "Homeostatic2"))
sseur <- subset(seur, idents=c("Homeostatic1", "Homeostatic2"))

dir.create("stats")
setwd("stats")
##########################
# stats on Fezf2 dataset
for (gt in levels(sseur$Treat)) {
print(gt)
obj <- subset(sseur, Treat==gt)
obj <- subset(obj, downsample=min(table(Idents(obj))))

df = obj@meta.data[, c('Hom1Final','Hom2Final','LV1_module','LV4_module','LV11_module','CellType','Replicate','sample','Layer')]
df$CellType <- as.factor(as.character(df$CellType))
df$sample <- as.factor(as.character(df$sample))
df$Replicate <- as.factor(as.character(df$Replicate))

module_names = colnames(df)[1:5]

ret=c()
for (i in 1:5) {
colnames(df)[i] = "col"
res=lmer(col ~ CellType + (1|Replicate) + (1|sample), REML=F, data=df)
colnames(df)[i] = module_names[i]
res2 = summary(res)
#print(res2)
saveRDS(res2, paste0(gt, "_", module_names[i], "_lmerTest_summary.rds"))
ret = c(ret, res2$coefficients[2,5])
}
ret.padj = p.adjust(ret, method='BH')
names(ret.padj) <- module_names
print(ret.padj)
write.table(ret.padj,file=paste0(gt,"_all_module_lmerTest_pvals.tsv"),sep="\t",quote=F,row.names=T,col.names=F)

# sig by layer
ret=c()
for (i in 1:5) {
colnames(df)[i] = "col"
#res=lmer(col ~ Layer + (1|Replicate)+ (1|sample),data=df)
res=lmer(col ~ Layer + (1|Replicate),REML=F,data=df)
res2 = lmer(col ~ (1|Replicate),REML=F,data=df)
colnames(df)[i] = module_names[i]
anov=anova(res,res2)
#print(anov)
ret=c(ret,anov$"Pr(>Chisq)"[2])
}
names(ret)=module_names
ret.padj = p.adjust(ret, method="BH")
print(ret.padj)
write.table(ret.padj,file=paste0(gt,"_all_module_by_layer_lmm_anova_pvals.tsv"),sep="\t",quote=F,row.names=T)

}
############################
# layer module score diff between the genotype
Idents(sseur) <- 'Treat'
obj <- subset(sseur, downsample=min(table(Idents(sseur))))
df = obj@meta.data[, c('Hom1Final','Hom2Final','LV1_module','LV4_module','LV11_module','CellType','Replicate','sample','Layer', 'Treat')]
df$CellType <- as.factor(as.character(df$CellType))
module_names = colnames(df)[1:5]
ret=c()
for (i in 1:5) {
colnames(df)[i] = "col"
res=lmer(col ~ Layer + Treat + (1|Replicate),REML=F,data=df)
res2 = lmer(col ~ Layer + (1|Replicate),REML=F,data=df)
colnames(df)[i] = module_names[i]
anov=anova(res,res2)
ret=c(ret,anov$"Pr(>Chisq)"[2])
}
names(ret)=module_names
ret.padj = p.adjust(ret, method="BH")
print(ret.padj)
write.table(ret.padj,file="Fezf2WTKO_layer_module_score_diff_by_Treat_lmm_anova_pvals.tsv",sep="\t",quote=F,row.names=T)

# celltype module score diff between the genotype
ret=c()
for (i in 1:5) {
colnames(df)[i] = "col"
res=lmer(col ~ CellType + Treat + (1|Replicate),REML=F,data=df)
res2 = lmer(col ~ CellType + (1|Replicate),REML=F,data=df)
colnames(df)[i] = module_names[i]
anov=anova(res,res2)
ret=c(ret,anov$"Pr(>Chisq)"[2])
}
names(ret)=module_names
ret.padj = p.adjust(ret, method="BH")
print(ret.padj)
write.table(ret.padj,file="Fezf2WTKO_CellType_module_score_diff_by_Treat_lmm_anova_pvals.tsv",sep="\t",quote=F,row.names=T)


##################################
# REF data module score stats
obj <- subset(sref, downsample=min(table(Idents(sref))))

df = obj@meta.data[, c('Hom1Final','Hom2Final','LV1_module','LV4_module','LV11_module','CellType','Replicate','sample','Layer')]
df$CellType <- as.factor(as.character(df$CellType))
df$sample <- as.factor(df$sample)
df$Replicate <- as.factor(df$Replicate)

module_names = colnames(df)[1:5]

# 
ret=c()
for (i in 1:5) {
colnames(df)[i] = "col"
res=lmer(col ~ CellType + (1|sample), REML=F, data=df)
colnames(df)[i] = module_names[i]
res2 = summary(res)
#print(res2)
saveRDS(res2, paste0("REF_", module_names[i], "_lmerTest_summary.rds"))
ret = c(ret, res2$coefficients[2,5])
}
ret.padj = p.adjust(ret, method='BH')
names(ret.padj) <- module_names
print(ret.padj)
write.table(ret.padj,file="REF_all_module_lmerTest_pvals.tsv",sep="\t",quote=F,row.names=T,col.names=F)

# sig by layer
ret=c()
for (i in 1:5) {
colnames(df)[i] = "col"
res=lmer(col ~ Layer + (1|Replicate),REML=F,data=df)
res2 = lmer(col ~ (1|Replicate),REML=F,data=df)
colnames(df)[i] = module_names[i]
anov=anova(res,res2)
#print(anov)
ret=c(ret,anov$"Pr(>Chisq)"[2])
}
names(ret)=module_names
ret.padj = p.adjust(ret, method="BH")
print(ret.padj)
write.table(ret.padj,file="REF_all_module_by_layer_lmm_anova_pvals.tsv",sep="\t",quote=F,row.names=T)



print("DONE")



