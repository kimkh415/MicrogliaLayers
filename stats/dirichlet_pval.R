#dirichlet regression on percent cells per celltype
#based on code from Chris Smillie et. al., UC paper

#library(Seurat)
library(DirichletReg)
library(reshape2)

setwd("/stanley/levin_dr_storage/kwanho/jeff_microglia/fezf2_wt_ko/composition_pval/dirichlet")
seur <- readRDS("../../ml_old/CompleteLayers_P14_Microglia_Cleanest_Reference.rds")
seur$anno <- Idents(seur)

counts <- table(seur$Identity, seur$anno)
percent <- counts/rowSums(counts)

my.data = as.data.frame.matrix(percent)
my.data$layer = c('Layer_1-4', 'Layer_1-4', 'Layer_5', 'Layer_5', 'Layer_6', 'Layer_6')
# Replace extreme value with very low occurrance to avoid error
percent <- percent + 1e-6  # percent + 0.0001
# Prep matrix
my.data$Y = DR_data(percent)

saveRDS(my.data, "prepared_data_for_DirichReg.rds")

# Calculate regression
fit1 = DirichReg(Y ~ layer, my.data, model='common')

# Get p-values
u = summary(fit1)
pvals = u$coef.mat[grep('Intercept', rownames(u$coef.mat), invert=T), 4]
v = names(pvals)
pvals = matrix(pvals, ncol=length(u$varnames))
rownames(pvals) = gsub('condition', '', v[1:nrow(pvals)])
colnames(pvals) = u$varnames
fit1$pvals = pvals
fit1$pvals

saveRDS(fit1, "WT_Mg_DirichReg_fit.rds")

tab <- melt(fit1$pvals)
tab = tab[order(tab$Var1),]
colnames(tab) <- c("Layer", "Cell Type", "pval")
write.table(fit1$pvals, "p_vals.tsv", sep='\t', quote=F, row.names=F, col.names=T)

# Adjust p values using Benjamini Hochberg method
tab$padj <- p.adjust(tab$pval, method='BH')
write.table(tab, "p_adjust.tsv", sep='\t', quote=F, row.names=F, col.names=T)

