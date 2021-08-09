#Identify Homeostatic microglia (Mg) marker genes considering sample to sample variability
# This code was run on R version 4.0.3

library(nebula)
library(Seurat)
library(ggplot2)
library(dplyr)
library(tibble)


# Path to the Seurat object that has 'Homeostatic1' and 'Homeostatic2' Mg identified and lableled under 'CellType' column and individual sample labeled under 'sample' column of the metadata.
args = commandArgs(trailingOnly=T)

seur.all <- readRDS(args[1])
Idents(seur.all) <- 'CellType'

seur <- subset(seur.all, idents=c('Homeostatic1','Homeostatic2'))

counts = seur@assays$RNA@counts

meta = seur@meta.data[,c('CellType','sample')]
meta$CellType = as.character(meta$CellType)

df = model.matrix(~CellType, data=meta)
# sample as a random effect
# uses negative binomial gamma mixed model (NBGMM) by default
re = nebula(counts, as.character(meta$sample), pred=df)
saveRDS(re, "nebula_nbgmm_res.rds")

res = re$summary
res <- res %>% filter(p_CellTypeHomeostatic2<0.05)
res <- res[order(res$p_CellTypeHomeostatic2),]
# remove mitochondrial and ribosomal genes from the result
res <- res[grep("^mt-|^Rp[sl]", res$gene, invert=T),]

m1 = res %>% subset(logFC_CellTypeHomeostatic2<0)
m2 = res %>% subset(logFC_CellTypeHomeostatic2>0)

outfname = "markers_nebula_Homs.xlsx"
write.xlsx(m1, file=outfname, sheetName="Homeostatic1")
write.xlsx(m2, file=outfname, sheetName="Homeostatic2", append=T)

