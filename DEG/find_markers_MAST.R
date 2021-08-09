library(Seurat)
library(dplyr)
library(xlsx)


# Finds MAST DEG for each cluster identified in 'seur'
# Input 'seur' is a Seurat object that has percentage mitochondrial gene expression under 'percent.mt',
# number of UMIs under 'nCount_RNA' and replicate information (numeric) under 'replicate' column of the metadata.
RunMAST <- function(seur) {

markers <- FindAllMarkers(seur, test.use='MAST', latent.vars=c('percent.mt','nCount_RNA','replicate'), only.pos=T)
m <- markers %>% filter(p_val_adj<=0.05)

saveRDS(m, "markers_MAST_all.rds")
outfname = "markers_MAST_all.xlsx"

for (i in levels(m$cluster)) {
  cur_markers = subset(m, cluster==i)
  if (nrow(cur_markers)>0) {
  if (i==levels(m$cluster)[2]) {write.xlsx(cur_markers[order(-cur_markers$avg_logFC),], file=outfname, sheetName=paste0(i,"_markers"))}
  else {write.xlsx(cur_markers[order(-cur_markers$avg_logFC),], file=outfname, sheetName=paste0(i,"_markers"), append=T)}
  }
}

return(m)
}


# Finds MAST DEG between 'Homeostatic1' and 'Homeostatic2' microglia
# Input same as above with the addition of fold-change cutoff.
RunMAST_hom <- function(seur, cutoff=.15) {

hom.markers <- FindMarkers(seur, ident.1='Homeostatic1', ident.2='Homeostatic2', test.use='MAST', latent.vars=c('percent.mt','nCount_RNA','replicate'), only.pos=F)
hm = hom.markers %>% filter(p_val_adj<=0.05)

saveRDS(hm, "markers_MAST_hom1_v_hom2.rds")

m1 = hm %>% subset(avg_logFC>cutoff)
m2 = hm %>% subset(avg_logFC<cutoff)
m2$avg_logFC = m2$avg_logFC*-1

outfname = "markers_MAST_hom1_v_hom2.xlsx"
write.xlsx(m1[order(-m1$avg_logFC),], file=outfname, sheetName="Homeostatic1_markers")
write.xlsx(m2[order(-m2$avg_logFC),], file=outfname, sheetName="Homeostatic2_markers", append=T)

}


# Finds MAST DEG for each cluster identified in 'seur'
# This function assums that the input Seurat object is annotated and contains the following clusters:
# "Homeostatic1", "Homeostatic2", "Apoe+" and "Ccr1+"
RunMAST2 <- function(seur) {
hom1.markers <- FindMarkers(seur, ident.1="Homeostatic1", ident.2="Homeostatic2", test.use='MAST', latent.vars=c('percent.mt','nCount_RNA','replicate'), only.pos=T)
hom2.markers <- FindMarkers(seur, ident.1="Homeostatic2", ident.2="Homeostatic1", test.use='MAST', latent.vars=c('percent.mt','nCount_RNA','replicate'), only.pos=T)
apoe.markers <- FindMarkers(seur, ident.1="Apoe+", ident.2=c("Homeostatic1","Homeostatic2"), test.use='MAST', latent.vars=c('percent.mt','nCount_RNA','replicate'), only.pos=T)
ccr1.markers <- FindMarkers(seur, ident.1="Ccr1+", ident.2=c("Homeostatic1","Homeostatic2"), test.use='MAST', latent.vars=c('percent.mt','nCount_RNA','replicate'), only.pos=T)

hom1.markers <- hom1.markers %>% filter(p_val_adj<=0.05)
hom2.markers <- hom2.markers %>% filter(p_val_adj<=0.05)
apoe.markers <- apoe.markers %>% filter(p_val_adj<=0.05)
ccr1.markers <- ccr1.markers %>% filter(p_val_adj<=0.05)

filename = "markers_MAST_P14_Cx3cr1_homs_apoe_ccr1.xlsx"
write.xlsx(hom1.markers[order(-hom1.markers$avg_logFC),], file=filename, sheetName="Homeostatic1")
write.xlsx(hom2.markers[order(-hom2.markers$avg_logFC),], file=filename, sheetName="Homeostatic2", append=T)
write.xlsx(apoe.markers[order(-apoe.markers$avg_logFC),], file=filename, sheetName="Apoe+", append=T)
write.xlsx(ccr1.markers[order(-ccr1.markers$avg_logFC),], file=filename, sheetName="Ccr1+", append=T)
}



