library(Seurat)
library(stringr)
library(dplyr)

options(stringsAsfactors=F)

setwd('/stanley/levin_dr_storage/kwanho/jeff_microglia/data/MerFISH')

print("Load data!")
cb = read.csv("VZG123_codebook.csv")
gene_names = as.character(cb$name)[1:115]

flist = list.files(pattern="^Run")
print(flist)
samples = substr(flist, 1, 15)
tab_list = lapply(flist, FUN=read.csv, header=F)

counts_list = c()
z = c()
area = c()
third_col = c()
spatial_x = c()
spatial_y = c()

dapi_tab = read.table("/home/unix/kwanho/microglia/data/MerFISH/dapi_mosaic/image_size.tsv", sep='\t', header=F)
dapi_tab$V1 = gsub('Dapi.tif','',dapi_tab$V1)
colnames(dapi_tab) <- c('sample','nRow','nCol')
rownames(dapi_tab) <- dapi_tab$sample
dapi_tab <- dapi_tab %>% select(-sample)

ind2sub <- function(ind, nrow, ncol) {
coord.x = ((ind-1)%%nrow)+1
coord.y = floor((ind-1)/nrow) + 1
return(c(coord.x, coord.y))
}

for (i in seq(1,length(flist)-1,2)) {
  print(samples[i])
  
  t1 = tab_list[[i]]
  t2 = tab_list[[i+1]]
  tab = rbind(t1, t2)
  colnames(tab) <- c('Z', 'area', 'thirdCol', gene_names)
  rownames(tab) <- paste0(samples[i], '_', rownames(tab))

  z = c(z, as.numeric(tab$Z))
  area = c(area, as.numeric(tab$area))
  pos_ind = as.numeric(tab$thirdCol)
  third_col = c(third_col, pos_ind)
  
  n_pixels_row = dapi_tab[samples[i], 1]
  n_pixels_col = dapi_tab[samples[i], 2]
  for (j in 1:length(pos_ind)) {
    coords <- ind2sub(pos_ind[j], n_pixels_row, n_pixels_col)
    spatial_x <- c(spatial_x, as.numeric(coords[2]))
    spatial_y <- c(spatial_y, as.numeric(coords[1]))
  }

  tab = t(tab[,4:ncol(tab)])
  print(dim(tab))
  counts_list[[samples[i]]] = tab
}

# combine into one count matrix
comb_counts = do.call(cbind, counts_list)
print(dim(comb_counts))
print(length(area))
meta = data.frame(z=z, area=area, location=third_col, spatial_x=spatial_x, spatial_y=spatial_y)
rownames(meta) = colnames(comb_counts)

seur <- CreateSeuratObject(comb_counts, meta.data=meta)
seur$run <- substr(colnames(seur), 1, 4)
seur$slice <- substr(colnames(seur), 5, 10)
seur$side <- substr(colnames(seur), 11, 15)
seur$sample <- str_split(colnames(seur), pattern="_", simplify=T)[,1]

emb <- matrix(cbind(seur$spatial_x, seur$spatial_y*-1), ncol=2, dimnames=list(colnames(seur), c('spatial_1','spatial_2')))
seur[["spatial"]] <- CreateDimReducObject(embeddings = emb, key = "spatial_", assay = DefaultAssay(seur))

saveRDS(seur, "analysis/final/seur_init.rds")
