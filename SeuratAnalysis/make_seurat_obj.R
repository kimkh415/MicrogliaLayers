library(Seurat)
library(stringr)
library(dplyr)
library(pracma)
library(Hmisc)
library(plyr)
options(future.globals.maxSize = 80 * 1024^3)


# takes as input 
# 1. Path to the 10X CellRanger counts output directory
# 2. Output directory
args = commandArgs(trailingOnly=T)


setwd(args[1])

# Locate counts
samples = list.files()
print(samples)
dirlist_10x = paste0(samples, "/filtered_feature_bc_matrix/")
print(dirlist_10x)
stopifnot(all(sapply(dirlist_10x, FUN=dir.exists)))

# Read counts mtx
mtx.list = lapply(dirlist_10x, FUN=Read10X)
for (i in 1:length(mtx.list)) {
colnames(mtx.list[[i]]) = paste0(samples[i], '_', colnames(mtx.list[[i]]))
}

# Create Seurat objects
print("creating seurate objects")
seur.list = lapply(mtx.list, FUN=CreateSeuratObject, min.cells=3, min.features=200)
names(seur.list) <- samples
print(seur.list)

for (i in 1:length(seur.list)) {
seur.list[[i]] <- seur.list[[i]] %>%
        PercentageFeatureSet(pattern="^mt-", col.name="percent.mt") %>%
        PercentageFeatureSet(pattern="^Rp[sl]", col.name="percent.rb")%>%
        NormalizeData() %>%
        FindVariableFeatures() %>%
        ScaleData(vars.to.regress=c('percent.mt', 'percent.rb'))
}

# Add cell cycle score
s.genes <- upFirst(tolower(cc.genes$s.genes))
g2m.genes <- upFirst(tolower(cc.genes$g2m.genes))
for (i in 1:length(seur.list)) {
seur.list[[i]] <- CellCycleScoring(seur.list[[i]], s.features = s.genes, g2m.features = g2m.genes, set.ident=F)
}

# Filtering cells
for (i in 1:length(seur.list)) {
seur.list[[i]] <- subset(seur.list[[i]], percent.mt<10 & nFeature_RNA>500)  # 4 cells with >1e5 nCount
}

# Merge datasets
seur = merge(seur.list[[1]], y=seur.list[2:length(seur.list)])


setwd(args[2])

# Noramlize RNA assay
seur <- NormalizeData(seur)

saveRDS(seur, "seur_merged_logNorm.rds")


