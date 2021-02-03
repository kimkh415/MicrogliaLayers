library(Seurat)
library(stringr)
options(stringsAsFactors=FALSE)


# Clean up cluster names by replacing " " and ":" with "_"
CleanClusterNames <- function(seur) {
a <- levels(seur)
newA <- str_replace_all(a, "[ :]", "_")
names(newA) <- a
return(RenameIdents(seur, newA))
}


# Orthologous genes from Ensembl Biomart
MouseToHuman<-function(mGenes)
{
ref <- read.table("/stanley/levin_dr/kwanho/ref/mouse_human_orthologous_genes.tsv", sep='\t', header=TRUE)
idx <- match(mGenes, ref$Gene.name)
res <- as.character(ref$Human.gene.stable.ID[idx])
names(res) <- ref$Human.gene.name[idx]
return(res)
}

args = commandArgs(trailingOnly=T)
outdir <- args[1]

obj.mg <- readRDS(args[2])  # Cx3cr1 egfp microglia dataset
obj.n <- readRDS(args[3])  # dataset with neurons
obj.mgko <- readRDS(args[4])  # microglia from Fezf2 KO dataset

obj.mg <- CleanClusterNames(obj.mg)
obj.n <- CleanClusterNames(obj.n)
obj.mgko <- CleanClusterNames(obj.mgko)

# Distinguish Upper (1-4) and Lower (5-6) layer microglia from the Fezf2 KO dataset as Mg_Hom_1 and Mg_Hom_2 respectively
obj.mgko <- subset(obj.mgko, idents=c("General_Microglia_1", "General_Microglia_2"))
Idents(obj.mgko) <- obj.mgko$Layer
a <- c("Mg_Hom_1_KO", "Mg_Hom_2_KO", "Mg_Hom_2_KO")
names(a) <- c("Layer 1-4", "Layer 5", "Layer 6")
obj.mgko <- RenameIdents(obj.mgko, a)

obj.mg <- subset(obj.mg, idents=c("Upper_Layer", "Lower_Layer"))
a <- c("Mg_Hom_1", "Mg_Hom_2")
names(a) <- c("Upper_Layer", "Lower_Layer")
obj.mg <- RenameIdents(obj.mg, a)


# merge datasets
print("merging data")
temp_obj <- merge(obj.mg, y=c(obj.n, obj.mgko))

# extract raw count data
count_raw <- temp_obj[['RNA']]@counts

print("combined matrix dimension")
print(dim(count_raw))

# normalize
count_norm <- apply(count_raw, 2, function(x)(x/sum(x))*10000)

# convert mouse gene names to human gene IDs (CPDB only supports human interaction database)
mgenes <- rownames(count_norm)
hgenes <- MouseToHuman(mgenes)
rownames(count_norm) <- hgenes

# drop genes without orthologs
count_norm <- count_norm[!(is.na(rownames(count_norm))), ]

# create output directory
system(paste0("mkdir -p ", outdir))

# output count matrix
out <- cbind(rownames(count_norm), count_norm)
colnames(out)[1] <- "Gene"
write.table(out, file=paste0(outdir, "/combined_count.tsv"), quote=FALSE, sep='\t', row.names=FALSE)

meta_data <- as.matrix(Idents(temp_obj))
out2 <- cbind(rownames(meta_data), meta_data)
colnames(out2) <- c("Cell", "cell_type")
write.table(out2, file=paste0(outdir, "/combined_meta.tsv"), quote=FALSE, sep='\t', row.names=FALSE)

print("DONE")
