library(Seurat)
library(stringr)
options(stringsAsFactors=FALSE)


CleanClusterNames <- function(seur) {
a <- levels(seur)
newA <- str_replace_all(a, "[ :/]", "_")
names(newA) <- a
return(RenameIdents(seur, newA))
}


mouse2human <- function(mgenes, getID=F)
{
require(biomaRt)
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
genesV2 = getLDS(attributes = c("mgi_symbol"), filters = c("mgi_symbol"), values = mgenes , mart = mouse, attributesL = c("hgnc_symbol", "ensembl_gene_id"), martL = human, uniqueRows=T)
hgenes <- as.character(genesV2[, 2])
if (getID) {
hgenes <- as.character(genesV2[, 3])
}
names(hgenes) <- as.character(genesV2[, 1])
#hgenes <- hgenes[which(duplicated(hgenes)==F)]
return(hgenes)
}


outdir <- "/stanley/levin_dr_storage/kwanho/jeff_microglia/analysis_071321/cpdb"

print("Import data!")
mg <- readRDS("~/microglia/paper/seur_P14_final_030221.rds")
mgwtko <- readRDS("/stanley/levin_dr_storage/kwanho/jeff_microglia/analysis_071321/Fezf2_WTKO_scRNA/remove_mito_ribo_genes/seur_final_Fezf2_WTKO_SCT_annotated.rds")
mgko <- subset(mgwtko, Treat=='KO')
mgwt <- subset(mgwtko, Treat=='WT')
obj.n <- readRDS("/stanley/levin_dr_storage/kwanho/jeff_microglia/data/fezf2_wt_ko/ex_neuron/subset_ExN/seur_sct_ExN_subset_clustered.rds")
Idents(obj.n) <- readRDS("/stanley/levin_dr_storage/kwanho/jeff_microglia/data/fezf2_wt_ko/ex_neuron/subset_ExN/metadata_ExNs_FinalSubtypeAssignments.rds")

# Select cell types to run interaction analysis with
print("Prepping data!")
mg <- subset(mg, idents=c("Homeostatic1","Homeostatic2"))
mgwt <- subset(mgwt, idents=c("Homeostatic1","Homeostatic2"))
mgko <- subset(mgko, idents=c("Homeostatic1","Homeostatic2"))

nwt = subset(obj.n, Genotype=="WT")
nko = subset(obj.n, Genotype=="KO")
nwt = subset(nwt, idents=levels(nwt)[grep("KO", levels(nwt), invert=T)])  # exclude KO Mismatch 1 and 2
nko = subset(nko, idents=setdiff(levels(nko), c('L5 CPN 2','L5 CStrPN','L6b Subplate')))  # exclude these with 14, 5 and 1 cells respectively

nwt <- CleanClusterNames(nwt)
nko <- CleanClusterNames(nko)

# merge datasets
ref = merge(mg, y=c(nwt))
wt = merge(mgwt, y=c(nwt))
ko = merge(mgko, y=c(nko))

# extract raw counts
ref.counts = ref@assays$RNA@counts
print(paste0("REF dim: ", dim(ref.counts)))
wt.counts = wt@assays$RNA@counts
print(paste0("WT dim: ", dim(wt.counts)))
ko.counts = ko@assays$RNA@counts
print(paste0("KO dim: ", dim(ko.counts)))

# normalize
print("Normalize!")
ref.norm <- apply(ref.counts, 2, function(x)(x/sum(x))*10000)
wt.norm <- apply(wt.counts, 2, function(x)(x/sum(x))*10000)
ko.norm <- apply(ko.counts, 2, function(x)(x/sum(x))*10000)

# convert mouse gene names to human gene IDs
print("Convert mouse gene symbol to human Ensembl ID")
ref.mgenes <- rownames(ref.norm)
ref.hgenes <- mouse2human(ref.mgenes, getID=T)
rownames(ref.norm) <- ref.hgenes[ref.mgenes]
wt.mgenes <- rownames(wt.norm)
wt.hgenes <- mouse2human(wt.mgenes, getID=T)
rownames(wt.norm) <- wt.hgenes[wt.mgenes]
ko.mgenes <- rownames(ko.norm)
ko.hgenes <- mouse2human(ko.mgenes, getID=T)
rownames(ko.norm) <- ko.hgenes[ko.mgenes]

# drop genes without orthologs
print("Check genes without ortholog!")
#wt.norm <- wt.norm[!(is.na(rownames(wt.norm))), ]
#ko.norm <- ko.norm[!(is.na(rownames(ko.norm))), ]
ref.naGenes = is.na(rownames(ref.norm))
if (length(ref.naGenes)>0) {
cat(paste0("REF: removing ", sum(ref.naGenes), " genes without orthologs...\n"))
ref.norm <- ref.norm[!ref.naGenes, ]
print(dim(ref.norm))
}
wt.naGenes = is.na(rownames(wt.norm))
if (length(wt.naGenes)>0) {
cat(paste0("WT: removing ", sum(wt.naGenes), " genes without orthologs...\n"))
wt.norm <- wt.norm[!wt.naGenes, ]
print(dim(wt.norm))
}
ko.naGenes = is.na(rownames(ko.norm))
if (length(ko.naGenes)>0) {
cat(paste0("WT: removing ", sum(ko.naGenes), " genes without orthologs...\n"))
ko.norm <- ko.norm[!ko.naGenes, ]
print(dim(ko.norm))
}

# Prep output
print("Prep output!")
ref.out = cbind(rownames(ref.norm), ref.norm)
colnames(ref.out)[1] = 'Gene'
ref.meta = as.matrix(Idents(ref))
ref.out2 = cbind(rownames(ref.meta), ref.meta)
colnames(ref.out2) <- c("Cell", "cell_type")
wt.out = cbind(rownames(wt.norm), wt.norm)
colnames(wt.out)[1] = 'Gene'
wt.meta = as.matrix(Idents(wt))
wt.out2 = cbind(rownames(wt.meta), wt.meta)
colnames(wt.out2) <- c("Cell", "cell_type")
ko.out = cbind(rownames(ko.norm), ko.norm)
colnames(ko.out)[1] = 'Gene'
ko.meta = as.matrix(Idents(ko))
ko.out2 = cbind(rownames(ko.meta), ko.meta)
colnames(ko.out2) <- c("Cell", "cell_type")

# Save output
print("Saving!")
if (!dir.exists(outdir)) { dir.create(outdir) }
write.table(ref.out, file=paste0(outdir, "/REF_count.tsv"), quote=FALSE, sep='\t', row.names=FALSE)
write.table(ref.out2, file=paste0(outdir, "/REF_meta.tsv"), quote=FALSE, sep='\t', row.names=FALSE)
write.table(wt.out, file=paste0(outdir, "/WT_count.tsv"), quote=FALSE, sep='\t', row.names=FALSE)
write.table(wt.out2, file=paste0(outdir, "/WT_meta.tsv"), quote=FALSE, sep='\t', row.names=FALSE)
write.table(ko.out, file=paste0(outdir, "/KO_count.tsv"), quote=FALSE, sep='\t', row.names=FALSE)
write.table(ko.out2, file=paste0(outdir, "/KO_meta.tsv"), quote=FALSE, sep='\t', row.names=FALSE)

print("DONE")
