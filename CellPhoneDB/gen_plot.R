library(viridis)
source("cpdb_my_dotplot.R")

# input argument:
# Path to a directory that contains combined tables of adjusted pvalues and mean exp of interacting pairs
args = commandArgs(trailingOnly=T)
res_path = args[1]

pval_file = file.path(res_path, "comb_pval.txt")
means_file = file.path(res_path, "comb_mean.txt")
if (!all(sapply(c(pval_file, means_file), FUN=file.exists))) {
print("check input file names.")
q()
}

pval = read.table(pval_file, sep='\t', header=T)
means = read.table(means_file, sep='\t', header=T)

df <- pval[,c(12:ncol(pval))]

rows.use = apply(df, 1, min) < 0.05

# resulting matrices
pval = pval[rows.use,]
means = means[rows.use,]

write.table(pval, file.path(res_path, "final_combined_pvalues.txt"), sep='\t', quote=F, row.names=F, col.names=T)
write.table(means, file.path(res_path, "final_combined_means.txt"), sep='\t', quote=F, row.names=F, col.names=T)

write.table(as.character(pval$interacting_pair), "final_sig_pairs.tsv", sep='\t', quote=F, row.names=F, col.names=F)


# Final dot plot
pval_file = file.path(res_path, "final_combined_pvalues.txt")
means_file = file.path(res_path, "final_combined_means.txt")
pval = read.table(pval_file, sep='\t', header=T)
means = read.table(means_file, sep='\t', header=T)

source("ordered_rows_updated.R")

row.order = unlist(row.ord)

# all pairs together
dot_plot(selected_rows=row.order, means_path=file.path(res_path, "final_combined_means.txt"), pvalues_path=file.path(res_path, "final_combined_pvalues.txt"), filename="cpdb_res_WTKO_all_interactions_new.pdf", height=0.5*length(row.order), width=27, my_palette=cols, sep.x=30.5)

# height of plots
a = list(General=7, KO=1.1, L14=3, L56=0.7, MgDriven=1.8, L5=6.2, L6=2.1, L15=1.5)

# plot by interaction class and WTKO separately
for(nam in names(row.ord)) {
cur.pairs = row.ord[[nam]]
dot_plot(selected_columns=grep("WT", colnames(pval), value=T), selected_rows=cur.pairs, means_path=file.path(res_path, "final_combined_means.txt"), pvalues_path=file.path(res_path, "final_combined_pvalues.txt"), filename=paste0("cpdb_res_WT_", nam, ".pdf"), height=a[[nam]], width=12, my_palette=cols)
dot_plot(selected_columns=grep("KO", colnames(pval), value=T), selected_rows=cur.pairs, means_path=file.path(res_path, "final_combined_means.txt"), pvalues_path=file.path(res_path, "final_combined_pvalues.txt"), filename=paste0("cpdb_res_KO_", nam, ".pdf"), height=a[[nam]], width=11, my_palette=cols)
}


library(stringr)
# heatmap
cnames = colnames(pval)[12:ncol(pval)]
a = str_split(cnames, '\\.', simplify=T)
mgs = c("WT_Hom1","WT_Hom2","KO_Hom1","KO_Hom2")
ns = names(table(a[,2]))
ns = ns[c(5:14,16:19,15,1:4)]  # adjust ordering of neuronal subtypes

mat=matrix(0, nrow=length(mgs), ncol=length(ns), dimnames=list(mgs, ns))

for (cn in cnames) {
cn.split = str_split(cn, '\\.', simplify=T)[1,]
mgtype = cn.split[1]
ntype = cn.split[2]
arr = pval[,cn]
count = sum(arr<=0.05)
mat[mgtype, ntype] = count
}

saveRDS(mat, "mat_data_for_heatmap.rds")

library(pheatmap)
library(RColorBrewer)
heatmap.cols = rev(colorRampPalette(brewer.pal(n = 7, name = "RdBu"))(50))
pdf("heatmap_num_interactions.pdf", height=3, width=10)
pheatmap(mat, cluster_rows=F, cluster_cols=F, scale='none', color=heatmap.cols)
dev.off()

