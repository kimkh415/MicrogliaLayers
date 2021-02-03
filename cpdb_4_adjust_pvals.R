library(reshape2)

args = commandArgs(trailingOnly=T)

# load separated CPDB output matrices
pval = read.table(args[1], sep='\t', header=T)
pval_r = read.table(args[2], sep='\t', header=T)
means = read.table(args[3], sep='\t', header=T)
means_r = read.table(args[4], sep='\t', header=T)

# combine them
comb_pval = cbind(pval, pval_r[,12:ncol(pval_r)])

# rows 238 and 239 have identical names for the 'interacting_pair' column
# both of these rows have pvalues of 1 across all columns
comb_pval = comb_pval[-239,]

df <- comb_pval[,c(2,12:ncol(comb_pval))]
dfm <- melt(df)
colnames(dfm) <- c("interacting_pair", "cell_types", "p.val")
dfm$p.adj <- p.adjust(dfm$p.val, method="BH")  # using Benjamini Hochberg method
write.table(dfm, "dfm_pval_padj.tsv", sep='\t', quote=F, row.names=F, col.names=T)

dfm2 <- dfm[,-3]
adj_pvals = dcast(dfm2, interacting_pair ~ cell_types)

saveRDS(adj_pvals, "adjusted_pvalues_BH.rds")

padj = adj_pvals[,1:77]
padj_r = adj_pvals[,c(1,78:ncol(adj_pvals))]

# match the rows in the means output for further filtering and plotting
means=means[-239,]
means_r = means_r[-239,]

# confirm the row order
ridx = match(means[,2], padj[,1])
#print(which(is.na(ridx)))
padj=padj[ridx,-1]
padj_r=padj_r[ridx,-1]
padj=cbind(means[,1:11], padj)
padj_r=cbind(means_r[,1:11], padj_r)

write.table(padj, "filtered_padj.txt", sep='\t', quote=F, row.names=F, col.names=T)
write.table(padj_r, "filtered_padj_r.txt", sep='\t', quote=F, row.names=F, col.names=T)
write.table(means, "filtered_means.txt", sep='\t', quote=F, row.names=F, col.names=T)
write.table(means_r, "filtered_means_r.txt", sep='\t', quote=F, row.names=F, col.names=T)

print("DONE")
