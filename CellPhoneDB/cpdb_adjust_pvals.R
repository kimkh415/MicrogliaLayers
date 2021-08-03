library(reshape2)
library(stringr)

args = commandArgs(trailingOnly=T)
res_path = args[1]

if (endsWith(res_path, '/')) {
res_path = sub('/$', '', res_path)
}

pval_file = file.path(res_path, "pvalues.txt")
means_file = file.path(res_path, "means.txt")
if (!all(sapply(c(pval_file, means_file), FUN=file.exists))) {
print("check input file names.")
q()
}


# Load CPDB output
pval = read.table(pval_file, sep='\t', header=T)
means = read.table(means_file, sep='\t', header=T)

# Adjust p-values
df <- pval[,c(1,12:ncol(pval))]  # for cpdb v2.x
dfm <- melt(df)
colnames(dfm) <- c("id_cp_interaction", "cell_types", "p.val")
dfm$p.adj <- p.adjust(dfm$p.val, method="BH")  # using Benjamini Hochberg method
dfm2 <- dfm[,-3]
df2 = dcast(dfm2, id_cp_interaction ~ cell_types)

# confirm the row order
ridx = match(means[,1], df2[,1])  # for cpdb v2.x
#print(which(is.na(ridx)))
padj=df2[ridx,-1]
padj=cbind(pval[,1:11], padj)

write.table(padj, file.path(res_path, "adjust_pvalues.txt"), sep='\t', quote=F, row.names=F, col.names=T)

print("DONE")
