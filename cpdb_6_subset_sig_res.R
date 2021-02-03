pval = read.table("filtered_padj.txt", sep='\t', header=T)
means = read.table("filtered_means.txt", sep='\t', header=T)
rows = as.character(read.table("rows.txt", sep='\t', header=F)$V1)

ridx=which(pval$interacting_pair %in% rows)

pval.sub=pval[ridx,]
means.sub=means[ridx,]

write.table(pval.sub, "subset_padj.txt", sep='\t', quote=F, row.names=F, col.names=T)
write.table(means.sub, "subset_means.txt", sep='\t', quote=F, row.names=F, col.names=T)

pval_r = read.table("filtered_padj_r.txt", sep='\t', header=T)
means_r = read.table("filtered_means_r.txt", sep='\t', header=T)
rows_r = as.character(read.table("rows_r.txt", sep='\t', header=F)$V1)

ridx_r=which(pval_r$interacting_pair %in% rows_r)

pval_r.sub=pval_r[ridx_r,]
means_r.sub=means_r[ridx_r,]

write.table(pval_r.sub, "subset_padj_r.txt", sep='\t', quote=F, row.names=F, col.names=T)
write.table(means_r.sub, "subset_means_r.txt", sep='\t', quote=F, row.names=F, col.names=T)

