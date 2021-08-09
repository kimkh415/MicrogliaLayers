# Filter WT Cellphonedb result
# Need to adjust p values before this step
library(stringr)


# path to a directory containing CellphoneDB result
args = commandArgs(trailingOnly=T)
res_path = args[1]

pval_file = file.path(res_path, "adjust_pvalues.txt")
means_file = file.path(res_path, "means.txt")
if (!all(sapply(c(pval_file, means_file), FUN=file.exists))) {
print("check input file names.")
q()
}

pval = read.table(pval_file, sep='\t', header=T)
means = read.table(means_file, sep='\t', header=T)

PN = c('L2_3_CPN_1','L2_3_CPN_2','L2_3_CPN_3','L2_3_CPN_4','L4_Stellate','L5_CPN_1','L5_CPN_2','L5_CStrPN','L5_NP','L5_PT','L6_CPN_1','L6_CPN_2','L6_CThPN_1','L6_CThPN_2','L6b_Subplate')
MG = c('Homeostatic1', 'Homeostatic2')

tab.pairs = expand.grid(MG, PN)
pn.mg.pairs = c(paste0(tab.pairs[,1], '.', tab.pairs[,2]), paste0(tab.pairs[,2], '.', tab.pairs[,1]))
stopifnot(all(sapply(pn.mg.pairs, FUN=function(x){x %in% colnames(pval)})))

info.cols = colnames(pval)[1:11]

pval = pval[,c(info.cols, pn.mg.pairs)]
means = means[,c(info.cols, pn.mg.pairs)]

df <- pval[,c(12:ncol(pval))]

# filering interactions by adjusted p-values
rows.use = apply(df, 1, min) < 0.05

# resulting matrices
pval = pval[rows.use,]
means = means[rows.use,]

# save!
write.table(pval, file.path(res_path, "mysig_pvalues.txt"), sep='\t', quote=F, row.names=F, col.names=T)
write.table(means, file.path(res_path, "mysig_means.txt"), sep='\t', quote=F, row.names=F, col.names=T)

write.table(as.character(pval$interacting_pair), "mysig_pairs_WT.tsv", sep='\t', quote=F, row.names=F, col.names=F)

source("cpdb_my_dotplot.R")
my.sep = (length(pn.mg.pairs)+1)/2
dot_plot(means_path=file.path(res_path, "mysig_means.txt"), pvalues_path=file.path(res_path, "mysig_pvalues.txt"), filename=file.path(res_path, "cpdb_WT_res.pdf"), height=32, width=20, sep.x=my.sep)

print("DONE!")



