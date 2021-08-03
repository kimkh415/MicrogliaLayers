library(stringr)

# Filter WT Cellphonedb result
res_path = "/stanley/levin_dr_storage/kwanho/jeff_microglia/analysis_071321/cpdb/KO/out"
prefix='KO'

#PN = c('L2_3_CPN_1','L2_3_CPN_2','L2_3_CPN_3','L2_3_CPN_4','L4_Stellate','L5_CPN_1','L5_CPN_2','L5_CStrPN','L5_NP','L5_PT','L6_CPN_1','L6_CPN_2','L6_CThPN_1','L6_CThPN_2','L6b_Subplate')
PN = c('L2_3_CPN_1','L2_3_CPN_2','L2_3_CPN_3','L2_3_CPN_4','L4_Stellate','L5_CPN_1','L6_CPN_1','L6_CPN_2','KO_Mismatch_1','KO_Mismatch_2','KO_Mismatch_3','KO_Mismatch_4')
MG = c('Homeostatic1', 'Homeostatic2')

pval_file = file.path(res_path, "adjust_pvalues.txt")
means_file = file.path(res_path, "means.txt")
if (!all(sapply(c(pval_file, means_file), FUN=file.exists))) {
print("check input file names.")
q()
}

pval = read.table(pval_file, sep='\t', header=T)
means = read.table(means_file, sep='\t', header=T)

tab.pairs = expand.grid(MG, PN)
pn.mg.pairs = c(paste0(tab.pairs[,1], '.', tab.pairs[,2]), paste0(tab.pairs[,2], '.', tab.pairs[,1]))
stopifnot(all(sapply(pn.mg.pairs, FUN=function(x){x %in% colnames(pval)})))

info.cols = colnames(pval)[1:11]

pval = pval[,c(info.cols, pn.mg.pairs)]
means = means[,c(info.cols, pn.mg.pairs)]

comb.pair = readRDS("combined_sig_pairs.rds")
rows.use = match(comb.pair, pval$interacting_pair)

# resulting matrices
pval = pval[rows.use,]
means = means[rows.use,]

# save!
write.table(pval, paste0(prefix, "_pvalues_for_comb.txt"), sep='\t', quote=F, row.names=F, col.names=T)
write.table(means, paste0(prefix, "_means_for_comb.txt"), sep='\t', quote=F, row.names=F, col.names=T)

print("DONE!")

