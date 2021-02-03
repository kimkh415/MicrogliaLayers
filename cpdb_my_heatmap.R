# Adopted from CPDB heatmap to change the order of clusters

library(pheatmap)
library(stringr)
heatmaps_plot = function(meta_file, pvalues_file, count_filename="heatmap_count.pdf", log_filename="heatmap_log_count.pdf", count_network_filename="count_network.txt", interaction_count_filename="interaction_count.txt", count_network_separator="\t", interaction_count_separator="\t", show_rownames = T, show_colnames = T,
                         scale="none", cluster_cols = F,border_color='white', cluster_rows = F, fontsize_row=11,
                         fontsize_col = 11, main = '',treeheight_row=0, family='Arial', treeheight_col = 0,
                         col1 = "dodgerblue4", col2 = 'peachpuff', col3 = 'deeppink4', meta_sep='\t', pvalues_sep='\t', pvalue=0.05){
  #######   Network

  meta = read.csv(meta_file, comment.char = '', sep=meta_sep)

  all_intr = read.table(pvalues_file, header=T, stringsAsFactors = F, sep=pvalues_sep, comment.char = '', check.names = F)
  intr_pairs = all_intr$interacting_pair
  all_intr = all_intr[,-c(1:11)]


  split_sep = '\\|'
  join_sep = '|'

  pairs1_all = unique(meta[,2])
  pairs1_all = pairs1_all[order(pairs1_all)]
  pairs1_all <- levels(pairs1_all)[c(20,22, 7:12, 14:16, 13, 18,19,17,21,23, 1:6)]

  pairs1 = c()
  for (i in 1:length(pairs1_all))
    for (j in 1:length(pairs1_all))
        pairs1 = c(pairs1,paste(pairs1_all[i],pairs1_all[j],sep=join_sep))

  all_count = matrix(ncol=3)
  colnames(all_count) = c('SOURCE','TARGET','count')
  count1 = c()
  for(i in 1:length(pairs1))
  {
    p1 = strsplit(pairs1[i], split_sep)[[1]][1]
    p2 = strsplit(pairs1[i], split_sep)[[1]][2]

    n1 = intr_pairs[which(all_intr[,pairs1[i]]<=pvalue)]
    pairs_rev = paste(p2, p1, sep=join_sep)
    n2 = intr_pairs[which(all_intr[,pairs_rev]<=pvalue)]

    if(p1!=p2)
      count1 = length(unique(n1))+length(unique(n2))
    else
      count1 = length(unique(n1))

    new_count = c(p1,p2,count1)
    names(new_count) = c('SOURCE','TARGET','count')
    all_count = rbind(all_count, new_count)
  }

  all_count = all_count[-1,]
  write.table(all_count, paste0(cluster_rows, "_", count_network_filename), sep=count_network_separator, quote=F, row.names = F)

  #######   count interactions
  print(length(pairs1))
  count1 = c()
  for(i in 1:length(pairs1))
  {
    p1 = strsplit(pairs1[i], split_sep)[[1]][1]
    p2 = strsplit(pairs1[i], split_sep)[[1]][2]

    n1 = intr_pairs[which(all_intr[,pairs1[i]]<=pvalue)]

    pairs_rev = paste(p2, p1, sep=join_sep)
    n2 = intr_pairs[which(all_intr[,pairs_rev]<=pvalue)]
    if(p1!=p2)
      count1 = c(count1,length(unique(n1))+length(unique(n2)))
    else
      count1 = c(count1,length(unique(n1)))

  }
  # print(count1)
  names(count1) <- pairs1
  write.table(count1, paste0(cluster_rows, "_", interaction_count_filename), sep=count_network_separator, quote=F, col.names = F) 
  if (any(count1)>0)
  {
  count_matrix <- matrix(count1, nrow=23, ncol=23)
  labels <- str_split(names(count1), '\\|', simplify=T)[,2]
  rownames(count_matrix) <- labels
  colnames(count_matrix) <- labels
  saveRDS(count_matrix, "count_mtx.rds")
  pheatmap(count_matrix, show_rownames = show_rownames, show_colnames = show_colnames, scale=scale, cluster_cols = cluster_cols,
             border_color=border_color, cluster_rows = cluster_rows, fontsize_row = fontsize_row, fontsize_col = fontsize_col,
             main = main, treeheight_row = treeheight_row, family = family,color = col.heatmap, treeheight_col = treeheight_col, filename = count_filename)
 col.heatmap <- colorRampPalette(c(col1,col2,col3 ))( 1000 )
pheatmap(count_matrix, show_rownames = show_rownames, show_colnames = show_colnames, scale=scale, cluster_cols = cluster_cols,
             border_color=border_color, cluster_rows = cluster_rows, fontsize_row = fontsize_row, fontsize_col = fontsize_col,
             main = main, treeheight_row = treeheight_row, family = family,color = col.heatmap, treeheight_col = treeheight_col, filename = count_filename)
 pheatmap(log(count_matrix+1), show_rownames = show_rownames, show_colnames = show_colnames, scale=scale, cluster_cols = cluster_cols,
             border_color=border_color, cluster_rows = cluster_rows, fontsize_row = fontsize_row, fontsize_col = fontsize_col,
             main = main, treeheight_row = treeheight_row, family = family,color = col.heatmap, treeheight_col = treeheight_col, filename = log_filename) 
  } else {
    stop("There are no significant results using p-value of: ", pvalue, call.=FALSE)
  }
}
