# Adopted from dotplot function from CPDB

library(ggplot2)
dot_plot = function(selected_rows = NULL,
                    selected_columns = NULL,
                    filename = 'plot.pdf',
		    pval_threshold = 1e-5,
                    width = 8,
                    height = 10,
		    sep.y=NULL,
		    sep.x=NULL,
		    my_palette=NULL,
		    colorbar_max=6,
                    means_path = './means.txt',
                    pvalues_path = './pvalues.txt',
                    means_separator = '\t',
                    pvalues_separator = '\t',
                    output_extension = '.pdf'
){

  all_pval = read.table(pvalues_path, header=T, stringsAsFactors = F, sep=means_separator, comment.char = '', check.names=F)
  all_means = read.table(means_path, header=T, stringsAsFactors = F, sep=pvalues_separator, comment.char = '', check.names=F)

  intr_pairs = all_pval$interacting_pair
  all_pval = all_pval[,-c(1:11)]
  all_means = all_means[,-c(1:11)]

  if(is.null(selected_rows)){
    selected_rows = intr_pairs
  }

  if(is.null(selected_columns)){
    selected_columns = colnames(all_pval)
  }

  sel_pval = all_pval[match(selected_rows, intr_pairs), selected_columns]
  sel_means = all_means[match(selected_rows, intr_pairs), selected_columns]

  df_names = expand.grid(selected_rows, selected_columns)
  pval = unlist(sel_pval)
  pval[pval<pval_threshold] = pval_threshold
  plot.data = cbind(df_names,pval)
  pr = unlist(as.data.frame(sel_means))
  plot.data = cbind(plot.data,log2(pr+1))
  colnames(plot.data) = c('pair', 'clusters', 'pvalue', 'mean')

  if (is.null(my_palette)) {
    my_palette <- colorRampPalette(c("black", "blue", "yellow", "red"), alpha=TRUE)(n=399)
  }

  ggplot(plot.data,aes(x=clusters,y=pair)) +
  geom_point(aes(size=-log10(pvalue),color=mean)) +
  scale_color_gradientn('Log2 mean (Molecule 1, Molecule 2)', colors=my_palette, limits=c(0,colorbar_max)) +
  geom_hline(yintercept=sep.y) +
  geom_vline(xintercept=sep.x) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text=element_text(size=14, colour = "black"),
        #axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5, family = 'Arial'),
	axis.text.x = element_blank(),
        axis.text.y = element_text(size=12, colour = "black", family = 'Arial'),
        axis.title=element_blank(),
        text = element_text('Arial'),
        panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"))

  if (output_extension == '.pdf') {
      ggsave(filename, width = width, height = height, device = cairo_pdf, limitsize=F)
  }
  else {
      ggsave(filename, width = width, height = height, limitsize=F)
  }
}
