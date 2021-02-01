library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(gplots)
library(dplyr)
library(randomForest)
library(ROCR)
library(pheatmap)


fitMethod <- function(dat)
{
y=factor(dat[,"CellType"])
print(levels(y))
print(table(y))
x <- dat %>% select(-CellType)
ret=tuneRF(x,y,doBest=T)
dev.off()
return(ret)
}


valMethod <- function(model, dat, cols, model_name="RF")
{
pdf(paste0(model_name, "_model_ROC.pdf"))
pred_probs <- predict(model, dat, type='prob')
classes <- levels(dat$CellType)
auc_vals <- c()
for (i in 1:length(classes))
{
true_val <- ifelse(dat$CellType==classes[i], 1, 0)
pred <- prediction(pred_probs[,i], true_val)
perf <- performance(pred, "tpr", "fpr")
if(i==1)
{
plot(perf,main="ROC Curve", col=cols[classes[i]])
}
else
{
plot(perf,main="ROC Curve", col=cols[classes[i]],add=TRUE)
}
auc.perf <- performance(pred, measure = "auc")
print(auc.perf@y.values)
auc_vals <- c(auc_vals, auc.perf@y.values[[1]])
}
legend("bottomright", legend=names(cols), col=cols, lty=1)
dev.off()

names(auc_vals) <- classes
write.table(auc_vals, paste0(model_name, "_model_AUC.tsv"))
}


testMethod <- function(model, dat)
{
y=predict(model,dat)
y_prob = predict(model, dat, type='prob')

return(list(y, y_prob))
}


plotRes <- function(out, prob, cols, prefix)
{
ord <- order(out)
out <- out[ord]
prob <- prob[ord,levels(out)]
anno = as.character(out)
names(anno) <- names(out)
anno <- as.data.frame(anno)
colnames(anno) <- "Predicted Label"
anno_cols <- list("Predicted Label"=cols)

pdf(paste0(prefix, "_prob_heatmap_rf.pdf"), height=8, width=12)
pheatmap(t(prob), legend=T, scale='column', show_colnames=F, main="Prediction Probabilities", treeheight_row = 0, treeheight_col = 0, cluster_rows=F, cluster_cols=F, annotation_col = anno, annotation_colors=anno_cols)
pheatmap(t(prob), legend=T, scale='none', show_colnames=F, main="Prediction Probabilities", treeheight_row = 0, treeheight_col = 0, cluster_rows=F, cluster_cols=F, annotation_col = anno, annotation_colors=anno_cols)
dev.off()

}



human2mouse <- function(hgenes)
{
require(biomaRt)

human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
 
genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = hgenes , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
 
mgenes <- as.character(genesV2[, 1])
names(mgenes) <- as.character(genesV2[, 2])

mgenes <- mgenes[which(duplicated(names(mgenes))==F)]

return(mgenes)
}


mouse2human <- function(mgenes)
{
require(biomaRt)

human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")

genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = mgenes , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)

hgenes <- as.character(genesV2[, 1])
names(hgenes) <- as.character(genesV2[, 2])

hgenes <- hgenes[which(duplicated(names(hgenes))==F)]

return(hgenes)
}

