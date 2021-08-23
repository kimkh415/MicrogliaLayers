library(Seurat)
library(ggrepel)
library(tidyverse)
library(RColorBrewer)
library(viridis)
library(clusterProfiler)
library(org.Mm.eg.db)


#### Cell Type Signature Scoring ###
MgIdent <- c("P2ry12", "Fcrls", "Trem2", "Tmem119", "Cx3cr1", "Hexb", "Tgfbr1", "Sparc", "P2ry13", "Olfml3", "Adgrg1", "C1qa", "C1qb", "C1qc", "Csf1r", "Fcgr3", "Ly86", "Laptm5") #Genes from Marsh et al. Biorxiv, 2020
listMgIdent <- list(MgIdent)
P14_Cx3cr1 <- AddModuleScore(object = P14_Cx3cr1, features = listMgIdent, name = "MgIdentity_Marsh")

BAM <- c("Mrc1", "Ms4a7", "Pf4", "Stab1", "Cbr2", "Cd163", "Lyve1") #Genes from Jordao et al. Science, 2019
listBAM <- list(BAM)
P14_Cx3cr1 <- AddModuleScore(object = P14_Cx3cr1, features = listBAM, name = "BAMIdentity_Jordao")


#### ex vivo Activation Signature Scoring ####
CNSExVivo <- c("Fos", "Junb", "Zfp36", "Jun", "Hspa1a", "Socs3", "Rgs1", "Egr1", "Btg2", "Fosb", "Hist1h1d", "Ier5", "1500015O10Rik", "Atf3", "Hist1h2ac", "Dusp1", "Hist1h1e", "Folr1", "Serpine1") #Genes from Marsh et al. Biorxiv, 2020
listCNSExVivo <- list(CNSExVivo)
P14_Cx3cr1 <- AddModuleScore(object = P14_Cx3cr1, features = listCNSExVivo, name = "CNS_ExVivoAct_Marsh")


#### Mg State Signature Scoring ####
Homeostatic1_Sig <- as.charcter(c("Fth1", "Cd81", "Actg1", "Ckb", "Selenop", "Ldhb", "Gnb2", "Fcer1g", "Itm2b", "Ctss")
listHomeostatic1 <- list(Homeostatic1_Sig)
Homeostatic2_Sig <- as.character(c("Gm42418", "AY036118", "Cx3cr1", "Tpt1", "Slco2b1", "Zfhx3", "Adap2os", "P2ry13", "Pld1", "Hpgds", "Maf", "Xist", "Ivns1abp", "Cd164", "Srsf5", "Jmjd1c", "Pag1", "Dip2b", "Irf2bp2", "Gfm2", "Il6st", "Mertk", "Kctd12", "Rnaset2a", "9930111J21Rik2", "Ints6l", "Kif21b", "Tnrc6b", "Eif4a2", "Dst", "Macf1", "Kcnq1ot1"))
listHomeostatic2 <- list(Homeostatic2_Sig)
Innate_Sig <- as.character(c("Ifi27l2a", "Ifitm3", "Rtp4", "Cdkn1a", "Bst2", "Slfn2", "Lgals3bp", "Isg15", "Oasl2", "Ifit3", "H2-D1", "H2-K1", "Ccl12", "Stat1", "Trim30a"))
listInnate_Sig <- list(Innate_Sig)
Inflam_Sig <- as.character(c("Mt1", "Cd63", "Cd9", "Ftl1", "Mif", "Ccl3", "Ccl4", "C3ar1", "Abcg1", "Ctsb", "Cd83", "Spp1", "Ctsz", "Cstb", "Prdx1"))
listInflam_Sig <- list(Inflam_Sig)
Homeostatic3_Sig <- as.character(c("Gm26870", "Gm10800", "Gm11168", "Gm10717", "Gm10801"))
listHomeostatic3 <- list(Homeostatic3_Sig)
Apoe_Sig <- as.character(c("Apoe", "Lyz2", "Fth1"))
listApoe_Sig <- list(Apoe_Sig)
Ccr1_Sig <- as.character(c("Ccr1", "Tmem176a", "Tmem52", "Tmem176b", "Ramp1"))
listCcr1_Sig <- list(Ccr1_Sig)
Proliferative_Sig <- as.character(c("Stmn1", "Hist1h1b", "Hmgb2", "Top2a", "Mki67", "Mcm6", "Hells", "Atad2", "Topbp1", "Lig1"))
listProliferative_Sig <- list(Proliferative_Sig)

P14_Cx3cr1 <- AddModuleScore(object = P14_Cx3cr1, features = listHomeostatic1, name = "Homeostatic1_Sig") #Example creation of Add Module Score to Seurat Object

DoHeatmap(subset(P14_Cx3cr1, downsample = 100), features = Homeostatic1_Sig, size = 3, angle = 0) #Example DoHeatmap for Mg State signature genes

VlnPlot(P14_Cx3cr1, features = "Homeostatic1_Sig", pt.size = 0, cols = c("#FFCC66","#CC0000"), idents = c("Homeostatic1", "Homeostatic2"))  + stat_summary(fun.data = mean_sdl, geom = "point", size = 10, color = "black", shape = 95) #Example Violin Plot for Mg State signature genes


#### Cluster Proportion Calculations and Plotting ####
counts = table(P14_Cx3cr1$dataset, Idents(P14_Cx3cr1))
probs = as.data.frame.matrix(counts/rowSums(counts))
probs$layer = str_split(rownames(counts), '_', simplify=T)[,2]
probs$replicate = str_split(rownames(counts), '_', simplify=T)[,3]
tab <- probs %>% gather(CellType, Fraction, -layer, -replicate) %>% group_by(CellType, layer) %>% summarise(avg=mean(Fraction), se=sd(Fraction)/sqrt(length(Fraction)))

tab_Homs <- tab[7:12,]
tab_Homs$CellType <- factor(tab_Homs$CellType, levels = c("Homeostatic1","Homeostatic2"))
tab_Homs$layer <- factor(tab_Homs$layer, levels = c("Upper","Middle","Lower"))
tab_nonHoms <- tab[c(1:6,13:21),]
tab_nonHoms$CellType <- factor(tab_nonHoms$CellType, levels = c("Homeostatic3","ApoeHigh","Innate Immune","Ccr1High","Inflammatory"))
tab_nonHoms$layer <- factor(tab_nonHoms$layer, levels = c("Upper","Middle","Lower"))

ggplot(data = tab_Homs, aes(x=CellType, y = avg, fill = layer)) +
  geom_bar(width = 0.75, stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin=avg-se, ymax=avg+se), width=.2, position = position_dodge(.75)) +
  theme_bw() +
  labs(y="P14 Homeostatic Mg Cluster Proportion", x="") +
  scale_y_continuous(limits = c(0,1))

ggplot(data = tab_nonHoms, aes(x=CellType, y = avg, fill = layer)) +
  geom_bar(width = 0.75, stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin=avg-se, ymax=avg+se), width=.2, position = position_dodge(.75)) +
  theme_bw() +
  labs(y="P14 nonHomeostatic Mg Cluster Proportion", x="") +
  scale_y_continuous(limits = c(0,.2))


#### Comparing of Homeostatic1 and Homeostatic2 by Volcano Plot ####
P14_Homs_complete <- read.csv(file = "/MAST_DE/P14_MAST_markers_complete.csv")
P14_Homs_complete <- P14_Homs_complete %>%
  mutate(threshold = p_val_adj < 0.05 & abs(avg_logFC) > 0.15)
P14_Homs_complete <-  P14_Homs_complete %>%
  mutate(Signature = "")
P14_Homs_complete <- P14_Homs_complete[order(P14_Homs_complete$threshold, decreasing = TRUE),]
P14_Homs_complete$Signature[c(1:593)] <- as.character(P14_Homs_complete$X[c(1:593)])
P14_Homs_complete <- P14_Homs_complete %>%
  mutate(color = factor(case_when(avg_logFC > 0.25 & p_val_adj < 0.05 ~ "Homeostatic1 Enriched", avg_logFC < -0.25 & p_val_adj < 0.05 ~ "Homeostatic2 Enriched", TRUE ~ "No Enrichment")))

ggplot(P14_Homs_complete) +
geom_point(aes(x = avg_logFC, y = -log10(p_val_adj), colour = color, stroke = 0)) +
geom_text_repel(aes(x = avg_logFC, y = -log10(p_val_adj), label = Signature)) +
ggtitle("P14 Cx3cr1 Homeostatic 1 vs Homeostatic 2") +
geom_vline(xintercept = c(-.15, .15), color = "black", alpha = 1.0, linetype = "dotted") +
geom_hline(yintercept = -log10(0.05), color = "black", alpha = 1.0, linetype = "dotted") +
coord_cartesian(xlim = c(-.8, .8), ylim = c(0,300)) +
xlab("Avg. log2 Fold Change") +
ylab("Adjusted p-value (-log10)") +
theme(legend.position = "none",
      plot.title = element_text(size = rel(1.5), hjust = 0.5),
      axis.title = element_text(size = rel(1.25)))+
theme_bw() +
scale_color_manual(name = "color",
                   values = c("Homeostatic1 Enriched" = "#FFCC66", "Homeostatic2 Enriched" = "#CC0000", "No Enrichment" = "grey"))


#### Gene Ontology Analysis
Homeostatic1 <- read_tsv(file = "MAST_DE/P14_Homeostatic1_MAST_markers.tsv")
Homeostatic1 <- Homeostatic1[c(1:223),] #logFC greater than 0.15
Hom1_genes <- Homeostatic1$X1
Background_Complete <- read.csv(file = "../KwanhoRun_Analyses/MAST_DE/P14_MAST_markers_complete.csv")
Background <- Background_Complete$X
Background <- as.character(Background)
ego_Hom1_P14 <- enrichGO(gene = Hom1_genes,
                       universe = Background,
                       keyType = "SYMBOL",
                       OrgDb = org.Mm.eg.db,
                       ont = "BP",
                       pAdjustMethod = "BH",
                       qvalueCutoff = 0.05,
                       readable = FALSE)
Hom1_GO_summary <- data.frame(ego_Hom1_P14)

dotplot(ego_Hom1_P14, showCategory = 10)


#### MERFISH Plotting ####
DimPlot(MerFISH, reduction = "tsne", pt.size = 0.1, cols = c("#66CD00", "#0000FF", "#00EEEE", "#CD7054", "#CDAD00", "#AB82FF"))
DoHeatmap(subset(MerFISH, downsample = 200), features = c("Neurod2", "Erbb4", "Tmem119", "Sox10", "Aldh1l1", "Slc6a20a"), group.colors = c("#0000FF", "#00EEEE", "#66CD00", "#CD7054", "#CDAD00", "#AB82FF"))
DimPlot(MerFISH, reduction = "spatial", pt.size = 0.1, cols = c("#0000FF", "#00EEEE", "#66CD00", "#CD7054", "#CDAD00", "#AB82FF"), split.by = "orig.ident", ncol = 3)
VlnPlot(MerFISH, features = c("Slc17a7", "Neurod2", "Erbb4", "Gabbr2", "Tmem119", "Fcrls", "Sox10", "Itpr2", "Aldh1l1", "Slc6a20a", "Slc6a13"), pt.size = 0, cols = c("#0000FF", "#00EEEE", "#66CD00", "#CD7054", "#CDAD00", "#AB82FF"))
FeaturePlot(MerFISH, reduction = "spatial", features = "F13a1", ncol = 3, pt.size = 0.1) + scale_colour_gradientn(colours = viridis(10, direction = -1))


#### MERFISH Ridge Plot of normalized cell distance from pia ####
NormDistances_Major_Agg$Type <- factor(NormDistances_Major_Agg$Type, levels = c("Pericyte","Astrocyte","Oligo","Mg","IN","PN"))
ggplot(NormDistances_Major_Agg, aes(x = value, y = Type, fill = Type)) + geom_density_ridges(scale = 1.5, rel_min_height = 0.01) +
  theme_bw() +
  scale_fill_manual(values = c("#AB82FF", "#CDAD00", "#CD7054", "#66CD00", "#00EEEE", "#0000FF")) +
  scale_x_continuous(limits = c(-0.05,1.05))  +
  geom_vline(xintercept = 0.09, linetype="dashed", color="black") +
  geom_vline(xintercept = 0.28, linetype="dashed", color="black") +
  geom_vline(xintercept = 0.45, linetype="dashed", color="black") +
  geom_vline(xintercept = 0.63, linetype="dashed", color="black") +
  geom_vline(xintercept = 0.95, linetype="dashed", color="black")


#### Setting cortical layer cut-offs based on PN subtype localization - use the same cut-offs for layer calling of cell types and Mg states ####
NormDistances_PN_region_long <- NormDistances_PN_region_long %>%
mutate(Layer = factor(case_when(value <= 1 & value > 0.63 ~ "Layer 6", value <= 0.63 & value > 0.45 ~ "Layer 5", value <= 0.45 & value > 0.28 ~ "Layer 4", value <= 0.28 & value > 0.09 ~ "Layer 2/3", value <= 0.09 & value >= 0 ~ "Layer 1")))


#### Counting and plotting cell types by established cortical layer ####
NormDistances_PN_region_long <- NormDistances_PN_region_long[NormDistances_PN_region_long$value!=0,]
Sum_PN_region_long <- NormDistances_PN_region_long %>%
  group_by(Class, Region, Layer, .drop = FALSE) %>%
summarize(count = n()) #data normalized in excel
ggplot(data = normPN_mean_bylayer, aes(x=Layer, y = mean, fill = Class)) +
geom_bar(width = 0.75, stat = "identity", position = position_dodge(0.9)) +
theme_bw() +
scale_fill_manual(values = c("#6495ED", "#458B00", "#FFB90F", "#FF7F00", "#7FFF00", "#EE2C2C", "#68228B", "#8B6508")) +
geom_errorbar(aes(ymin=mean-SEM, ymax=mean+SEM), width=0.2, position = position_dodge(0.9)) +
geom_point(data = normPN_raw_bylayer, aes(x=Layer, y=NormPN), position = position_dodge(width = 0.9), size = .6) +
theme(legend.position = "none") +
labs(y="", x="") +
scale_y_continuous(limits = c(0,1.25))


#### Mg State QC Calling ####
Mg_Forced <- read.csv(file = "QC/mg_labels_Zscoring_Forced/Combined_mynewscores_Forced.csv")
Mg_Forced <- Mg_Forced %>%
  mutate(MgForcedCall = MgForcedCall)
Mg_Forced_AllScores <- Mg_Forced[,c(2:8,21)]
moltenForced <- melt(Mg_Forced_AllScores)
moltenForced$variable <- factor(moltenForced$variable, levels = c("Hom1_Pre", "Hom2_Pre", "Apoe_Pre", "Innate_Pre", "Inflam_Pre", "Ccr1_Pre", "Dividing_Pre"))
moltenForced$MgForcedCall <- factor(moltenForced$MgForcedCall, levels = c("Hom1_Fin", "Hom2_Fin", "Apoe_Fin", "Innate_Fin", "Inflam_Fin", "Ccr1_Fin", "Dividing_Fin"))
ggplot(data = moltenForced, aes(x=variable, y=value, fill=variable))+
geom_boxplot(outlier.size = 0)+
scale_fill_manual(values = c("#FFCC66", "#CC0000", "#1f78b4", "#a6cee3", "#66a61e", "#7570b3", "#f781bf"))+
facet_grid(MgForcedCall~.)+
theme_bw() +
scale_y_continuous(limits = c(-.75,1.15))


#### Mg Neighborhood PN composition ####
PNneighbors_Hom1_raw <- read.csv(file = "CellularNeighborhoodComposition/neighbors_FINAL/PNneighbors_Homeostatic1_raw.csv")
PNneighbors_Hom1_raw <- PNneighbors_Hom1_raw[,-1]
PNneighbors_Hom1_raw_Long <- gather(PNneighbors_Hom1_raw, factor_key=FALSE)
PNneighbors_Hom1_raw_Long$State <- as.character(rep("Homeostatic1", length(80)))
PNneighbors_Hom1_raw_Long$key <- factor(PNneighbors_Hom1_raw_Long$key, levels = c("L23_CPN", "L4_Stellate", "L5_PT", "L5_NP", "L5_CStrPN", "L6_CPN", "L6_CThPN", "L6b_Subplate"))
PNneighbors_Hom2_raw <- read.csv(file = "CellularNeighborhoodComposition/neighbors_FINAL/PNneighbors_Homeostatic2_raw.csv")
PNneighbors_Hom2_raw <- PNneighbors_Hom2_raw[,-1]
PNneighbors_Hom2_raw_Long <- gather(PNneighbors_Hom2_raw, factor_key=FALSE)
PNneighbors_Hom2_raw_Long$State <- as.character(rep("Homeostatic2", length(80)))
PNneighbors_Hom2_raw_Long$key <- factor(PNneighbors_Hom2_raw_Long$key, levels = c("L23_CPN", "L4_Stellate", "L5_PT", "L5_NP", "L5_CStrPN", "L6_CPN", "L6_CThPN", "L6b_Subplate"))
PNneighbors_HomeosFinal_raw <- rbind(PNneighbors_Hom1_raw_Long, PNneighbors_Hom2_raw_Long)
PNneighbors_HomeosFinal_raw$key <- as.character(PNneighbors_HomeosFinal_raw$key)
colnames(PNneighbors_HomeosFinal_raw) <- c("PN", "value", "State") #summarize data in excel
PNneighbors_HomeosSummary <- read.csv(file = "CellularNeighborhoodComposition/neighbors_FINAL/PNneighborhood_Homeostatics_meanSEM.csv")
ggplot(data = PNneighbors_HomeosSummary, aes(x=PN, y = Mean, fill = State)) +
geom_bar(width = 0.75, stat = "identity", position = position_dodge(0.9)) +
theme_bw() +
labs(y="Mg Neighbor Composition", x="") +
scale_fill_manual(values = c("#FFCC66","#CC0000")) +
ylim(0,.6) +
geom_errorbar(aes(ymin=Mean-SEM, ymax=Mean+SEM), width=0.2, position = position_dodge(0.9)) +
geom_point(data = PNneighbors_HomeosFinal_raw, aes(x=PN, y=value), position = position_dodge(width = 0.9), size = 0.6)


#### Ligand Receptor Plotting in MERFISH data ####
mg.cells = WhichCells(MgMerFISH, expression=Plxna4>0)  #Thresholds: Plxna4>0, Nrp2>0
n.cells = WhichCells(PNMerFISH, expression=Sema3a>.01)   #Thresholds: Sema3a>0.01, Sema3c>0.01

DimPlot(MerFISH, reduction='spatial', split.by='sample', group.by="CellType", cells.highlight= list(mg.cells, n.cells), cols.highlight=c("steelblue2", "red3"), cols= "grey", ncol=3, pt.size=.1, sizes.highlight = 0.1)
