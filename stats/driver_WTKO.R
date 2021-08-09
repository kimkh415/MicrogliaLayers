.libPaths(.libPaths()[2])

library(lme4)
library(dplyr)
library(lmerTest)

## prep data
#flist = c("MicrogliaDensity_Summary_P7_Fezf2_RawDataBins.csv","MicrogliaDensity_SummaryStatistics_P14_Fezf2_RawDataBins.csv","MicrogliaDensity_RawData_P60_Fezf2_GroupedBins.csv")
#dat <- read.csv(flist[1], header=T)
#dat <- rbind(dat, read.csv(flist[2], header=T))
#dat <- rbind(dat, read.csv(flist[3], header=T))
#colnames(dat) <- c('layer', 'image','genotype','age','y')
#tmp <- as.character(dat$image)
#dat$mouse <- as.factor(gsub("_Image[0-9]$|_[0-9]$", "", tmp))
#dat$layer <- as.factor(gsub(" \\(Bins 1-5\\)| \\(Bins 6-10\\)| \\(Bins 6-11\\)", "", dat$layer))
#dat$image <- as.factor(dat$image)
#dat$genotype <- as.factor(dat$genotype)
#dat$age <- as.factor(dat$age)
#saveRDS(dat, "mg_density_WTKO.rds")

# read prepped data
dat <- readRDS("mg_density_WTKO.rds")

# mg density difference between WT and KO
# in individual layer (upper and lower separately)
# and in individual age
mg.genotype = c()
for (a in levels(dat$age)) {
sdat = dat %>% filter(age==a)
for (l in levels(dat$layer)) {
ssdat = sdat %>% filter(layer==l)
print(table(ssdat$age))
print(table(ssdat$layer))
fit1=lmer(y~genotype+(1|mouse),data=ssdat)
#fit2=lmer(y~(1|mouse),data=ssdat)
lname = gsub(" |-", "_", l)
res = summary(fit1)
mg.genotype = c(mg.genotype, res$coefficients[2,"Pr(>|t|)"])
saveRDS(res, paste0("lmerTest_summary_lmer_mgDensity_diff_by_genotype_", a, lname, ".rds"))
}
}
names(mg.genotype) <- rep(levels(dat$age), each=nlevels(dat$layer))
names(mg.genotype) <- paste0(names(mg.genotype), rep(levels(dat$layer), times=nlevels(dat$age)))





# Difference between upper and lower layer
# in individual genotype
mg.layer = c()
for (a in levels(dat$age)) {
sdat = dat %>% filter(age==a)
for (gt in levels(dat$genotype)) {
ssdat = sdat %>% filter(genotype==gt)
print(table(ssdat$age))
print(table(ssdat$genotype))
#fit1=lmer(y~layer+(1|mouse)+(1|image),data=ssdat)
fit1=lmer(y~layer+(1|image),data=ssdat)
#fit1=lm(y~layer+mouse+image,data=ssdat)
gtname = gsub(" |-", "_", gt)
res = summary(fit1)
#f = res$fstatistic
#pval = pf(f[1],f[2],f[3],lower.tail=F)
pval = res$coefficients[2,5]
mg.layer = c(mg.layer, pval)
saveRDS(res, paste0("lmerTest_summary_lmer_noMouse_mgDensity_diff_by_layer_", a, gtname, ".rds"))
}
}
names(mg.layer) <- rep(levels(dat$age), each=nlevels(dat$genotype))
names(mg.layer) <- paste0(names(mg.layer), rep(levels(dat$genotype), times=nlevels(dat$age)))


# Difference between the genotype in each age
dat$age <- factor(dat$age, levels=c("P7","P14","P60"))

ret = c()
for (a in levels(dat$age)) {
print(a)
sdat = dat %>% filter(age==a)
fit1 = lmer(y~layer+genotype+(1|mouse), data=sdat)
fit2 = lmer(y~layer+(1|mouse), data=sdat)
anov=anova(fit1, fit2)
ret=c(ret,anov$"Pr(>Chisq)"[2])
}
names(ret)=levels(dat$age)
ret.padj = p.adjust(ret, method="BH")
print(ret.padj)
write.table(ret.padj,file="mg_density_in_layers_by_genotype_lmm_anova_pvals.tsv",sep="\t",quote=F,row.names=T)

