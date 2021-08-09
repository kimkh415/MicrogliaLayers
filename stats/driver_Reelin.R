.libPaths(.libPaths()[2])

library(lme4)
library(dplyr)
library(lmerTest)

## prep data
#dat <- read.csv("ReelinWTKO_P14_Satb2MicrogliaDensityAnalysis.csv", header=T)
#colnames(dat) <- c('image','bin','genotype','mg_density','satb2_density','bin_geno','ratio')
#tmp <- as.character(dat$image)
#dat$mouse <- as.factor(gsub("_Image[0-9]$", "", tmp))
#dat$image <- as.factor(dat$image)
#dat$genotype <- as.factor(dat$genotype)
#dat$bin_geno <- as.factor(dat$bin_geno)
#dat$bin <- as.factor(dat$bin)
#saveRDS(dat, "mg_density_Reelin.rds")

# read prepped data
dat <- readRDS("mg_density_Reelin.rds")




#########
# 1. Try to fit a lmm considering all the random effects (image and mouse)
# 2. If the model is too complex (output singular model), fit a lm model instead
# pvalue obtained from lmerTest package
# when comparing multiple groups (e.g. ratio difference among bins in Reelin WT),
# another model was fit without a fixed-effect then performed anova to test whether
# the two models are significantly different.
########




# Ratio difference among the bins in individual genotype
ratio.bin=c()
for (gt in levels(dat$genotype)) {
sdat = dat %>% filter(genotype==gt)
fit1=lmer(ratio~bin+(1|image)+(1|mouse), data=sdat, REML=F)
fit2=lmer(ratio~(1|image)+(1|mouse), data=sdat, REML=F)
gtname = gsub(" ", "", gt)
res=anova(fit1,fit2,test="LRT")
print(res)
ratio.bin=c(ratio.bin, res[2,"Pr(>Chisq)"])
saveRDS(res, paste0("anova_lmer_res_ratio_diff_among_bins_", gtname, ".rds"))
}
names(ratio.bin) = levels(dat$genotype)

# Ratio difference between the genotype in each bin
ratio.genotype=c()
for (b in levels(dat$bin)) {
sdat = dat %>% filter(bin==b)
fit1=lmer(ratio~genotype+(1|mouse), data=sdat)
res = summary(fit1)
ratio.genotype = c(ratio.genotype, res$coefficients[2,5])
saveRDS(res, paste0("lmerTest_summary_lmer_ratio_diff_bet_genotype_", b, ".rds"))
}
names(ratio.genotype) = levels(dat$bin)

mg.bin=c()
# Mg density difference among the bins in individual genotype
for (gt in levels(dat$genotype)) {
sdat = dat %>% filter(genotype==gt)
fit1=lmer(mg_density~bin+(1|image)+(1|mouse), data=sdat, REML=F)
fit2=lmer(mg_density~(1|image)+(1|mouse), data=sdat, REML=F)
gtname = gsub(" ", "", gt)
res=anova(fit1,fit2,test="LRT")
print(res)
mg.bin=c(mg.bin, res[2,"Pr(>Chisq)"])
saveRDS(res, paste0("anova_lmer_res_mgDensity_diff_among_bins_", gtname, ".rds"))
}
names(mg.bin) = levels(dat$genotype)

# Mg density difference between the genotype in each bin
# lmer fit produces singular model for some
mg.genotype=c()
for (b in levels(dat$bin)) {
sdat = dat %>% filter(bin==b)
#fit1=lmer(mg_density~genotype+(1|mouse), data=sdat)
#fit1=lm(mg_density~genotype, data=sdat)
#fit1=lm(mg_density~genotype+mouse, data=sdat)
fit1=lm(mg_density~genotype, data=sdat)
res=summary(fit1)
#print(res)
mg.genotype=c(mg.genotype,res$coefficients[2,4])
#saveRDS(res, paste0("lmerTest_summary_res_lm_mgDensity_diff_bet_genotype_", b, ".rds"))
}
names(mg.genotype) = levels(dat$bin)
print(mg.genotype)

# CPN density difference among the bins in individual genotype
# lmer fit produces singular model for some
cpn.bin=c()
for (gt in levels(dat$genotype)) {
sdat = dat %>% filter(genotype==gt)
#fit1=lmer(satb2_density~bin+(1|image)+(1|mouse), data=sdat)
#fit2=lmer(satb2_density~(1|image)+(1|mouse), data=sdat)
fit1=lm(satb2_density~bin+image+mouse, data=sdat)
fit2=lm(satb2_density~image+mouse, data=sdat)
gtname = gsub(" ", "", gt)
res=anova(fit1,fit2,test="LRT")
#print(res)
cpn.bin=c(cpn.bin, res[2,"Pr(>Chi)"])
saveRDS(res, paste0("anova_lm_res_cpnDensity_diff_among_bins_", gtname, ".rds"))
}
names(cpn.bin) = levels(dat$genotype)

# CPN density difference between the genotype in each bin
cpn.genotype=c()
for (b in levels(dat$bin)) {
sdat = dat %>% filter(bin==b)
fit1=lmer(satb2_density~genotype+(1|mouse), data=sdat)
res=summary(fit1)
cpn.genotype=c(cpn.genotype,res$coefficients[2,5])
#saveRDS(res, paste0("lmerTest_summary_res_lmer_cpnDensity_diff_bet_genotype_", b, ".rds"))
}
names(cpn.genotype) = levels(dat$bin)


print(ratio.bin)
print(ratio.genotype)
print(mg.bin)
print(mg.genotype)
print(cpn.bin)
print(cpn.genotype)



