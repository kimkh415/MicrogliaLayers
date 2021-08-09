# This code was run on R v4.0.3

library(lme4)
library(dplyr)
library(lmerTest)

## prep data
#dat <- read.csv("MicrogliaDensity_Summary_P7_P14_P60_WT_RawDataBins.csv", header=T)
#dat$Genotype <- "Fezf2 WT"  # this dataset includes WT only
#colnames(dat) <- c('layer', 'image','genotype','age','y')
#tmp <- as.character(dat$image)
#dat$mouse <- as.factor(gsub("_Image[0-9]$|_[0-9]$", "", tmp))
#saveRDS(dat, "mg_density_WTonly.rds")

# read prepped data
dat <- readRDS("mg_density_WTonly.rds")

# test difference by layer(Region)
fit1=lmer(y~layer+age+(1|mouse),data=dat)
fit2=lmer(y~ age+(1|mouse),data=dat)
anov=anova(fit1,fit2,test="LRT")
saveRDS(anov, "anova_lmer_res_WTonly_diff_bet_layer.rds")

# Within each time point, test layer as fixed effect
mg.layer=c()
for (a in levels(dat$age)) {
sdat = dat %>% filter(age==a)
print(table(sdat$age))
fit1=lmer(y~layer+(1|image),data=sdat)
res = summary(fit1)
pval = res$coefficients[2,5]
mg.layer=c(mg.layer, pval)
saveRDS(res, paste0("lmerTest_summary_lmer_WTonly_mgDensity_diff_bet_layer_", a, ".rds"))
}
names(mg.layer)=levels(dat$age)
mg.layer

# within each layer, test age as fixed effect
for(a in levels(dat$layer)) {
sdat = dat %>% filter(layer==a)
fit1 = lmer(y~age+(1|mouse), data=sdat)
fit2 = lmer(y~(1|mouse), data=sdat)
anov=anova(fit1, fit2, test='LRT')
print(anov)
}


