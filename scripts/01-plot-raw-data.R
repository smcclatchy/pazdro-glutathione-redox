# initial plots to determine whether or not data are normal before
# considering transforming them

# load tidyverse library for access to ggplot2
library(tidyverse)

# load data
load("../data/Heart-GSSG-enviornment.RData")
# the RData file just loaded contains the following:
# control - a cross2 object for qtl mapping
# kinship_lmm and kinship_loco - kinship matrices for qtl mapping in the DO
# Pazdro_GigaMUGA_genoprobs_qced_8state_sorted - allele probabilities for the
#                                               8 DO founder strains A thru G
# pheno - GSH and GSSG phenotypes with animal IDs, sex, generation
# probs - appears to be a duplicate of Pazdro_GigaMUGA_genoprobs_qced_8state_sorted
# qtlscan_HeartGSSG - scan object including sex + generation (sexgen)
# sex - values for animal sex as 0 or 1
# query_genes - function for querying a genomic database
# query_genes_mgi - function for querying MGI
# query_variants - function for variant querying
# rankZ - function for rankZ transformation


# plot the raw data from pheno

pdf(file = "../results/raw-plots.pdf")
hist(pheno$Heart_GSH)
hist(pheno$Heart_GSSG)
pheno %>% 
  ggplot(aes(x = Heart_GSH, y = Heart_GSSG)) + 
  geom_point() +
  stat_smooth(method = "lm")

# separate linear models by sex
pheno %>% 
  ggplot(aes(x = Heart_GSH, y = Heart_GSSG, color = as.factor(sex))) +
  geom_point() +
  stat_smooth(method = "lm")

# separate linear models by generation
pheno %>% 
  ggplot(aes(x = Heart_GSH, y = Heart_GSSG, color = as.factor(generation))) +
  geom_point() +
  stat_smooth(method = "lm")

qqnorm(pheno$Heart_GSH, main = "Heart GSH Q-Q Plot")
qqline(pheno$Heart_GSH)
qqnorm(pheno$Heart_GSSG, main = "Heart GSSG Q-Q Plot")
qqline(pheno$Heart_GSSG)

# boxplots by sex and generation
pheno %>% 
  ggplot(aes(x=as.factor(sex), y=Heart_GSH)) +
  geom_boxplot() +
  geom_point(position=position_jitter(width=0.25, height=0)) +
  xlab("sex") 

pheno %>% 
  ggplot(aes(x=as.factor(sex), y=Heart_GSSG)) +
  geom_boxplot() +
  geom_point(position=position_jitter(width=0.25, height=0)) +
  xlab("sex")  

pheno %>% 
  ggplot(aes(x=as.factor(generation), y=Heart_GSH)) +
  geom_boxplot() +
  geom_point(aes(color = as.factor(sex)), position=position_jitter(width=0.25, height=0)) +
  xlab("generation")  

pheno %>% 
  ggplot(aes(x=as.factor(generation), y=Heart_GSSG)) +
  geom_boxplot() +
  geom_point(aes(color = as.factor(sex)), position=position_jitter(width=0.25, height=0)) +
  xlab("generation")  

# linear model of GSSG as response to GSH
GSHlm <- lm(pheno$Heart_GSSG ~ pheno$Heart_GSH)
summary(GSHlm)
intercept <- GSHlm$coefficients[[1]]
slope <- GSHlm$coefficients[[2]]
# equation for GSSG as a response to GSH
# GSSG = 2.36 + -0.075*GSH
plot(pheno$Heart_GSH, pheno$Heart_GSSG)
abline(a = intercept, b = slope)

# look at residuals vs fitted values
resid <- GSHlm$residuals
hist(resid)
plot(fitted(GSHlm), residuals(GSHlm),
     main = "Residuals vs Fitted values")
cor(resid, fitted.values(GSHlm))

# look at residuals vs predictor
plot(pheno$Heart_GSH, residuals(GSHlm),
     main = "Residuals vs GSH values")
cor(resid, pheno$Heart_GSH) # 1.183569e-16

# look at residuals vs response
plot(pheno$Heart_GSSG, residuals(GSHlm),
     main = "Residuals vs GSSG values")
cor(resid, pheno$Heart_GSSG) # 0.9482969

# look at GSH/GSSG ratio and residuals
hist(pheno$Heart_GSHGSSGRatio)
plot(pheno$Heart_GSHGSSGRatio, resid,  
     main = "Residuals vs GSH/GSSG ratio values\nfor lm(GSSG ~ GSH)")
cor(resid, pheno$Heart_GSHGSSGRatio) # -0.2577074

# remember to turn the PDF device off!
dev.off()


