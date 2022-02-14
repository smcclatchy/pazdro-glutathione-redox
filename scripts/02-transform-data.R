# transform phenotypes in the cross object (called control)
# use both log10 and rankZ transformations
# plot transformed variables

# load tidyverse including ggplot2
library(tidyverse)

# load data only if you don't already have it in your RStudio Environment
# remove the hashtag and run the code below if you don't have it

# load("../data/Heart-GSSG-enviornment.RData")

# log base 10 transform
log10(control$pheno[, 2:4])

# log transform each individual phenotype
logHeart_GSH <- log10(control$pheno[, "Heart_GSH"])
logHeart_GSSG <- log10(control$pheno[, "Heart_GSSG"])
logHeart_TotalGSH <- log10(control$pheno[, "Heart_TotalGSH"])

# put them back into the phenotypes in the cross object
control$pheno <- cbind(control$pheno, 
               logHeart_GSH, 
               logHeart_GSSG,
               logHeart_TotalGSH)

head(control$pheno)

# now plot them as for raw phenotypes
# output as pdf
pdf(file = "../results/log-plots.pdf")
as.data.frame(control$pheno) %>% 
  ggplot(aes(x= `logHeart_GSH`, `logHeart_GSSG`)) + 
  geom_point() +
  scale_x_log10() +
  scale_y_log10()

as.data.frame(control$pheno) %>%
  ggplot(aes(x = `logHeart_GSH`)) + 
  geom_histogram() +
  scale_x_log10() 

as.data.frame(control$pheno) %>% 
  ggplot(aes(`logHeart_GSSG`)) + 
  geom_histogram() +
  scale_x_log10() 

qqnorm(logHeart_GSH)
qqline(logHeart_GSH)
qqnorm(logHeart_GSSG)
qqline(logHeart_GSSG)

# remember to turn the PDF device off!
dev.off()

# rank Z transform variable
zHeart_GSH <- rankZ(control$pheno[, "Heart_GSH"])
zHeart_GSSG <- rankZ(control$pheno[, "Heart_GSSG"])
zHeart_TotalGSH <- rankZ(control$pheno[, "Heart_TotalGSH"])

# put the rankZ-transformed variables back into the phenotypes
control$pheno <- cbind(control$pheno, 
      zHeart_GSH, 
      zHeart_GSSG,
      zHeart_TotalGSH)
head(control$pheno)

# plot rankZ transformed variables
# output as pdf

pdf(file = "../results/rankZ-plots.pdf")
hist(zHeart_GSH)
hist(zHeart_GSSG)
plot(zHeart_GSH, zHeart_GSSG)
qqnorm(zHeart_GSH)
qqline(zHeart_GSH)
qqnorm(zHeart_GSSG)
qqline(zHeart_GSSG)

# remember to turn the PDF device off!
dev.off()
# much better!