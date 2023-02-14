library(tidyverse)

## functions for scatterplot matrix
## put (absolute) correlations on the upper panels,
## with size proportional to the correlations.
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y, use = "complete.obs"))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}

## put histograms on the diagonal
panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...)
}

# load heart, liver and kidney phenotype data
load("../data/Heart-GSSG-enviornment.RData")
liver_pheno <- read_csv("../../beccas_qtl2/data/R01_GSH_DO_pheno.csv")
kidney_pheno <- read_csv(file = "../data/kidney_pheno_covar.csv")

# change variables from character to numeric
liver_pheno <- liver_pheno %>% mutate_if(is.character, as.numeric)
kidney_pheno <- kidney_pheno %>% mutate_if(is.character, as.numeric)

# place animal ID into liver and kidney data
liver_pheno <- liver_pheno %>% mutate("id" = rownames(liver_pheno))
kidney_pheno <- kidney_pheno %>% mutate("id" = rownames(kidney_pheno))

# calculate redox potential for each tissue
pheno <- pheno %>% mutate("heart_redox_potential" = (-264 + 31 * log10(Heart_GSSG/Heart_GSH^2)))
liver_pheno <- liver_pheno %>% mutate("liver_redox_potential" = (-264 + 31 * log10(Liver_Unadj_GSSG/Liver_GSH^2)))
kidney_pheno <- kidney_pheno %>% mutate("kidney_redox_potential" = (-264 + 31 * log10(Kidney_GSSG/Kidney_GSH^2)))


# combine all into a single data frame
pheno$id <- str_remove(pheno$id, "DO-")
kidneyLiver <- kidney_pheno %>% left_join(liver_pheno, by = "id")
kidneyLiverHeart <- kidneyLiver %>% left_join(pheno, by = "id")


# scatterplot matrix of GSH, GSSG and total GSH 
pairs(kidneyLiverHeart[, c(2,3,11,28,32,33)], upper.panel = panel.cor, diag.panel = panel.hist)

# explore tissue by tissue
# kidney and liver
pairs(kidneyLiverHeart[, c(2:3, 6, 11, 28, 31)], upper.panel = panel.cor, diag.panel = panel.hist)

# kidney and heart
pairs(kidneyLiverHeart[, c(2:3, 6, 10, 32, 33)], upper.panel = panel.cor, diag.panel = panel.hist)

# liver and heart
pairs(kidneyLiverHeart[, c(10:12, 30:32)], upper.panel = panel.cor, diag.panel = panel.hist)

# explore linear models of GSH in kidney and heart
# rankZ transform both
kidneyLiverHeart <- kidneyLiverHeart %>% 
  mutate("zkidney_GSH" = rankZ(Kidney_GSH), "zHeart_GSH" = rankZ(Heart_GSH))

summary(lm(kidneyLiverHeart$Heart_GSH ~ kidneyLiverHeart$Kidney_GSH))
gshModel <- lm(kidneyLiverHeart$Heart_GSH ~ kidneyLiverHeart$Kidney_GSH)
hist(gshModel$residuals) # residuals not normally distributed

# with rankZ transformed data
summary(lm(kidneyLiverHeart$zHeart_GSH ~ kidneyLiverHeart$zkidney_GSH))

# now with no intercept
summary(lm(kidneyLiverHeart$Heart_GSH ~ kidneyLiverHeart$Kidney_GSH - 1))
summary(lm(kidneyLiverHeart$zHeart_GSH ~ kidneyLiverHeart$zkidney_GSH - 1))

# scatterplot and overlay linear model
kidneyLiverHeart %>% ggplot(aes(Kidney_GSH, Heart_GSH)) + 
  geom_point() + geom_smooth(method = "lm")
# now with no intercept
kidneyLiverHeart %>% ggplot(aes(Kidney_GSH, Heart_GSH)) + 
  geom_point() + geom_smooth(method = "lm", formula = 'y ~ x - 1')

# rankZ transformed
kidneyLiverHeart %>% ggplot(aes(zkidney_GSH, zHeart_GSH)) + 
  geom_point() + geom_smooth(method = "lm")
kidneyLiverHeart %>% ggplot(aes(zkidney_GSH, zHeart_GSH)) + 
  geom_point() + geom_smooth(method = "lm", formula = 'y ~ x - 1')






