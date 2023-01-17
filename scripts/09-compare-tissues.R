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
liver_pheno <- read_csv("R01_GSH_DO_pheno.csv")
kidney_pheno <- read_csv(file = "kidney_pheno_covar.csv")

# change variables from character to numeric
liver_pheno <- liver_pheno %>% mutate_if(is.character, as.numeric)
kidney_pheno <- kidney_pheno %>% mutate_if(is.character, as.numeric)

# place animal ID into liver and kidney data
liver_pheno <- liver_pheno %>% mutate("id" = rownames(liver_pheno))
kidney_pheno <- kidney_pheno %>% mutate("id" = rownames(kidney_pheno))

# combine all into a single data frame
kidneyLiver <- kidney_pheno %>% left_join(liver_pheno, by = "id")
kidneyLiverHeart <- kidneyLiver %>% left_join(pheno, by = "id")

# scatterplot matrix of GSH, GSSG and total GSH 
pairs(kidneyLiverHeart[, c(1:3, 10:12, 30:32)], upper.panel = panel.cor, diag.panel = panel.hist)

# explore tissue by tissue
# kidney and liver
pairs(kidneyLiverHeart[, c(1:3, 10:12)], upper.panel = panel.cor, diag.panel = panel.hist)

# kidney and heart
pairs(kidneyLiverHeart[, c(1:3, 30:32)], upper.panel = panel.cor, diag.panel = panel.hist)

# liver and heart
pairs(kidneyLiverHeart[, c(10:12, 30:32)], upper.panel = panel.cor, diag.panel = panel.hist)
