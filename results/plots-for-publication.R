# scatterplots and boxplots for publication

# load tidyverse library for access to ggplot2
library(tidyverse)
library(patchwork)

# load data
load("../data/Heart-GSSG-enviornment.RData")
covariates <- read.csv("../data/R01_GSH_DO_covar.csv")


pdf(file = "../results/pub-quality-raw-plots.pdf")

control$pheno %>%
  ggplot(aes(x = Heart_GSH)) +
  geom_histogram(binwidth = .5) +
  labs(x = "GSH",
       title = "Histogram of Heart GSH Measurements")

control$pheno %>%
  ggplot(aes(x = Heart_GSSG)) +
  geom_histogram() +
  labs(x = "GSSG",
       title = "Histogram of Heart GSSG Measurements")

control$pheno %>%
  ggplot(aes(sample = Heart_GSH)) +
  stat_qq_line() +
  stat_qq() +
  labs(x = "theoretical quantiles",
       y = "sample quantiles",
       title = "Quantile-quantile Plot of Heart GSH Measurements")

control$pheno %>%
  ggplot(aes(sample = Heart_GSSG)) +
  stat_qq_line() +
  stat_qq() +
  labs(x = "theoretical quantiles",
       y = "sample quantiles",
       title = "Quantile-quantile Plot of Heart GSSG Measurements")

control$pheno <- control$pheno %>% left_join(covariates, by = "id")

control$pheno %>%
  ggplot(aes(x = sex.y, y = Heart_GSH, color = sex.y)) +
  geom_boxplot(outlier.shape = NA)  +
  labs(x = "",
       y = "GSH",
       title = "Boxplots of Reduced Heart Glutathione (GSH) Measurements by Sex") +
  geom_jitter(width = .15) +
  guides(color=guide_legend(title="Sex"))

control$pheno %>%
  ggplot(aes(x = sex.y, y = Heart_GSSG, color = sex.y)) +
  geom_boxplot(outlier.shape = NA)  +
  labs(x = "",
       y = "GSSG",
       title = "Boxplots of Heart Glutathione Disulfide (GSSG) by Sex") +
  geom_jitter(width = .15) +
  guides(color=guide_legend(title="Sex"))

control$pheno %>%
  filter(generation.y == 30 | generation.y == 32) %>%
  ggplot(aes(x=as.factor(generation.y), y=Heart_GSH)) +
  geom_boxplot(outlier.shape = NA)  +
  geom_point(aes(color = sex.y), position=position_jitter(width=0.15, height=0)) +
  labs(x = "generation",
       y = "GSH",
       title = "Boxplots of Heart Glutathione (GSH) by Sex and Generation") +
  guides(color=guide_legend(title="Sex"))

control$pheno %>%
  filter(generation.y == 30 | generation.y == 32) %>%
  ggplot(aes(x=as.factor(generation.y), y=Heart_GSSG)) +
  geom_boxplot(outlier.shape = NA)  +
  geom_point(aes(color = as.factor(sex.y)), 
             position=position_jitter(width=0.15, height=0)) +
  labs(x = "generation",
       y = "GSSG",
       title = "Boxplots of Heart Glutathione Disulfide (GSSG) by Sex and Generation") +
  guides(color=guide_legend(title="Sex"))

# separate linear models by sex
control$pheno %>%
  ggplot(aes(x = Heart_GSH, y = Heart_GSSG, color = sex.y)) +
  geom_point() +
  stat_smooth(method = "lm") +
  labs(x = "GSH",
       y = "GSSG",
       title = "GSH vs GSSG by Sex") +
  guides(color=guide_legend(title="Sex"))

# separate linear models by generation
control$pheno %>%
  filter(generation.y == 30 | generation.y == 32) %>%
  ggplot(aes(x = Heart_GSH, y = Heart_GSSG, color = as.factor(generation.y))) +
  geom_point() +
  stat_smooth(method = "lm")+
  labs(x = "GSH",
       y = "GSSG",
       title = "GSH vs GSSG by Generation") +
  guides(color=guide_legend(title="Generation"))

# remember to turn the PDF device off!
dev.off()

# create a two-panel plot for the original phenotypes and a 3-panel plot for the
# derivations
histGSH <- control$pheno %>%
  ggplot(aes(x = Heart_GSH)) +
  geom_histogram(binwidth = .5) +
  labs(x = "GSH (nmol/mg)")

histGSSG <- control$pheno %>%
  ggplot(aes(x = Heart_GSSG)) +
  geom_histogram() +
  labs(x = "GSSG (nmol/mg)")

histTotal <- control$pheno %>%
  ggplot(aes(x = Heart_TotalGSH)) +
  geom_histogram(binwidth = .5) +
  labs(x = "Total GSH (nmol/mg)")

histRatio <- control$pheno %>%
  ggplot(aes(x = Heart_GSHGSSGRatio)) +
  geom_histogram() +
  labs(x = "GSSG/GSH ratio (nmol/mg)")

histRedoxpot <- control$pheno %>%
  ggplot(aes(x = redoxPot)) +
  geom_histogram(binwidth) +
  labs(x = "Redox potential")

histGSH + histGSSG + histTotal + histRatio + histRedoxpot +
  plot_annotation(title = 'Glutatathione (GSH) and reduced glutathione (GSSG) data distribution')

# set graphical parameters to 3 rows, 2 columns
par(mfrow=c(3,2))

scanPlotGSH <- plot_scan1(zScanSex, map = control$pmap,
           lodcolumn = "zHeart_GSH", main = "GSH",
           sub = "with sex as covariate")

scanPlotGSG <- plot_scan1(zScanSex, map = control$pmap,
                          lodcolumn = "zHeart_GSSG", main = "GSSG",
                          sub = "with sex as covariate")

scanPlotTotal <- plot_scan1(zScanSex, map = control$pmap,
                          lodcolumn = "zHeart_TotalGSH", main = "Total GSH",
                          sub = "with sex as covariate")

scanPlotRatio <- plot_scan1(zScanSex, map = control$pmap,
                          lodcolumn = "zGSHGSSG_Ratio", main = "GSH/GSSG Ratio",
                          sub = "with sex as covariate")

scanPlotRedoxPot <- plot_scan1(zScanSex, map = control$pmap,
                          lodcolumn = "zredoxPot", main = "Redox potential",
                          sub = "with sex as covariate")
