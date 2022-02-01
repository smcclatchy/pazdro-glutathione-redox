library(tidyverse)
library(qtl2)

probs <- readRDS("~/Projects/Pazdro/beccas_qtl2/data/Pazdro_GigaMUGA_genoprobs_qced_8state.rds")
pheno <- read.csv("~/Projects/Pazdro/beccas_qtl2/data/R01_GSH_DO_pheno.csv")
pheno <- pheno %>% left_join(covar, by = "id")
covar <- read.csv("~/Projects/Pazdro/beccas_qtl2/data/R01_GSH_DO_covar.csv")
covar$sex
covar <- covar %>% mutate(sex = ifelse(sex == "M", 1, 0))
gmap <- read_csv("~/Projects/Pazdro/beccas_qtl2/data/R01_GSH_DO_gmap_sorted.csv")
pmap <- read_csv("~/Projects/Pazdro/beccas_qtl2/data/R01_GSH_DO_pmap_sorted.csv")

for(i in c(2,3,4,5,7)){
  qtl = scan1(genoprobs = probs, pheno = pheno[i], addcovar = pheno$sex)
  plot_scan1(x = qtl, map = gmap, main = paste0("Genomescan for ", names(pheno[i]), " [scan1]"))
}
