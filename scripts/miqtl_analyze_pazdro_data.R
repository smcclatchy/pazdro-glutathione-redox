options(stringsAsFactors = FALSE)
library(tidyverse)

## qtl2
library(qtl2)
# devtools::install_github("gkeele/miqtl")
library(miqtl)

genoprobs36 <- readRDS("data/Pazdro_GigaMUGA_genoprobs_qced_36state_sorted.rds")
## Convert to miQTL
physical_map_df <- read.csv("data/physical_map_sorted.csv", header = TRUE) 
physical_map <- physical_map_df %>%
  qtl2convert::map_df_to_list(pos_column = "pos")
genetic_map_df <- read.csv("data/genetic_map_sorted.csv", header = TRUE)
genetic_map <- genetic_map_df %>%
  qtl2convert::map_df_to_list(pos_column = "pos")
holder_cross_object <- list(gmap = genetic_map,
                            pmap = physical_map)
## Make genomecache for miQTL
## 13 GB (basically the 36-state and 8-state genoprobs combined)
convert.qtl2.to.HAPPY(qtl2.object = genoprobs36, 
                      cross.object = holder_cross_object,
                      HAPPY.output.path = "data/Pazdro_genomecache", 
                      diplotype.order = "qtl2")

## Mapping in qtl2
## Compare with miQTL
genoprobs8 <- readRDS("data/Pazdro_GigaMUGA_genoprobs_qced_8state_sorted.rds")
pheno_dat <- read.csv("data/R01_GSH_DO_pheno.csv", header = TRUE)
pheno <- pheno_dat %>%
  mutate(zLiverGSH = rankZ(Liver_GSH), # add whichever phenotypes you want here
         zLiverAdjGSSG = rankZ(Liver_Adj_GSSG),
         zLiverTotGSH = rankZ(Liver_Adj_Total_GSH),
         zLiverNADP = rankZ(Liver_NADP),
         zLiverNADPH = rankZ(Liver_NADPH),
         zLiverNADPHRatio = rankZ(Liver_NADP_NADPH_Ratio),
         zLiverGSHRatio = rankZ(Liver_Adj_GSH_GSSG_Ratio)
         ) %>%
  column_to_rownames("id") %>%
  as.matrix
covar_dat <- read.csv("data/R01_GSH_DO_covar.csv", header = TRUE) %>%
  mutate(generation = as.factor(generation)) %>%
  column_to_rownames("id")
covar_mat <- cbind(model.matrix(~sex, data = covar_dat)[,-1, drop = FALSE],
                   model.matrix(~generation, data = covar_dat)[,-1, drop = FALSE])
rownames(covar_mat) <- rownames(covar_dat)

K <- calc_kinship(genoprobs8, type = "loco")
map_df <- read.csv("data/physical_map_sorted.csv", header = TRUE) 
map <- map_df %>%
  qtl2convert::map_df_to_list(pos_column = "pos")

qtl2_scan <- scan1(genoprobs = genoprobs8, 
                   pheno = pheno[,21:27], # change numbers to reflect phenotypes you want to scan
                   addcovar = covar_mat[,1,drop=FALSE],
                   kinship = K)

plot(qtl2_scan, map)
find_peaks(qtl2_scan, map = map, threshold = 6, sort_by = "pos")

for (i in 1:length(dimnames(qtl2_scan)[[2]])) {
  plot(qtl2_scan, map, chr = 2, # change chromosome number as needed
       main = dimnames(qtl2_scan)[[2]][i],
       lodcolumn = i)
  }

for (i in 1:length(dimnames(qtl2_scan)[[2]])) {
  plot(qtl2_scan, map, chr = 14, # change chromosome number as needed
       main = dimnames(qtl2_scan)[[2]][i],
       lodcolumn = i)
}

for (i in 1:length(dimnames(qtl2_scan)[[2]])) {
  plot(qtl2_scan, map, chr = 16, # change chromosome number as needed
       main = dimnames(qtl2_scan)[[2]][i],
       lodcolumn = i)
}


## Mapping in miqtl
library(miqtl)

miqtl_pheno_dat <- pheno %>%
  as.data.frame %>%
  rownames_to_column("SUBJECT.NAME") %>%
  left_join(covar_mat %>%
              as.data.frame %>%
              rownames_to_column("SUBJECT.NAME"))

# include generation as a covariate
# miqtl_pheno_dat <- pheno %>% 
#   as.data.frame %>%
#   rownames_to_column("SUBJECT.NAME") %>%
#   left_join(covar_mat %>%
#               as.data.frame %>%
#               rownames_to_column("SUBJECT.NAME") %>%
#               mutate(generation = as.factor(generation)))

miqtl_scan <- scan.h2lmm(genomecache = "data/Pazdro_genomecache/", 
                         data = miqtl_pheno_dat, 
                         formula = rankZ(Liver_GSH) ~ 1 + sexM, # change phenotype here
                         K = K[["2"]], # change chromosome here
                         model = "additive", 
                         use.multi.impute = FALSE, 
                         chr = 2) # change chromosome here too, and below in genome plotter
genome.plotter.chr(miqtl_scan, chr = 2, use.lod = TRUE, main = "miQTL (run like qtl2)")

miqtl_mi11_scan <- scan.h2lmm(genomecache = "data/Pazdro_genomecache/", 
                         data = miqtl_pheno_dat, 
                         formula = rankZ(Liver_GSH) ~ 1 + sexM, # change phenotype here
                         K = K[["2"]],  # change chromosome here
                         model = "additive", 
                         use.multi.impute = TRUE, 
                         num.imp = 11,
                         chr = 2) # change chromosome here too, and below in genome plotter
genome.plotter.chr(miqtl_mi11_scan, chr = 2, use.lod = TRUE, main = "miQTL (11 imputations)", y.max.manual = 10)
# change y.max.manual to maximum y value shown on previous plot (might be only 7, 9, etc)



