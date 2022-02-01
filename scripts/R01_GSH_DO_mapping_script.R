# R01 GSH DO Mapping Code 
# Updated July 2020
# Becca Gould 

#Make a folder under users (for me, my user is "becca") and title it based on your project. For mine, it's R01_GSH_DO_mapping. Then make a data, results, scripts, and docs folder.

  setwd("/users/becca/R01_GSH_DO_mapping/data")

#confirm where you are through getwd()
  getwd() 

#I can use "~" instead of "/users/becca/" everytime as it represents my home base

#load the command line tools - see https://github.com/Rdatatable/data.table/wiki/Installation for more information - must do every time you open up the Rproject!
library(qtl2)
library (tidyverse)
library (readxl)
library(yaml)
library(devtools)
library(RSQLite)
library(jsonlite)
library (data.table)
library (RcppEigen)
library (RSQLite)
  
  
#For all plots, use these parameters to keep everything consistent. You can adjust them accordingly.
  #par(mar=c(4.1, 4.1, 2.6, 2.6))


####################################################
## Helpful Sites
####################################################

#https://kbroman.org/qtl2/assets/vignettes/user_guide.html#QTL_analysis_in_Diversity_Outbred_mice
#https://smcclatchy.github.io/mapping/ 


####################################################
## Editing the phenotype file to make it R/qtl2-friendly
####################################################
#to have the phenotype file for reference - can be used when plotting the data to see if it needs to be transformed
  pheno <- read.csv(file = "~/R01_GSH_DO_mapping/data/R01_GSH_DO_pheno_covar.csv", header = TRUE)
#make row names the ID of each sample - did not end up using for genome scans
  rownames(pheno) <- pheno$id
#checking pheno file
  pheno[1:10,]

#change sexes to numeric variables
  pheno$sex[pheno$sex == "M"] <- 1
  pheno$sex[pheno$sex == "F"] <- 0
#check pheno file
  pheno[1:10,]


####################################################
## Read in the control file (gm.json)
####################################################

#Load in the control file to tell R/qtl2 to read the data and title it according to your project. For mine, it's R01_GSH_DO_QTLdata. 
#For this to work, all of the files in the control file need to be in the folder. Ex: if calling for the "genoprobs" file, it needs to actually be in the folder.
  R01_GSH_DO_QTLdata <- read_cross2(file = "~/R01_GSH_DO_mapping/data/R01_GSH_DO_control.json")

####################################################
## Genotype probabilities and allele probabilities - provided by Belinda and Vivek
####################################################

#read in the genoprobs file that is sorted by chromosomes in numerical order - the 8state.rds is the allele probabilities, the 32state.rds is the genotype probabilities
#^this is actually the ALlELE probabilities, but for simplicity, we will call it "probs"
  probs <- readRDS("~/R01_GSH_DO_mapping/data/Pazdro_GigaMUGA_genoprobs_qced_8state_sorted.rds")

  nrow(data.frame(R01_GSH_DO_QTLdata$gmap[1]))
  #should be 10415
  dim(probs[[1]])
  #should be 347 individuals, 8 alleles, 10415 markers


####################################################
## Variant files
####################################################

#Will need these for the final lesson episodes on SNP association mapping and QTL analysis in Diversity Outbred mice. Make sure they are the most updated versions!
  query_variants <- create_variant_query_func("~/R01_GSH_DO_mapping/data/cc_variants.sqlite")
  query_genes_mgi <- create_gene_query_func("~/R01_GSH_DO_mapping/data/mouse_genes_mgi.sqlite")
  query_genes <- create_gene_query_func("~/R01_GSH_DO_mapping/data/mouse_genes.sqlite")

####################################################
## Calculating kinship
####################################################

#calculate the kinship loco
#you can increase the cores amount if you have more cores in your computer. For mine, I have 18 cores available so to speed it up, I'll use 10 of them.
  kinship_loco <- calc_kinship(probs = probs, "loco", use_allele_probs = TRUE, cores = 10)

#Create the r plot of the kinship matrix
#pdf('output/rplot_kinship.pdf')
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  image(1:nrow(kinship_loco[[1]]), 1:ncol(kinship_loco[[1]]), kinship_loco[[1]][,ncol(kinship_loco[[1]]):1], xlab = "Samples", 
        ylab = "Samples", yaxt = "n", main = "Kinship between samples", 
        breaks = 0:100/100, col = heat.colors(length(0:100) - 1))


####################################################
## add covariates
####################################################

#adding sex only as covariate
  sex = model.matrix(~sex, data = pheno)[,-1]

#adding sex and generation as covariates
  sexgen = model.matrix(~sex + generation, data = pheno)[,-1]

#adding generation only as covariates
  gen = model.matrix(~generation, data = pheno)[,-1]


####################################################
## Add X covariates
## did not do because accounting for sex as an added covariate
####################################################
Xcovar <- get_x_covar(R01_GSH_DO_QTLdata)

####################################################
## checking if data should be transformed
####################################################

#gives you the names of each phenotype
  R01_GSH_DO_QTLdata[["pheno"]]
  
#setting the parameters for the plots
  par(mar=c(4.1, 4.1, 2.6, 2.6))

#For Liver GSH -- Conclusion: use log transformation
  hist(x = pheno[,"Liver_GSH"], main = "Liver GSH Histogram")
  pheno$logLiverGSH = log(R01_GSH_DO_QTLdata$pheno[,"Liver_GSH"])
  hist(x = pheno[,"logLiverGSH"], main = "Log - Liver GSH Histogram")

#For Liver Adj GSSG -- Conclusion: Use log transformation
  hist(x = pheno[,"Liver_Adj_GSSG"], main = "Liver Adj GSSG Histogram")
  pheno$logLiverAdjGSSG = log(R01_GSH_DO_QTLdata$pheno[,"Liver_Adj_GSSG"])
  hist(x = pheno[,"logLiverAdjGSSG"], main = "Log - Liver Adj GSSG Histogram")

#For Liver Adj Total GSH -- Conclusion: Use log transformation
  hist(x = pheno[,"Liver_Adj_Total_GSH"], main = "Liver Adj Total GSH Histogram")
  pheno$logLiverAdjTotalGSH = log(R01_GSH_DO_QTLdata$pheno[,"Liver_Adj_Total_GSH"])
  hist(x = pheno[,"logLiverAdjTotalGSH"], main = "Log - Liver Adj Total GSH Histogram")

#For Liver Adj GSH.GSSG -- Conclusion: Use log transformation
  hist(x = pheno[,"Liver_Adj_GSH_GSSG_Ratio"], main = "Liver Adj GSH/GSSG Histogram")
  pheno$logLiverAdjGSH_GSSGRatio = log(R01_GSH_DO_QTLdata$pheno[,"Liver_Adj_GSH_GSSG_Ratio"])
  hist(x = pheno[,"logLiverAdjGSH_GSSGRatio"], main = "Log - Liver Adj GSH/GSSG Histogram")

#For Liver NADH -- Conclusion: Use log transformation
  hist(x = pheno[,"Liver_NADH"], main = "Liver NADH Histogram")
  pheno$logLiverNADH = log(R01_GSH_DO_QTLdata$pheno[,"Liver_NADH"])
  hist(x = pheno[,"logLiverNADH"], main = "Log - Liver NADH Histogram")

#For Liver NADP -- Conclusion: Use log transformation
  hist(x = pheno[,"Liver_NADP"], main = "Liver NADP Histogram")
  pheno$logLiverNADP = log(R01_GSH_DO_QTLdata$pheno[,"Liver_NADP"])
  hist(x = pheno[,"logLiverNADP"], main = "Log - Liver NADP Histogram")

#For Liver NADPH -- Conclusion: Use log transformation
  hist(x = pheno[,"Liver_NADPH"], main = "Liver NADPH Histogram")
  pheno$logLiverNADPH = log(R01_GSH_DO_QTLdata$pheno[,"Liver_NADPH"])
  hist(x = pheno[,"logLiverNADPH"], main = "Log - Liver NADPH Histogram")

#For Liver NADP/NADPH -- Conclusion: Use log transformation
  hist(x = pheno[,"Liver_NADP_NADPH_Ratio"], main = "Liver NADP/NADPH Histogram")
  pheno$logLiverNADP_NADPHRatio = log(R01_GSH_DO_QTLdata$pheno[,"Liver_NADP_NADPH_Ratio"])
  hist(x = pheno[,"logLiverNADP_NADPHRatio"], main = "Log - Liver NADP/NADPH Histogram")

#For AST -- Conclusion: Use log transformation
  hist(x = pheno[,"AST"], main = "AST Histogram")
  pheno$logAST = log(R01_GSH_DO_QTLdata$pheno[,"AST"])
  hist(x = pheno[,"logAST"], main = "Log - AST Histogram")

#For ALT -- Conclusion: Use log transformation
  hist(x = pheno[,"ALT"], main = "ALT Histogram")
  pheno$logALT = log(R01_GSH_DO_QTLdata$pheno[,"ALT"])
  hist(x = pheno[,"logALT"], main = "Log - ALT Histogram")

#For AST/ALT -- Conclusion: code didn't work
  hist(x = pheno[,"AST.ALT"], main = "AST/ALT Histogram")
  pheno$logAST_ALT = log(R01_GSH_DO_QTLdata$pheno[,"AST.ALT"])
  hist(x = pheno[,"logAST.ALT"], main = "Log - AST/ALT Histogram")

#For BUN -- Conclusion: Use log transformation
  hist(x = pheno[,"BUN"], main = "BUN Histogram")
  pheno$logBUN = log(R01_GSH_DO_QTLdata$pheno[,"BUN"])
  hist(x = pheno[,"logBUN"], main = "Log - BUN Histogram")

#For Glucose -- Conclusion: Use log transformation
  hist(x = pheno[,"Glucose"], main = "Glucose Histogram")
  pheno$logGlucose = log(R01_GSH_DO_QTLdata$pheno[,"Glucose"])
  hist(x = pheno[,"logGlucose"], main = "Log - Glucose Histogram")

#For Liver Weight -- Conclusion: Use log transformation
  hist(x = pheno[,"Liver_Weight"], main = "Liver Weight Histogram")
  pheno$logLiver_Weight = log(R01_GSH_DO_QTLdata$pheno[,"Liver_Weight"])
  hist(x = pheno[,"logLiver_Weight"], main = "Log - Liver Weight Histogram")

#For Initial Weight -- Conclusion: Use log transformation
  hist(x = pheno[,"Initial_Weight"], main = "Initial Weight Histogram")
  pheno$logInitial_Weight = log(R01_GSH_DO_QTLdata$pheno[,"Initial_Weight"])
  hist(x = pheno[,"logInitial_Weight"], main = "Log - Initial Weight Histogram")

#For Final Weight -- Conclusion: Do not transform
  hist(x = pheno[,"Final_Weight"], main = "Final Weight Histogram")
  pheno$logFinal_Weight = log(R01_GSH_DO_QTLdata$pheno[,"Final_Weight"])
  hist(x = pheno[,"logFinal_Weight"], main = "Log - Final Weight Histogram")

#For Liver Weight/Final Weight -- Conclusion: Use log transformation
  hist(x = pheno[,"Liver_Weight_Final_Weight_Ratio"], main = "Liver Wt/Final Wt Histogram")
  pheno$logLiverWeight_FinalWeight = log(R01_GSH_DO_QTLdata$pheno[,"Liver_Weight_Final_Weight_Ratio"])
  hist(x = pheno[,"logLiverWeight_FinalWeight"], main = "Log - Liver Wt/Final Wt Histogram")

#For Liver Unadj GSSG -- Conclusion: Use log transformation
  hist(x = pheno[,"Liver_Unadj_GSSG"], main = "Liver Unadj GSSG Histogram")
  pheno$logLiverUnadjGSSG = log(R01_GSH_DO_QTLdata$pheno[,"Liver_Unadj_GSSG"])
  hist(x = pheno[,"logLiverUnadjGSSG"], main = "Log - Liver Unadj GSSG Histogram")

#For Liver Unadj Total GSH -- Conclusion: Use log transformation
  hist(x = pheno[,"Liver_Unadj_Total_GSH"], main = "Liver Unadj Total GSH Histogram")
  pheno$logLiverUnadjTotalGSH = log(R01_GSH_DO_QTLdata$pheno[,"Liver_Unadj_Total_GSH"])
  hist(x = pheno[,"logLiverUnadjTotalGSH"], main = "Log - Liver Unadj Total GSH Histogram")

#For Liver Unadj GSH.GSSG -- Conclusion: Use log transformation
  hist(x = pheno[,"Liver_Unadj_GSH_GSSG_Ratio"], main = "Liver Unadj GSH/GSSG Histogram")
  pheno$logLiverUnadjGSH_GSSGRatio = log(R01_GSH_DO_QTLdata$pheno[,"Liver_Unadj_GSH_GSSG_Ratio"])
  hist(x = pheno[,"logLiverUnadjGSH_GSSGRatio"], main = "Log - Liver Unadj GSH/GSSG Histogram")

####################################################
## Perform the genome scans with permutation tests
####################################################

#for all pheno - did not end up using because my pheno file has sex + generation in it - do not use!
#qtlscan <- scan1(genoprobs = probs, pheno = R01_GSH_DO_QTLdata$pheno, kinship = kinship_loco, addcovar = sexgen, Xcovar = Xcovar, cores=10)

#if not using transformed data, could call from control file
#qtlscan_pheno <- scan1(genoprobs = probs, pheno = R01_GSH_DO_QTLdata$pheno[,"____"], kinship = kinship_loco, addcovar = sexgen, cores=10)

#for Liver GSH -- using transformed data
  qtlscan_LiverGSH_sexgen <- scan1(genoprobs = probs, pheno = pheno["logLiverGSH"], kinship = kinship_loco, addcovar = sexgen, cores=10)
  qtlscan_LiverGSH_sex <- scan1(genoprobs = probs, pheno = pheno["logLiverGSH"], kinship = kinship_loco, addcovar = sex, cores=10)
  qtlscan_LiverGSH_gen <- scan1(genoprobs = probs, pheno = pheno["logLiverGSH"], kinship = kinship_loco, addcovar = gen, cores=10)
  
  perm_LiverGSH_sexgen <- scan1perm(genoprobs = probs, pheno = pheno["logLiverGSH"], addcovar = sexgen, n_perm = 1000, cores=10)
  perm_LiverGSH_sex <- scan1perm(genoprobs = probs, pheno = pheno["logLiverGSH"], addcovar = sex, n_perm = 1000, cores=10)
  perm_LiverGSH_gen <- scan1perm(genoprobs = probs, pheno = pheno["logLiverGSH"], addcovar = gen, n_perm = 1000, cores=10)
  
#for Liver Adj GSSG -- using transformed data
  qtlscan_LiverAdjGSSG_sexgen <- scan1(genoprobs = probs, pheno = pheno["logLiverAdjGSSG"], kinship = kinship_loco, addcovar = sexgen, cores=10)
  qtlscan_LiverAdjGSSG_sex <- scan1(genoprobs = probs, pheno = pheno["logLiverAdjGSSG"], kinship = kinship_loco, addcovar = sex, cores=10)
  qtlscan_LiverAdjGSSG_gen <- scan1(genoprobs = probs, pheno = pheno["logLiverAdjGSSG"], kinship = kinship_loco, addcovar = gen, cores=10)
  
  perm_LiverAdjGSSG_sexgen <- scan1perm(genoprobs = probs, pheno = pheno["logLiverAdjGSSG"], addcovar = sexgen, n_perm = 1000, cores=10)
  perm_LiverAdjGSSG_sex <- scan1perm(genoprobs = probs, pheno = pheno["logLiverAdjGSSG"], addcovar = sex, n_perm = 1000, cores=10)
  perm_LiverAdjGSSG_gen <- scan1perm(genoprobs = probs, pheno = pheno["logLiverAdjGSSG"], addcovar = gen, n_perm = 1000, cores=10)
  
#for Liver Adj Total GSH -- using transformed data
  qtlscan_LiverAdjTotalGSH_sexgen <- scan1(genoprobs = probs, pheno = pheno["logLiverAdjTotalGSH"], kinship = kinship_loco, addcovar = sexgen, cores=10)
  qtlscan_LiverAdjTotalGSH_sex <- scan1(genoprobs = probs, pheno = pheno["logLiverAdjTotalGSH"], kinship = kinship_loco, addcovar = sex, cores=10)
  qtlscan_LiverAdjTotalGSH_gen <- scan1(genoprobs = probs, pheno = pheno["logLiverAdjTotalGSH"], kinship = kinship_loco, addcovar = gen, cores=10)
  
  perm_LiverAdjTotalGSH_sexgen <- scan1perm(genoprobs = probs, pheno = pheno["logLiverAdjTotalGSH"], addcovar = sexgen, n_perm = 1000, cores=10)
  perm_LiverAdjTotalGSH_sex <- scan1perm(genoprobs = probs, pheno = pheno["logLiverAdjTotalGSH"], addcovar = sex, n_perm = 1000, cores=10)
  perm_LiverAdjTotalGSH_gen <- scan1perm(genoprobs = probs, pheno = pheno["logLiverAdjTotalGSH"], addcovar = gen, n_perm = 1000, cores=10)
  
#for Liver Adj Total GSH/GSSG Ratio -- using transformed data
  qtlscan_LiverAdjGSH_GSSG_sexgen <- scan1(genoprobs = probs, pheno = pheno["logLiverAdjGSH_GSSGRatio"], kinship = kinship_loco, addcovar = sexgen, cores=10)
  qtlscan_LiverAdjGSH_GSSG_sex <- scan1(genoprobs = probs, pheno = pheno["logLiverAdjGSH_GSSGRatio"], kinship = kinship_loco, addcovar = sex, cores=10)
  qtlscan_LiverAdjGSH_GSSG_gen <- scan1(genoprobs = probs, pheno = pheno["logLiverAdjGSH_GSSGRatio"], kinship = kinship_loco, addcovar = gen, cores=10)
  
  perm_LiverAdjGSH_GSSG_sexgen <- scan1perm(genoprobs = probs, pheno = pheno["logLiverAdjGSH_GSSGRatio"], addcovar = sexgen, n_perm = 1000, cores=10)
  perm_LiverAdjGSH_GSSG_sex <- scan1perm(genoprobs = probs, pheno = pheno["logLiverAdjGSH_GSSGRatio"], addcovar = sex, n_perm = 1000, cores=10)
  perm_LiverAdjGSH_GSSG_gen <- scan1perm(genoprobs = probs, pheno = pheno["logLiverAdjGSH_GSSGRatio"], addcovar = gen, n_perm = 1000, cores=10)
  
#for Liver NADH -- using transformed data
  qtlscan_LiverNADH_sexgen <- scan1(genoprobs = probs, pheno = pheno["logLiverNADH"], kinship = kinship_loco, addcovar = sexgen, cores=10)
  qtlscan_LiverNADH_sex <- scan1(genoprobs = probs, pheno = pheno["logLiverNADH"], kinship = kinship_loco, addcovar = sex, cores=10)
  qtlscan_LiverNADH_gen <- scan1(genoprobs = probs, pheno = pheno["logLiverNADH"], kinship = kinship_loco, addcovar = gen, cores=10)
  
  perm_LiverNADH_sexgen <- scan1perm(genoprobs = probs, pheno = pheno["logLiverNADH"], addcovar = sexgen, n_perm = 100, cores=10)
  perm_LiverNADH_sex <- scan1perm(genoprobs = probs, pheno = pheno["logLiverNADH"], addcovar = sex, n_perm = 1000, cores=10)
  perm_LiverNADH_gen <- scan1perm(genoprobs = probs, pheno = pheno["logLiverNADH"], addcovar = gen, n_perm = 100, cores=10)  
  
#for Liver NADP -- using transformed data
  qtlscan_LiverNADP_sexgen <- scan1(genoprobs = probs, pheno = pheno["logLiverNADP"], kinship = kinship_loco, addcovar = sexgen, cores=10)
  qtlscan_LiverNADP_sex <- scan1(genoprobs = probs, pheno = pheno["logLiverNADP"], kinship = kinship_loco, addcovar = sex, cores=10)
  qtlscan_LiverNADP_gen <- scan1(genoprobs = probs, pheno = pheno["logLiverNADP"], kinship = kinship_loco, addcovar = gen, cores=10)
  
  perm_LiverNADP_sexgen <- scan1perm(genoprobs = probs, pheno = pheno["logLiverNADP"], addcovar = sexgen, n_perm = 1000, cores=10)
  perm_LiverNADP_sex <- scan1perm(genoprobs = probs, pheno = pheno["logLiverNADP"], addcovar = sex, n_perm = 1000, cores=10)
  perm_LiverNADP_gen <- scan1perm(genoprobs = probs, pheno = pheno["logLiverNADP"], addcovar = gen, n_perm = 1000, cores=10)  
  
#for Liver NADPH -- using transformed data
  qtlscan_LiverNADPH_sexgen <- scan1(genoprobs = probs, pheno = pheno["logLiverNADPH"], kinship = kinship_loco, addcovar = sexgen, cores=10)
  qtlscan_LiverNADPH_sex <- scan1(genoprobs = probs, pheno = pheno["logLiverNADPH"], kinship = kinship_loco, addcovar = sex, cores=10)
  qtlscan_LiverNADPH_gen <- scan1(genoprobs = probs, pheno = pheno["logLiverNADPH"], kinship = kinship_loco, addcovar = gen, cores=10)
  
  perm_LiverNADPH_sexgen <- scan1perm(genoprobs = probs, pheno = pheno["logLiverNADPH"], addcovar = sexgen, n_perm = 1000, cores=10)
  perm_LiverNADPH_sex <- scan1perm(genoprobs = probs, pheno = pheno["logLiverNADPH"], addcovar = sex, n_perm = 1000, cores=10)
  perm_LiverNADPH_gen <- scan1perm(genoprobs = probs, pheno = pheno["logLiverNADPH"], addcovar = gen, n_perm = 1000, cores=10)  
  
#for Liver NADP/NADPH Ratio -- using transformed data
  qtlscan_LiverNADP_NADPH_sexgen <- scan1(genoprobs = probs, pheno = pheno["logLiverNADP_NADPHRatio"], kinship = kinship_loco, addcovar = sexgen, cores=10)
  qtlscan_LiverNADP_NADPH_sex <- scan1(genoprobs = probs, pheno = pheno["logLiverNADP_NADPHRatio"], kinship = kinship_loco, addcovar = sex, cores=10)
  qtlscan_LiverNADP_NADPH_gen <- scan1(genoprobs = probs, pheno = pheno["logLiverNADP_NADPHRatio"], kinship = kinship_loco, addcovar = gen, cores=10)
  
  perm_LiverNADP_NADPH_sexgen <- scan1perm(genoprobs = probs, pheno = pheno["logLiverNADP_NADPHRatio"], addcovar = sexgen, n_perm = 1000, cores=10)
  perm_LiverNADP_NADPH_sex <- scan1perm(genoprobs = probs, pheno = pheno["logLiverNADP_NADPHRatio"], addcovar = sex, n_perm = 1000, cores=10)
  perm_LiverNADP_NADPH_gen <- scan1perm(genoprobs = probs, pheno = pheno["logLiverNADP_NADPHRatio"], addcovar = gen, n_perm = 1000, cores=10) 
  
#for Liver Unadj GSSG -- using transformed data
  qtlscan_LiverUnadjGSSG_sexgen <- scan1(genoprobs = probs, pheno = pheno["logLiverUnadjGSSG"], kinship = kinship_loco, addcovar = sexgen, cores=10)
  qtlscan_LiverUnadjGSSG_sex <- scan1(genoprobs = probs, pheno = pheno["logLiverUnadjGSSG"], kinship = kinship_loco, addcovar = sex, cores=10)
  qtlscan_LiverUnadjGSSG_gen <- scan1(genoprobs = probs, pheno = pheno["logLiverUnadjGSSG"], kinship = kinship_loco, addcovar = gen, cores=10)
  
  perm_LiverUnadjGSSG_sexgen <- scan1perm(genoprobs = probs, pheno = pheno["logLiverUnadjGSSG"], addcovar = sexgen, n_perm = 1000, cores=10)
  perm_LiverUnadjGSSG_sex <- scan1perm(genoprobs = probs, pheno = pheno["logLiverUnadjGSSG"], addcovar = sex, n_perm = 1000, cores=10)
  perm_LiverUnadjGSSG_gen <- scan1perm(genoprobs = probs, pheno = pheno["logLiverUnadjGSSG"], addcovar = gen, n_perm = 1000, cores=10)
  
#for Liver Unadj Total GSH -- using transformed data
  qtlscan_LiverUnadjTotalGSH_sexgen <- scan1(genoprobs = probs, pheno = pheno["logLiverUnadjTotalGSH"], kinship = kinship_loco, addcovar = sexgen, cores=10)
  qtlscan_LiverUnadjTotalGSH_sex <- scan1(genoprobs = probs, pheno = pheno["logLiverUnadjTotalGSH"], kinship = kinship_loco, addcovar = sex, cores=10)
  qtlscan_LiverUnadjTotalGSH_gen <- scan1(genoprobs = probs, pheno = pheno["logLiverUnadjTotalGSH"], kinship = kinship_loco, addcovar = gen, cores=10)
  
  perm_LiverUnadjTotalGSH_sexgen <- scan1perm(genoprobs = probs, pheno = pheno["logLiverUnadjTotalGSH"], addcovar = sexgen, n_perm = 100, cores=10)
  perm_LiverUnadjTotalGSH_sex <- scan1perm(genoprobs = probs, pheno = pheno["logLiverUnadjTotalGSH"], addcovar = sex, n_perm = 100, cores=10)
  perm_LiverUnadjTotalGSH_gen <- scan1perm(genoprobs = probs, pheno = pheno["logLiverUnadjTotalGSH"], addcovar = gen, n_perm = 100, cores=10)
  
  #for Liver Unadj Total GSH/GSSG Ratio -- using transformed data
  qtlscan_LiverUnadjGSH_GSSG_sexgen <- scan1(genoprobs = probs, pheno = pheno["logLiverUnadjGSH_GSSGRatio"], kinship = kinship_loco, addcovar = sexgen, cores=10)
  qtlscan_LiverUnadjGSH_GSSG_sex <- scan1(genoprobs = probs, pheno = pheno["logLiverUnadjGSH_GSSGRatio"], kinship = kinship_loco, addcovar = sex, cores=10)
  qtlscan_LiverUnadjGSH_GSSG_gen <- scan1(genoprobs = probs, pheno = pheno["logLiverUnadjGSH_GSSGRatio"], kinship = kinship_loco, addcovar = gen, cores=10)
  
  perm_LiverUnadjGSH_GSSG_sexgen <- scan1perm(genoprobs = probs, pheno = pheno["logLiverUnadjGSH_GSSGRatio"], addcovar = sexgen, n_perm = 1000, cores=10)
  perm_LiverUnadjGSH_GSSG_sex <- scan1perm(genoprobs = probs, pheno = pheno["logLiverUnadjGSH_GSSGRatio"], addcovar = sex, n_perm = 1000, cores=10)
  perm_LiverUnadjGSH_GSSG_gen <- scan1perm(genoprobs = probs, pheno = pheno["logLiverUnadjGSH_GSSGRatio"], addcovar = gen, n_perm = 1000, cores=10)  
  
  
####################################################

#Running a permutation test for the X chromosome specifically
  #chr_lengths <- R01_GSH_DO_QTLdata$gmap
  #Xperm_LiverGSH_sexgen <- scan1perm(genoprobs = probs, pheno = pheno["logLiverGSH"], addcovar = sexgen, Xcovar = Xcovar, n_perm = 10, perm_Xsp = TRUE, chr_lengths = R01_GSH_DO_QTLdata$gmap, cores = 10)
  #Error in sum(lengths[!is_x_chr]) : invalid 'type' (list) of argument

####################################################
## Plotting using permutations
## to get just the plot without permutations, don't add the "abline" section, just use plot_scan1
####################################################

#set the parameters for the plot
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  ylim <- c(0, maxlod(qtlscan_LiverGSH_sex)*1.02)
  #remove ylim from plot_scan1 to get default max y-axis LOD score

#For Liver GSH
  #sex and gen covariates
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  threshold_LiverGSH_sexgen = summary(perm_LiverGSH_sexgen, alpha = c(0.1, 0.05, 0.01))
  plot_scan1(x = qtlscan_LiverGSH_sexgen, map = R01_GSH_DO_QTLdata$gmap,  main = "Genome Scan for Liver GSH - sex and gen", ylim = ylim)
  abline(h = threshold_LiverGSH_sexgen, col = c("red", "blue", "green"), lwd = 2)
  find_peaks(scan1_output = qtlscan_LiverGSH_sexgen, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_LiverGSH_sexgen, alpha = 0.1), peakdrop = 1.8, prob = 0.95, expand2markers = FALSE)
  peaks_LiverGSH_sexgen <- find_peaks(scan1_output = qtlscan_LiverGSH_sexgen, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_LiverGSH_sexgen, alpha = 0.1), peakdrop = 1.8, prob = 0.95, expand2markers = FALSE)
  
  #sex covariate
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  threshold_LiverGSH_sex = summary(perm_LiverGSH_sex, alpha = c(0.1, 0.05, 0.01))
  plot_scan1(x = qtlscan_LiverGSH_sex, map = R01_GSH_DO_QTLdata$gmap,  main = "Genome Scan for Liver GSH - sex", ylim = ylim)
  abline(h = threshold_LiverGSH_sex, col = c("red", "blue", "green"), lwd = 2)
  find_peaks(scan1_output = qtlscan_LiverGSH_sex, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_LiverGSH_sex, alpha = 0.1), peakdrop = 1.8, prob = 0.95, expand2markers = FALSE)
  peaks_LiverGSH_sex <- find_peaks(scan1_output = qtlscan_LiverGSH_sex, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_LiverGSH_sex, alpha = 0.1), peakdrop = 1.8, prob = 0.95, expand2markers = FALSE)
  
  #gen covariate
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  threshold_LiverGSH_gen = summary(perm_LiverGSH_gen, alpha = c(0.1, 0.05, 0.01))
  plot_scan1(x = qtlscan_LiverGSH_gen, map = R01_GSH_DO_QTLdata$gmap,  main = "Genome Scan for Liver GSH - gen", ylim = ylim)
  abline(h = threshold_LiverGSH_gen, col = c("red", "blue", "green"), lwd = 2)
  find_peaks(scan1_output = qtlscan_LiverGSH_gen, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_LiverGSH_gen, alpha = 0.1), peakdrop = 1.8, prob = 0.95, expand2markers = FALSE)
  peaks_LiverGSH_gen <- find_peaks(scan1_output = qtlscan_LiverGSH_gen, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_LiverGSH_gen, alpha = 0.1), peakdrop = 1.8, prob = 0.95, expand2markers = FALSE)  
  
#For Liver Adj GSSG
  #sex and gen covariates
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  threshold_LiverAdjGSSG_sexgen = summary(perm_LiverAdjGSSG_sexgen, alpha = c(0.1, 0.05, 0.01))
  plot_scan1(x = qtlscan_LiverAdjGSSG_sexgen, map = R01_GSH_DO_QTLdata$gmap,  main = "Genome Scan for Liver Adj GSSG - sex and gen", ylim = ylim)
  abline(h = threshold_LiverAdjGSSG_sexgen, col = c("red", "blue", "green"), lwd = 2)
  find_peaks(scan1_output = qtlscan_LiverAdjGSSG_sexgen, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_LiverAdjGSSG_sexgen, alpha = 0.1), peakdrop = 1.8, prob = 0.95, expand2markers = FALSE)
  peaks_LiverAdjGSSG_sexgen <- find_peaks(scan1_output = qtlscan_LiverAdjGSSG_sexgen, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_LiverAdjGSSG_sexgen, alpha = 0.1), peakdrop = 1.8, prob = 0.95, expand2markers = FALSE)

  #sex covariate
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  threshold_LiverAdjGSSG_sex = summary(perm_LiverAdjGSSG_sex, alpha = c(0.1, 0.05, 0.01))
  plot_scan1(x = qtlscan_LiverAdjGSSG_sex, map = R01_GSH_DO_QTLdata$gmap,  main = "Genome Scan for Liver Adj GSSG - sex", ylim = ylim)
  abline(h = threshold_LiverAdjGSSG_sex, col = c("red", "blue", "green"), lwd = 2)
  find_peaks(scan1_output = qtlscan_LiverAdjGSSG_sex, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_LiverAdjGSSG_sex, alpha = 0.1), peakdrop = 1.8, prob = 0.95, expand2markers = FALSE)
  peaks_LiverAdjGSSG_sex <- find_peaks(scan1_output = qtlscan_LiverAdjGSSG_sex, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_LiverAdjGSSG_sex, alpha = 0.1), peakdrop = 1.8, prob = 0.95, expand2markers = FALSE)
  
  #gen covariate
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  threshold_LiverAdjGSSG_gen = summary(perm_LiverAdjGSSG_gen, alpha = c(0.1, 0.05, 0.01))
  plot_scan1(x = qtlscan_LiverAdjGSSG_gen, map = R01_GSH_DO_QTLdata$gmap,  main = "Genome Scan for Liver Adj GSSG - gen", ylim = ylim)
  abline(h = threshold_LiverAdjGSSG_gen, col = c("red", "blue", "green"), lwd = 2)
  find_peaks(scan1_output = qtlscan_LiverAdjGSSG_gen, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_LiverAdjGSSG_gen, alpha = 0.1), peakdrop = 1.8, prob = 0.95, expand2markers = FALSE)
  peaks_LiverAdjGSSG_gen <- find_peaks(scan1_output = qtlscan_LiverAdjGSSG_gen, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_LiverAdjGSSG_gen, alpha = 0.1), peakdrop = 1.8, prob = 0.95, expand2markers = FALSE)

#For Liver Adj Total GSH
  #sex and gen covariates
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  threshold_LiverAdjTotalGSH_sexgen = summary(perm_LiverAdjTotalGSH_sexgen, alpha = c(0.1, 0.05, 0.01))
  plot_scan1(x = qtlscan_LiverAdjTotalGSH_sexgen, map = R01_GSH_DO_QTLdata$gmap,  main = "Genome Scan for Liver Adj Total GSH - sex and gen", ylim = ylim)
  abline(h = threshold_LiverAdjTotalGSH_sexgen, col = c("red", "blue", "green"), lwd = 2)
  find_peaks(scan1_output = qtlscan_LiverAdjTotalGSH_sexgen, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_LiverAdjTotalGSH_sexgen, alpha = 0.1), peakdrop = 1.8, prob = 0.95, expand2markers = FALSE)
  peaks_LiverAdjTotalGSH_sexgen <- find_peaks(scan1_output = qtlscan_LiverAdjTotalGSH_sexgen, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_LiverAdjTotalGSH_sexgen, alpha = 0.1), peakdrop = 1.8, prob = 0.95, expand2markers = FALSE)
  
  #sex covariate
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  threshold_LiverAdjTotalGSH_sex = summary(perm_LiverAdjTotalGSH_sex, alpha = c(0.1, 0.05, 0.01))
  plot_scan1(x = qtlscan_LiverAdjTotalGSH_sex, map = R01_GSH_DO_QTLdata$gmap,  main = "Genome Scan for Liver Adj Total GSH - sex", ylim = ylim)
  abline(h = threshold_LiverAdjTotalGSH_sex, col = c("red", "blue", "green"), lwd = 2)
  find_peaks(scan1_output = qtlscan_LiverAdjTotalGSH_sex, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_LiverAdjTotalGSH_sex, alpha = 0.1), peakdrop = 1.8, prob = 0.95, expand2markers = FALSE)
  peaks_LiverAdjTotalGSH_sex <- find_peaks(scan1_output = qtlscan_LiverAdjTotalGSH_sex, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_LiverAdjTotalGSH_sex, alpha = 0.1), peakdrop = 1.8, prob = 0.95, expand2markers = FALSE)
  
  #gen covariate
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  threshold_LiverAdjTotalGSH_gen = summary(perm_LiverAdjTotalGSH_gen, alpha = c(0.1, 0.05, 0.01))
  plot_scan1(x = qtlscan_LiverAdjTotalGSH_gen, map = R01_GSH_DO_QTLdata$gmap,  main = "Genome Scan for Liver Adj Total GSH - gen", ylim = ylim)
  abline(h = threshold_LiverAdjTotalGSH_gen, col = c("red", "blue", "green"), lwd = 2)
  find_peaks(scan1_output = qtlscan_LiverAdjTotalGSH_gen, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_LiverAdjTotalGSH_gen, alpha = 0.1), peakdrop = 1.8, prob = 0.95, expand2markers = FALSE)
  peaks_LiverAdjTotalGSH_gen <- find_peaks(scan1_output = qtlscan_LiverAdjTotalGSH_gen, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_LiverAdjTotalGSH_gen, alpha = 0.1), peakdrop = 1.8, prob = 0.95, expand2markers = FALSE)
  
#For Liver Adj GSH/GSSG Ratio
  #sex and gen covariates
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  threshold_LiverAdjGSH_GSSG_sexgen = summary(perm_LiverAdjGSH_GSSG_sexgen, alpha = c(0.1, 0.05, 0.01))
  plot_scan1(x = qtlscan_LiverAdjGSH_GSSG_sexgen, map = R01_GSH_DO_QTLdata$gmap,  main = "Genome Scan for Liver Adj GSH/GSSG Ratio - sex and gen", ylim = ylim)
  abline(h = threshold_LiverAdjGSH_GSSG_sexgen, col = c("red", "blue", "green"), lwd = 2)
  find_peaks(scan1_output = qtlscan_LiverAdjGSH_GSSG_sexgen, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_LiverAdjGSH_GSSG_sexgen, alpha = 0.1), peakdrop = 1.8, prob = 0.95, expand2markers = FALSE)
  peaks_LiverAdjGSH_GSSG_sexgen <- find_peaks(scan1_output = qtlscan_LiverAdjGSH_GSSG_sexgen, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_LiverAdjGSH_GSSG_sexgen, alpha = 0.1), peakdrop = 1.8, prob = 0.95, expand2markers = FALSE)
  
  #sex covariate
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  threshold_LiverAdjGSH_GSSG_sex = summary(perm_LiverAdjGSH_GSSG_sex, alpha = c(0.1, 0.05, 0.01))
  plot_scan1(x = qtlscan_LiverAdjGSH_GSSG_sex, map = R01_GSH_DO_QTLdata$gmap,  main = "Genome Scan for Liver Adj GSH/GSSG - sex", ylim = ylim)
  abline(h = threshold_LiverAdjGSH_GSSG_sex, col = c("red", "blue", "green"), lwd = 2)
  find_peaks(scan1_output = qtlscan_LiverAdjGSH_GSSG_sex, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_LiverAdjGSH_GSSG_sex, alpha = 0.1), peakdrop = 1.8, prob = 0.95, expand2markers = FALSE)
  peaks_LiverAdjGSH_GSSG_sex <- find_peaks(scan1_output = qtlscan_LiverAdjGSH_GSSG_sex, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_LiverAdjGSH_GSSG_sex, alpha = 0.1), peakdrop = 1.8, prob = 0.95, expand2markers = FALSE)
  
  #gen covariate
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  threshold_LiverAdjGSH_GSSG_gen = summary(perm_LiverAdjGSH_GSSG_gen, alpha = c(0.1, 0.05, 0.01))
  plot_scan1(x = qtlscan_LiverAdjGSH_GSSG_gen, map = R01_GSH_DO_QTLdata$gmap,  main = "Genome Scan for Liver Adj GSH/GSSG Ratio - gen", ylim = ylim)
  abline(h = threshold_LiverAdjGSH_GSSG_gen, col = c("red", "blue", "green"), lwd = 2)
  find_peaks(scan1_output = qtlscan_LiverAdjGSH_GSSG_gen, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_LiverAdjGSH_GSSG_gen, alpha = 0.1), peakdrop = 1.8, prob = 0.95, expand2markers = FALSE)
  peaks_LiverAdjGSH_GSSG_gen <- find_peaks(scan1_output = qtlscan_LiverAdjGSH_GSSG_gen, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_LiverAdjGSH_GSSG_gen, alpha = 0.1), peakdrop = 1.8, prob = 0.95, expand2markers = FALSE)  

#For Liver NADH
  #sex and gen covariates
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  threshold_LiverNADH_sexgen = summary(perm_LiverNADH_sexgen, alpha = c(0.1, 0.05, 0.01))
  plot_scan1(x = qtlscan_LiverNADH_sexgen, map = R01_GSH_DO_QTLdata$gmap,  main = "Genome Scan for Liver NADH - sex and gen (100 perm)", ylim = ylim)
  abline(h = threshold_LiverNADH_sexgen, col = c("red", "blue", "green"), lwd = 2)
  find_peaks(scan1_output = qtlscan_LiverNADH_sexgen, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_LiverNADH_sexgen, alpha = 0.1), peakdrop = 1.8, prob = 0.95, expand2markers = FALSE)
  peaks_LiverNADH_sexgen <- find_peaks(scan1_output = qtlscan_LiverNADH_sexgen, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_LiverNADH_sexgen, alpha = 0.1), peakdrop = 1.8, prob = 0.95, expand2markers = FALSE)
  
  #sex covariate
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  threshold_LiverNADH_sex = summary(perm_LiverNADH_sex, alpha = c(0.1, 0.05, 0.01))
  plot_scan1(x = qtlscan_LiverNADH_sex, map = R01_GSH_DO_QTLdata$gmap,  main = "Genome Scan for Liver NADH - sex", ylim = ylim)
  abline(h = threshold_LiverNADH_sex, col = c("red", "blue", "green"), lwd = 2)
  find_peaks(scan1_output = qtlscan_LiverNADH_sex, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_LiverNADH_sex, alpha = 0.1), peakdrop = 1.8, prob = 0.95, expand2markers = FALSE)
  peaks_LiverNADH_sex <- find_peaks(scan1_output = qtlscan_LiverNADH_sex, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_LiverNADH_sex, alpha = 0.1), peakdrop = 1.8, prob = 0.95, expand2markers = FALSE)
  
  #gen covariate
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  threshold_LiverNADH_gen = summary(perm_LiverNADH_gen, alpha = c(0.1, 0.05, 0.01))
  plot_scan1(x = qtlscan_LiverNADH_gen, map = R01_GSH_DO_QTLdata$gmap,  main = "Genome Scan for Liver NADH - gen (100 perm)", ylim = ylim)
  abline(h = threshold_LiverNADH_gen, col = c("red", "blue", "green"), lwd = 2)
  find_peaks(scan1_output = qtlscan_LiverNADH_gen, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_LiverNADH_gen, alpha = 0.1), peakdrop = 1.8, prob = 0.95, expand2markers = FALSE)
  peaks_LiverNADH_gen <- find_peaks(scan1_output = qtlscan_LiverNADH_gen, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_LiverNADH_gen, alpha = 0.1), peakdrop = 1.8, prob = 0.95, expand2markers = FALSE)    

#For Liver NADP
  #sex and gen covariates
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  threshold_LiverNADP_sexgen = summary(perm_LiverNADP_sexgen, alpha = c(0.1, 0.05, 0.01))
  plot_scan1(x = qtlscan_LiverNADP_sexgen, map = R01_GSH_DO_QTLdata$gmap,  main = "Genome Scan for Liver NADP - sex and gen", ylim = ylim)
  abline(h = threshold_LiverNADP_sexgen, col = c("red", "blue", "green"), lwd = 2)
  find_peaks(scan1_output = qtlscan_LiverNADP_sexgen, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_LiverNADP_sexgen, alpha = 0.1), peakdrop = 1.8, prob = 0.95, expand2markers = FALSE)
  peaks_LiverNADP_sexgen <- find_peaks(scan1_output = qtlscan_LiverNADP_sexgen, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_LiverNADP_sexgen, alpha = 0.1), peakdrop = 1.8, prob = 0.95, expand2markers = FALSE)
  
  #sex covariate
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  threshold_LiverNADP_sex = summary(perm_LiverNADP_sex, alpha = c(0.1, 0.05, 0.01))
  plot_scan1(x = qtlscan_LiverNADP_sex, map = R01_GSH_DO_QTLdata$gmap,  main = "Genome Scan for Liver NADP - sex", ylim = ylim)
  abline(h = threshold_LiverNADP_sex, col = c("red", "blue", "green"), lwd = 2)
  find_peaks(scan1_output = qtlscan_LiverNADP_sex, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_LiverNADP_sex, alpha = 0.1), peakdrop = 1.8, prob = 0.95, expand2markers = FALSE)
  peaks_LiverNADP_sex <- find_peaks(scan1_output = qtlscan_LiverNADP_sex, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_LiverNADP_sex, alpha = 0.1), peakdrop = 1.8, prob = 0.95, expand2markers = FALSE)
  
  #gen covariate
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  threshold_LiverNADP_gen = summary(perm_LiverNADP_gen, alpha = c(0.1, 0.05, 0.01))
  plot_scan1(x = qtlscan_LiverNADP_gen, map = R01_GSH_DO_QTLdata$gmap,  main = "Genome Scan for Liver NADP - gen", ylim = ylim)
  abline(h = threshold_LiverNADP_gen, col = c("red", "blue", "green"), lwd = 2)
  find_peaks(scan1_output = qtlscan_LiverNADP_gen, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_LiverNADP_gen, alpha = 0.1), peakdrop = 1.8, prob = 0.95, expand2markers = FALSE)
  peaks_LiverNADP_gen <- find_peaks(scan1_output = qtlscan_LiverNADP_gen, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_LiverNADP_gen, alpha = 0.1), peakdrop = 1.8, prob = 0.95, expand2markers = FALSE)      

#For Liver NADPH
  #sex and gen covariates
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  threshold_LiverNADPH_sexgen = summary(perm_LiverNADPH_sexgen, alpha = c(0.1, 0.05, 0.01))
  plot_scan1(x = qtlscan_LiverNADPH_sexgen, map = R01_GSH_DO_QTLdata$gmap,  main = "Genome Scan for Liver NADPH - sex and gen", ylim = ylim)
  abline(h = threshold_LiverNADPH_sexgen, col = c("red", "blue", "green"), lwd = 2)
  find_peaks(scan1_output = qtlscan_LiverNADPH_sexgen, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_LiverNADPH_sexgen, alpha = 0.1), peakdrop = 1.8, prob = 0.95, expand2markers = FALSE)
  peaks_LiverNADPH_sexgen <- find_peaks(scan1_output = qtlscan_LiverNADPH_sexgen, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_LiverNADPH_sexgen, alpha = 0.1), peakdrop = 1.8, prob = 0.95, expand2markers = FALSE)
  
  #sex covariate
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  threshold_LiverNADPH_sex = summary(perm_LiverNADPH_sex, alpha = c(0.1, 0.05, 0.01))
  plot_scan1(x = qtlscan_LiverNADPH_sex, map = R01_GSH_DO_QTLdata$gmap,  main = "Genome Scan for Liver NADPH - sex", ylim = ylim)
  abline(h = threshold_LiverNADPH_sex, col = c("red", "blue", "green"), lwd = 2)
  find_peaks(scan1_output = qtlscan_LiverNADPH_sex, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_LiverNADPH_sex, alpha = 0.1), peakdrop = 1.8, prob = 0.95, expand2markers = FALSE)
  peaks_LiverNADPH_sex <- find_peaks(scan1_output = qtlscan_LiverNADPH_sex, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_LiverNADPH_sex, alpha = 0.1), peakdrop = 1.8, prob = 0.95, expand2markers = FALSE)
  
  #gen covariate
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  threshold_LiverNADPH_gen = summary(perm_LiverNADPH_gen, alpha = c(0.1, 0.05, 0.01))
  plot_scan1(x = qtlscan_LiverNADPH_gen, map = R01_GSH_DO_QTLdata$gmap,  main = "Genome Scan for Liver NADPH - gen", ylim = ylim)
  abline(h = threshold_LiverNADPH_gen, col = c("red", "blue", "green"), lwd = 2)
  find_peaks(scan1_output = qtlscan_LiverNADPH_gen, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_LiverNADPH_gen, alpha = 0.1), peakdrop = 1.8, prob = 0.95, expand2markers = FALSE)
  peaks_LiverNADPH_gen <- find_peaks(scan1_output = qtlscan_LiverNADPH_gen, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_LiverNADPH_gen, alpha = 0.1), peakdrop = 1.8, prob = 0.95, expand2markers = FALSE)        
  
#For Liver NADP/NADPH
  #sex and gen covariates
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  threshold_LiverNADP_NADPH_sexgen = summary(perm_LiverNADP_NADPH_sexgen, alpha = c(0.1, 0.05, 0.01))
  plot_scan1(x = qtlscan_LiverNADP_NADPH_sexgen, map = R01_GSH_DO_QTLdata$gmap,  main = "Genome Scan for Liver NADP/NADPH Ratio - sex and gen", ylim = ylim)
  abline(h = threshold_LiverNADP_NADPH_sexgen, col = c("red", "blue", "green"), lwd = 2)
  find_peaks(scan1_output = qtlscan_LiverNADP_NADPH_sexgen, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_LiverNADP_NADPH_sexgen, alpha = 0.1), peakdrop = 1.8, prob = 0.95, expand2markers = FALSE)
  peaks_LiverNADP_NADPH_sexgen <- find_peaks(scan1_output = qtlscan_LiverNADP_NADPH_sexgen, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_LiverNADP_NADPH_sexgen, alpha = 0.1), peakdrop = 1.8, prob = 0.95, expand2markers = FALSE)
  
  #sex covariate
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  threshold_LiverNADP_NADPH_sex = summary(perm_LiverNADP_NADPH_sex, alpha = c(0.1, 0.05, 0.01))
  plot_scan1(x = qtlscan_LiverNADP_NADPH_sex, map = R01_GSH_DO_QTLdata$gmap,  main = "Genome Scan for Liver NADP/NADPH Ratio - sex", ylim = ylim)
  abline(h = threshold_LiverNADP_NADPH_sex, col = c("red", "blue", "green"), lwd = 2)
  find_peaks(scan1_output = qtlscan_LiverNADP_NADPH_sex, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_LiverNADP_NADPH_sex, alpha = 0.1), peakdrop = 1.8, prob = 0.95, expand2markers = FALSE)
  peaks_LiverNADP_NADPH_sex <- find_peaks(scan1_output = qtlscan_LiverNADP_NADPH_sex, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_LiverNADP_NADPH_sex, alpha = 0.1), peakdrop = 1.8, prob = 0.95, expand2markers = FALSE)
  
  #gen covariate
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  threshold_LiverNADP_NADPH_gen = summary(perm_LiverNADP_NADPH_gen, alpha = c(0.1, 0.05, 0.01))
  plot_scan1(x = qtlscan_LiverNADP_NADPH_gen, map = R01_GSH_DO_QTLdata$gmap,  main = "Genome Scan for Liver NADP/NADPH Ratio - gen", ylim = ylim)
  abline(h = threshold_LiverNADP_NADPH_gen, col = c("red", "blue", "green"), lwd = 2)
  find_peaks(scan1_output = qtlscan_LiverNADP_NADPH_gen, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_LiverNADP_NADPH_gen, alpha = 0.1), peakdrop = 1.8, prob = 0.95, expand2markers = FALSE)
  peaks_LiverNADP_NADPH_gen <- find_peaks(scan1_output = qtlscan_LiverNADP_NADPH_gen, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_LiverNADP_NADPH_gen, alpha = 0.1), peakdrop = 1.8, prob = 0.95, expand2markers = FALSE)          
  
#For Liver Unadj GSSG
  #sex and gen covariates
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  threshold_LiverUnadjGSSG_sexgen = summary(perm_LiverUnadjGSSG_sexgen, alpha = c(0.1, 0.05, 0.01))
  plot_scan1(x = qtlscan_LiverUnadjGSSG_sexgen, map = R01_GSH_DO_QTLdata$gmap,  main = "Genome Scan for Liver Unadj GSSG - sex and gen", ylim = ylim)
  abline(h = threshold_LiverUnadjGSSG_sexgen, col = c("red", "blue", "green"), lwd = 2)
  find_peaks(scan1_output = qtlscan_LiverUnadjGSSG_sexgen, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_LiverUnadjGSSG_sexgen, alpha = 0.1), peakdrop = 1.8, prob = 0.95, expand2markers = FALSE)
  peaks_LiverUnadjGSSG_sexgen <- find_peaks(scan1_output = qtlscan_LiverUnadjGSSG_sexgen, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_LiverUnadjGSSG_sexgen, alpha = 0.1), peakdrop = 1.8, prob = 0.95, expand2markers = FALSE)
  
  #sex covariate
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  threshold_LiverUnadjGSSG_sex = summary(perm_LiverUnadjGSSG_sex, alpha = c(0.1, 0.05, 0.01))
  plot_scan1(x = qtlscan_LiverUnadjGSSG_sex, map = R01_GSH_DO_QTLdata$gmap,  main = "Genome Scan for Liver Unadj GSSG - sex", ylim = ylim)
  abline(h = threshold_LiverUnadjGSSG_sex, col = c("red", "blue", "green"), lwd = 2)
  find_peaks(scan1_output = qtlscan_LiverUnadjGSSG_sex, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_LiverUnadjGSSG_sex, alpha = 0.1), peakdrop = 1.8, prob = 0.95, expand2markers = FALSE)
  peaks_LiverUnadjGSSG_sex <- find_peaks(scan1_output = qtlscan_LiverUnadjGSSG_sex, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_LiverUnadjGSSG_sex, alpha = 0.1), peakdrop = 1.8, prob = 0.95, expand2markers = FALSE)
  
  #gen covariate
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  threshold_LiverUnadjGSSG_gen = summary(perm_LiverUnadjGSSG_gen, alpha = c(0.1, 0.05, 0.01))
  plot_scan1(x = qtlscan_LiverUnadjGSSG_gen, map = R01_GSH_DO_QTLdata$gmap,  main = "Genome Scan for Liver Unadj GSSG - gen (100 perm)", ylim = ylim)
  abline(h = threshold_LiverUnadjGSSG_gen, col = c("red", "blue", "green"), lwd = 2)
  find_peaks(scan1_output = qtlscan_LiverUnadjGSSG_gen, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_LiverUnadjGSSG_gen, alpha = 0.1), peakdrop = 1.8, prob = 0.95, expand2markers = FALSE)
  peaks_LiverUnadjGSSG_gen <- find_peaks(scan1_output = qtlscan_LiverUnadjGSSG_gen, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_LiverUnadjGSSG_gen, alpha = 0.1), peakdrop = 1.8, prob = 0.95, expand2markers = FALSE)
  
#For Liver Unadj Total GSH
  #sex and gen covariates
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  threshold_LiverUnadjTotalGSH_sexgen = summary(perm_LiverUnadjTotalGSH_sexgen, alpha = c(0.1, 0.05, 0.01))
  plot_scan1(x = qtlscan_LiverUnadjTotalGSH_sexgen, map = R01_GSH_DO_QTLdata$gmap,  main = "Genome Scan for Liver Unadj Total GSH - sex and gen (100 perm)", ylim = ylim)
  abline(h = threshold_LiverUnadjTotalGSH_sexgen, col = c("red", "blue", "green"), lwd = 2)
  find_peaks(scan1_output = qtlscan_LiverUnadjTotalGSH_sexgen, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_LiverUnadjTotalGSH_sexgen, alpha = 0.1), peakdrop = 1.8, prob = 0.95, expand2markers = FALSE)
  peaks_LiverUnadjTotalGSH_sexgen <- find_peaks(scan1_output = qtlscan_LiverUnadjTotalGSH_sexgen, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_LiverUnadjTotalGSH_sexgen, alpha = 0.1), peakdrop = 1.8, prob = 0.95, expand2markers = FALSE)
  
  #sex covariate
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  threshold_LiverUnadjTotalGSH_sex = summary(perm_LiverUnadjTotalGSH_sex, alpha = c(0.1, 0.05, 0.01))
  plot_scan1(x = qtlscan_LiverUnadjTotalGSH_sex, map = R01_GSH_DO_QTLdata$gmap,  main = "Genome Scan for Liver Unadj Total GSH - sex (100 perm)", ylim = ylim)
  abline(h = threshold_LiverUnadjTotalGSH_sex, col = c("red", "blue", "green"), lwd = 2)
  find_peaks(scan1_output = qtlscan_LiverUnadjTotalGSH_sex, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_LiverUnadjTotalGSH_sex, alpha = 0.1), peakdrop = 1.8, prob = 0.95, expand2markers = FALSE)
  peaks_LiverUnadjTotalGSH_sex <- find_peaks(scan1_output = qtlscan_LiverUnadjTotalGSH_sex, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_LiverUnadjTotalGSH_sex, alpha = 0.1), peakdrop = 1.8, prob = 0.95, expand2markers = FALSE)
  
  #gen covariate
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  threshold_LiverUnadjTotalGSH_gen = summary(perm_LiverUnadjTotalGSH_gen, alpha = c(0.1, 0.05, 0.01))
  plot_scan1(x = qtlscan_LiverUnadjTotalGSH_gen, map = R01_GSH_DO_QTLdata$gmap,  main = "Genome Scan for Liver Unadj Total GSH - gen (100 perm)", ylim = ylim)
  abline(h = threshold_LiverUnadjTotalGSH_gen, col = c("red", "blue", "green"), lwd = 2)
  find_peaks(scan1_output = qtlscan_LiverUnadjTotalGSH_gen, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_LiverUnadjTotalGSH_gen, alpha = 0.1), peakdrop = 1.8, prob = 0.95, expand2markers = FALSE)
  peaks_LiverUnadjTotalGSH_gen <- find_peaks(scan1_output = qtlscan_LiverUnadjTotalGSH_gen, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_LiverUnadjTotalGSH_gen, alpha = 0.1), peakdrop = 1.8, prob = 0.95, expand2markers = FALSE)
  
#For Liver Unadj GSH/GSSG Ratio
  #sex and gen covariates
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  threshold_LiverUnadjGSH_GSSG_sexgen = summary(perm_LiverUnadjGSH_GSSG_sexgen, alpha = c(0.1, 0.05, 0.01))
  plot_scan1(x = qtlscan_LiverUnadjGSH_GSSG_sexgen, map = R01_GSH_DO_QTLdata$gmap,  main = "Genome Scan for Liver Unadj GSH/GSSG Ratio - sex and gen", ylim = ylim)
  abline(h = threshold_LiverUnadjGSH_GSSG_sexgen, col = c("red", "blue", "green"), lwd = 2)
  find_peaks(scan1_output = qtlscan_LiverUnadjGSH_GSSG_sexgen, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_LiverUnadjGSH_GSSG_sexgen, alpha = 0.1), peakdrop = 1.8, prob = 0.95, expand2markers = FALSE)
  peaks_LiverUnadjGSH_GSSG_sexgen <- find_peaks(scan1_output = qtlscan_LiverUnadjGSH_GSSG_sexgen, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_LiverUnadjGSH_GSSG_sexgen, alpha = 0.1), peakdrop = 1.8, prob = 0.95, expand2markers = FALSE)
  
  #sex covariate
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  threshold_LiverUnadjGSH_GSSG_sex = summary(perm_LiverUnadjGSH_GSSG_sex, alpha = c(0.1, 0.05, 0.01))
  plot_scan1(x = qtlscan_LiverUnadjGSH_GSSG_sex, map = R01_GSH_DO_QTLdata$gmap,  main = "Genome Scan for Liver Unadj GSH/GSSG - sex", ylim = ylim)
  abline(h = threshold_LiverUnadjGSH_GSSG_sex, col = c("red", "blue", "green"), lwd = 2)
  find_peaks(scan1_output = qtlscan_LiverUnadjGSH_GSSG_sex, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_LiverUnadjGSH_GSSG_sex, alpha = 0.1), peakdrop = 1.8, prob = 0.95, expand2markers = FALSE)
  peaks_LiverUnadjGSH_GSSG_sex <- find_peaks(scan1_output = qtlscan_LiverUnadjGSH_GSSG_sex, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_LiverUnadjGSH_GSSG_sex, alpha = 0.1), peakdrop = 1.8, prob = 0.95, expand2markers = FALSE)
  
  #gen covariate
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  threshold_LiverUnadjGSH_GSSG_gen = summary(perm_LiverUnadjGSH_GSSG_gen, alpha = c(0.1, 0.05, 0.01))
  plot_scan1(x = qtlscan_LiverUnadjGSH_GSSG_gen, map = R01_GSH_DO_QTLdata$gmap,  main = "Genome Scan for Liver Unadj GSH/GSSG Ratio - gen", ylim = ylim)
  abline(h = threshold_LiverUnadjGSH_GSSG_gen, col = c("red", "blue", "green"), lwd = 2)
  find_peaks(scan1_output = qtlscan_LiverUnadjGSH_GSSG_gen, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_LiverUnadjGSH_GSSG_gen, alpha = 0.1), peakdrop = 1.8, prob = 0.95, expand2markers = FALSE)
  peaks_LiverUnadjGSH_GSSG_gen <- find_peaks(scan1_output = qtlscan_LiverUnadjGSH_GSSG_gen, map = R01_GSH_DO_QTLdata$gmap, threshold = summary(perm_LiverUnadjGSH_GSSG_gen, alpha = 0.1), peakdrop = 1.8, prob = 0.95, expand2markers = FALSE)   
  
####################################################
## Identify significant LOD peaks by estimating QTL effects (coefficients)
## to find peaks - must specify a specific threshold level
## plot one chromosome at a time - scan1blup function best for DO, scan1coef used as well
## set xlim for interval that contains the significant peak once you've plotted the entire chromosome
####################################################

#For Liver GSH -- sex covariate ---- Chromosome 2
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  chr = 2
  coef_blup_LiverGSH_chr2 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["logLiverGSH"], kinship = kinship_loco[[chr]], addcovar = sex, cores = 10)
  plot_coefCC(x = coef_blup_LiverGSH_chr2, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSH_sex, main = "Genome Scan for Liver GSH - sex", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(40,60)
  plot_coefCC(x = coef_blup_LiverGSH_chr2, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverGSH_sex, main = "Genome Scan for Liver GSH - sex", legend = "bottomleft", bgcolor="gray95", xlim = xlim)

#For Liver Adj GSSG -- sex covariate ---- Chromosome 2
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  chr = 2
  coef_blup_LiverAdjGSSG_chr2 <- scan1blup(genoprobs = probs[,chr], pheno = pheno["logLiverAdjGSSG"], kinship = kinship_loco[[chr]], addcovar = sex, cores = 10)
  plot_coefCC(x = coef_blup_LiverAdjGSSG_chr2, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverAdjGSSG_sex, main = "Genome Scan for Liver Adj GSSG - sex", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(40,60)
  plot_coefCC(x = coef_blup_LiverAdjGSSG_chr2, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverAdjGSSG_sex, main = "Genome Scan for Liver Adj GSSG - sex", legend = "bottomleft", bgcolor="gray95", xlim = xlim)

#For Liver Adj GSSG -- sex covariate ---- Chromosome 5
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  chr = 5
  coef_blup_LiverAdjGSSG_chr5 <- scan1blup(genoprobs = probs[,chr], pheno = pheno["logLiverAdjGSSG"], kinship = kinship_loco[[chr]], addcovar = sex, cores = 10)
  plot_coefCC(x = coef_blup_LiverAdjGSSG_chr5, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverAdjGSSG_sex, main = "Genome Scan for Liver Adj GSSG - sex", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(40, 70)
  plot_coefCC(x = coef_blup_LiverAdjGSSG_chr5, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverAdjGSSG_sex, main = "Genome Scan for Liver Adj GSSG - sex", legend = "bottomleft", bgcolor="gray95", xlim = xlim)

#For Liver NADP -- sex covariate ---- Chromosome 2
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  chr = 2
  coef_blup_LiverNADP_chr2 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["logLiverNADP"], kinship = kinship_loco[[chr]], addcovar = sex, cores = 10)
  plot_coefCC(x = coef_blup_LiverNADP_chr2, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverNADP_sex, main = "Genome Scan for Liver NADP - sex", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(40,60)
  plot_coefCC(x = coef_blup_LiverNADP_chr2, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverNADP_sex, main = "Genome Scan for Liver NADP - sex", legend = "bottomleft", bgcolor="gray95", xlim = xlim)

#For Liver NADPH -- sex covariate ---- Chromosome 14
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  chr = 14
  coef_blup_LiverNADPH_chr14 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["logLiverNADPH"], kinship = kinship_loco[[chr]], addcovar = sex, cores = 10)
  plot_coefCC(x = coef_blup_LiverNADPH_chr14, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverNADPH_sex, main = "Genome Scan for Liver NADPH - sex", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(25,55)
  plot_coefCC(x = coef_blup_LiverNADPH_chr14, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverNADPH_sex, main = "Genome Scan for Liver NADPH - sex", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
  
#For Liver NADP/NADPH Ratio -- sex covariate ---- Chromosome 2
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  chr = 2
  coef_blup_LiverNADP_NADPH_chr2 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["logLiverNADP_NADPHRatio"], kinship = kinship_loco[[chr]], addcovar = sex, cores = 10)
  plot_coefCC(x = coef_blup_LiverNADP_NADPH_chr2, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverNADP_NADPH_sex, main = "Genome Scan for Liver NADP/NADPH Ratio - sex", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(40,70)
  plot_coefCC(x = coef_blup_LiverNADP_NADPH_chr2, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverNADP_NADPH_sex, main = "Genome Scan for Liver NADP/NADPH Ratio - sex", legend = "bottomleft", bgcolor="gray95", xlim = xlim)

#For Liver NADP/NADPH Ratio -- sex covariate ---- Chromosome 14
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  chr = 14
  coef_blup_LiverNADP_NADPH_chr14 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["logLiverNADP_NADPHRatio"], kinship = kinship_loco[[chr]], addcovar = sex, cores = 10)
  plot_coefCC(x = coef_blup_LiverNADP_NADPH_chr14, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverNADP_NADPH_sex, main = "Genome Scan for Liver NADP/NADPH Ratio - sex", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(10,45)
  plot_coefCC(x = coef_blup_LiverNADP_NADPH_chr14, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverNADP_NADPH_sex, main = "Genome Scan for Liver NADP/NADPH Ratio - sex", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
  
#For Liver NADP/NADPH Ratio -- sexgen covariate ---- Chromosome 7
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  chr = 7
  coef_blup_LiverNADP_NADPH_chr7_sexgen <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["logLiverNADP_NADPHRatio"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 10)
  plot_coefCC(x = coef_blup_LiverNADP_NADPH_chr7_sexgen, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverNADP_NADPH_sexgen, main = "Genome Scan for Liver NADP/NADPH Ratio - sexgen", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(0,90)
  plot_coefCC(x = coef_blup_LiverNADP_NADPH_chr7_sexgen, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverNADP_NADPH_sexgen, main = "Genome Scan for Liver NADP/NADPH Ratio - sexgen", legend = "bottomleft", bgcolor="gray95", xlim = xlim)  
 
#For Liver NADP/NADPH Ratio -- gen covariate ---- Chromosome 7
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  chr = 7
  coef_blup_LiverNADP_NADPH_chr7_gen <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["logLiverNADP_NADPHRatio"], kinship = kinship_loco[[chr]], addcovar = gen, cores = 10)
  plot_coefCC(x = coef_blup_LiverNADP_NADPH_chr7_gen, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverNADP_NADPH_gen, main = "Genome Scan for Liver NADP/NADPH Ratio - gen", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(0,90)
  plot_coefCC(x = coef_blup_LiverNADP_NADPH_chr7_gen, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverNADP_NADPH_gen, main = "Genome Scan for Liver NADP/NADPH Ratio - gen", legend = "bottomleft", bgcolor="gray95", xlim = xlim)    
  
#For Liver Adj Total GSH -- sex covariate ---- Chromosome 2
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  chr = 2
  coef_blup_LiverAdjTotalGSH_chr2 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["logLiverAdjTotalGSH"], kinship = kinship_loco[[chr]], addcovar = sex, cores = 10)
  plot_coefCC(x = coef_blup_LiverAdjTotalGSH_chr2, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverAdjTotalGSH_sex, main = "Genome Scan for Liver Adj Total GSH - sex", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(40,60)
  plot_coefCC(x = coef_blup_LiverAdjTotalGSH_chr2, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverAdjTotalGSH_sex, main = "Genome Scan for Liver Adj Total GSH - sex", legend = "bottomleft", bgcolor="gray95", xlim = xlim)

#For Liver Unadj GSSG -- sex covariate ---- Chromosome 2
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  chr = 2
  coef_blup_LiverUnadjGSSG_chr2 <- scan1blup(genoprobs = probs[,chr], pheno = pheno["logLiverAdjGSSG"], kinship = kinship_loco[[chr]], addcovar = sex, cores = 10)
  plot_coefCC(x = coef_blup_LiverUnadjGSSG_chr2, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverUnadjGSSG_sex, main = "Genome Scan for Liver Unadj GSSG - sex", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(40, 60)
  plot_coefCC(x = coef_blup_LiverUnadjGSSG_chr2, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverUnadjGSSG_sex, main = "Genome Scan for Liver Unadj GSSG - sex", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
  
#For Liver Unadj Total GSH -- sex covariate ---- Chromosome 2
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  chr = 2
  coef_blup_LiverUnadjTotalGSH_chr2 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["logLiverUnadjTotalGSH"], kinship = kinship_loco[[chr]], addcovar = sex, cores = 10)
  plot_coefCC(x = coef_blup_LiverUnadjTotalGSH_chr2, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverUnadjTotalGSH_sex, main = "Genome Scan for Liver Unadj Total GSH - sex", legend = "bottomleft", bgcolor="gray95")
  xlim <- c(40,60)
  plot_coefCC(x = coef_blup_LiverUnadjTotalGSH_chr2, map = R01_GSH_DO_QTLdata$gmap, scan1_output = qtlscan_LiverUnadjTotalGSH_sex, main = "Genome Scan for Liver Unadj Total GSH - sex", legend = "bottomleft", bgcolor="gray95", xlim = xlim)
  
    
####################################################
## Connecting to SNP and Gene Databases
## SNP Associations
## To get a table of the SNPs with the largest LOD scores, use the function top_snps(). This will show all SNPs with LOD score within some amount (the default is 1.5) of the maximum SNP LOD score. 
####################################################

#For Liver GSH -- sex covariate -- Chromosome 2
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  
  peak_Mbp <- max(scan1_output = qtlscan_LiverGSH_sex, map = R01_GSH_DO_QTLdata$pmap)$pos
  variants_LiverGSH_chr2 <- query_variants(2, peak_Mbp - 1, peak_Mbp + 1)
  out_snps_LiverGSH_chr2 <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["logLiverGSH"], kinship = kinship_loco[["2"]], addcovar = sex, query_func = query_variants,
                          chr = 2, start=peak_Mbp-1, end=peak_Mbp+1, keep_all_snps = TRUE)
  plot_snpasso(out_snps_LiverGSH_chr2$lod, out_snps_LiverGSH_chr2$snpinfo, main = "Liver GSH SNPs - sex ")
  top <- top_snps(out_snps_LiverGSH_chr2$lod, out_snps_LiverGSH_chr2$snpinfo)
  print(top[,c(1, 8:15, 20)], row.names=FALSE)
  
  LiverGSH_Chr2_Genes_MGI <- query_genes_mgi(2, peak_Mbp-1, peak_Mbp+1)
  plot(out_snps_LiverGSH_chr2$lod, out_snps_LiverGSH_chr2$snpinfo, drop_hilit=1.5, genes = LiverGSH_Chr2_Genes_MGI, main = "Liver GSH Genes MGI - sex")
  
  LiverGSH_Chr2_Genes <- query_genes(2, peak_Mbp-1, peak_Mbp+1)
  plot(out_snps_LiverGSH_chr2$lod, out_snps_LiverGSH_chr2$snpinfo, drop_hilit=1.5, genes = LiverGSH_Chr2_Genes, main = "Liver GSH Genes - sex")
  
#For Liver Adj GSSG -- sex covar -- chromosome 2
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  
  peak_Mbp <- max(scan1_output = qtlscan_LiverAdjGSSG_sex, map = R01_GSH_DO_QTLdata$pmap)$pos
  variants_LiverAdjGSSG_chr2 <- query_variants(2, peak_Mbp - 1, peak_Mbp + 1)
  out_snps_LiverAdjGSSG_chr2 <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["logLiverAdjGSSG"], kinship = kinship_loco[["2"]], addcovar = sex, query_func = query_variants,
                                      chr = 2, start=peak_Mbp-1, end=peak_Mbp+1, keep_all_snps = TRUE)
  plot_snpasso(out_snps_LiverAdjGSSG_chr2$lod, out_snps_LiverAdjGSSG_chr2$snpinfo, main = "Liver Adj GSSG SNPs - sex")
  top <- top_snps(out_snps_LiverAdjGSSG_chr2$lod, out_snps_LiverAdjGSSG_chr2$snpinfo)
  print(top[,c(1, 8:15, 20)], row.names=FALSE)
  
  LiverAdjGSSG_Chr2_Genes_MGI <- query_genes_mgi(2, peak_Mbp-1, peak_Mbp+1)
  plot(out_snps_LiverAdjGSSG_chr2$lod, out_snps_LiverAdjGSSG_chr2$snpinfo, drop_hilit=1.5, genes = LiverAdjGSSG_Chr2_Genes_MGI, main = "Liver Adj GSSG Genes MGI - sex")
  
  LiverAdjGSSG_Chr2_Genes <- query_genes(2, peak_Mbp-1, peak_Mbp+1)
  plot(out_snps_LiverAdjGSSG_chr2$lod, out_snps_LiverAdjGSSG_chr2$snpinfo, drop_hilit=1.5, genes = LiverAdjGSSG_Chr2_Genes, main = "Liver Adj GSSG Genes - sex")
  
#For Liver Adj GSSG -- sex covar -- chromosome 5
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  
  peak_Mbp <- max(scan1_output = qtlscan_LiverAdjGSSG_sex, map = R01_GSH_DO_QTLdata$pmap)$pos
  variants_LiverAdjGSSG_chr5 <- query_variants(5, peak_Mbp - 1, peak_Mbp + 1)
  out_snps_LiverAdjGSSG_chr5 <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["logLiverAdjGSSG"], kinship = kinship_loco[["5"]], addcovar = sex, query_func = query_variants,
                                          chr = 5, start=peak_Mbp-1, end=peak_Mbp+1, keep_all_snps = TRUE)
  plot_snpasso(out_snps_LiverAdjGSSG_chr5$lod, out_snps_LiverAdjGSSG_chr5$snpinfo, main = "Liver Adj GSSG SNPs - sex ")
  top <- top_snps(out_snps_LiverAdjGSSG_chr5$lod, out_snps_LiverAdjGSSG_chr5$snpinfo)
  print(top[,c(1, 8:15, 20)], row.names=FALSE)
  
  LiverAdjGSSG_Chr5_Genes_MGI <- query_genes_mgi(5, peak_Mbp-1, peak_Mbp+1)
  plot(out_snps_LiverAdjGSSG_chr5$lod, out_snps_LiverAdjGSSG_chr5$snpinfo, drop_hilit=1.5, genes = LiverAdjGSSG_Chr5_Genes_MGI, main = "Liver Adj GSSG Genes MGI")
  
  LiverAdjGSSG_Chr5_Genes <- query_genes(5, peak_Mbp-1, peak_Mbp+1)
  plot(out_snps_LiverAdjGSSG_chr5$lod, out_snps_LiverAdjGSSG_chr5$snpinfo, drop_hilit=1.5, genes = LiverAdjGSSG_Chr5_Genes, main = "Liver Adj GSSG Genes - sex")
  
#For Liver NADP -- sex covariate -- Chromosome 2
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  
  peak_Mbp <- max(scan1_output = qtlscan_LiverNADP_sex, map = R01_GSH_DO_QTLdata$pmap)$pos
  variants_LiverNADP_chr2 <- query_variants(2, peak_Mbp - 1, peak_Mbp + 1)
  out_snps_LiverNADP_chr2 <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["logLiverNADP"], kinship = kinship_loco[["2"]], addcovar = sex, query_func = query_variants,
                                      chr = 2, start=peak_Mbp-1, end=peak_Mbp+1, keep_all_snps = TRUE)
  plot_snpasso(out_snps_LiverNADP_chr2$lod, out_snps_LiverNADP_chr2$snpinfo, main = "Liver NADP SNPs - sex ")
  top <- top_snps(out_snps_LiverNADP_chr2$lod, out_snps_LiverNADP_chr2$snpinfo)
  print(top[,c(1, 8:15, 20)], row.names=FALSE)
  
  LiverNADP_Chr2_Genes_MGI <- query_genes_mgi(2, peak_Mbp-1, peak_Mbp+1)
  plot(out_snps_LiverNADP_chr2$lod, out_snps_LiverNADP_chr2$snpinfo, drop_hilit=1.5, genes = LiverNADP_Chr2_Genes_MGI, main = "Liver NADP Genes MGI - sex")
  
  LiverNADP_Chr2_Genes <- query_genes(2, peak_Mbp-1, peak_Mbp+1)
  plot(out_snps_LiverNADP_chr2$lod, out_snps_LiverNADP_chr2$snpinfo, drop_hilit=1.5, genes = LiverNADP_Chr2_Genes, main = "Liver NADP Genes - sex")
  
#For Liver NADPH -- sex covariate -- Chromosome 14
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  
  peak_Mbp <- max(scan1_output = qtlscan_LiverNADPH_sex, map = R01_GSH_DO_QTLdata$pmap)$pos
  variants_LiverNADPH_chr14 <- query_variants(14, peak_Mbp - 1, peak_Mbp + 1)
  out_snps_LiverNADPH_chr14 <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["logLiverNADPH"], kinship = kinship_loco[["14"]], addcovar = sex, query_func = query_variants,
                                       chr = 14, start=peak_Mbp-1, end=peak_Mbp+1, keep_all_snps = TRUE)
  plot_snpasso(out_snps_LiverNADPH_chr14$lod, out_snps_LiverNADPH_chr14$snpinfo, main = "Liver NADPH SNPs - sex ")
  top <- top_snps(out_snps_LiverNADPH_chr14$lod, out_snps_LiverNADPH_chr14$snpinfo)
  print(top[,c(1, 8:15, 20)], row.names=FALSE)
  
  LiverNADPH_Chr14_Genes_MGI <- query_genes_mgi(14, peak_Mbp-1, peak_Mbp+1)
  plot(out_snps_LiverNADPH_chr14$lod, out_snps_LiverNADPH_chr14$snpinfo, drop_hilit=1.5, genes = LiverNADPH_Chr14_Genes_MGI, main = "Liver NADPH Genes MGI - sex")
  
  LiverNADPH_Chr14_Genes <- query_genes(14, peak_Mbp-1, peak_Mbp+1)
  plot(out_snps_LiverNADPH_chr14$lod, out_snps_LiverNADPH_chr14$snpinfo, drop_hilit=1.5, genes = LiverNADPH_Chr14_Genes, main = "Liver NADPH Genes - sex")  

#For Liver NADP/NADPH Ratio -- sex covariate -- Chromosome 2
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  
  peak_Mbp <- max(scan1_output = qtlscan_LiverNADP_NADPH_sex, map = R01_GSH_DO_QTLdata$pmap)$pos
  variants_LiverNADP_NADPH_chr2 <- query_variants(2, peak_Mbp - 1, peak_Mbp + 1)
  out_snps_LiverNADP_NADPH_chr2 <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["logLiverNADP_NADPHRatio"], kinship = kinship_loco[["2"]], addcovar = sex, query_func = query_variants,
                                         chr = 2, start=peak_Mbp-1, end=peak_Mbp+1, keep_all_snps = TRUE)
  plot_snpasso(out_snps_LiverNADP_NADPH_chr2$lod, out_snps_LiverNADP_NADPH_chr2$snpinfo, main = "Liver NADP/NADPH Ratio SNPs - sex ")
  top <- top_snps(out_snps_LiverNADP_NADPH_chr2$lod, out_snps_LiverNADP_NADPH_chr2$snpinfo)
  print(top[,c(1, 8:15, 20)], row.names=FALSE)
  
  LiverNADP_NADPH_Chr2_Genes_MGI <- query_genes_mgi(2, peak_Mbp-1, peak_Mbp+1)
  plot(out_snps_LiverNADP_NADPH_chr2$lod, out_snps_LiverNADP_NADPH_chr2$snpinfo, drop_hilit=1.5, genes = LiverNADP_NADPH_Chr2_Genes_MGI, main = "Liver NADP/NADPH Genes MGI - sex")
  
  LiverNADP_NADPH_Chr2_Genes <- query_genes(2, peak_Mbp-1, peak_Mbp+1)
  plot(out_snps_LiverNADP_NADPH_chr2$lod, out_snps_LiverNADP_NADPH_chr2$snpinfo, drop_hilit=1.5, genes = LiverNADP_NADPH_Chr2_Genes, main = "Liver NADP/NADPH Genes - sex")  

#For Liver NADP/NADPH Ratio -- sex covariate -- Chromosome 14
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  
  peak_Mbp <- max(scan1_output = qtlscan_LiverNADP_NADPH_sex, map = R01_GSH_DO_QTLdata$pmap)$pos
  variants_LiverNADP_NADPH_chr14 <- query_variants(14, peak_Mbp - 1, peak_Mbp + 1)
  out_snps_LiverNADP_NADPH_chr14 <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["logLiverNADP_NADPHRatio"], kinship = kinship_loco[["14"]], addcovar = sex, query_func = query_variants,
                                         chr = 14, start=peak_Mbp-1, end=peak_Mbp+1, keep_all_snps = TRUE)
  plot_snpasso(out_snps_LiverNADP_NADPH_chr14$lod, out_snps_LiverNADP_NADPH_chr14$snpinfo, main = "Liver NADP/NADPH Ratio SNPs - sex ")
  top <- top_snps(out_snps_LiverNADP_NADPH_chr14$lod, out_snps_LiverNADP_NADPH_chr14$snpinfo)
  print(top[,c(1, 8:15, 20)], row.names=FALSE)
  
  LiverNADP_NADPH_Chr14_Genes_MGI <- query_genes_mgi(14, peak_Mbp-1, peak_Mbp+1)
  plot(out_snps_LiverNADP_NADPH_chr14$lod, out_snps_LiverNADP_NADPH_chr14$snpinfo, drop_hilit=1.5, genes = LiverNADP_NADPH_Chr14_Genes_MGI, main = "Liver NADP/NADPH Ratio Genes MGI - sex")
  
  LiverNADP_NADPH_Chr14_Genes <- query_genes(14, peak_Mbp-1, peak_Mbp+1)
  plot(out_snps_LiverNADP_NADPH_chr14$lod, out_snps_LiverNADP_NADPH_chr14$snpinfo, drop_hilit=1.5, genes = LiverNADP_NADPH_Chr14_Genes, main = "Liver NADP/NADPH Ratio Genes - sex")  

#For Liver NADP/NADPH Ratio -- sexgen covariate -- Chromosome 7
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  
  peak_Mbp <- max(scan1_output = qtlscan_LiverNADP_NADPH_sexgen, map = R01_GSH_DO_QTLdata$pmap)$pos
  variants_LiverNADP_NADPH_chr7_sexgen <- query_variants(7, peak_Mbp - 1, peak_Mbp + 1)
  out_snps_LiverNADP_NADPH_chr7_sexgen <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["logLiverNADP_NADPHRatio"], kinship = kinship_loco[["7"]], addcovar = sexgen, query_func = query_variants,
                                              chr = 7, start=peak_Mbp-1, end=peak_Mbp+1, keep_all_snps = TRUE)
  plot_snpasso(out_snps_LiverNADP_NADPH_chr7_sexgen$lod, out_snps_LiverNADP_NADPH_chr7_sexgen$snpinfo, main = "Liver NADP/NADPH Ratio SNPs - sexgen ")
  top <- top_snps(out_snps_LiverNADP_NADPH_chr7_sexgen$lod, out_snps_LiverNADP_NADPH_chr7_sexgen$snpinfo)
  print(top[,c(1, 8:15, 20)], row.names=FALSE)
  
  LiverNADP_NADPH_Chr7_Genes_MGI_sexgen <- query_genes_mgi(7, peak_Mbp-1, peak_Mbp+1)
  plot(out_snps_LiverNADP_NADPH_chr7_sexgen$lod, out_snps_LiverNADP_NADPH_chr7_sexgen$snpinfo, drop_hilit=1.5, genes = LiverNADP_NADPH_Chr7_Genes_MGI_sexgen, main = "Liver NADP/NADPH Ratio Genes MGI - sexgen")
  
  LiverNADP_NADPH_Chr7_Genes_sexgen <- query_genes(7, peak_Mbp-1, peak_Mbp+1)
  plot(out_snps_LiverNADP_NADPH_chr7_sexgen$lod, out_snps_LiverNADP_NADPH_chr7_sexgen$snpinfo, drop_hilit=1.5, genes = LiverNADP_NADPH_Chr7_Genes_sexgen, main = "Liver NADP/NADPH Ratio Genes - sexgen")  
  
#For Liver NADP/NADPH Ratio -- gen covariate -- Chromosome 7
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  
  peak_Mbp <- max(scan1_output = qtlscan_LiverNADP_NADPH_gen, map = R01_GSH_DO_QTLdata$pmap)$pos
  variants_LiverNADP_NADPH_chr7_gen <- query_variants(7, peak_Mbp - 1, peak_Mbp + 1)
  out_snps_LiverNADP_NADPH_chr7_gen <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["logLiverNADP_NADPHRatio"], kinship = kinship_loco[["7"]], addcovar = gen, query_func = query_variants,
                                             chr = 7, start=peak_Mbp-1, end=peak_Mbp+1, keep_all_snps = TRUE)
  plot_snpasso(out_snps_LiverNADP_NADPH_chr7_gen$lod, out_snps_LiverNADP_NADPH_chr7_gen$snpinfo, main = "Liver NADP/NADPH Ratio SNPs - gen ")
  top <- top_snps(out_snps_LiverNADP_NADPH_chr7_gen$lod, out_snps_LiverNADP_NADPH_chr7_gen$snpinfo)
  print(top[,c(1, 8:15, 20)], row.names=FALSE)
  
  LiverNADP_NADPH_Chr7_Genes_MGI_gen <- query_genes_mgi(7, peak_Mbp-1, peak_Mbp+1)
  plot(out_snps_LiverNADP_NADPH_chr7_gen$lod, out_snps_LiverNADP_NADPH_chr7_gen$snpinfo, drop_hilit=1.5, genes = LiverNADP_NADPH_Chr7_Genes_MGI_gen, main = "Liver NADP/NADPH Ratio Genes MGI - gen")
  
  LiverNADP_NADPH_Chr7_Genes_gen <- query_genes(7, peak_Mbp-1, peak_Mbp+1)
  plot(out_snps_LiverNADP_NADPH_chr7_gen$lod, out_snps_LiverNADP_NADPH_chr7_gen$snpinfo, drop_hilit=1.5, genes = LiverNADP_NADPH_Chr7_Genes_gen, main = "Liver NADP/NADPH Ratio Genes - gen")  
  
#For Liver Unadj GSSG -- sex covar -- chromosome 2
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  
  peak_Mbp <- max(scan1_output = qtlscan_LiverUnadjGSSG_sex, map = R01_GSH_DO_QTLdata$pmap)$pos
  variants_LiverUnadjGSSG_chr2 <- query_variants(2, peak_Mbp - 1, peak_Mbp + 1)
  out_snps_LiverUnadjGSSG_chr2 <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["logLiverUnadjGSSG"], kinship = kinship_loco[["2"]], addcovar = sex, query_func = query_variants,
                                          chr = 2, start=peak_Mbp-1, end=peak_Mbp+1, keep_all_snps = TRUE)
  plot_snpasso(out_snps_LiverUnadjGSSG_chr2$lod, out_snps_LiverUnadjGSSG_chr2$snpinfo, main = "Liver Unadj GSSG SNPs - sex")
  top <- top_snps(out_snps_LiverUnadjGSSG_chr2$lod, out_snps_LiverUnadjGSSG_chr2$snpinfo)
  print(top[,c(1, 8:15, 20)], row.names=FALSE)
  
  LiverUnadjGSSG_Chr2_Genes_MGI <- query_genes_mgi(2, peak_Mbp-1, peak_Mbp+1)
  plot(out_snps_LiverUnadjGSSG_chr2$lod, out_snps_LiverUnadjGSSG_chr2$snpinfo, drop_hilit=1.5, genes = LiverUnadjGSSG_Chr2_Genes_MGI, main = "Liver Unadj GSSG Genes MGI - sex")
  
  LiverUnadjGSSG_Chr2_Genes <- query_genes(2, peak_Mbp-1, peak_Mbp+1)
  plot(out_snps_LiverUnadjGSSG_chr2$lod, out_snps_LiverUnadjGSSG_chr2$snpinfo, drop_hilit=1.5, genes = LiverUnadjGSSG_Chr2_Genes, main = "Liver Unadj GSSG Genes - sex")
  
#For Liver Unadj Total GSH -- sex covariate -- Chromosome 2
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  
  peak_Mbp <- max(scan1_output = qtlscan_LiverUnadjTotalGSH_sex, map = R01_GSH_DO_QTLdata$pmap)$pos
  variants_LiverUnadjTotalGSH_chr2 <- query_variants(2, peak_Mbp - 1, peak_Mbp + 1)
  out_snps_LiverUnadjTotalGSH_chr2 <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["logLiverUnadjTotalGSH"], kinship = kinship_loco[["2"]], addcovar = sex, query_func = query_variants,
                                      chr = 2, start=peak_Mbp-1, end=peak_Mbp+1, keep_all_snps = TRUE)
  plot_snpasso(out_snps_LiverUnadjTotalGSH_chr2$lod, out_snps_LiverUnadjTotalGSH_chr2$snpinfo, main = "Liver Unadj Total GSH SNPs - sex ")
  top <- top_snps(out_snps_LiverUnadjTotalGSH_chr2$lod, out_snps_LiverUnadjTotalGSH_chr2$snpinfo)
  print(top[,c(1, 8:15, 20)], row.names=FALSE)
  
  LiverUnadjTotalGSH_Chr2_Genes_MGI <- query_genes_mgi(2, peak_Mbp-1, peak_Mbp+1)
  plot(out_snps_LiverUnadjTotalGSH_chr2$lod, out_snps_LiverUnadjTotalGSH_chr2$snpinfo, drop_hilit=1.5, genes = LiverUnadjTotalGSH_Chr2_Genes_MGI, main = "Liver Unadj Total GSH Genes MGI - sex")
  
  LiverUnadjTotalGSH_Chr2_Genes <- query_genes(2, peak_Mbp-1, peak_Mbp+1)
  plot(out_snps_LiverUnadjTotalGSH_chr2$lod, out_snps_LiverUnadjTotalGSH_chr2$snpinfo, drop_hilit=1.5, genes = LiverUnadjTotalGSH_Chr2_Genes, main = "Liver Unadj Total GSH Genes - sex")  
  
  
####################################################
## GWAS SNP Association Scan
## Make a Manhattan plot of the results; use altcol to define a color alternate for chromosomes and gap=0 to have no gap between chromosomes
####################################################

#set parameters for the plot
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  
#For Liver GSH
  #sex and gen covariates
  out_gwas_LiverGSH_sexgen <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["logLiverGSH"], kinship = kinship_loco, addcovar = sexgen, query_func=query_variants, cores=10)
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  plot(out_gwas_LiverGSH_sexgen$lod, out_gwas_LiverGSH_sexgen$snpinfo, altcol="green4", gap=0, main = " Liver GSH GWAS - sex and gen")
  
  #sex covariate
  out_gwas_LiverGSH_sex <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["logLiverGSH"], kinship = kinship_loco, addcovar = sex, query_func=query_variants, cores=10)
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  plot(out_gwas_LiverGSH_sex$lod, out_gwas_LiverGSH_sex$snpinfo, altcol="green4", gap=0, main = " Liver GSH GWAS - sex")
  
  #gen covariate
  out_gwas_LiverGSH_gen <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["logLiverGSH"], kinship = kinship_loco, addcovar = gen, query_func=query_variants, cores=10)
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  plot(out_gwas_LiverGSH_gen$lod, out_gwas_LiverGSH_gen$snpinfo, altcol="green4", gap=0, main = " Liver GSH GWAS - gen")
  
#For Liver Adj GSSG
  #sex and gen covariates
  out_gwas_LiverAdjGSSG_sexgen <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["logLiverAdjGSSG"], kinship = kinship_loco, addcovar = sexgen, query_func=query_variants, cores=10)
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  plot(out_gwas_LiverAdjGSSG_sexgen$lod, out_gwas_LiverAdjGSSG_sexgen$snpinfo, altcol="green4", gap=0)
  
  #sex covariate
  out_gwas_LiverAdjGSSG_sex <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["logLiverAdjGSSG"], kinship = kinship_loco, addcovar = sex, query_func=query_variants, cores=10)
  plot(out_gwas_LiverAdjGSSG_sex$lod, out_gwas_LiverAdjGSSG_sex$snpinfo, altcol="green4", gap=0)
  
  #gen covariate
  out_gwas_LiverAdjGSSG_gen <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["logLiverAdjGSSG"], kinship = kinship_loco, addcovar = gen, query_func=query_variants, cores=10)
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  plot(out_gwas_LiverAdjGSSG_gen$lod, out_gwas_LiverAdjGSSG_gen$snpinfo, altcol="green4", gap=0)
  
#For Liver Adj Total GSH
  #sex and gen covariates
  out_gwas_LiverAdjTotalGSH_sexgen <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["logLiverAdjTotalGSH"], kinship = kinship_loco, addcovar = sexgen, query_func=query_variants, cores=10)
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  plot(out_gwas_LiverAdjTotalGSH_sexgen$lod, out_gwas_LiverAdjTotalGSH_sexgen$snpinfo, altcol="green4", gap=0)
  
  #sex covariate
  out_gwas_LiverAdjTotalGSH_sex <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["logLiverAdjTotalGSH"], kinship = kinship_loco, addcovar = sex, query_func=query_variants, cores=10)
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  plot(out_gwas_LiverAdjTotalGSH_sex$lod, out_gwas_LiverAdjTotalGSH_sex$snpinfo, altcol="green4", gap=0)
  
  #gen covariate
  out_gwas_LiverAdjTotalGSH_gen <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["logLiverAdjTotalGSH"], kinship = kinship_loco, addcovar = gen, query_func=query_variants, cores=10)
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  plot(out_gwas_LiverAdjTotalGSH_gen$lod, out_gwas_LiverAdjTotalGSH_gen$snpinfo, altcol="green4", gap=0)

#For Liver Adj GSH/GSSG Ratio
  #sex and gen covariates
  out_gwas_LiverAdjGSH_GSSG_sexgen <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["logLiverAdjGSH_GSSGRatio"], kinship = kinship_loco, addcovar = sexgen, query_func=query_variants, cores=10)
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  plot(out_gwas_LiverAdjGSH_GSSG_sexgen$lod, out_gwas_LiverAdjGSH_GSSG_sexgen$snpinfo, altcol="green4", gap=0)
  
  #sex covariate
  out_gwas_LiverAdjGSH_GSSG_sex <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["logLiverAdjGSH_GSSGRatio"], kinship = kinship_loco, addcovar = sex, query_func=query_variants, cores=10)
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  plot(out_gwas_LiverAdjGSH_GSSG_sex$lod, out_gwas_LiverAdjGSH_GSSG_sex$snpinfo, altcol="green4", gap=0)
  
  #gen covariate
  out_gwas_LiverAdjGSH_GSSG_gen <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["logLiverAdjGSH_GSSGRatio"], kinship = kinship_loco, addcovar = gen, query_func=query_variants, cores=10)
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  plot(out_gwas_LiverAdjGSH_GSSG_gen$lod, out_gwas_LiverAdjGSH_GSSG_gen$snpinfo, altcol="green4", gap=0)

#For Liver NADH
  #sex and gen covariates
  out_gwas_LiverNADH_sexgen <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["logLiverNADH"], kinship = kinship_loco, addcovar = sexgen, query_func=query_variants, cores=10)
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  plot(out_gwas_LiverNADH_sexgen$lod, out_gwas_LiverNADH_sexgen$snpinfo, altcol="green4", gap=0)
  
  #sex covariate
  out_gwas_LiverNADH_sex <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["logLiverNADH"], kinship = kinship_loco, addcovar = sex, query_func=query_variants, cores=10)
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  plot(out_gwas_LiverNADH_sex$lod, out_gwas_LiverNADH_sex$snpinfo, altcol="green4", gap=0)
  
  #gen covariate
  out_gwas_LiverNADH_gen <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["logLiverNADH"], kinship = kinship_loco, addcovar = gen, query_func=query_variants, cores=10)
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  plot(out_gwas_LiverNADH_gen$lod, out_gwas_LiverNADH_gen$snpinfo, altcol="green4", gap=0)  
  
#For Liver NADP
  #sex and gen covariates
  out_gwas_LiverNADP_sexgen <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["logLiverNADP"], kinship = kinship_loco, addcovar = sexgen, query_func=query_variants, cores=10)
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  plot(out_gwas_LiverNADP_sexgen$lod, out_gwas_LiverNADP_sexgen$snpinfo, altcol="green4", gap=0)
  
  #sex covariate
  out_gwas_LiverNADP_sex <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["logLiverNADP"], kinship = kinship_loco, addcovar = sex, query_func=query_variants, cores=10)
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  plot(out_gwas_LiverNADP_sex$lod, out_gwas_LiverNADP_sex$snpinfo, altcol="green4", gap=0)
  
  #gen covariate
  out_gwas_LiverNADP_gen <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["logLiverNADP"], kinship = kinship_loco, addcovar = gen, query_func=query_variants, cores=10)
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  plot(out_gwas_LiverNADP_gen$lod, out_gwas_LiverNADP_gen$snpinfo, altcol="green4", gap=0)    
  
#For Liver NADPH
  #sex and gen covariates
  out_gwas_LiverNADPH_sexgen <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["logLiverNADPH"], kinship = kinship_loco, addcovar = sexgen, query_func=query_variants, cores=10)
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  plot(out_gwas_LiverNADPH_sexgen$lod, out_gwas_LiverNADPH_sexgen$snpinfo, altcol="green4", gap=0)
  
  #sex covariate
  out_gwas_LiverNADPH_sex <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["logLiverNADPH"], kinship = kinship_loco, addcovar = sex, query_func=query_variants, cores=10)
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  plot(out_gwas_LiverNADPH_sex$lod, out_gwas_LiverNADPH_sex$snpinfo, altcol="green4", gap=0)
  
  #gen covariate
  out_gwas_LiverNADPH_gen <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["logLiverNADPH"], kinship = kinship_loco, addcovar = gen, query_func=query_variants, cores=10)
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  plot(out_gwas_LiverNADPH_gen$lod, out_gwas_LiverNADPH_gen$snpinfo, altcol="green4", gap=0)    
  
#For Liver NADP/NADPH Ratio
  #sex and gen covariates
  out_gwas_LiverNADP_NADPH_sexgen <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["logLiverNADP_NADPHRatio"], kinship = kinship_loco, addcovar = sexgen, query_func=query_variants, cores=10)
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  plot(out_gwas_LiverNADP_NADPH_sexgen$lod, out_gwas_LiverNADP_NADPH_sexgen$snpinfo, altcol="green4", gap=0)
  
  #sex covariate
  out_gwas_LiverNADP_NADPH_sex <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["logLiverNADP_NADPHRatio"], kinship = kinship_loco, addcovar = sex, query_func=query_variants, cores=10)
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  plot(out_gwas_LiverNADP_NADPH_sex$lod, out_gwas_LiverNADP_NADPH_sex$snpinfo, altcol="green4", gap=0)
  
  #gen covariate
  out_gwas_LiverNADP_NADPH_gen <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["logLiverNADP_NADPHRatio"], kinship = kinship_loco, addcovar = gen, query_func=query_variants, cores=10)
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  plot(out_gwas_LiverNADP_NADPH_gen$lod, out_gwas_LiverNADP_NADPH_gen$snpinfo, altcol="green4", gap=0)      

#For Liver Unadj GSSG
  #sex and gen covariates
  out_gwas_LiverUnadjGSSG_sexgen <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["logLiverUnadjGSSG"], kinship = kinship_loco, addcovar = sexgen, query_func=query_variants, cores=10)
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  plot(out_gwas_LiverUnadjGSSG_sexgen$lod, out_gwas_LiverUnadjGSSG_sexgen$snpinfo, altcol="green4", gap=0)
  
  #sex covariate
  out_gwas_LiverUnadjGSSG_sex <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["logLiverUnadjGSSG"], kinship = kinship_loco, addcovar = sex, query_func=query_variants, cores=10)
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  plot(out_gwas_LiverUnadjGSSG_sex$lod, out_gwas_LiverUnadjGSSG_sex$snpinfo, altcol="green4", gap=0)
  
  #gen covariate
  out_gwas_LiverUnadjGSSG_gen <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["logLiverUnadjGSSG"], kinship = kinship_loco, addcovar = gen, query_func=query_variants, cores=10)
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  plot(out_gwas_LiverUnadjGSSG_gen$lod, out_gwas_LiverUnadjGSSG_gen$snpinfo, altcol="green4", gap=0)
  
#For Liver Unadj Total GSH
  #sex and gen covariates
  out_gwas_LiverUnadjTotalGSH_sexgen <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["logLiverUnadjTotalGSH"], kinship = kinship_loco, addcovar = sexgen, query_func=query_variants, cores=10)
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  plot(out_gwas_LiverUnadjTotalGSH_sexgen$lod, out_gwas_LiverUnadjTotalGSH_sexgen$snpinfo, altcol="green4", gap=0)
  
  #sex covariate
  out_gwas_LiverUnadjTotalGSH_sex <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["logLiverUnadjTotalGSH"], kinship = kinship_loco, addcovar = sex, query_func=query_variants, cores=10)
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  plot(out_gwas_LiverUnadjTotalGSH_sex$lod, out_gwas_LiverUnadjTotalGSH_sex$snpinfo, altcol="green4", gap=0)
  
  #gen covariate
  out_gwas_LiverUnadjTotalGSH_gen <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["logLiverUnadjTotalGSH"], kinship = kinship_loco, addcovar = gen, query_func=query_variants, cores=10)
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  plot(out_gwas_LiverUnadjTotalGSH_gen$lod, out_gwas_LiverUnadjTotalGSH_gen$snpinfo, altcol="green4", gap=0)
  
#For Liver Unadj GSH/GSSG Ratio
  #sex and gen covariates
  out_gwas_LiverUnadjGSH_GSSG_sexgen <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["logLiverUnadjGSH_GSSGRatio"], kinship = kinship_loco, addcovar = sexgen, query_func=query_variants, cores=10)
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  plot(out_gwas_LiverUnadjGSH_GSSG_sexgen$lod, out_gwas_LiverUnadjGSH_GSSG_sexgen$snpinfo, altcol="green4", gap=0)
  
  #sex covariate
  out_gwas_LiverUnadjGSH_GSSG_sex <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["logLiverUnadjGSH_GSSGRatio"], kinship = kinship_loco, addcovar = sex, query_func=query_variants, cores=10)
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  plot(out_gwas_LiverUnadjGSH_GSSG_sex$lod, out_gwas_LiverUnadjGSH_GSSG_sex$snpinfo, altcol="green4", gap=0)
  
  #gen covariate
  out_gwas_LiverUnadjGSH_GSSG_gen <- scan1snps(genoprobs = probs, map = R01_GSH_DO_QTLdata$pmap, pheno = pheno["logLiverUnadjGSH_GSSGRatio"], kinship = kinship_loco, addcovar = gen, query_func=query_variants, cores=10)
  par(mar=c(4.1, 4.1, 2.6, 2.6))
  plot(out_gwas_LiverUnadjGSH_GSSG_gen$lod, out_gwas_LiverUnadjGSH_GSSG_gen$snpinfo, altcol="green4", gap=0)  
  
####################################################
## Heritability calculation
## Duy calculated heritability using the qtl2 function est_herit to estimate heritability with linear mixed model
## This calculates the ratio of genetic variance to total variance
####################################################

#make kinship function using linear mixed model, not loco
#default type of kinship is "overall" aka "linear mixed model" -- did not need to specify a type
kinship_lmm <- calc_kinship(probs = probs, use_allele_probs = TRUE, cores = 10)
  
#Liver GSH
  #sex and gen covariates
  herit_LiverGSH_sexgen <- est_herit(pheno["logLiverGSH"], kinship_lmm, sexgen, cores = 10)
  #sex covariate
  herit_LiverGSH_sex <- est_herit(pheno["logLiverGSH"], kinship_lmm, sex, cores = 10)
  #gen covariate
  herit_LiverGSH_gen <- est_herit(pheno["logLiverGSH"], kinship_lmm, gen, cores = 10)
 
#Liver Adj GSSG
  #sex and gen covariates
  herit_LiverAdjGSSG_sexgen <- est_herit(pheno["logLiverAdjGSSG"], kinship_lmm, sexgen, cores = 10)
  #sex covariate
  herit_LiverAdjGSSG_sex <- est_herit(pheno["logLiverAdjGSSG"], kinship_lmm, sex, cores = 10)
  #gen covariate
  herit_LiverAdjGSSG_gen <- est_herit(pheno["logLiverAdjGSSG"], kinship_lmm, gen, cores = 10)
  
#Liver Adj Total GSH
  #sex and gen covariates
  herit_LiverAdjTotalGSH_sexgen <- est_herit(pheno["logLiverAdjTotalGSH"], kinship_lmm, sexgen, cores = 10)
  #sex covariate
  herit_LiverAdjTotalGSH_sex <- est_herit(pheno["logLiverAdjTotalGSH"], kinship_lmm, sex, cores = 10)
  #gen covariate
  herit_LiverAdjTotalGSH_gen <- est_herit(pheno["logLiverAdjTotalGSH"], kinship_lmm, gen, cores = 10)
  
#Liver Adj GSH/GSSG Ratio
  #sex and gen covariates
  herit_LiverAdjGSH_GSSG_sexgen <- est_herit(pheno["logLiverAdjGSH_GSSGRatio"], kinship_lmm, sexgen, cores = 10)
  #sex covariate
  herit_LiverAdjGSH_GSSG_sex <- est_herit(pheno["logLiverAdjGSH_GSSGRatio"], kinship_lmm, sex, cores = 10)
  #gen covariate
  herit_LiverAdjGSH_GSSG_gen <- est_herit(pheno["logLiverAdjGSH_GSSGRatio"], kinship_lmm, gen, cores = 10)
  
#Liver NADH
  #sex and gen covariates
  herit_LiverNADH_sexgen <- est_herit(pheno["logLiverNADH"], kinship_lmm, sexgen, cores = 10)
  #sex covariate
  herit_LiverNADH_sex <- est_herit(pheno["logLiverNADH"], kinship_lmm, sex, cores = 10)
  #gen covariate
  herit_LiverNADH_gen <- est_herit(pheno["logLiverNADH"], kinship_lmm, gen, cores = 10)  

#Liver NADP
  #sex and gen covariates
  herit_LiverNADP_sexgen <- est_herit(pheno["logLiverNADP"], kinship_lmm, sexgen, cores = 10)
  #sex covariate
  herit_LiverNADP_sex <- est_herit(pheno["logLiverNADP"], kinship_lmm, sex, cores = 10)
  #gen covariate
  herit_LiverNADP_gen <- est_herit(pheno["logLiverNADP"], kinship_lmm, gen, cores = 10)  
  
#Liver NADPH
  #sex and gen covariates
  herit_LiverNADPH_sexgen <- est_herit(pheno["logLiverNADPH"], kinship_lmm, sexgen, cores = 10)
  #sex covariate
  herit_LiverNADPH_sex <- est_herit(pheno["logLiverNADPH"], kinship_lmm, sex, cores = 10)
  #gen covariate
  herit_LiverNADPH_gen <- est_herit(pheno["logLiverNADPH"], kinship_lmm, gen, cores = 10)  
  
#Liver NADP/NADH Ratio
  #sex and gen covariates
  herit_LiverNADP_NADPH_sexgen <- est_herit(pheno["logLiverNADP_NADPHRatio"], kinship_lmm, sexgen, cores = 10)
  #sex covariate
  herit_LiverNADP_NADPH_sex <- est_herit(pheno["logLiverNADP_NADPHRatio"], kinship_lmm, sex, cores = 10)
  #gen covariate
  herit_LiverNADP_NADPH_gen <- est_herit(pheno["logLiverNADP_NADPHRatio"], kinship_lmm, gen, cores = 10)  
  
#AST
  #sex and gen covariates
  herit_AST_sexgen <- est_herit(pheno["logAST"], kinship_lmm, sexgen, cores = 10)
  #sex covariate
  herit_AST_sex <- est_herit(pheno["logAST"], kinship_lmm, sex, cores = 10)
  #gen covariate
  herit_AST_gen <- est_herit(pheno["logAST"], kinship_lmm, gen, cores = 10)    
  
#ALT
  #sex and gen covariates
  herit_ALT_sexgen <- est_herit(pheno["logALT"], kinship_lmm, sexgen, cores = 10)
  #sex covariate
  herit_ALT_sex <- est_herit(pheno["logALT"], kinship_lmm, sex, cores = 10)
  #gen covariate
  herit_ALT_gen <- est_herit(pheno["logALT"], kinship_lmm, gen, cores = 10)    
  
#BUN
  #sex and gen covariates
  herit_BUN_sexgen <- est_herit(pheno["logBUN"], kinship_lmm, sexgen, cores = 10)
  #sex covariate
  herit_BUN_sex <- est_herit(pheno["logBUN"], kinship_lmm, sex, cores = 10)
  #gen covariate
  herit_BUN_gen <- est_herit(pheno["logBUN"], kinship_lmm, gen, cores = 10)   

#Glucose
  #sex and gen covariates
  herit_Glucose_sexgen <- est_herit(pheno["logGlucose"], kinship_lmm, sexgen, cores = 10)
  #sex covariate
  herit_Glucose_sex <- est_herit(pheno["logGlucose"], kinship_lmm, sex, cores = 10)
  #gen covariate
  herit_Glucose_gen <- est_herit(pheno["logGlucose"], kinship_lmm, gen, cores = 10)       

#Liver Weight
  #sex and gen covariates
  herit_LiverWeight_sexgen <- est_herit(pheno["logLiver_Weight"], kinship_lmm, sexgen, cores = 10)
  #sex covariate
  herit_LiverWeight_sex <- est_herit(pheno["logLiver_Weight"], kinship_lmm, sex, cores = 10)
  #gen covariate
  herit_LiverWeight_gen <- est_herit(pheno["logLiver_Weight"], kinship_lmm, gen, cores = 10)       
  
#Initial Weight
  #sex and gen covariates
  herit_InitialWeight_sexgen <- est_herit(pheno["logInitial_Weight"], kinship_lmm, sexgen, cores = 10)
  #sex covariate
  herit_InitialWeight_sex <- est_herit(pheno["logInitial_Weight"], kinship_lmm, sex, cores = 10)
  #gen covariate
  herit_InitialWeight_gen <- est_herit(pheno["logInitial_Weight"], kinship_lmm, gen, cores = 10)   
  
#Final Weight
  #sex and gen covariates
  herit_FinalWeight_sexgen <- est_herit(pheno["logFinal_Weight"], kinship_lmm, sexgen, cores = 10)
  #sex covariate
  herit_FinalWeight_sex <- est_herit(pheno["logFinal_Weight"], kinship_lmm, sex, cores = 10)
  #gen covariate
  herit_FinalWeight_gen <- est_herit(pheno["logFinal_Weight"], kinship_lmm, gen, cores = 10)   
  
#Final Weight
  #sex and gen covariates
  herit_LiverWeight_FinalWeight_sexgen <- est_herit(pheno["logLiverWeight_FinalWeight"], kinship_lmm, sexgen, cores = 10)
  #sex covariate
  herit_LiverWeight_FinalWeight_sex <- est_herit(pheno["logLiverWeight_FinalWeight"], kinship_lmm, sex, cores = 10)
  #gen covariate
  herit_LiverWeight_FinalWeight_gen <- est_herit(pheno["logLiverWeight_FinalWeight"], kinship_lmm, gen, cores = 10)    
  
#Liver Unadj GSSG
  #sex and gen covariates
  herit_LiverUnadjGSSG_sexgen <- est_herit(pheno["logLiverUnadjGSSG"], kinship_lmm, sexgen, cores = 10)
  #sex covariate
  herit_LiverUnadjGSSG_sex <- est_herit(pheno["logLiverUnadjGSSG"], kinship_lmm, sex, cores = 10)
  #gen covariate
  herit_LiverUnadjGSSG_gen <- est_herit(pheno["logLiverUnadjGSSG"], kinship_lmm, gen, cores = 10)
  
#Liver Unadj Total GSH
  #sex and gen covariates
  herit_LiverUnadjTotalGSH_sexgen <- est_herit(pheno["logLiverUnadjTotalGSH"], kinship_lmm, sexgen, cores = 10)
  #sex covariate
  herit_LiverUnadjTotalGSH_sex <- est_herit(pheno["logLiverUnadjTotalGSH"], kinship_lmm, sex, cores = 10)
  #gen covariate
  herit_LiverUnadjTotalGSH_gen <- est_herit(pheno["logLiverUnadjTotalGSH"], kinship_lmm, gen, cores = 10)
  
#Liver Unadj GSH/GSSG Ratio
  #sex and gen covariates
  herit_LiverUnadjGSH_GSSG_sexgen <- est_herit(pheno["logLiverUnadjGSH_GSSGRatio"], kinship_lmm, sexgen, cores = 10)
  #sex covariate
  herit_LiverUnadjGSH_GSSG_sex <- est_herit(pheno["logLiverUnadjGSH_GSSGRatio"], kinship_lmm, sex, cores = 10)
  #gen covariate
  herit_LiverUnadjGSH_GSSG_gen <- est_herit(pheno["logLiverUnadjGSH_GSSGRatio"], kinship_lmm, gen, cores = 10)  
  
  