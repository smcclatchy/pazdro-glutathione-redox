library(qtl2)

# load data only if you don't already have it in your RStudio Environment
# remove the hashtag and run the code below if you don't have it

# load("../data/Heart-GSSG-enviornment.RData")

# run a genome scan on the rankZ transformed phenotypes
# use sex and generation as covariates

pheno$sex <- pheno$sex
pheno$generation <- pheno$generation
# create a data frame with sex and generation
addcovarSexGen <- model.matrix(~ sex + generation, 
                               data = pheno)[,-1]

zScanSex <- scan1(genoprobs = probs, 
               pheno = control$pheno[, c("zHeart_GSH", 
                                         "zHeart_GSSG", 
                                         "zHeart_TotalGSH", 
                                         "zGSHGSSG_Ratio")],
               addcovar = addcovarSexGen[, "sex"])

# run separate LOD score thresholds for sex and generation for comparison
# with one another and with combined sex + generation
Sexperms <- scan1perm(genoprobs = probs, 
                      pheno = control$pheno[, c("zHeart_GSH", 
                                                "zHeart_GSSG", 
                                                "zHeart_TotalGSH", 
                                                "zGSHGSSG_Ratio")], 
                      addcovar = addcovarSexGen[, "sex"],
                      n_perm = 10)
summary(Sexperms)
Sexthresholds <- summary(Sexperms)

zScanGen <- scan1(genoprobs = probs, 
                  pheno = control$pheno[, c("zHeart_GSH", 
                                            "zHeart_GSSG", 
                                            "zHeart_TotalGSH", 
                                            "zGSHGSSG_Ratio")],
                  addcovar = addcovarSexGen[, "generation"])

Genperms <- scan1perm(genoprobs = probs, 
                      pheno = control$pheno[, c("zHeart_GSH", 
                                                "zHeart_GSSG", 
                                                "zHeart_TotalGSH", 
                                                "zGSHGSSG_Ratio")], 
                      addcovar = addcovarSexGen[, "generation"],
                      n_perm = 10)

summary(Genperms)
Genthresholds <- summary(Genperms)

zScanSexGen <- scan1(genoprobs = probs, 
                  pheno = control$pheno[, c("zHeart_GSH", 
                                            "zHeart_GSSG", 
                                            "zHeart_TotalGSH", 
                                            "zGSHGSSG_Ratio")],
                  addcovar = addcovarSexGen)


# LOD score significance thresholds for combined sex + generation
covarperms <- scan1perm(genoprobs = probs, 
                   pheno = control$pheno[, c("zHeart_GSH", 
                                             "zHeart_GSSG", 
                                             "zHeart_TotalGSH", 
                                             "zGSHGSSG_Ratio")], 
                   addcovar = addcovarSexGen,
                   n_perm = 10) # using only 10 permutations for expediency;
                                # replace with 1000 before going to press
summary(covarperms)
thresholds <- summary(covarperms)

Sexthresholds
Genthresholds
thresholds


# use the thresholds to find LOD peaks that exceed the thresholds
find_peaks(scan1_output = zScanSex, 
           map = control$gmap, 
           threshold = Sexthresholds, 
           prob = 0.95, 
           expand2markers = FALSE)

find_peaks(scan1_output = zScanGen, 
           map = control$gmap, 
           threshold = Genthresholds, 
           prob = 0.95, 
           expand2markers = FALSE)

find_peaks(scan1_output = zScanSexGen, 
           map = control$gmap, 
           threshold = thresholds, 
           prob = 0.95, 
           expand2markers = FALSE)

pdf(file = "../results/CovariaterankZscans.pdf")
# plot scans for each phenotype
plot_scan1(zScanSex, map = control$gmap, 
           lodcolumn = "zHeart_GSH", main = "GSH with sex as covariate")
abline(h = Sexthresholds[1], col = "red", lwd = 2)

plot_scan1(zScanSex, map = control$gmap, 
           lodcolumn = "zHeart_GSSG", main = "GSSG with sex as covariate")
abline(h = Sexthresholds[2], col = "red", lwd = 2)

plot_scan1(zScanSex, map = control$gmap, 
           lodcolumn = "zHeart_TotalGSH", main = "Total GSH with sex as covariate")
abline(h = Sexthresholds[3], col = "red", lwd = 2)

plot_scan1(zScanSex, map = control$gmap, 
           lodcolumn = "zGSHGSSG_Ratio", main = "GSH/GSSG Ratio with sex as covariate")
abline(h = Sexthresholds[4], col = "red", lwd = 2)

plot_scan1(zScanGen, map = control$gmap, 
           lodcolumn = "zHeart_GSH", main = "GSH with generation as covariate")
abline(h = Genthresholds[1], col = "red", lwd = 2)

plot_scan1(zScanGen, map = control$gmap, 
           lodcolumn = "zHeart_GSSG", main = "GSSG with generation as covariate")
abline(h = Genthresholds[2], col = "red", lwd = 2)

plot_scan1(zScanGen, map = control$gmap, 
           lodcolumn = "zHeart_TotalGSH", main = "Total GSH with generation as covariate")
abline(h = Genthresholds[3], col = "red", lwd = 2)

plot_scan1(zScanGen, map = control$gmap, 
           lodcolumn = "zGSHGSSG_Ratio", main = "GSH/GSSG Ratio with generation as covariate")
abline(h = Genthresholds[4], col = "red", lwd = 2)

plot_scan1(zScanSexGen, map = control$gmap, 
           lodcolumn = "zHeart_GSH", main = "GSH with sex and generation as covariate")
abline(h = thresholds[1], col = "red", lwd = 2)

plot_scan1(zScanSexGen, map = control$gmap, 
           lodcolumn = "zHeart_GSSG", main = "GSSG with sex and generation as covariate")
abline(h = thresholds[2], col = "red", lwd = 2)

plot_scan1(zScanSexGen, map = control$gmap, 
           lodcolumn = "zHeart_TotalGSH", main = "Total GSH with sex and generation as covariate")
abline(h = thresholds[3], col = "red", lwd = 2)

plot_scan1(zScanSexGen, map = control$gmap, 
           lodcolumn = "zGSHGSSG_Ratio", main = "GSH/GSSG Ratio with sex and generation as covariate")
abline(h = thresholds[4], col = "red", lwd = 2)

dev.off()
