library(qtl2)

# load data only if you don't already have it in your RStudio Environment
# remove the hashtag and run the code below if you don't have it

# load("../data/Heart-GSSG-enviornment.RData")

# run a genome scan on the rankZ transformed phenotypes
# use sex and generation as covariates

# create a data frame with sex and generation
addcovarSexGenBatch <- model.matrix(~ sex + generation + Set.No,
                               data = control$pheno)[,-1]

# genome scans with sex, generation and batch as individual covariates
# and all 3 combined as covariates
zScanSex <- scan1(genoprobs = probs,
               pheno = control$pheno[, c("zHeart_GSH",
                                         "zHeart_GSSG",
                                         "zHeart_TotalGSH",
                                         "zGSHGSSG_Ratio")],
               kinship = kinship_loco,
               addcovar = addcovarSexGenBatch[, "sex"])

zScanGen <- scan1(genoprobs = probs,
                  pheno = control$pheno[, c("zHeart_GSH",
                                            "zHeart_GSSG",
                                            "zHeart_TotalGSH",
                                            "zGSHGSSG_Ratio")],
                  kinship = kinship_loco,
                  addcovar = addcovarSexGenBatch[, "generation"])

zScanBatch <- scan1(genoprobs = probs,
                    pheno = control$pheno[, c("zHeart_GSH",
                                              "zHeart_GSSG",
                                              "zHeart_TotalGSH",
                                              "zGSHGSSG_Ratio")],
                    kinship = kinship_loco,
                    addcovar = addcovarSexGenBatch[, "Set.No"])

zScanSexGenBatch <- scan1(genoprobs = probs,
                          pheno = control$pheno[, c("zHeart_GSH",
                                                    "zHeart_GSSG",
                                                    "zHeart_TotalGSH",
                                                    "zGSHGSSG_Ratio")],
                          kinship = kinship_loco,
                          addcovar = addcovarSexGenBatch)

# run separate LOD score thresholds for sex and generation and batch for
# comparison with one another and with combined sex + generation + batch
Sexperms <- scan1perm(genoprobs = probs,
                      pheno = control$pheno[, c("zHeart_GSH",
                                                "zHeart_GSSG",
                                                "zHeart_TotalGSH",
                                                "zGSHGSSG_Ratio")],
                      addcovar = addcovarSexGenBatch[, "sex"],
                      n_perm = 10)
summary(Sexperms)
Sexthresholds <- summary(Sexperms)

Genperms <- scan1perm(genoprobs = probs,
                      pheno = control$pheno[, c("zHeart_GSH",
                                                "zHeart_GSSG",
                                                "zHeart_TotalGSH",
                                                "zGSHGSSG_Ratio")],
                      addcovar = addcovarSexGenBatch[, "generation"],
                      n_perm = 10)

summary(Genperms)
Genthresholds <- summary(Genperms)

Batchperms <- scan1perm(genoprobs = probs,
                      pheno = control$pheno[, c("zHeart_GSH",
                                                "zHeart_GSSG",
                                                "zHeart_TotalGSH",
                                                "zGSHGSSG_Ratio")],
                      addcovar = addcovarSexGenBatch[, "Set.No"],
                      n_perm = 10)
summary(Batchperms)
Batchthresholds <- summary(Batchperms)

# LOD score significance thresholds for combined sex + generation + batch
covarperms <- scan1perm(genoprobs = probs,
                   pheno = control$pheno[, c("zHeart_GSH",
                                             "zHeart_GSSG",
                                             "zHeart_TotalGSH",
                                             "zGSHGSSG_Ratio")],
                   addcovar = addcovarSexGenBatch,
                   n_perm = 10) # using only 10 permutations for expediency;
                                # replace with 1000 before going to press
summary(covarperms)
thresholds <- summary(covarperms)

Sexthresholds
Genthresholds
Batchthresholds
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

find_peaks(scan1_output = zScanBatch,
           map = control$gmap,
           threshold = Batchthresholds,
           prob = 0.95,
           expand2markers = FALSE)

find_peaks(scan1_output = zScanSexGenBatch,
           map = control$gmap,
           threshold = thresholds,
           prob = 0.95,
           expand2markers = FALSE)

# plot scans for each phenotype with different covariates side-by-side
# set plots on a 2x2 grid

par(mfcol = c(2,2))

# GSH scans
pdf(file = "../results/CovariateGSHscans.pdf")

plot_scan1(zScanSex, map = control$gmap,
           lodcolumn = "zHeart_GSH", main = "GSH",
           sub = "with sex as covariate")
abline(h = Sexthresholds[1], col = "red", lwd = 2)

plot_scan1(zScanGen, map = control$gmap,
           lodcolumn = "zHeart_GSH", main = "GSH",
           sub = "with generation as covariate")
abline(h = Genthresholds[1], col = "red", lwd = 2)

plot_scan1(zScanBatch, map = control$gmap,
           lodcolumn = "zHeart_GSH", main = "GSH",
           sub = "with batch as covariate")
abline(h = Batchthresholds[1], col = "red", lwd = 2)

plot_scan1(zScanSexGenBatch, map = control$gmap,
           lodcolumn = "zHeart_GSH", main = "GSH",
           sub = "with sex, generation and batch as covariate")
abline(h = thresholds[1], col = "red", lwd = 2)

dev.off()

# GSSG scans
pdf(file = "../results/CovariateGSSGscans.pdf")

plot_scan1(zScanSex, map = control$gmap,
           lodcolumn = "zHeart_GSSG", main = "GSSG",
           sub = "with sex as covariate")
abline(h = Sexthresholds[2], col = "red", lwd = 2)

plot_scan1(zScanGen, map = control$gmap,
           lodcolumn = "zHeart_GSSG", main = "GSSG",
           sub = "with generation as covariate")
abline(h = Genthresholds[2], col = "red", lwd = 2)

plot_scan1(zScanBatch, map = control$gmap,
           lodcolumn = "zHeart_GSSG", main = "GSSG",
           sub = "with batch as covariate")
abline(h = Batchthresholds[2], col = "red", lwd = 2)

plot_scan1(zScanSexGenBatch, map = control$gmap,
           lodcolumn = "zHeart_GSSG", main = "GSSG",
           sub = "with sex, generation and batch as covariate")
abline(h = thresholds[2], col = "red", lwd = 2)

dev.off()

# Total GSH scans
pdf(file = "../results/CovariateTotalGSHscans.pdf")

plot_scan1(zScanSex, map = control$gmap,
           lodcolumn = "zHeart_TotalGSH", main = "Total GSH",
           sub = "with sex as covariate")
abline(h = Sexthresholds[3], col = "red", lwd = 2)

plot_scan1(zScanGen, map = control$gmap,
           lodcolumn = "zHeart_TotalGSH", main = "Total GSH",
           sub = "with generation as covariate")
abline(h = Genthresholds[3], col = "red", lwd = 2)

plot_scan1(zScanBatch, map = control$gmap,
           lodcolumn = "zHeart_TotalGSH", main = "Total GSH",
           sub = "with batch as covariate")
abline(h = Batchthresholds[3], col = "red", lwd = 2)

plot_scan1(zScanSexGenBatch, map = control$gmap,
           lodcolumn = "zHeart_TotalGSH", main = "Total GSH",
           sub = "with sex, generation and batch as covariate")
abline(h = thresholds[3], col = "red", lwd = 2)

dev.off()

# GSH/GSSG ratio scans
pdf(file = "../results/CovariateRatioscans.pdf")

plot_scan1(zScanGen, map = control$gmap,
           lodcolumn = "zGSHGSSG_Ratio", main = "GSH/GSSG Ratio",
           sub = "with generation as covariate")
abline(h = Genthresholds[4], col = "red", lwd = 2)

plot_scan1(zScanSex, map = control$gmap,
           lodcolumn = "zGSHGSSG_Ratio", main = "GSH/GSSG Ratio",
           sub = "with sex as covariate")
abline(h = Sexthresholds[4], col = "red", lwd = 2)

plot_scan1(zScanBatch, map = control$gmap,
           lodcolumn = "zGSHGSSG_Ratio", main = "GSH/GSSG Ratio",
           sub = "with batch as covariate")
abline(h = Batchthresholds[4], col = "red", lwd = 2)

plot_scan1(zScanSexGenBatch, map = control$gmap,
           lodcolumn = "zGSHGSSG_Ratio", xlab = NULL, main = "GSH/GSSG Ratio",
           sub = "with sex, generation and batch as covariate")
abline(h = thresholds[4], col = "red", lwd = 2)

dev.off()
