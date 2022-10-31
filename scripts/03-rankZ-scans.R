library(qtl2)

# load data only if you don't already have it in your RStudio Environment
# remove the hashtag and run the code below if you don't have it

# load("../data/Heart-GSSG-enviornment.RData")

# run a genome scan on the rankZ transformed phenotypes
zScan <- scan1(genoprobs = probs, 
             pheno = control$pheno[, c("zHeart_GSH", 
                                       "zHeart_GSSG", 
                                       "zHeart_TotalGSH", 
                                       "zGSHGSSG_Ratio")])

# look at the first several rows of the scan
head(zScan)

# LOD score significance thresholds
perms <- scan1perm(genoprobs = probs, 
                   pheno = control$pheno[, c("zHeart_GSH", 
                                             "zHeart_GSSG", 
                                             "zHeart_TotalGSH", 
                                             "zGSHGSSG_Ratio")], 
                   n_perm = 10) # using only 10 permutations for expediency;
                                # replace with 1000 before going to press
summary(perms)
thresholds <- summary(perms)

# use the thresholds to find LOD peaks that exceed the thresholds
find_peaks(scan1_output = zScan, 
           map = control$gmap, 
           threshold = thresholds, 
           prob = 0.95, 
           expand2markers = FALSE)
# heart GSH has 3 peaks above the threshold (chr 14, 16 & 19)
# heart GSSG has 1 peak above the threshold on chr 10

pdf(file = "../results/rankZscans.pdf")
# plot scans for each phenotype
plot_scan1(zScan, map = control$gmap, 
           lodcolumn = "zHeart_GSH", main = "RankZ transformed GSH genome scan")
abline(h = thresholds[1], col = "red", lwd = 2)

plot_scan1(zScan, map = control$gmap, 
           lodcolumn = "zHeart_GSSG", main = "RankZ transformed GSSG genome scan")
abline(h = thresholds[2], col = "red", lwd = 2)

plot_scan1(zScan, map = control$gmap, 
           lodcolumn = "zHeart_TotalGSH", main = "RankZ transformed total GSH genome scan")
abline(h = thresholds[3], col = "red", lwd = 2)

plot_scan1(zScan, map = control$gmap, 
           lodcolumn = "zGSHGSSG_Ratio", main = "RankZ transformed GSH/GSSG Ratio")
abline(h = thresholds[3], col = "red", lwd = 2)

# for comparison repeat a genome scan on the untransformed phenotypes
rawScan <- scan1(genoprobs = probs,
                 pheno = control$pheno[, c("Heart_GSH",
                                           "Heart_GSSG",
                                           "Heart_TotalGSH",
                                           "Heart_GSHGSSGRatio")])

# LOD score significance thresholds
rawperms <- scan1perm(genoprobs = probs,
                      pheno = control$pheno[, c("Heart_GSH",
                                                "Heart_GSSG",
                                                "Heart_TotalGSH",
                                                "Heart_GSHGSSGRatio")],
                      n_perm = 10) # using only 10 permutations for expediency;
rawthresholds <- summary(rawperms)
rawthresholds
find_peaks(scan1_output = rawScan,
           map = control$gmap,
           threshold = rawthresholds,
           prob = 0.95,
           expand2markers = FALSE)

plot_scan1(rawScan, map = control$gmap, 
           lodcolumn = "Heart_GSHGSSGRatio", 
           main = "Untransformed GSH/GSSG Ratio")
abline(h = rawthresholds[4], col = "red", lwd = 2)

plot_scan1(rawScan, map = control$gmap, 
           lodcolumn = "Heart_GSH", 
           main = "Untransformed Heart GSH genome scan")
abline(h = rawthresholds[1], col = "red", lwd = 2)

plot_scan1(rawScan, map = control$gmap, 
           lodcolumn = "Heart_GSSG", 
           main = "Untransformed Heart GSSG genome scan")
abline(h = rawthresholds[2], col = "red", lwd = 2)

plot_scan1(rawScan, map = control$gmap, 
           lodcolumn = "Heart_TotalGSH", 
           main = "Untransformed Total Heart GSH genome scan")
abline(h = rawthresholds[3], col = "red", lwd = 2)

dev.off()

