library(qtl2)

# load data only if you don't already have it in your RStudio Environment
# remove the hashtag and run the code below if you don't have it

# load("../data/Heart-GSSG-enviornment.RData")

# run a genome scan on the rankZ transformed phenotypes
zScan <- scan1(genoprobs = probs, 
             pheno = control$pheno[, c("zHeart_GSH", "zHeart_GSSG", "zHeart_TotalGSH")])

# look at the first several rows of the scan
head(zScan)

# LOD score significance thresholds
perms <- scan1perm(genoprobs = probs, 
                   pheno = control$pheno[, c("zHeart_GSH", "zHeart_GSSG", "zHeart_TotalGSH")], 
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

# plot scans for each phenotype
plot_scan1(zScan, map = control$gmap, 
           lodcolumn = "zHeart_GSH", main = "RankZ transformed GSH genome scan")
plot_scan1(zScan, map = control$gmap, 
           lodcolumn = "zHeart_GSSG", main = "RankZ transformed GSSG genome scan")
plot_scan1(zScan, map = control$gmap, 
           lodcolumn = "zHeart_TotalGSH", main = "RankZ transformed total GSH genome scan")



