library(devtools)
devtools::install_github("gkeele/musppr")

# The maxlod argument should be the permutation output from qtl2, so a vector of
# maximum LOD scores from the permutation scans. # For threshold estimation, the
# fwer argument is the multiple-testing corrected error rate. So for 95%
# significance, that's a 5% (0.05) error rate. -Greg Keele  
thresh <- musppr::pull_fwer_thresh(maxlod = Sexperms, fwer = 0.05)

# For calculating the genome-wide p-value at a QTL, you need the lod_obs 
# argument, which is the LOD score of your peak of interest, as you would pull 
# with find_peaks in qtl2. -Greg Keele 
qtl.lod.peak <- find_peaks(zScanSex, map = control$pmap, threshold = 7)[1, 5]
pval <- musppr::pull_fwer_pval(lod_obs = qtl.lod.peak, maxlod = Sexperms)

# Both functions allow you to calculate it either empirically (use_gev = FALSE)
# or with an extreme value distribution (use_gev = TRUE). They should be fairly
# similar. If they differ a lot, I'd probably trust the GEV more, or possibly
# run more permutations if you want to use an empirical estimate.
