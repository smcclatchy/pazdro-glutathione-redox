#estimating QTL effects and connecting to SNP and Gene Datbases
#Heart GSH QTL 2022

library(qtl2)

# load genome scans and thresholds from previous script
# sex is the only covariate
load("../data/genome-scans-thresholds-sex.RData")
Sexthresholds
plot_scan1(zScanSex, map = control$gmap,
           lodcolumn = "zHeart_GSH", main = "GSH",
           sub = "with sex as covariate")
abline(h = Sexthresholds[1], col = "red", lwd = 2)
find_peaks(scan1_output = zScanSex,
           map = control$pmap,
           threshold = Sexthresholds,
           prob = 0.95)

# create database query functions
query_variants <- create_variant_query_func("../data/cc_variants.sqlite")
query_genes <- create_gene_query_func("../data/mouse_genes_mgi.sqlite")

####################################################
## Estimate QTL Effects (Coefficients) 
####################################################

#For Heart GSH --- Chromosome 14
par(mar=c(4.1, 4.1, 2.6, 2.6))

# estimate QTL effects by founder strain
# using gmap (cM)
chr = 14
coef_blup_HeartGSH_chr14 <- scan1blup(genoprobs =  probs[,chr], 
                                      pheno = control$pheno["zHeart_GSH"], 
                                      kinship = kinship_loco[[chr]], 
                                      addcovar = sex)

plot_coefCC(x = coef_blup_HeartGSH_chr14, map = control$pmap,
            scan1_output = zScanSex,
            main = "Chromosome 14 Heart GSH BLUPs", 
            sub="with sex as covariate",
            legend = "topright", legend_ncol = 1, bgcolor="gray95")

# repeat for chromosome 16 & 19 peaks
# chr 16
chr = 16
coef_blup_HeartGSH_chr16 <- scan1blup(genoprobs =  probs[,chr], 
                                      pheno = control$pheno["zHeart_GSH"], 
                                      kinship = kinship_loco[[chr]], 
                                      addcovar = sex)

plot_coefCC(x = coef_blup_HeartGSH_chr16, map = control$pmap,
            scan1_output = zScanSex,
            main = "Chromosome 16 Heart GSH BLUPs", 
            legend = "bottom", legend_ncol = 4,
            sub="with sex as covariate", bgcolor="gray95")

# chr 19
chr = 19
coef_blup_HeartGSH_chr19 <- scan1blup(genoprobs =  probs[,chr], 
                                      pheno = control$pheno["zHeart_GSH"], 
                                      kinship = kinship_loco[[chr]], 
                                      addcovar = sex)

plot_coefCC(x = coef_blup_HeartGSH_chr19, map = control$pmap,
            scan1_output = zScanSex,
            main = "Chromosome 19 Heart GSH BLUPs", 
            legend = "bottom", legend_ncol = 4,
            sub="with sex as covariate", bgcolor="gray95")

#could use ci_lo or ci_hi, but in this case, I want a specific chromosome 14 peak
#start = pmap_peaksHeartGSH[LODzsexp$chr ==  chr,"ci_lo"]
#end = pmap_peaksHeartGSH[LODzsexp$chr == chr, "ci_hi"] 

pander(LODzsexp)
#based on pmap_peaksHeartGSH, peak of interest is ~52.79 Mbp
#Becca typically does +/- of the QTL interval
variants_HeartGSH_chr14 <- query_variants(chr, 45, 60)
out_snps_HeartGSH_chr14 <- scan1snps(genoprobs = probs, map = control$pmap, 
                                     pheno = control$pheno["zHeart_GSH"], 
                                     kinship = kinship_loco[[chr]], 
                                     addcovar = sexgen, 
                                     query_func = query_variants,
                                     chr = chr, 
                                     keep_all_snps = TRUE)
plot_snpasso(out_snps_HeartGSH_chr14$lod, out_snps_HeartGSH_chr14$snpinfo, 
             main = "Heart GSH SNPs")

HeartGSH_Genes_MGI_chr14 <- query_genes(chr = chr, start = 45, end = 60)
plot(out_snps_HeartGSH_chr14$lod, out_snps_HeartGSH_chr14$snpinfo, drop_hilit=1.5, genes = HeartGSH_Genes_MGI_chr14, main = "Heart GSH Genes MGI")

#For Heart GSH --- Chromosome 16
par(mar=c(4.1, 4.1, 2.6, 2.6))

#estimate QTL effects by founder strain
#using gmap (cM)
chr = 16
coef_blup_HeartGSH_chr16 <- scan1blup(genoprobs =  probs[,chr], 
                                      pheno = control$pheno["zHeart_GSH"], 
                                      kinship = kinship_loco[[chr]], 
                                      addcovar = sexgen)
plot_coefCC(x = coef_blup_HeartGSH_chr16, map = control$gmap, 
            scan1_output = zScanGen, 
            main = "Heart GSH BLUPs plotted with CC Founders - generation", 
            legend = "bottomleft", bgcolor="gray95")

plot_coefCC(x = coef_blup_HeartGSH_chr16, map = control$pmap, 
            scan1_output = zScanSex, 
            main = "Heart GSH BLUPs plotted with CC Founders - sex", 
            legend = "bottomleft", bgcolor="gray95",
            xlim = c(90, max(control$pmap$`16`)))

xlim <- c(30,40)
plot_coefCC(x = coef_blup_HeartGSH_chr16, map = control$gmap, scan1_output = zScanGen, main = "Heart GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)

#using pmap (Mbp)
chr = 16
#could use ci_lo or ci_hi, but in this case, I want a specific chromosome 16 peak
#start = pmap_peaksHeartGSH[LODzsexp$chr ==  chr,"ci_lo"]
#end = pmap_peaksHeartGSH[LODzsexp$chr == chr, "ci_hi"] 

pander(LODzsexp)
#based on pmap_peaksHeartGSH, peak of interest is ~96.99176 Mbp
#Becca typically does +/- of the QTL interval
variants_HeartGSH_chr16 <- query_variants(chr, 90, max(control$pmap$`16`))
out_snps_HeartGSH_chr16 <- scan1snps(genoprobs = probs, map = control$pmap, 
                                     pheno = control$pheno["zHeart_GSH"], 
                                     kinship = kinship_loco[[chr]], 
                                     addcovar = sexgen, 
                                     query_func = query_variants,
                                     chr = chr, 
                                     start = 90, end = max(control$pmap$`16`), 
                                     keep_all_snps = TRUE)
plot_snpasso(out_snps_HeartGSH_chr16$lod, 
             out_snps_HeartGSH_chr16$snpinfo, main = "Heart GSH SNPs")

top <- top_snps(out_snps_HeartGSH_chr16$lod, out_snps_HeartGSH_chr16$snpinfo)
print(top, row.names=FALSE)

HeartGSH_Genes_MGI_chr16 <- query_genes(chr = chr, start = 95, end = 99)
plot(out_snps_HeartGSH_chr16$lod, out_snps_HeartGSH_chr16$snpinfo, drop_hilit=1.5, genes = HeartGSH_Genes_MGI_chr16, main = "Heart GSH Genes MGI")


#For Heart GSH --- Chromosome 19
par(mar=c(4.1, 4.1, 2.6, 2.6))

#estimate QTL effects by founder strain
#using gmap (cM)
chr = 19
coef_blup_HeartGSH_chr19 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zHeart_GSH"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 10)
plot_coefCC(x = coef_blup_HeartGSH_chr19, map = control$gmap, scan1_output = zScanGen, main = "Heart GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
xlim <- c(50,59)
plot_coefCC(x = coef_blup_HeartGSH_chr19, map = control$gmap, scan1_output = zScanGen, main = "Heart GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)

#using pmap (Mbp)
chr = 19
#could use ci_lo or ci_hi, but in this case, I want a specific chromosome 19 peak
#start = pmap_peaksHeartGSH[pmap_peaksHeartGSH$chr ==  chr,"ci_lo"]
#end = pmap_peaksHeartGSH[pmap_peaksHeartGSH$chr == chr, "ci_hi"] 

pander(pmap_peaksHeartGSH)
#based on pmap_peaksHeartGSH, peak of interest is ~57.231541 Mbp
#Becca typically does +/- of the QTL interval
variants_HeartGSH_chr19 <- query_variants(chr, 56, 58.5)
out_snps_HeartGSH_chr19 <- scan1snps(genoprobs = probs, map = control$pmap, pheno = pheno["zHeart_GSH"], kinship = kinship_loco[[chr]], addcovar = sexgen, query_func = query_variants,
                                     chr = chr, start = 56, end = 58.5, keep_all_snps = TRUE)
plot_snpasso(out_snps_HeartGSH_chr19$lod, out_snps_HeartGSH_chr19$snpinfo, main = "Heart GSH SNPs")

HeartGSH_Genes_MGI_chr19 <- query_genes(chr = chr, start = 56, end = 58.5)
plot(out_snps_HeartGSH_chr19$lod, out_snps_HeartGSH_chr19$snpinfo, drop_hilit=1.5, genes = HeartGSH_Genes_MGI_chr19, main = "Heart GSH Genes MGI")

#For Heart GSSG --- Chromosome 10
par(mar=c(4.1, 4.1, 2.6, 2.6))

#estimate QTL effects by founder strain
#using gmap (cM)
chr = 10
coef_blup_HeartGSSG_chr10 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zHeart_GSSG"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 10)
plot_coefCC(x = coef_blup_HeartGSSG_chr10, map = control$gmap, scan1_output = zScanGen, main = "Heart GSSG BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
xlim <- c(57,70)
plot_coefCC(x = coef_blup_HeartGSSG_chr10, map = control$gmap, scan1_output = zScanGen, main = "Heart GSSG BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)

#using pmap (Mbp)
chr = 10
#could use ci_lo or ci_hi, but in this case, I want a specific chromosome 10 peak
#start = pmap_peaksHeartGSH[pmap_peaksHeartGSH$chr ==  chr,"ci_lo"]
#end = pmap_peaksHeartGSH[pmap_peaksHeartGSH$chr == chr, "ci_hi"] 

pander(pmap_peaksHeartGSH)
#based on pmap_peaksHeartGSH, peak of interest is ~124.103853 Mbp
#Becca typically does +/- of the QTL interval
variants_HeartGSSG_chr10 <- query_variants(chr, 123, 126)
out_snps_HeartGSSG_chr10 <- scan1snps(genoprobs = probs, map = control$pmap, pheno = pheno["zHeart_GSSG"], kinship = kinship_loco[[chr]], addcovar = sexgen, query_func = query_variants,
                                     chr = chr, start = 123, end = 126, keep_all_snps = TRUE)
plot_snpasso(out_snps_HeartGSSG_chr10$lod, out_snps_HeartGSSG_chr10$snpinfo, main = "Heart GSSG SNPs")

HeartGSSG_Genes_MGI_chr10 <- query_genes(chr = chr, start = 123, end = 126)
plot(out_snps_HeartGSSG_chr10$lod, out_snps_HeartGSSG_chr10$snpinfo, drop_hilit=1.5, genes = HeartGSSG_Genes_MGI_chr10, main = "Heart GSSG Genes MGI")

###For Heart GSHGSSG_Ratio --- Chromosome 2
par(mar=c(4.1, 4.1, 2.6, 2.6))

#estimate QTL effects by founder strain
#using gmap (cM)
chr = 2
coef_blup_HeartGSHGSSG_Ratio_chr2 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zGSHGSSG_Ratio"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
plot_coefCC(x = coef_blup_HeartGSHGSSG_Ratio_chr2, map = control$gmap, scan1_output = zScanGen, main = "Heart GSH/GSSG Ratio BLUPs plotted with CC Founders", legend = "topright", bgcolor="gray95")
xlim <- c(0,20)
plot_coefCC(x = coef_blup_HeartGSHGSSG_Ratio_chr2, map = control$gmap, scan1_output = zScanGen, main = "Heart GSH/GSSG Ratio BLUPs plotted with CC Founders", legend = "topright", bgcolor="gray95", xlim = xlim)

#using pmap (Mbp)
chr = 2
#could use ci_lo or ci_hi, but in this case, I want a specific chromosome 2 peak
#start = pmap_peaksHeartGSH[pmap_peaksHeartGSH$chr ==  chr,"ci_lo"]
#end = pmap_peaksHeartGSH[pmap_peaksHeartGSH$chr == chr, "ci_hi"] 

pander(pmap_peaksHeartGSH)
#based on pmap_peaksHeartGSH, peak of interest is ~34.941792 Mbp
#Becca typically does +/- of the QTL interval
variants_HeartGSHGSSG_Ratio_chr2 <- query_variants(chr, 33.5, 35.5)
out_snps_HeartGSHGSSG_Ratio_chr2 <- scan1snps(genoprobs = probs, map = control$pmap, pheno = pheno["zGSHGSSG_Ratio"], kinship = kinship_loco[[chr]], addcovar = sexgen, query_func = query_variants,
                                      chr = chr, start = 33.5, end = 35.5, keep_all_snps = TRUE)
plot_snpasso(out_snps_HeartGSHGSSG_Ratio_chr2$lod, out_snps_HeartGSHGSSG_Ratio_chr2$snpinfo, main = "Heart GSHGSSG Ratio SNPs")

HeartGSHGSSG_Ratio_Genes_MGI_chr2 <- query_genes(chr = chr, start = 33.5, end = 35.5)
plot(out_snps_HeartGSHGSSG_Ratio_chr2$lod, out_snps_HeartGSHGSSG_Ratio_chr2$snpinfo, drop_hilit=1.5, genes = HeartGSHGSSG_Ratio_Genes_MGI_chr2, main = "Heart GSHGSSG Ratio Genes MGI")

###For Heart GSHGSSG_Ratio --- Chromosome 14
par(mar=c(4.1, 4.1, 2.6, 2.6))

#estimate QTL effects by founder strain
#using gmap (cM)
chr = 14
coef_blup_HeartGSHGSSG_Ratio_chr14 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zGSHGSSG_Ratio"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 14)
plot_coefCC(x = coef_blup_HeartGSHGSSG_Ratio_chr14, map = control$gmap, scan1_output = zScanGen, main = "Heart GSH/GSSG Ratio BLUPs plotted with CC Founders", legend = "topright", bgcolor="gray95")
xlim <- c(22,30)
plot_coefCC(x = coef_blup_HeartGSHGSSG_Ratio_chr14, map = control$gmap, scan1_output = zScanGen, main = "Heart GSH/GSSG Ratio BLUPs plotted with CC Founders", legend = "topright", bgcolor="gray95", xlim = xlim)

#using pmap (Mbp)
chr = 14
#could use ci_lo or ci_hi, but in this case, I want a specific chromosome 14 peak
#start = pmap_peaksHeartGSH[pmap_peaksHeartGSH$chr ==  chr,"ci_lo"]
#end = pmap_peaksHeartGSH[pmap_peaksHeartGSH$chr == chr, "ci_hi"] 

pander(pmap_peaksHeartGSH)
#based on pmap_peaksHeartGSH, peak of interest is ~54.143953 Mbp
#Becca typically does +/- of the QTL interval
variants_HeartGSHGSSG_Ratio_chr14 <- query_variants(chr, 53, 55)
out_snps_HeartGSHGSSG_Ratio_chr14 <- scan1snps(genoprobs = probs, map = control$pmap, pheno = pheno["zGSHGSSG_Ratio"], kinship = kinship_loco[[chr]], addcovar = sexgen, query_func = query_variants,
                                              chr = chr, start = 53, end = 55, keep_all_snps = TRUE)
plot_snpasso(out_snps_HeartGSHGSSG_Ratio_chr14$lod, out_snps_HeartGSHGSSG_Ratio_chr14$snpinfo, main = "Heart GSHGSSG Ratio SNPs")

HeartGSHGSSG_Ratio_Genes_MGI_chr14 <- query_genes(chr = chr, start = 53, end = 55)
plot(out_snps_HeartGSHGSSG_Ratio_chr14$lod, out_snps_HeartGSHGSSG_Ratio_chr14$snpinfo, drop_hilit=1.5, genes = HeartGSHGSSG_Ratio_Genes_MGI_chr14, main = "Heart GSHGSSG Ratio Genes MGI")

###For Heart GSHGSSG_Ratio --- Chromosome 16
par(mar=c(4.1, 4.1, 2.6, 2.6))

#estimate QTL effects by founder strain
#using gmap (cM)
chr = 16
coef_blup_HeartGSHGSSG_Ratio_chr16 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zGSHGSSG_Ratio"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 16)
plot_coefCC(x = coef_blup_HeartGSHGSSG_Ratio_chr16, map = control$gmap, scan1_output = zScanGen, main = "Heart GSH/GSSG Ratio BLUPs plotted with CC Founders", legend = "topright", bgcolor="gray95")
xlim <- c(30,40)
plot_coefCC(x = coef_blup_HeartGSHGSSG_Ratio_chr16, map = control$gmap, scan1_output = zScanGen, main = "Heart GSH/GSSG Ratio BLUPs plotted with CC Founders", legend = "topright", bgcolor="gray95", xlim = xlim)

#using pmap (Mbp)
chr = 16
#could use ci_lo or ci_hi, but in this case, I want a specific chromosome 16 peak
#start = pmap_peaksHeartGSH[pmap_peaksHeartGSH$chr ==  chr,"ci_lo"]
#end = pmap_peaksHeartGSH[pmap_peaksHeartGSH$chr == chr, "ci_hi"] 

pander(pmap_peaksHeartGSH)
#based on pmap_peaksHeartGSH, peak of interest is ~96.735776 Mbp
#Becca typically does +/- of the QTL interval
variants_HeartGSHGSSG_Ratio_chr16 <- query_variants(chr, 95, 97.5)
out_snps_HeartGSHGSSG_Ratio_chr16 <- scan1snps(genoprobs = probs, map = control$pmap, pheno = pheno["zGSHGSSG_Ratio"], kinship = kinship_loco[[chr]], addcovar = sexgen, query_func = query_variants,
                                               chr = chr, start = 95, end = 97.5, keep_all_snps = TRUE)
plot_snpasso(out_snps_HeartGSHGSSG_Ratio_chr16$lod, out_snps_HeartGSHGSSG_Ratio_chr16$snpinfo, main = "Heart GSHGSSG Ratio SNPs")

HeartGSHGSSG_Ratio_Genes_MGI_chr16 <- query_genes(chr = chr, start = 95, end = 97.5)
plot(out_snps_HeartGSHGSSG_Ratio_chr16$lod, out_snps_HeartGSHGSSG_Ratio_chr16$snpinfo, drop_hilit=1.5, genes = HeartGSHGSSG_Ratio_Genes_MGI_chr16, main = "Heart GSHGSSG Ratio Genes MGI")

###For Heart GSHGSSG_Ratio --- Chromosome 19
par(mar=c(4.1, 4.1, 2.6, 2.6))

#estimate QTL effects by founder strain
#using gmap (cM)
chr = 19
coef_blup_HeartGSHGSSG_Ratio_chr19 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zGSHGSSG_Ratio"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 10)
plot_coefCC(x = coef_blup_HeartGSHGSSG_Ratio_chr19, map = control$gmap, scan1_output = zScanGen, main = "Heart GSH/GSSG Ratio BLUPs plotted with CC Founders", legend = "topright", bgcolor="gray95")
xlim <- c(50,57)
plot_coefCC(x = coef_blup_HeartGSHGSSG_Ratio_chr19, map = control$gmap, scan1_output = zScanGen, main = "Heart GSH/GSSG Ratio BLUPs plotted with CC Founders", legend = "topright", bgcolor="gray95", xlim = xlim)

#using pmap (Mbp)
chr = 19
#could use ci_lo or ci_hi, but in this case, I want a specific chromosome 19 peak
#start = pmap_peaksHeartGSH[pmap_peaksHeartGSH$chr ==  chr,"ci_lo"]
#end = pmap_peaksHeartGSH[pmap_peaksHeartGSH$chr == chr, "ci_hi"] 

pander(pmap_peaksHeartGSH)
#based on pmap_peaksHeartGSH, peak of interest is ~57.511685 Mbp
#Becca typically does +/- of the QTL interval
variants_HeartGSHGSSG_Ratio_chr19 <- query_variants(chr, 56, 58)
out_snps_HeartGSHGSSG_Ratio_chr19 <- scan1snps(genoprobs = probs, map = control$pmap, pheno = pheno["zGSHGSSG_Ratio"], kinship = kinship_loco[[chr]], addcovar = sexgen, query_func = query_variants,
                                               chr = chr, start = 56, end = 58, keep_all_snps = TRUE)
plot_snpasso(out_snps_HeartGSHGSSG_Ratio_chr19$lod, out_snps_HeartGSHGSSG_Ratio_chr19$snpinfo, main = "Heart GSHGSSG Ratio SNPs")

HeartGSHGSSG_Ratio_Genes_MGI_chr19 <- query_genes(chr = chr, start = 56, end = 58)
plot(out_snps_HeartGSHGSSG_Ratio_chr19$lod, out_snps_HeartGSHGSSG_Ratio_chr19$snpinfo, drop_hilit=1.5, genes = HeartGSHGSSG_Ratio_Genes_MGI_chr19, main = "Heart GSHGSSG Ratio Genes MGI")

###For Heart heart_redoxPot --- Chromosome 2
par(mar=c(4.1, 4.1, 2.6, 2.6))

#estimate QTL effects by founder strain
#using gmap (cM)
chr = 2
coef_blup_heart_redoxPot_chr2 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zheart_redoxPot"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
plot_coefCC(x = coef_blup_heart_redoxPot_chr2, map = control$gmap, scan1_output = zScanGen, main = "Heart Redox Potential BLUPs plotted with CC Founders", legend = "bottomright", bgcolor="gray95")
xlim <- c(2,5)
plot_coefCC(x = coef_blup_heart_redoxPot_chr2, map = control$gmap, scan1_output = zScanGen, main = "Heart Redox Potential BLUPs plotted with CC Founders", legend = "bottomright", bgcolor="gray95", xlim = xlim)

#using pmap (Mbp)
chr = 2
#could use ci_lo or ci_hi, but in this case, I want a specific chromosome 2 peak
#start = pmap_peaksHeartGSH[pmap_peaksHeartGSH$chr ==  chr,"ci_lo"]
#end = pmap_peaksHeartGSH[pmap_peaksHeartGSH$chr == chr, "ci_hi"] 

pander(pmap_peaksHeartGSH)
#based on pmap_peaksHeartGSH, peak of interest is ~34.941792 Mbp
#Becca typically does +/- of the QTL interval
variants_heart_redoxPot_chr2 <- query_variants(chr, 34, 35.5)
out_snps_heart_redoxPot_chr2 <- scan1snps(genoprobs = probs, map = control$pmap, pheno = pheno["zheart_redoxPot"], kinship = kinship_loco[[chr]], addcovar = sexgen, query_func = query_variants,
                                               chr = chr, start = 34, end = 35.5, keep_all_snps = TRUE)
plot_snpasso(out_snps_heart_redoxPot_chr2$lod, out_snps_heart_redoxPot_chr2$snpinfo, main = "Heart Redox Potential SNPs")

heart_redoxPot_Genes_MGI_chr2 <- query_genes(chr = chr, start = 34, end = 35.5)
plot(out_snps_heart_redoxPot_chr2$lod, out_snps_heart_redoxPot_chr2$snpinfo, drop_hilit=1.5, genes = heart_redoxPot_Genes_MGI_chr2, main = "Heart Redox Potential Genes MGI")

###For Heart heart_redoxPot --- Chromosome 5
par(mar=c(4.1, 4.1, 2.6, 2.6))

#estimate QTL effects by founder strain
#using gmap (cM)
chr = 5
coef_blup_heart_redoxPot_chr5 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zheart_redoxPot"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 5)
plot_coefCC(x = coef_blup_heart_redoxPot_chr5, map = control$gmap, scan1_output = zScanGen, main = "Heart Redox Potential BLUPs plotted with CC Founders", legend = "topright", bgcolor="gray95")
xlim <- c(50,58)
plot_coefCC(x = coef_blup_heart_redoxPot_chr5, map = control$gmap, scan1_output = zScanGen, main = "Heart Redox Potential BLUPs plotted with CC Founders", legend = "topright", bgcolor="gray95", xlim = xlim)

#using pmap (Mbp)
chr = 5
#could use ci_lo or ci_hi, but in this case, I want a specific chromosome 5 peak
#start = pmap_peaksHeartGSH[pmap_peaksHeartGSH$chr ==  chr,"ci_lo"]
#end = pmap_peaksHeartGSH[pmap_peaksHeartGSH$chr == chr, "ci_hi"] 

pander(pmap_peaksHeartGSH)
#based on pmap_peaksHeartGSH, peak of interest is ~121.65796 Mbp
#Becca typically does +/- of the QTL interval
variants_heart_redoxPot_chr5 <- query_variants(chr, 120.5, 122.5)
out_snps_heart_redoxPot_chr5 <- scan1snps(genoprobs = probs, map = control$pmap, pheno = pheno["zheart_redoxPot"], kinship = kinship_loco[[chr]], addcovar = sexgen, query_func = query_variants,
                                          chr = chr, start = 120.5, end = 122.5, keep_all_snps = TRUE)
plot_snpasso(out_snps_heart_redoxPot_chr5$lod, out_snps_heart_redoxPot_chr5$snpinfo, main = "Heart Redox Potential SNPs")

heart_redoxPot_Genes_MGI_chr5 <- query_genes(chr = chr, start = 120.5, end = 122.5)
plot(out_snps_heart_redoxPot_chr5$lod, out_snps_heart_redoxPot_chr5$snpinfo, drop_hilit=1.5, genes = heart_redoxPot_Genes_MGI_chr5, main = "Heart Redox Potential Genes MGI")

###For Heart heart_redoxPot --- Chromosome 16
par(mar=c(4.1, 4.1, 2.6, 2.6))

#estimate QTL effects by founder strain
#using gmap (cM)
chr = 16
coef_blup_heart_redoxPot_chr16 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zheart_redoxPot"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 10)
plot_coefCC(x = coef_blup_heart_redoxPot_chr16, map = control$gmap, scan1_output = zScanGen, main = "Heart Redox Potential BLUPs plotted with CC Founders", legend = "topleft", bgcolor="gray95")
xlim <- c(30,40)
plot_coefCC(x = coef_blup_heart_redoxPot_chr16, map = control$gmap, scan1_output = zScanGen, main = "Heart Redox Potential BLUPs plotted with CC Founders", legend = "topright", bgcolor="gray95", xlim = xlim)

#using pmap (Mbp)
chr = 16
#could use ci_lo or ci_hi, but in this case, I want a specific chromosome 16 peak
#start = pmap_peaksHeartGSH[pmap_peaksHeartGSH$chr ==  chr,"ci_lo"]
#end = pmap_peaksHeartGSH[pmap_peaksHeartGSH$chr == chr, "ci_hi"] 

pander(pmap_peaksHeartGSH)
#based on pmap_peaksHeartGSH, peak of interest is ~96.735776 Mbp
#Becca typically does +/- of the QTL interval
variants_heart_redoxPot_chr16 <- query_variants(chr, 95.5, 97.5)
out_snps_heart_redoxPot_chr16 <- scan1snps(genoprobs = probs, map = control$pmap, pheno = pheno["zheart_redoxPot"], kinship = kinship_loco[[chr]], addcovar = sexgen, query_func = query_variants,
                                          chr = chr, start = 95.5, end = 97.5, keep_all_snps = TRUE)
plot_snpasso(out_snps_heart_redoxPot_chr16$lod, out_snps_heart_redoxPot_chr16$snpinfo, main = "Heart Redox Potential SNPs")

heart_redoxPot_Genes_MGI_chr16 <- query_genes(chr = chr, start = 95.5, end = 97.5)
plot(out_snps_heart_redoxPot_chr16$lod, out_snps_heart_redoxPot_chr16$snpinfo, drop_hilit=1.16, genes = heart_redoxPot_Genes_MGI_chr16, main = "Heart Redox Potential Genes MGI")


####################################################
## GWAS SNP Association Scan
## Make a Manhattan plot of the results; use altcol to define a color alternate for chromosomes and gap=0 to have no gap between chromosomes
####################################################

out_gwas_HeartGSH <- scan1snps(genoprobs = probs, map = control$pmap, pheno = pheno["zHeart_GSH"], kinship = kinship_loco, addcovar = sexgen, query_func=query_variants, cores=10)
par(mar=c(4.1, 4.1, 2.6, 2.6))
plot(out_gwas_HeartGSH$lod, out_gwas_HeartGSH$snpinfo, altcol="green4", gap=0, main = "Heart GSH GWAS", ylim = c(0,6))

out_gwas_HeartGSSG <- scan1snps(genoprobs = probs, map = control$pmap, pheno = pheno["zHeart_GSSG"], kinship = kinship_loco, addcovar = sexgen, query_func=query_variants, cores=10)
par(mar=c(4.1, 4.1, 2.6, 2.6))
plot(out_gwas_HeartGSSG$lod, out_gwas_HeartGSSG$snpinfo, altcol="green4", gap=0, main = "Heart GSSG GWAS", ylim = c(0,6))

out_gwas_HeartGSHGSSGratio <- scan1snps(genoprobs = probs, map = control$pmap, pheno = pheno["zGSHGSSG_Ratio"], kinship = kinship_loco, addcovar = sexgen, query_func=query_variants, cores=10)
par(mar=c(4.1, 4.1, 2.6, 2.6))
plot(out_gwas_HeartGSHGSSGratio$lod, out_gwas_HeartGSHGSSGratio$snpinfo, altcol="green4", gap=0, main = "Heart GSH/GSSG Ratio GWAS", ylim = c(0,6))

out_gwas_HeartTotalGSH <- scan1snps(genoprobs = probs, map = control$pmap, pheno = pheno["zHeart_TotalGSH"], kinship = kinship_loco, addcovar = sexgen, query_func=query_variants, cores=10)
par(mar=c(4.1, 4.1, 2.6, 2.6))
plot(out_gwas_HeartTotalGSH$lod, out_gwas_HeartTotalGSH$snpinfo, altcol="green4", gap=0, main = "Heart Total GSH GWAS", ylim = c(0,6))

out_gwas_HeartRedoxPot <- scan1snps(genoprobs = probs, map = control$pmap, pheno = pheno["zheart_redoxPot"], kinship = kinship_loco, addcovar = sexgen, query_func=query_variants, cores=10)
par(mar=c(4.1, 4.1, 2.6, 2.6))
plot(out_gwas_HeartRedoxPot$lod, out_gwas_HeartRedoxPot$snpinfo, altcol="green4", gap=0, main = "Heart Redox Potential GWAS", ylim = c(0,6))

##################################################################
## Checking other glutathione genes 
##################################################################

#set working directory
pdf(file = "GSH Other Genes - QTL Results - RankZ sexgen.pdf")
##NOW SAVING ALL PLOTS AND TABLES ONTO A PDF##

#Glutathione Peroxidase (Gpx1) - Chr 9 59.24 cM
chr = 9
coef_blup_HeartGSH_chr9 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zHeartGSH"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
plot_coefCC(x = coef_blup_HeartGSH_chr9, map = control$gmap, scan1_output = qtlscan_HeartGSH, main = "Heart GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
xlim <- c(58.5,60)
plot_coefCC(x = coef_blup_HeartGSH_chr9, map = control$gmap, scan1_output = qtlscan_HeartGSH, main = "Gpx1 Position -- Heart GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)

#Glutathione cysteine ligase - catalytic subunit (Gclc) - Chr 9 43.36 cM
chr = 9
#coef_blup_HeartGSH_chr9 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zHeartGSH"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
#plot_coefCC(x = coef_blup_HeartGSH_chr9, map = control$gmap, scan1_output = qtlscan_HeartGSH, main = "Heart GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
xlim <- c(42.5,44)
plot_coefCC(x = coef_blup_HeartGSH_chr9, map = control$gmap, scan1_output = qtlscan_HeartGSH, main = "Gclc Position -- Heart GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)

#Glutamate-cysteine ligase – modifier subunit (Gclm) - Chr 3 52.94 cM
chr = 3
coef_blup_HeartGSH_chr3 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zHeartGSH"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
plot_coefCC(x = coef_blup_HeartGSH_chr3, map = control$gmap, scan1_output = qtlscan_HeartGSH, main = "Heart GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
xlim <- c(51.5,53.5)
plot_coefCC(x = coef_blup_HeartGSH_chr3, map = control$gmap, scan1_output = qtlscan_HeartGSH, main = "Gclm Position -- Heart GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)

#Glutathione synthetase (Gss) - Chr 2 77.26 cM
chr = 2
coef_blup_HeartGSH_chr2 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zHeartGSH"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
plot_coefCC(x = coef_blup_HeartGSH_chr2, map = control$gmap, scan1_output = qtlscan_HeartGSH, main = "Heart GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
xlim <- c(76.5,78)
plot_coefCC(x = coef_blup_HeartGSH_chr2, map = control$gmap, scan1_output = qtlscan_HeartGSH, main = "Gss Position -- Heart GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)


#Glutathione reductase  (Gsr) - Chr 8 20.69 cM
chr = 8
coef_blup_HeartGSH_chr8 <- scan1blup(genoprobs =  probs[,chr], pheno = pheno["zHeartGSH"], kinship = kinship_loco[[chr]], addcovar = sexgen, cores = 2)
plot_coefCC(x = coef_blup_HeartGSH_chr8, map = control$gmap, scan1_output = qtlscan_HeartGSH, main = "Heart GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95")
xlim <- c(10,22.5)
plot_coefCC(x = coef_blup_HeartGSH_chr8, map = control$gmap, scan1_output = qtlscan_HeartGSH, main = "Gsr Position -- Heart GSH BLUPs plotted with CC Founders", legend = "bottomleft", bgcolor="gray95", xlim = xlim)

dev.off()