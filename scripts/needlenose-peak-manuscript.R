library(qtl2)
library(purrr)

# load rankZ function
source("rankZ.R")

####################################################
## Create the probability plotting function (made by Greg Keele)
####################################################

prob_plot <- function(pheno_vec,
                      pheno_name = NULL,
                      genoprobs,
                      qtl_chr,
                      qtl_marker,
                      cols = gray(10000:1/10000),
                      label_col = as.character(qtl2::CCcolors),
                      founders = c("AJ", "B6", "129", "NOD", "NZO", "CAST", 
                                   "PWK", "WSB"),
                      main = "") { 
  sorted_pheno <- sort(pheno_vec)
  image(genoprobs[[qtl_chr]][names(sorted_pheno), rev(LETTERS[1:8]), qtl_marker] * 2, 
        yaxt = "n", xaxt = "n", col = cols)
  axis(2, at = seq(0, 8, 1 + 1/8)/8, labels = FALSE,
       lty = 0, srt = 90, las = 2)
  mtext(text = main, side = 3, padj = -1, cex = 1.25)
  mtext(text = rev(founders), side = 2, col = rev(label_col), 
        at = seq(0, 1, length.out = 8), 
        las = 1, cex = 1.25, adj = 1.25, font = 2)
  mtext(text = paste("lower", "<--", 
                     ifelse(is.null(pheno_name), 
                            "phenotype", pheno_name), "-->", "higher"), 
        side = 1, padj = 1.25, cex = 1.25)
  }

# redo figures for draft needlenose peak manuscript
liver_DO_data <- read_cross2(file = "../data/R01_GSH_DO_control.json")
probs <- readRDS("../data/Pazdro_GigaMUGA_genoprobs_qced_8state_sorted.rds")
kinship_loco <- calc_kinship(probs = probs, "loco", use_allele_probs = TRUE)

# convert sex to numeric 0 (F) or 1 (M)
liver_DO_data$covar$sex <- match(liver_DO_data$covar$sex, c("F", "M")) - 1
sex = model.matrix(~sex, data = liver_DO_data$covar)[,-1]

# rank z transform phenotypes
zLiver_GSH <- rankZ(liver_DO_data$pheno[, "Liver_GSH"])
zLiver_NADPH <- rankZ(liver_DO_data$pheno[, "Liver_NADPH"])

# put the rankZ-transformed variables back into the phenotypes
liver_DO_data$pheno <- cbind(liver_DO_data$pheno,
                       zLiver_GSH,
                       zLiver_NADPH)

save(liver_DO_data, file = "../data/transformed-data.Rdata")

zScanSex <- scan1(genoprobs = probs,
                  pheno = liver_DO_data$pheno[, c("zLiver_GSH", "zLiver_NADPH")],
                  kinship = kinship_loco,
                  addcovar = sex)

# Sexperms <- scan1perm(genoprobs = probs,
#                       pheno = liver_DO_data$pheno[, c("zLiver_GSH", "zLiver_NADPH")],
#                       addcovar = sex,
#                       n_perm = 100)
# summary(Sexperms)
# save(zScanSex, Sexperms, file = "../data/genome-scans-thresholds-sex.RData")

# panel labels A-F for chromosomes 2 and 14 figures
panel_labels <-  LETTERS[1:6]

# set y-axis limits to 10 for LOD scores
ylim <- c(0, 10)


# chromosome 2 GSH peak coefficients 
coef_c2 <- scan1coef(probs[, "2"], liver_DO_data$pheno[,"zLiver_GSH"],
                     kinship_loco$`2`, addcovar = sex)
blup_c2 <- scan1blup(probs[, "2"], liver_DO_data$pheno[,"zLiver_GSH"],
                     kinship_loco$`2`, addcovar = sex)

# figure 1 for GSH on chr 2
# if the 'plot.new has not been called yet' error appears, make sure there 
# isn't a dangling dev.off() command (close the output png graphics device) and
# run plot.new()
plot.new()
png(filename = "../results/GSH-scan.png")
plot_scan1(zScanSex, map = liver_DO_data$gmap, lodcolumn = "zLiver_GSH", 
           ylim = ylim)
mtext(panel_labels[1], cex = 2, adj = 0, line = 2)
dev.off()

png(filename = "../results/chr2-GSH-scan.png")
plot_scan1(zScanSex, map = liver_DO_data$gmap, lodcolumn = "zLiver_GSH", 
           chr = 2, ylim = ylim) 
mtext(panel_labels[2], cex = 2, adj = 0, line = 2)
dev.off()

png(filename = "../results/chr2-GSH-blup-plot.png")
xlim <- c(45, 65)
plot_coefCC(blup_c2, liver_DO_data$gmap["2"], scan1_output=zScanSex, 
            bgcolor="gray95", legend="bottomleft", xlim=xlim, yaxs = "r")
mtext(panel_labels[3], cex = 2, adj = 0, line = 2)
dev.off()

png(filename = "../results/chr2-GSH-coef-plot.png")
plot_coefCC(coef_c2, liver_DO_data$gmap["2"], scan1_output=zScanSex, 
            bgcolor="gray95", legend="bottomleft", yaxs = "r")
mtext(panel_labels[4], cex = 2, adj = 0, line = 2)
dev.off()



# chromosome 14 NADPH peak coefficients 
coef_c14 <- scan1coef(probs[, "14"], 
                      liver_DO_data$pheno[,"zLiver_NADPH"],
                      kinship_loco$`14`, 
                      addcovar = sex)
blup_c14 <- scan1blup(probs[, "14"], 
                      liver_DO_data$pheno[,"zLiver_NADPH"],
                      kinship_loco$`14`, 
                      addcovar = sex)

# chromosome 14 NADPH peak for figure 2
png(filename = "../results/NADPH-scan.png")
plot_scan1(zScanSex, map = liver_DO_data$gmap, lodcolumn = "zLiver_NADPH",
                               ylim = ylim)
mtext(panel_labels[1], cex = 2, adj = 0, line = 2)
dev.off()

png(filename = "../results/chr14-NADPH-scan.png")
plot_scan1(zScanSex, map = liver_DO_data$gmap,  lodcolumn = "zLiver_NADPH",
                               chr = 14, ylim = ylim)
mtext(panel_labels[2], cex = 2, adj = 0, line = 2)
dev.off()

png(filename = "../results/chr14-NADPH-blup-plot.png")
xlim <- c(35, 50)
plot_coefCC(blup_c14, liver_DO_data$gmap["14"], scan1_output=zScanSex, 
            bgcolor="gray95", legend="bottomleft", xlim=xlim)
mtext(panel_labels[3], cex = 2, adj = 0, line = 2)
dev.off()

png(filename = "../results/chr14-NADPH-coef-plot.png")
plot_coefCC(coef_c14, liver_DO_data$gmap["14"], scan1_output=zScanSex, 
            bgcolor="gray95", legend="bottomleft")
mtext(panel_labels[4], cex = 2, adj = 0, line = 2)
dev.off()

# prob plots for peak markers on chr 2 and 14
## Set threshold to 7 and use peakdrop to bring in secondary peak on chr2
## use both genetic and physical maps
qtl_gmap <- find_peaks(zScanSex, map = liver_DO_data$gmap, threshold = 7, 
                       peakdrop = 6.5)
qtl_pmap <- find_peaks(zScanSex, map = liver_DO_data$pmap, threshold = 7, 
                       peakdrop = 6.5)
# Add marker information
qtl_gmap$marker.id <- find_marker(map = liver_DO_data$gmap, chr = qtl_gmap$chr, 
                                  pos = qtl_gmap$pos)
qtl_pmap$marker.id <- find_marker(map = liver_DO_data$pmap, chr = qtl_pmap$chr, 
                                  pos = qtl_pmap$pos)

for (i in 1:dim(qtl_gmap)[1]) {
  png(filename = paste("../results/", qtl_gmap$chr[i], qtl_gmap$marker.id[i], 
                       "-probplot.png", sep = ""))
  prob_plot(pheno_vec = liver_DO_data$pheno[, qtl_gmap$lodcolumn[i]],
            genoprobs = probs,
            qtl_chr = qtl_gmap$chr[i],
            qtl_marker = qtl_gmap$marker.id[i],
            main = paste(qtl_gmap$marker.id[i], "\n", qtl_gmap$pos[i],
                         mtext(panel_labels[i+4], cex = 2, adj = 0, line = 2)
            )
            )
  dev.off()
  }
 