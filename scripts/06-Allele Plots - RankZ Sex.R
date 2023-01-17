# - Allele Plots

#Creates Heat Maps that shows the individual CC founder allele contribution at a specific SNP from the Giga MUGA


#load the command line tools 
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
library (writexl)
library (pander)


####################################################
## Finding and sorting through all of the QTL data
####################################################

#tells you all of the qtlscans that you have
  ls(pattern = "zScan")


#use cbind to combine all of the qtlscans + take those 2D tables and combining them with another table over and over again
## scans is an R object containing your genome scans from scan1() that are loaded in to the R environment 
#scans <- cbind(qtlscan_LiverGSH, qtlscan_LiverGSSG, qtlscan_LiverTotalGSH, qtlscan_LiverGSH_GSSGRatio, qtlscan_LiverGSH_GSSGcovar, qtlscan_LiverNADH, qtlscan_LiverNADP, qtlscan_LiverNADPH, qtlscan_LiverNADP_NADPHRatio, qtlscan_ALT, qtlscan_AST, qtlscan_BUN)

#thresholds <- cbind(threshold_LiverGSH, threshold_LiverGSSG, threshold_LiverTotalGSH, threshold_LiverGSH_GSSGRatio, threshold_LiverGSH_GSSGcovar, threshold_LiverNADH, threshold_LiverNADP, threshold_LiverNADPH, threshold_LiverNADP_NADPHRatio, threshold_ALT, threshold_AST, threshold_BUN)
#head(thresholds)

#scans variable created from "Exporting QTL Results" -- for RankZ transformed pheno 
#head(scans)
  

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
                      founders = c("AJ", "B6", "129", "NOD", "NZO", "CAST", "PWK", "WSB"),
                      main = "") {
  sorted_pheno <- sort(pheno_vec)
  image(genoprobs[[qtl_chr]][names(sorted_pheno), rev(LETTERS[1:8]), qtl_marker] * 2,
        yaxt = "n", xaxt = "n", col = cols)
  axis(2, at = seq(0, 8, 1 + 1/8)/8, labels = FALSE,
       lty = 0, srt = 90, las = 2)
  mtext(text = main, side = 3, padj = -1, cex = 1.25)
  mtext(text = rev(founders), side = 2, col = rev(label_col), at = seq(0, 1, length.out = 8),
        las = 1, cex = 1.25, adj = 1.25, font = 2)
  mtext(text = paste("lower", "<--", ifelse(is.null(pheno_name), "phenotype", pheno_name), "-->", "higher"), side = 1, padj = 1.25, cex = 1.25)
}

    
####################################################
## Alter the pheno file accordingly 
## pheno file needs to be a data frame for QTL analysis, but for these allele probability plots, it must be a matrix
####################################################

  #need to make the pheno file a matrix so that it runs in the code (currently a data frame)
  #first need to identify what specifically to make part of the phenotype matrix (only need transformed data!)
    names(pheno)
    #from this, I've identified I only need columns 23-27
    #pheno_mat is the matrix of outcomes (phenotypes)
    pheno_mat <- as.matrix(pheno[c(23:27)])
    
  #check rownames to make sure they are already set as the write row names (they are)
    rownames(pheno[c(23:27)])


####################################################
## Review and print the QTL peaks from all of the QTL scans 
####################################################

## I just set threshold to 6 (tells you all of the important qtl peaks with a LOD score > 6)
## map is the qtl2 map you want to use (gmap or pmap)
  qtl_gmap <- find_peaks(zScanSex, map = control$gmap, threshold = 6, peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  qtl_gmap
  
  qtl_pmap <- find_peaks(zScanSex, map = control$pmap, threshold = 6, peakdrop = 1.8, drop = 1.5, expand2markers = FALSE)
  qtl_pmap
  
#Add marker information
  qtl_gmap$marker.id <- find_marker(map = control$gmap, chr = qtl_gmap$chr, pos = qtl_gmap$pos)
  qtl_gmap$marker.id
  qtl_gmap
  
  qtl_pmap$marker.id <- find_marker(map = control$pmap, chr = qtl_pmap$chr, pos = qtl_pmap$pos)
  qtl_pmap$marker.id
  qtl_pmap
  
#set wd 
write_xlsx(list("QTL List RankZ Sex - cM" = qtl_gmap,
                "QTL List RankZ Sex - Mbp" = qtl_pmap),
           "QTL List - RankZ sex.xlsx")
#gives print out of all LOD peaks > 6
#later edited by Becca --> "Final QTL results - RankZ"

####################################################
## ALLELE PLOTS CODE - LOOP
####################################################
  
#set working directory to store the plots
  pdf(file = "allele-plots_cM - RankZ sex.pdf") # create a file called allele-plots.pdf
  # loop through all qtl_gmap above lod threshold of 6 and create an individual plot
  for (i in 1:dim(qtl_gmap)[1]) {
    prob_plot(pheno_vec = pheno_mat[,qtl_gmap$lodcolumn[i]],
              genoprobs = probs,
              qtl_chr = qtl_gmap$chr[i],
              qtl_marker = qtl_gmap$marker.id[i],
              main = paste("lodindex", qtl_gmap$lodindex[i], "Chr", qtl_gmap$chr[i], qtl_gmap$marker.id[i], qtl_gmap$pos[i], qtl_gmap$lodcolumn[i]))
  }
  # be sure to turn the graphics output off at the end!
  dev.off()  
  
#set working directory to store the plots
  pdf(file = "allele-plots_Mbp - RankZ sex.pdf") # create a file called allele-plots.pdf
  # loop through all qtl_gmap above lod threshold of 6 and create an individual plot
  for (i in 1:dim(qtl_pmap)[1]) {
    prob_plot(pheno_vec = pheno_mat[,qtl_pmap$lodcolumn[i]],
              genoprobs = probs,
              qtl_chr = qtl_pmap$chr[i],
              qtl_marker = qtl_pmap$marker.id[i],
              main = paste(qtl_pmap$lodindex[i], "Chr", qtl_pmap$chr[i], qtl_pmap$marker.id[i], qtl_pmap$pos[i], qtl_pmap$lodcolumn[i]))
  }
  # be sure to turn the graphics output off at the end!
  dev.off()  
  

  


  






