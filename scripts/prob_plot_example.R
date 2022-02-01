## scans is an R object containing your genome scans from scan1()

## plotting function
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

## Grab QTL, I just set threshold to 6
# map is the qtl2 map
qtl <- find_peaks(scans, map, threshold = 6)
# Add marker information
qtl$marker.id <- find_marker(map, chr = qtl$chr, pos = qtl$pos)

## pheno_mat is the matrix of outcomes
prob_plot(pheno_vec = pheno_mat[,qtl$lodcolumn[1]],
          genoprobs = genoprobs,
          qtl_chr = qtl$chr[1],
          qtl_marker = qtl$marker.id[1])

