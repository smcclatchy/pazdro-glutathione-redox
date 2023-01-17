library(qtl2)

# convert genotype probabilities allele probabilities without inserting 
# pseudomarkers.

allele_prob <- genoprob_to_alleleprob(Pazdro_GigaMUGA_genoprobs_qced_8state_sorted)

# calculate kinship matrices for allele probabilities
allele_kinship <- calc_kinship(allele_prob, "loco")

# create a numeric covariate for sex; be sure to include the individual IDs as 
# names.
sex <- setNames(control$pheno$sex, rownames(control$pheno))

# run a genome scan of rankZ transformed heart phenotypes using allele 
# probabilities with sex as additive covariate.
allele_scan <- scan1(allele_prob, control$pheno[, 5:8], allele_kinship, sex)

# plot results of allele scan
plot(allele_scan, control$gmap, lodcolumn = "zHeart_GSH", main = "GSH")
plot(allele_scan, control$gmap, lodcolumn = "zHeart_GSSG", main = "GSSG")
plot(allele_scan, control$gmap, lodcolumn = "zHeart_TotalGSH", 
     main = "Total GSH")
plot(allele_scan, control$gmap, lodcolumn = "zGSHGSSG_Ratio", 
     main = "GSH/GSSG Ratio")

# look at the QTL effects for chr14 GSH peak.
coef_c14 <- scan1coef(allele_prob[,"14"], 
                      control$pheno["zHeart_GSH"], 
                      allele_kinship[["14"]], 
                      sex)
plot_coefCC(coef_c14, control$gmap["14"], bgcolor="gray95", 
            legend="bottomleft")

c14effects <- scan1coef(probs[,"14"], control$pheno["zHeart_GSH"])

# download SQLite databases for mouse genome build 38 on Figshare
download.file("https://ndownloader.figshare.com/files/18533342", 
              "cc_variants.sqlite")
download.file("https://ndownloader.figshare.com/files/17609252", 
              "mouse_genes_mgi.sqlite")

# create a function for querying the CC variants
query_variants <- create_variant_query_func("../data/cc_variants.sqlite")

# grab the variants in the interval 49-55 Mbp on chromosome 14
variants_14_49_55 <- query_variants(14, 49, 55)

# create a function for querying the MGI mouse gene annotations
query_genes <- create_gene_query_func("../data/mouse_genes_mgi.sqlite")

# grab the genes overlapping the interval 49-55 Mbp on chromosome 14
genes_14_49_55 <- query_genes(14, 49, 55)

# find the location of the inferred QTL. The calculations were only performed at 
# the marker positions, and so we can use max(), giving both the scan1() output 
# and the physical map, and then pull out the position from the results.
peak_Mbp <- max(zScanSex, control$pmap)$pos

# pull out the variants in the confidence interval on chr 14 using the query 
# function we defined above:
# confidence interval is 48.90916-55.3026
# peak is 54.23953
variants <- query_variants(14, peak_Mbp - 5.33, peak_Mbp + 1.07)

#  All of these steps are combined into a single function scan1snps(), which takes the genotype probabilities, a physical map of those locations, the phenotype, the kinship matrix, covariates, the query function for grabbing the SNPs in a given interval, and the chromosome, start, and end positions for the interval. (If insert_pseudomarkers() was used to insert pseudomarkers, you will need to use interp_map() to get interpolated Mbp positions for the pseudomarkers.)
  
scan_snps <- scan1snps(probs, control$pmap, 
                       control$pheno[,c("zHeart_GSH", "zHeart_GSSG")],
                       kinship_loco[["14"]], 
                       sex, 
                       query_func=query_variants,
                        chr=14, 
                       start=peak_Mbp-5.33, 
                       end=peak_Mbp+1.07, 
                       keep_all_snps=TRUE)
  
#  The output is a list with two components: lod is a matrix of LOD scores (with a single column, since we’re using just one phenotype), and snpinfo is a data frame with SNP information. With the argument keep_all_snps=TRUE, the snpinfo data frame contains information about all of the variants in the region with an index column indicating the equivalence classes.
  
#  The function plot_snpasso() can be used to plot the results, with points at each of the SNPs. The default is to plot all SNPs. In this case, there are 27737 variants in the region, but only 150 distinct ones.
  
plot(scan_snps$lod, scan_snps$snpinfo)

#  We can use our query_genes() function to identify the genes in the region, and plot_genes() to plot their locations. But also plot_snpasso() can take the gene locations with the argument genes and then display them below the SNP association results. Here, we are also highlighting the top SNPs in the SNP association plot using the drop_hilit argument. SNPs with LOD score within drop_hilit of the maximum are shown in pink.
  
genes <- query_genes(14, peak_Mbp - 5.33, peak_Mbp + 1.07)
plot(scan_snps$lod, scan_snps$snpinfo, drop_hilit=1.5, genes=genes)
  
#  To get a table of the SNPs with the largest LOD scores, use the function top_snps(). This will show all SNPs with LOD score within some amount (the default is 1.5) of the maximum SNP LOD score. We’re going to display just a subset of the columns.
  
top <- top_snps(scan_snps$lod, scan_snps$snpinfo)
print(top[,c(1, 8:15, 20)], row.names=FALSE)
  
 # You can include a visualization of the strain distribution pattern of the SNPs in the association plot with the argument sdp_panel=TRUE.
  
plot(scan_snps$lod, scan_snps$snpinfo, drop_hilit=1.5, genes=genes, sdp_panel=TRUE)
  
#  The scan1snps() function can also be used to perform a genome-wide SNP association scan, by providing a variant query function but leaving chr, start, and end unspecified. In this case it’s maybe best to use keep_all_snps=FALSE (the default) and only save the index SNPs.
  
out_gwas <- scan1snps(probs, control$pmap, 
                      control$pheno[,c("zHeart_GSH", "zHeart_GSSG")],
                      kinship_loco, sex, 
                      query_func=query_variants)
  
# We can make a Manhattan plot of the results as follows. We use altcol to define a color for alternate chromosomes and gap=0 to have no gap between chromosomes.
  
plot(out_gwas$lod, out_gwas$snpinfo, altcol="green4", gap=0)
