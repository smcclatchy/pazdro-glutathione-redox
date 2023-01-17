# retrieve heart expression data from EBI Expression Atlas
# https://www.ebi.ac.uk/gxa
# cross-reference this data with genes from qtl peak intervals
# 
# install packages from Bioconductor
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("SummarizedExperiment")
BiocManager::install("biomaRt")

library(SummarizedExperiment)
library(biomaRt)
library(tidyverse)

# example expression data - see EBI URL above for others
# replace URL below with your URL of interest, then
# update the output filename at the end of this script
# experiments-content/E-MTAB-599 E-MTAB-4644 E-MTAB-3725
load(url("https://www.ebi.ac.uk/gxa/experiments-content/E-MTAB-599/static/E-MTAB-599-atlasExperimentSummary.Rdata"))

# Get counts.
counts = assays(experiment_summary$rnaseq)$counts

# Get gene annotation. (unfortunately, this data set doesn't have gene annotation).
annot = rowRanges(experiment_summary$rnaseq)

# Get sample info.
samples = colData(experiment_summary$rnaseq)
names(samples)

# what tissues or organs for expression data? which stages?
unique(samples$organism_part)
unique(samples$developmental_stage)

# find samples with expression in heart
samples[which(samples$organism_part == "heart"),]

# get sample IDs for these
rownames(samples) # all sample IDs
heart_sample_IDs <- rownames(samples[which(samples$organism_part == "heart"), ]) # only those for heart

# find count data for these sample IDs
colnames(counts)
colnames(counts) == heart_sample_IDs
which(colnames(counts) == heart_sample_IDs)
heart_counts <- counts[, which(colnames(counts) == heart_sample_IDs)]

# change heart counts to data frame to make it easier to deal with
# make row names one of the variables so ensembl IDs are in the data frame
heart_counts <- as.data.frame(heart_counts)

# retrieve Ensembl gene IDs for qtl peak intervals using Biomart
# select mouse genome database
ensembl <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl")
# filters are just that - how we want to filter the data 
# (e.g. by chromosome & region)
filters = listFilters(ensembl)
filters[1:6,]

# attributes are the data we want to retrieve from a query (e.g. Ensembl IDs)
attributes = listAttributes(ensembl)
attributes[1:6,]

chr14_conf_int <- list(14, 49239713, 54817541)

# query chromosome 14 from 49.239713	to 54.817541 Mbp
chr14_genes <- getBM(attributes = c("ensembl_gene_id", "description", "chromosome_name", 
                     "start_position", "end_position", "gene_biotype",
                     "phenotype_description", "name_1006", 
                     "strain_name", "strain_gender"),
      filters = c("chromosome_name", "start", "end"),
      values = chr14_conf_int,
      mart = ensembl)

names(chr14_genes)

# intersect Ensembl gene IDs from interval with expression results
# Ensembl IDs are row names of heart expression data
# unless you only have one sample! then just
#  use names(heart_counts) instead of rownames(heart_counts)
#  no need to unlist
# first column of biomart query are Ensembl IDs
rownames(heart_counts)
(chr14_genes)$ensembl_gene_id

# number of matching IDs between the two
sum(chr14_genes$ensembl_gene_id %in% unlist(rownames(heart_counts)))

sum(chr14_genes$ensembl_gene_id %in% unlist(rownames(heart_counts)))
heart_counts <- cbind(heart_counts, "ensembl_gene_id"=rownames(heart_counts))

# now combine counts and biomart data into one
chr14_genes <- chr14_genes %>% left_join(heart_counts, 
                                         by = "ensembl_gene_id")
head(chr14_genes)

# update the output filename to reflect the expression source data
write.table(chr14_genes, file = "E-MTAB-599_chr14_genes.tsv", sep = "\t",  
            quote = FALSE, row.names = FALSE)
