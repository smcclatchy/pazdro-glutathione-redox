# this script produces a 5-panel figure of histograms for GSH, GSSH, the total,
# ratio and redox potential

# load libraries
library(tidyverse)


# calculate and add in redox potential
redoxPot <- (-264 + 31 * log10(control$pheno[,"Heart_GSSG"]/control$pheno[,"Heart_GSH"]^2))
control$pheno <- cbind(control$pheno, redoxPot)

# add in sex and generation from pheno
control$pheno <- cbind(control$pheno,
                       pheno[match(rownames(control$pheno),
                                   pheno$id), c("id", "sex", "generation")])

# add in sample runs
control$pheno <- cbind(control$pheno,
                       sampleRuns[match(rownames(control$pheno),
                                        sampleRuns$ID),])

# make sure both sets of IDs match up, then remove the column with the ids and 
# leave them as row names

# sum of rownames not equal to id/ID should be 0
sum(rownames(control$pheno) != control$pheno$id)
sum(rownames(control$pheno) != control$pheno$ID)

# what columns are id/ID in? Remove those columns
which(names(control$pheno) %in% c("id", "ID"))
head(control$pheno[, -c(which(names(control$pheno) %in% c("id", "ID")))])
control$pheno <- control$pheno[, -c(which(names(control$pheno) %in% c("id", "ID")))]

# rename the variables that will be plotted
control$pheno2 <- control$pheno %>% rename(c("Glutathione (GSH)"="Heart_GSH", 
                        "Reduced glutathione (GSSG)"="Heart_GSSG", 
                        "Total GSH"="Heart_TotalGSH", 
                        "Ratio GSH:GSSG"="Heart_GSHGSSGRatio",
                        "Redox potential (mV)"="redoxPot")
                        ) 


control$pheno2 <- control$pheno %>%
  select("Glutathione (GSH)", "Reduced glutathione (GSSG)", "Total GSH", 
         "Ratio GSH:GSSG", "Redox potential (mV)") %>% 
  pivot_longer(everything(), names_to = "indicator", values_to = "value")

ggplot(control$pheno2, aes(value)) +
  geom_histogram() +
  facet_wrap(~indicator, nrow = 3, ncol = 2, scales = "free_x", dir = "v") +
  labs(x=NULL, caption = "Histograms of glutathione measurements and derivatives. ") +
  theme_classic()
