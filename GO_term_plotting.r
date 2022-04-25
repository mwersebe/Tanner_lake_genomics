# Matthew Wersebe
# University of Oklahoma
# Tanners Lake Genomics
# GO Term enrichment Visualization
# April 25, 2022
################################################################
# Load required libraries

library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)


# Set working directory:

setwd("/home/giovanni/Tanners_lake/Genes")

# Read in data:

terms <- read_tsv("GO.terms.tsv", col_names = TRUE)

ggplot(data = terms, aes(x = term, y = fold_enrich))+
  geom_bar(aes(fill = FDR_p_value), stat = "identity")+
  theme(panel.grid.major.x=element_blank(), 
        panel.grid.minor.x=element_blank(), 
        panel.grid.major.y=element_blank(), 
        panel.grid.minor.y=element_blank(),
        panel.background=element_rect(fill="white"),
        axis.text.x = element_text(angle=90,size = 9), 
        axis.text.y = element_text(size = 15))+
  xlab("GO Term for Molecular Function")+
  ylab("Fold Enrichment of Term")+
  labs(fill = "FDR\nAdjusted\nP-value")+
  ggtitle("GO Term Enrichment Analysis")+
  theme(text = element_text(size = 15))
