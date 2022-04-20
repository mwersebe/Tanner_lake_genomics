#Matthew Wersebe
#University of Oklahoma, Norman
#Tanners Lake Genomics Data Analysis:
#April 20, 2022
##############################################################################

# TAFT outputs:

# Load required packages:

library(tidyverse)
library(ggplot2)

##############################################################################
setwd("/home/giovanni/Tanners_lake/TAFT")

#Read in results:

TAFT_out <- read_tsv("run_ALL_taft_out.tsv", col_names = TRUE)

TAFT_out %>% dplyr::select(Waples_Test_P) %>% sapply(as.numeric) %>% as.vector ->waples_p

adjust_fdr = p.adjust(waples_p, method = "fdr", n = length(waples_p))

adjust_bonf = p.adjust(waples_p, method = "bonferroni", n = length(waples_p))

log_p_bonf <- -log(adjust_bonf)

log_p_fdr <- -log(adjust_fdr)

TAFT_out %>% mutate(adjust_fdr) %>% mutate(adjust_bonf) %>% mutate(log_p_bonf) %>% mutate(log_p_fdr)-> TAFT_out


ggplot(TAFT_out, aes(POS, log_p_fdr))+
  geom_point()+
  geom_hline(yintercept = -log(0.05))+
  facet_wrap(vars(`CHROM`))+
  ggtitle("Genome Wide TAFT Manhattan Plot- FDR Correction")+
  ylab("-log(FDR Adjusted Waples Test P-value)")+
  xlab("Genomic Position")+
  theme_bw()+
  theme(text = element_text(size = 15))+
  theme(axis.text.x=element_blank())

ggplot(TAFT_out, aes(POS, log_p_bonf))+
  geom_point()+
  geom_hline(yintercept = -log(0.05))+
  facet_wrap(vars(`CHROM`))+
  ggtitle("Genome Wide TAFT Manhattan Plot- Bonferroni Correction")+
  ylab("-log(Bonferroni Adjusted Waples Test P-value)")+
  xlab("Genomic Position")+
  theme_bw()+
  theme(text = element_text(size = 15))+
  theme(axis.text.x=element_blank())

subset(TAFT_out, TAFT_out$adjust_fdr <= 0.05)

subset(TAFT_out, TAFT_out$adjust_bonf <= 0.05)

#174 outliers FDR

#19 outliers Bonferroni




