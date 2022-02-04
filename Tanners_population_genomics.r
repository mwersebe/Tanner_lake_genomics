#Matthew Wersebe
#University of Oklahoma, Norman
#Tanners Lake Genomics Data Analysis:
#Jan 30, 2022
################################################################################

# load required libraries
library("adegenet")
library("tidyr")
library("ggplot2")
library("vcfR")

################################################################################
setwd("/home/giovanni/Tanners_lake")
############################################################################
# Read in LD pruned VCF file and population map, convert to gen light:

Tan_LDprune <- read.vcfR("Tanners_CHRall_filter_BiSnps_intron_LDprune.vcf.gz")

popmap <- read.table("pop.map", header = F)
head(popmap)

all(colnames(Tan_LDprune@gt)[-1] == popmap$V1)
 #should return TRUE

Tan_LDprune.genlight <- vcfR2genlight(Tan_LDprune)

ploidy(Tan_LDprune.genlight) <- 2

pop(Tan_LDprune.genlight) <- popmap$V2

###############################################################################
# Population Genetic Structure with PCA:


Tan_PCA <- glPca(Tan_LDprune.genlight, parallel = T, n.cores = 8)

# Selected nf = 3

#Add population tags:
Tan_PCA_scores <- as.data.frame(Tan_PCA$scores)
Tan_PCA_scores$pop <- pop(Tan_LDprune.genlight)


#Extract Percent explained variance:
var.exp <- (100*Tan_PCA$eig/sum(Tan_PCA$eig))
head(var.exp)

#PCA Bi plot:
PCA_plot <- ggplot(Tan_PCA_scores, aes(x=PC1, y=PC2, colour = pop))+
  geom_point(size=2)+
  stat_ellipse(level = 0.95, size = 1)+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  xlab("PC1 (7.48%)")+
  ylab("PC2 (4.34%)")+
  theme_bw()

PCA_plot

Tan_PCA_scores

#Tan 22-24 02A and 04B clustering with younger clones. Tan 16-18 10A loadings similar to younger clones.

##########################################################################################
library(dplyr)
library(readr)
library(tidyr)
library(purrr)
library(magrittr)
library(hablar)
#Tanners Temporal Allele Frequency Changes:

#Using the intronic LD pruned datat set, extracted Tan 22-24 'polymorphic' SNPs.
#Polymorphic being anything that had a frequency of >= doubleton (i.e., excluding all singletons)
#Extracted these SNPs from with AF and MAF from 10-12 and Lake Clones (0) subpops as well.

AlleleFreqs <- read_tsv("Tanners_polymorphic.tsv", col_names = T)
head(AlleleFreqs)

AlleleFreqs %>% pivot_longer(cols = c("NumberSamples_22": "MAF_0"), names_to = c("Measure","Population"), names_sep = "_", values_to = "Count") %>%
  pivot_wider(names_from = "Measure", values_from = "Count") -> AlleleFreqs
  

AlleleFreqs

# If the value of MAF = AlleleFreq, ALT is the minor allele (MA)
# If the value of MAF < AlleleFreq, REF is the minor allele (MA)
# Write code to standardize the Freq of the ALT allele. 

AlleleFreqs %>% mutate(ALT_freq = if_else(AlleleFreq > MAF, AlleleFreq, MAF)) -> AlleleFreqs
AlleleFreqs

# Group the data by starting ALT_freq 
# Bins: SAF >= 0.140 -0.25; >0.25 - 0.5, >0.5 -0.75, >0.75


Starting <- AlleleFreqs %>% filter(Population == "22") %>% select(ALT_freq) %>% sapply(as.numeric) %>% as.vector

Starting_Freq <- rep(Starting, each = 7)

AlleleFreqs <- AlleleFreqs %>% mutate(Starting_Freq)
AlleleFreqs <- AlleleFreqs %>% convert(num(Population))

LowFreqAlleles <- AlleleFreqs %>% group_by(Starting_Freq >= 0.14 & Starting_Freq < 0.25) %>% filter(`Starting_Freq >= 0.14 & Starting_Freq < 0.25` == "TRUE") %>% convert(num(Population))

MedFreqAlleles <- AlleleFreqs %>% group_by(Starting_Freq >= 0.25 & Starting_Freq < 0.5) %>% filter(`Starting_Freq >= 0.25 & Starting_Freq < 0.5` == "TRUE") %>% convert(num(Population))

MedHighFreqAlleles <- AlleleFreqs %>% group_by(Starting_Freq >= 0.5 & Starting_Freq < 0.75) %>% filter(`Starting_Freq >= 0.5 & Starting_Freq < 0.75` == "TRUE") %>% convert(num(Population))

HighFreqAlleles <- AlleleFreqs %>% group_by(Starting_Freq >= 0.75 & Starting_Freq <= 1.0) %>% filter(`Starting_Freq >= 0.75 & Starting_Freq <= 1` == "TRUE") %>% convert(num(Population))

# Plot ALT Allele Frequency over time usings spaghetti plots. Facet Wrap per chromosome.


ggplot(AlleleFreqs, aes(x = Population, y = ALT_freq, group = BasePair))+
  geom_line()+
  scale_x_reverse()+
  facet_wrap(vars(`#CHROM`))+
  ggtitle("Allele Frequency Change Over Time")

ggplot(LowFreqAlleles, aes(x = Population, y = ALT_freq, group = BasePair))+
  geom_line()+
  scale_x_reverse()+
  facet_wrap(vars(`#CHROM`))+
  ggtitle("Low Starting Frequency")

ggplot(MedFreqAlleles, aes(x = Population, y = ALT_freq, group = BasePair))+
  geom_line()+
  scale_x_reverse()+
  facet_wrap(vars(`#CHROM`))+
  ggtitle("Medium Starting Frequency")

ggplot(MedHighFreqAlleles, aes(x = Population, y = ALT_freq, group = BasePair))+
  geom_line()+
  scale_x_reverse()+
  facet_wrap(vars(`#CHROM`))+
  ggtitle("Medium-High Starting Frequency")

ggplot(HighFreqAlleles, aes(x = Population, y = ALT_freq, group = BasePair))+
  geom_line()+
  scale_x_reverse()+
  facet_wrap(vars(`#CHROM`))+
  ggtitle("High Starting Frequency")

###############################################################################

# Plot by allele type: Doubletons, tripletons, etc. 

Doubletons <- AlleleFreqs %>% group_by(Starting_Freq >= 0.14 & Starting_Freq < 0.17) %>% filter(`Starting_Freq >= 0.14 & Starting_Freq < 0.17` == "TRUE") %>% convert(num(Population))

ggplot(Doubletons, aes(x = Population, y = ALT_freq, group = BasePair))+
  geom_line()+
  scale_x_reverse()+
  facet_wrap(vars(`#CHROM`))+
  ggtitle("Allele Trajectory for Doubletons")


Tripletons <- AlleleFreqs %>% group_by(Starting_Freq >= 0.21 & Starting_Freq <= 0.25) %>% filter(`Starting_Freq >= 0.21 & Starting_Freq <= 0.25` == "TRUE")

ggplot(Tripletons, aes(x = Population, y = ALT_freq, group = BasePair))+
  geom_line()+
  scale_x_discrete(limits = rev)+
  facet_wrap(vars(`#CHROM`))+
  ggtitle("Allele Trajectory for Tripletons")

