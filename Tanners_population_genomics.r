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
library("poppr")
library("RColorBrewer")
library("ape")
library("pegas")

################################################################################
setwd("/home/giovanni/Tanners_lake")
############################################################################
# Read in LD pruned VCF file and population map, convert to gen light:

Tan_LDprune <- read.vcfR("Tanners_CHRall_filter_BiSnps_LDprune_nomiss.vcf.gz")

popmap <- read.table("pop.map", header = F)
head(popmap)

all(colnames(Tan_LDprune@gt)[-1] == popmap$V1)
 #should return TRUE

Tan_LDprune.genlight <- vcfR2genlight(Tan_LDprune)

ploidy(Tan_LDprune.genlight) <- 2

pop(Tan_LDprune.genlight) <- popmap$V2

###############################################################################
# Population Genetic Structure with PCA:
library(parallel)

Tan_PCA <- glPca(Tan_LDprune.genlight, parallel = T, n.cores = 8, nf=3)

# Selected nf = 3

#Add population tags:
Tan_PCA_scores <- as.data.frame(Tan_PCA$scores)
Tan_PCA_scores$pop <- pop(Tan_LDprune.genlight)


#Extract Percent explained variance:
var.exp <- (100*Tan_PCA$eig/sum(Tan_PCA$eig))
head(var.exp)
cols <- brewer.pal(n = nPop(Tan_LDprune.genlight), name = "Paired")

#PCA Bi plot:
PCA_plot <- ggplot(Tan_PCA_scores, aes(x=PC1, y=PC2, colour = pop))+
  geom_point(size=2)+
  stat_ellipse(level = 0.95, size = 1)+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  scale_color_manual(values = cols)+
  xlab("PC1 (7.6%)")+
  ylab("PC2 (4.64%)")+
  theme_bw()

PCA_plot

Tan_PCA_scores

#Tan 22-24 02A and 04B clustering with younger clones. Tan 16-18 10A loadings similar to younger clones.

# Genetic Distance Tree:

Tan.dist <- bitwise.dist(Tan_LDprune.genlight)

Tan.tree <- aboot(Tan_LDprune.genlight, tree = "upgma", distance = bitwise.dist, sample = 100, showtree = F, cutoff = 0, quiet = T)

write.tree(Tan.tree,file="Tan.tre")
tree<-read.tree("Tan.tre")
library(phytools)
plot(tree)




plot.phylo(Tan.tree, cex = 0.8, font = 2, adj = 0, tip.color = cols[pop(Tan_LDprune.genlight)])
  nodelabels(Tan.tree$node.label, adj = c(1.3, -0.5), frame = "n", cex = 0.8,font = 2, xpd = TRUE) 
  legend('topleft', legend = c("10-12","16-18","18-20", "2-4", "22-24", "6-8", "LC"), fill = cols, border = FALSE, bty = "n", cex = 1.0) 
  axis(side = 1) 
  title(xlab = "Genetic distance (proportion of loci that are different)")

##########################################################################################
# Genetic Summary Statistics:
library("hierfstat")
  
popmap <- read.table("popmap.layers", header = F)

Tan_LDprune.genid <- vcfR2genind(Tan_LDprune)

ploidy(Tan_LDprune.genid) <- 2

pop(Tan_LDprune.genid) <- popmap$V2

tan.bs <- basic.stats(Tan_LDprune.genid)

tan.bs$perloc[,7]

Observed_perloc_fst <- as.data.frame(tan.bs$perloc[,7])
names(Observed_perloc_fst)[1] <- "fst"

#mutate with outlier status: check Simulated for threshold:

Observed_perloc_fst <- Observed_perloc_fst %>% mutate(outlier = ifelse(fst > threshold, "outlier", "background"))

Observed_perloc_fst$index <- 1:nrow(Observed_perloc_fst)
outliers = subset(Observed_perloc_fst, Observed_perloc_fst$outlier == "outlier")

#Write the outlier SNPs to a file.
write.table(outliers, "Fst_outliers_snps.txt", append = F, sep = "\t", dec = ".", row.names = F, col.names = T)


#Per Locus Fst Plot:
PerLocFst <- ggplot(tan.bs$perloc, aes(x=Fst)) + 
  geom_histogram(color="blue", fill="lightblue", binwidth=0.02)+
  ylab("Number of Sites")+
  xlab("Site-wise Weir and Cockerham Fst")+
  theme_bw()+
  ggtitle("Observed Site-wise Fst")+
  geom_vline(xintercept = as.numeric(tan.bs$overall[7]), color = "red")

PerLocFst

#Per Locus Fis Plot:
PerLocFis <- ggplot(tan.bs$perloc, aes(x=Fis)) + 
  geom_histogram(color="blue", fill="lightblue", binwidth=0.02)+
  ylab("Number of Sites")+
  xlab("Site-wise Weir and Cockerham Fis")+
  theme_bw()+
  ggtitle("Observed Site-wise Fis")+
  geom_vline(xintercept = as.numeric(tan.bs$overall[9]), color = "red")

PerLocFis


#Overall Genetic Distance as table:
library("reshape2")
library("scales")

tan.gendist <- genet.dist(Tan_LDprune.genid)

tan.gendist <- melt(as.matrix(tan.gendist))

tan.gendist

gendist <- ggplot(tan.gendist, aes(Var1, Var2)) + # x and y axes => Var1 and Var2
  geom_tile(aes(fill = value)) + # background colours are mapped according to the value column
  geom_text(aes(fill = tan.gendist$value, label = round(tan.gendist$value, 3))) +
  scale_fill_gradient(low = "lightblue", high = "darkslategray") + 
  theme(panel.grid.major.x=element_blank(), 
        panel.grid.minor.x=element_blank(), 
        panel.grid.major.y=element_blank(), 
        panel.grid.minor.y=element_blank(),
        panel.background=element_rect(fill="white"), 
        axis.text.x = element_text(angle=90, hjust = 1,vjust=1,size = 12,face = "bold"), 
        plot.title = element_text(size=20,face="bold"),
        axis.text.y = element_text(size = 12,face = "bold")) + 
  ggtitle("Pairwise Genetic Distance") + theme(legend.title=element_text(face="bold", size=14)) + scale_y_discrete(name="") + 
  scale_x_discrete(name="") + labs(fill="Chord Dist.")

gendist



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

