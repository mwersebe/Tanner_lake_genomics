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
  theme_light()+
  theme(text = element_text(size = 15))+
  labs(color = " Temporal \n Subpopulation")

PCA_plot

Tan_PCA_scores

#Tan 22-24 02A and 04B clustering with younger clones. Tan 16-18 10A loadings similar to younger clones.

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

Observed_perloc_fst <- Observed_perloc_fst %>% mutate(outlier_99 = ifelse(fst > threshold, "outlier", "background"))

Observed_perloc_fst$index <- 1:nrow(Observed_perloc_fst)


# calculate P values for each Fst value in the observed data set

pvalue<- vector(mode = "numeric", length = length(Observed_perloc_fst[,1]))

for (i in 1:length(Observed_perloc_fst[,1])){
  
  pvalue[i] <- (sum(perloc_fst[,1] >= Observed_perloc_fst[i,1], na.rm =T)/115800)
  
}

Observed_perloc_fst <- cbind.data.frame(Observed_perloc_fst, pvalue)

head(Observed_perloc_fst)

#Fix the INF problem with multiple correction:
Observed_perloc_fst[7229,4] <- 0.00000001
Observed_perloc_fst[7783,4] <- 0.00000001
Observed_perloc_fst[8323,4] <- 0.00000001
Observed_perloc_fst[9235,4] <- 0.00000001
Observed_perloc_fst[9414,4] <- 0.00000001
Observed_perloc_fst[10806,4] <- 0.00000001
Observed_perloc_fst[12670,4] <- 0.00000001
Observed_perloc_fst[13445,4] <- 0.00000001
Observed_perloc_fst[14790,4] <- 0.00000001
Observed_perloc_fst[14799,4] <- 0.00000001
Observed_perloc_fst[18761,4] <- 0.00000001

# Adjust P values for multiple comparisons:

adjust_fdr = p.adjust(Observed_perloc_fst[,4], method = "fdr", n = length(Observed_perloc_fst[,1]))

adjust_bonf = p.adjust(Observed_perloc_fst[,4], method = "bonferroni", n = length(Observed_perloc_fst[,1]))

Observed_perloc_fst <- cbind.data.frame(Observed_perloc_fst, adjust_fdr, adjust_bonf)

head(Observed_perloc_fst)

subset(Observed_perloc_fst, Observed_perloc_fst$adjust_fdr <= 0.01)

library(readr)

write_tsv(Observed_perloc_fst, "Observed_perloc_fst.tsv", quote = "none")


#Write the outlier SNPs to a file. Not run.

#outliers = subset(Observed_perloc_fst, Observed_perloc_fst$adjust_fdr <= 0.01)
#write.table(outliers, "Fst_outliers_snps.txt", append = F, sep = "\t", dec = ".", row.names = F, col.names = T)


#Per Locus Fst Plot:
PerLocFst <- ggplot(tan.bs$perloc, aes(x=Fst)) + 
  geom_histogram(color="blue", fill="lightblue", binwidth=0.02)+
  ylab("Number of Sites")+
  xlab("Site-wise Weir and Cockerham Fst")+
  theme_bw()+
  ggtitle("Observed Site-wise Fst")+
  geom_vline(xintercept = as.numeric(tan.bs$overall[7]), color = "red")+
  theme(text = element_text(size = 15))

PerLocFst

# Site wise Fst Manhattan plot:

man_sites <- read_tsv("Observed_perloc_fst.tsv", col_names = T)

ggplot(man_sites, aes(POS, log_P_bonf))+
  geom_point()+
  geom_hline(yintercept = -log(0.05))+
  facet_wrap(vars(`CHROM`))+
  ggtitle("Genome Wide Fst Manhattan Plot")+
  ylab("-log(Bonferroni Adjusted P-value)")+
  xlab("Genomic Position")+
  theme_bw()+
  theme(text = element_text(size = 15))+
  theme(axis.text.x=element_blank())
  

subset(man_sites, man_sites$log_P_bonf > 0)        


#Overall Genetic Distance as table:
library("reshape2")
library("scales")

tan.gendist <- genet.dist(Tan_LDprune.genid)

tan.gendist <- melt(as.matrix(tan.gendist))

test.df <- as.data.frame(tan.gendist)

test.df$index<- c(5,4,6,2,1,3,8,7,9)

tan.gendist <- test.df[order(test.df$index),]

library(readr)
write_tsv(tan.gendist, "fixed.gendist.tsv")

fixed <- read.table("fixed.gendist.tsv", header = T)
fixed

ggplot(fixed, aes(Var2, Var1)) + # x and y axes => Var1 and Var2
  geom_tile(aes(fill = value)) + # background colors are mapped according to the value column
  geom_text(aes(fill = value, label = round(value, 3)))+
  scale_fill_gradient(low = "lightblue", high = "darkslategray") + 
  theme(panel.grid.major.x=element_blank(), 
        panel.grid.minor.x=element_blank(), 
        panel.grid.major.y=element_blank(), 
        panel.grid.minor.y=element_blank(),
        panel.background=element_rect(fill="white"), 
        axis.text.x = element_text(angle=90, hjust = 1,vjust=1,size = 12,face = "bold"), 
        plot.title = element_text(size=20,face="bold"),
        axis.text.y = element_text(size = 12,face = "bold")) + 
  ggtitle("Pairwise Genetic Distance") + theme(legend.title=element_text(face="bold", size=14))+ scale_y_discrete(name="") + 
  scale_x_discrete(name="") + labs(fill="Chord Dist.")

