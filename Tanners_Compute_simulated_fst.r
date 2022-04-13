#Matthew Wersebe
#University of Oklahoma, Norman
#Tanners Lake Genomics Data Analysis:
#March 09, 2022
################################################################################
# load required libraries
library("adegenet")
library("tidyr")
library("ggplot2")
library("vcfR")
library("hierfstat")
library("qpcR")
library("dplyr")



#set working Directory:
setwd("/home/giovanni/Tanners_demography/MaximumLHood/Simulated_Distribution/ThreePopContigunity_maxL")

#Vector containing the VCF files:

files <- list.files("/home/giovanni/Tanners_demography/MaximumLHood/Simulated_Distribution/ThreePopContigunity_maxL", pattern = "*.vcf")

#initialize list/vectors for the basic stats:

# Read in the pop map:
popmap <- read.table("pop.map", header = F)

# Lists for the basic stats of interest:

sim_overall_fst = vector(mode = "numeric", length = length(files))

sim_overall_fis = vector(mode = "numeric", length = length(files))

sim_perloc_fst = vector(mode = "list", length = length(files))

sim_perloc_fis = vector(mode = "list", length = length(files))


for (i in 1:100){
  #read in VCF; overwrite each time to save memory
  
  sim_vcf <- read.vcfR(files[i])
  
  #convert VCF to the genind object for; overwrite each time to save memory
  
  Sim.genid <- vcfR2genind(sim_vcf)
  
  #indicate ploidy and population assignments
  
  ploidy(Sim.genid) <- 2
  
  pop(Sim.genid) <- popmap$V2
  
  #run basic stats to get Fst and Fis:
  
  sim.bs <- basic.stats(Sim.genid)
  
  print("Done with basic stats")
  
  # Extract the overall stats values:
  
  sim_overall_fst[i] = as.numeric(sim.bs$overall[7])
  
  sim_overall_fis[i] = as.numeric(sim.bs$overall[9])
  
  # Extract Per Locus stats values:
  
  sim_perloc_fst[[i]] = as.numeric(sim.bs$perloc[,7])
  
  sim_perloc_fis[[i]] = as.numeric(sim.bs$perloc[,9])
    
  print("Done")
}

#sanitiy check
sim_overall_fst[100]

sim_perloc_fst[[95]]

# Convert these lists to a dataframe

perloc_fst <- as.data.frame(Reduce(qpcR:::cbind.na, sim_perloc_fst))

# Convert the multi-column data frame to a single column; update name.

perloc_fst <- data.frame(unlist(x = perloc_fst))
 names(perloc_fst)[1] <- "fst"
 
# Take a look and see the top values:
 
perloc_fst %>%
   arrange(desc(fst)) %>%
   slice(1:100)


# Annotate the status of the Fst measures; greater than 99th percentile:

threshold <- quantile(perloc_fst$fst,0.995, na.rm = T)
threshold
perloc_fst <- perloc_fst %>% mutate(outlier = ifelse(fst > threshold, "outlier", "background"))
 
head(perloc_fst)

 
# Plot Nice Figures:

ggplot(perloc_fst, aes(x=fst)) + 
  geom_histogram(color="blue", fill="lightblue", binwidth=0.02)+
  ylab("Number of Sites")+
  xlab("Site-wise Weir and Cockerham Fst")+
  theme_bw()+
  ggtitle("Simulated Site-wise Fst Distribution")+
  geom_vline(xintercept = mean(perloc_fst[,1], na.rm = T), color = "red")

ggplot(perloc, aes(x=Fis)) + 
  geom_histogram(color="blue", fill="lightblue", binwidth=0.02)+
  ylab("Number of Sites")+
  xlab("Site-wise Weir and Cockerham Fis")+
  theme_bw()+
  geom_vline(xintercept = as.numeric(sim.bs$overall[9]), color = "red")



