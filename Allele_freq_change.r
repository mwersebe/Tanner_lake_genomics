# Matthew Wersebe
# University of Oklahoma
# April 8, 2022 
# Biallelic Drift modeling
# Tanners Lake Genomics
################################################################################

library(learnPopGen)

# Model different trajectories for allele frequencies:

#Different starting freqs:

freqs <- seq(from = 0.05, to = 0.95, by = 0.05)

#List to store the data frames

sims <- vector(mode = "list", length = 19)

# Loop through the starting freqs and model drift for 50 generations:

for (i in 1:19){
  freqs[i]
  sim <- genetic.drift(p0 = freqs[i], Ne=3000, nrep = 100, time = 50, show = "p")
  
  sims[[i]] <- as.data.frame(sim[1:50,])
  
}

# Assess if starting freq has an effect on ending freq?

diff <- matrix(data = NULL, nrow = 100, ncol = 19)

for (i in 1:19){
#need nested for loop to work through each section of the list object:
  for (g in 1:100){
  diff[g,i] <- sims[[i]][50,g]- sims[[i]][1,g]
}
}

diff <- c(t(diff))
diff

freq <- rep(freqs,each = 100)

test <- cbind.data.frame(freq, diff)
test

boxplot(diff~freq)

###############################################################################
#Actual SNPs from Tanners lake:

#load libraries:
library(dplyr)
library(readr)
library(tidyr)
library(magrittr)
library(hablar)

################################################################################

setwd("/home/giovanni/Tanners_lake")

Allelefreqs <- read_tsv("All_allele_freqs.tsv", col_names = T)
head(Allelefreqs)

Allelefreqs %>% pivot_longer(cols = c("AF_100" : "MAF_0"), names_to = c("Measure","Population"), names_sep = "_", values_to = "Count") %>%
  pivot_wider(names_from = "Measure", values_from = "Count") %>% mutate(ALT_freq = if_else(AF > MAF, AF, MAF)) %>% mutate(REF_freq = 1-ALT_freq) -> Allelefreqs

Starting <- Allelefreqs %>% filter(Population == "100") %>% select(REF_freq) %>% sapply(as.numeric) %>% as.vector
Starting_Freq <- rep(Starting, each = 3)

Allelefreqs <- Allelefreqs %>% mutate(Starting_Freq)
#Allelefreqs <- Allelefreqs %>% convert(num(Population))

MID <- Allelefreqs %>% filter(Population == "50") %>% select(REF_freq) %>% sapply(as.numeric) %>% as.vector
MID_Freq <- rep(MID, each = 3)

Allelefreqs <- Allelefreqs %>% mutate(MID_Freq)


TOP <- Allelefreqs %>% filter(Population == "0") %>% select(REF_freq) %>% sapply(as.numeric) %>% as.vector
TOP_Freq <- rep(TOP, each = 3)

Allelefreqs <- Allelefreqs %>% mutate(TOP_Freq)

Allelefreqs %>% mutate(Delta_100_50 = Starting_Freq - MID_Freq) %>% mutate(Delta_50_0 = MID_Freq - TOP_Freq) -> Allelefreqs

Allelefreqs <- Allelefreqs %>% convert(num(Population))

BigChangeHighFreq <- 
  Allelefreqs %>% group_by(Starting_Freq >= 0.75 & Starting_Freq <= 1.0) %>% filter(`Starting_Freq >= 0.75 & Starting_Freq <= 1` == "TRUE") %>% filter(Delta_100_50 >= 0.35)

BigChangeHighFreq


BigChangeLowFreq <- 
  Allelefreqs %>% group_by(Starting_Freq >= 0 & Starting_Freq <= 0.25) %>% filter(`Starting_Freq >= 0 & Starting_Freq <= 0.25` == "TRUE") %>% filter(Delta_100_50 <= -0.35)

BigChangeLowFreq
  
###########################################################################  
  

ggplot(BigChangeHighFreq, aes(x = Population, y = REF_freq, group = POS))+
  geom_line()+
  scale_x_reverse()+
  facet_wrap(vars(`CHROM`))+
  ylab("Reference Allele Frequency")+
  xlab("Generation (Approx)")+
  ggtitle("Sites with large reduction in initial frequency")+
  theme_light()+
  theme(text = element_text(size = 15)) 


ggplot(BigChangeLowFreq, aes(x = Population, y = REF_freq, group = POS))+
  geom_line()+
  scale_x_reverse()+
  facet_wrap(vars(`CHROM`))+
  ylab("Reference Allele Frequency")+
  xlab("Generation (Approx)")+
  ggtitle("Sites with large increase in initial frequency")+
  theme_light()+
  theme(text = element_text(size = 15))