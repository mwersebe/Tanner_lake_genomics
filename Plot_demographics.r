# Matthew Wersebe
# University of Oklahoma
# April 8, 2022 
# Plotting 2D Site Frequency Sprectra
# Tanners Lake Genomics
################################################################################
# Set WD for plotting Observed SFS:

setwd("/home/giovanni/Tanners_demography/MaximumLHood/ThreePopContigunity")

# Load required libraries:

library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(ggpubr)

##############################################################################
## OBSERVED ## 
#Read in the Site Freq. Spectrum; as.matrix is critical for melt to work properly:

SFS_1_0 <- as.matrix(read.table("ThreePopContigunity_jointMAFpop1_0.obs", skip=1, stringsAsFactors = F, header=T, row.names = 1))

# Use melt to reformat:

SFS_1_0 <- melt(SFS_1_0)

# Rename the columns properly
colnames(SFS_1_0) <- c("MID", "TOP", "Count")

OBS_1_0 <- ggplot(SFS_1_0, aes(x = TOP, y = MID, fill = Count)) +
  geom_tile(color = "white")+
  coord_fixed()+
  theme(axis.text.x = element_blank(), axis.text.y = element_blank())+
  ggtitle("TOP & MID Observed")+
  scale_fill_gradient(low = "white", high = "darkred", guide = "colorbar")+
  guides(fill = guide_colourbar(barwidth = 0.5,
                                barheight = 8))

###############################################################################

SFS_2_0 <- as.matrix(read.table("ThreePopContigunity_jointMAFpop2_0.obs", skip=1, stringsAsFactors = F, header=T, row.names = 1))

# Use melt to reformat:

SFS_2_0 <- melt(SFS_2_0)
head(SFS_2_0)

# Rename the columns properly
colnames(SFS_2_0) <- c("DEEP", "TOP", "Count")

OBS_2_0 <-ggplot(SFS_2_0, aes(x = TOP, y = DEEP, fill = Count)) +
  geom_tile(color = "white")+
  coord_fixed()+
  theme(axis.text.x = element_blank(), axis.text.y = element_blank())+
  ggtitle("TOP & DEEP Observed")+
  scale_fill_gradient(low = "white", high = "darkred", guide = "colorbar")+
  guides(fill = guide_colourbar(barwidth = 0.5,
                                barheight = 8))

###############################################################################

SFS_2_1 <- as.matrix(read.table("ThreePopContigunity_jointMAFpop2_1.obs", skip=1, stringsAsFactors = F, header=T, row.names = 1))

# Use melt to reformat:

SFS_2_1 <- melt(SFS_2_1)
head(SFS_2_1)

# Rename the columns properly
colnames(SFS_2_1) <- c("DEEP", "MID", "Count")

OBS_2_1 <- ggplot(SFS_2_1, aes(x = MID, y = DEEP, fill = Count)) +
  geom_tile(color = "white")+
  coord_fixed()+
  theme(axis.text.x = element_blank(), axis.text.y = element_blank())+
  ggtitle("MID & DEEP Observed")+
  scale_fill_gradient(low = "white", high = "darkred", guide = "colorbar")+
  guides(fill = guide_colourbar(barwidth = 0.5,
                                barheight = 8))
################################################################################

## Simulated ## 
SFS_1_0_sim <- as.matrix(read.table("ThreePopContigunity_jointMAFpop1_0.txt", skip=1, stringsAsFactors = F, header=T, row.names = 1))

# Use melt to reformat:

SFS_1_0_sim <- melt(SFS_1_0_sim)

# Rename the columns properly
colnames(SFS_1_0_sim) <- c("MID", "TOP", "Count")

SIM_1_0 <- ggplot(SFS_1_0_sim, aes(x = TOP, y = MID, fill = Count)) +
  geom_tile(color = "white")+
  coord_fixed()+
  theme(axis.text.x = element_blank(), axis.text.y = element_blank())+
  ggtitle("TOP & MID Simulated")+
  scale_fill_gradient(low = "white", high = "darkred", guide = "colorbar")+
  guides(fill = guide_colourbar(barwidth = 0.5,
                                barheight = 8))

###############################################################################

SFS_2_0_sim <- as.matrix(read.table("ThreePopContigunity_jointMAFpop2_0.txt", skip=1, stringsAsFactors = F, header=T, row.names = 1))

# Use melt to reformat:

SFS_2_0_sim <- melt(SFS_2_0_sim)
head(SFS_2_0_sim)

# Rename the columns properly
colnames(SFS_2_0_sim) <- c("DEEP", "TOP", "Count")

SIM_2_0 <- ggplot(SFS_2_0_sim, aes(x = TOP, y = DEEP, fill = Count)) +
  geom_tile(color = "white")+
  coord_fixed()+
  theme(axis.text.x = element_blank(), axis.text.y = element_blank())+
  ggtitle("TOP & DEEP Simulated")+
  scale_fill_gradient(low = "white", high = "darkred", guide = "colorbar")+
  guides(fill = guide_colourbar(barwidth = 0.5,
                                barheight = 8))

###############################################################################

SFS_2_1_sim <- as.matrix(read.table("ThreePopContigunity_jointMAFpop2_1.txt", skip=1, stringsAsFactors = F, header=T, row.names = 1))
head(SFS_2_1_sim)
# Use melt to reformat:

SFS_2_1_sim <- melt(SFS_2_1_sim)
head(SFS_2_1_sim)

# Rename the columns properly
colnames(SFS_2_1_sim) <- c("DEEP", "MID", "Count")

SIM_2_1 <- ggplot(SFS_2_1_sim, aes(x = MID, y = DEEP, fill = Count)) +
  geom_tile(color = "white")+
  coord_fixed()+
  theme(axis.text.x = element_blank(), axis.text.y = element_blank())+
  ggtitle("MID & DEEP Simulated")+
  scale_fill_gradient(low = "white", high = "darkred", guide = "colorbar")+
  guides(fill = guide_colourbar(barwidth = 0.5,
                                barheight = 8))
################################################################################

# Make a pretty stacked plot:

OBS_1_0
OBS_2_0
OBS_2_1

SIM_1_0
SIM_2_0
SIM_2_1

ggarrange(OBS_1_0, OBS_2_0, OBS_2_1, SIM_1_0, SIM_2_0, SIM_2_1, nrow = 2, ncol = 3)
  