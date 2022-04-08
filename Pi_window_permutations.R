#Matthew Wersebe
#University of Oklahoma, Norman
#Tanners Lake Genomics Data Analysis:
#March 28, 2022
##############################################################################

# Pi tests

##############################################################################

# Load Required Libraries

library("ggplot2")
library("dplyr")

#############################################################################
# set working directory:

setwd("/home/giovanni/Tanners_lake/PopGenWindows")

##############################################################################
# Deep Population:

#Read in the Outlier Windows:
DEEP_outliers_pi <- read.table("DEEP_outlier_pi.txt", header = F)
head(DEEP_outliers_pi)

#Calculate the Mean for pi
Outlier_pi_DEEP <- mean(DEEP_outliers_pi$V5)

#Read in the Windows:
DEEP_windows_pi <- read.table("DEEP_pi_windows.txt", header =F)
head(DEEP_windows_pi)

dim(DEEP_windows_pi)
#set seed to get consistent results:
set.seed(123)

DEEP_perms <- vector(mode = "list", length = 1000)
DEEP_perms_mean_pi <- vector(mode = "numeric", length = 1000)

for (i in 1:1000){
  DEEP_perms[[i]] <-sample_n(DEEP_windows_pi, 1000, replace = F)
  DEEP_perms_mean_pi[i] <- mean(DEEP_perms[[i]]$V5, na.rm =T)
}

DEEP_perms_mean_pi <- as.data.frame(DEEP_perms_mean_pi)


DEEP <- ggplot(DEEP_perms_mean_pi, aes(x=DEEP_perms_mean_pi)) + 
  geom_histogram(color="blue", fill="lightblue", binwidth=0.0001)+
  ylab(" ")+
  xlab(" ")+
  theme_bw()+
  ggtitle("DEEP Population")+
  geom_vline(xintercept = Outlier_pi_DEEP, color = "red")
DEEP

################################################################################
# Mid Population:

#Read in the Outlier Windows:
MID_outliers_pi <- read.table("MID_outlier_pi.txt", header = F)
head(MID_outliers_pi)

#Calculate the Mean for pi
Outlier_pi_MID <- mean(MID_outliers_pi$V5)

#Read in the Windows:
MID_windows_pi <- read.table("MID_pi_windows.txt", header =F)
head(MID_windows_pi)

dim(MID_windows_pi)
#set seed to get consistent results:
set.seed(231)

MID_perms <- vector(mode = "list", length = 1000)
MID_perms_mean_pi <- vector(mode = "numeric", length = 1000)

for (i in 1:1000){
  MID_perms[[i]] <-sample_n(MID_windows_pi, 1000, replace = F)
  MID_perms_mean_pi[i] <- mean(MID_perms[[i]]$V5, na.rm =T)
}

MID_perms_mean_pi <- as.data.frame(MID_perms_mean_pi)


MID <- ggplot(MID_perms_mean_pi, aes(x=MID_perms_mean_pi)) + 
  geom_histogram(color="blue", fill="lightblue", binwidth=0.0001)+
  ylab(" ")+
  xlab("Window Estimate of Pi")+
  theme_bw()+
  ggtitle("MID Population")+
  geom_vline(xintercept = Outlier_pi_MID, color = "red")
MID

###############################################################################
# Top Population:

#Read in the Outlier Windows:
TOP_outliers_pi <- read.table("TOP_outlier_pi.txt", header = F)
head(TOP_outliers_pi)

#Calculate the Mean for pi
Outlier_pi_TOP <- mean(TOP_outliers_pi$V5)

#Read in the Windows:
TOP_windows_pi <- read.table("TOP_pi_windows.txt", header =F)
head(TOP_windows_pi)

dim(TOP_windows_pi)
#set seed to get consistent results:
set.seed(312)

TOP_perms <- vector(mode = "list", length = 1000)
TOP_perms_mean_pi <- vector(mode = "numeric", length = 1000)

for (i in 1:1000){
  TOP_perms[[i]] <-sample_n(TOP_windows_pi, 1000, replace = F)
  TOP_perms_mean_pi[i] <- mean(TOP_perms[[i]]$V5, na.rm =T)
}

TOP_perms_mean_pi <- as.data.frame(TOP_perms_mean_pi)


TOP <- ggplot(TOP_perms_mean_pi, aes(x=TOP_perms_mean_pi)) + 
  geom_histogram(color="blue", fill="lightblue", binwidth=0.0001)+
  ylab("Number of Windows")+
  xlab(" ")+
  theme_bw()+
  ggtitle("TOP Population")+
  geom_vline(xintercept = Outlier_pi_TOP, color = "red")
TOP
##############################################################################
#Pretty Figure:
library(ggpubr)

Stitched <- ggarrange(TOP, MID, DEEP, ncol = 3, nrow = 1, widths = c(1,1,1.1))
Stitched


###############################################################################
# Empirical P-value for MID Population:

Outlier_pi <- mean(TOP_outliers_pi$V5)
Outlier_pi

pvalue <- sum(MID_perms_mean_pi[,1] <= Outlier_pi, na.rm =T)/1000
 pvalue  

#Only 5 are smaller that this:

 
################################################################################
# Tajima D permutation test:
################################################################################
 
 #Read in the Outlier Windows:
 DEEP_outliers_TD <- read.table("DEEP_outlier_TajimaD.txt", header = T)
DEEP_outliers_TD
 
 #Calculate the Mean for Tajima D 
 Outlier_TD_DEEP <- mean(DEEP_outliers_TD$TajimaD)
 Outlier_TD_DEEP
 
 #Read in the Windows:
 DEEP_windows_TD <- read.table("Tanners_DEEP_tajimaD.tsv", header =T)
 head(DEEP_windows_TD)
 
 dim(DEEP_windows_TD)
 #set seed to get consistent results:
 set.seed(123)
 
 DEEP_perms_TD <- vector(mode = "list", length = 1000)
 DEEP_perms_mean_TD <- vector(mode = "numeric", length = 1000)
 
 for (i in 1:1000){
   DEEP_perms_TD[[i]] <-sample_n(DEEP_windows_TD, 1000, replace = F)
   DEEP_perms_mean_TD[i] <- mean(DEEP_perms_TD[[i]]$TajimaD, na.rm =T)
 }
 
 DEEP_perms_mean_TD <- as.data.frame(DEEP_perms_mean_TD)
 
 head(DEEP_perms_mean_TD)
 
 DEEP_D <- ggplot(DEEP_perms_mean_TD, aes(x=DEEP_perms_mean_TD)) + 
   geom_histogram(color="blue", fill="lightblue", binwidth=0.01)+
   ylab(" ")+
   xlab(" ")+
   theme_bw()+
   ggtitle("DEEP Population")+
   geom_vline(xintercept = Outlier_TD_DEEP, color = "red")
 DEEP_D
 
