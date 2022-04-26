# Matthew Wersebe
# University of Oklahoma
# Tanners Lake LC50 Analysis
# April 26 2022
#############################################################################

# Set required WD

setwd("/home/giovanni/Tanners_lake/LC50_data")

# Load Required Libaries 

library(MASS)
library(brglm)

# Get list of files in WD for automatically running the GLMs

files <- list.files("/home/giovanni/Tanners_lake/LC50_data", pattern = "*.txt")

# Initalize vectors for results:

LC_50s <- vector(mode = "numeric", length = length(files))

models <- vector(mode = "list", length = length(files))

# For loop to run on files vector consequtively:

for(i in 1:length(files)){
  
data <- read.table(files[i], header = T)

response <- cbind(data$Dead, data$Alive)

models[[i]] <- brglm(response~Concentration, family= binomial, data)

LC_50s[i] <-dose.p(models[[i]],p=0.5)

}

LC_50s

# Read in new data frame with pop and clone IDs:

frame <- read.table("data.frame", header = T, stringsAsFactors = T)

frame <- cbind.data.frame(frame, LC_50s)
frame

# Plot the data:

library(ggplot2)
library(ggpval)

ggplot(frame, aes(x = POP, y = LC_50s))+
  geom_boxplot(outlier.shape = NA, fill = "white", color = "#3366FF")+
  geom_jitter(width = 0.2)+
  ylab(bquote(''~LC[50]~' (mg/L'~Cl^-1~')'))+
  xlab("Temporal subpopulation")+
  ggtitle(bquote('Population Specific '~LC[50]~ 'Estimates'))+
  theme_light()

# Run anova:
#check Levels
levels(frame$POP)

library(dplyr)

group_by(frame, POP)%>%
  summarise(
    count = n(),
    mean = mean(LC_50s, na.rm = T),
    median = median(LC_50s, na.rm = T),
    sd = sd(LC_50s, na.rm = T),
    IQR = IQR(LC_50s, na.rm = T)
  )

#Run non-parametric Anova:

kruskal.test(LC_50s ~ POP, data = frame)

# Post Hoc Test:

pairwise.wilcox.test(frame$LC_50s, frame$POP, p.adjust.method = "none")


