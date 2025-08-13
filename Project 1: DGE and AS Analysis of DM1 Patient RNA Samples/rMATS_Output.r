#Set working directory to Alternative_Splicing_Analysis directory

#Install required package dplyr
BiocManager::install("dplyr")

#load package
library(dplyr)

#reads file in current directory and saves file to variable
SE.MATS.JCEC <- read.delim("./SE.MATS.JCEC.txt")

#first part of filtering file with following cutoffs
filtered_SE <- SE.MATS.JCEC %>%
  filter(abs(SE.MATS.JCEC$IncLevelDifference) > 0.1 & SE.MATS.JCEC$FDR < 0.05)
#filters IncLevelDifference column by absolute values > 0.1 and FDR column for values < 0.05

#creates csv  file
write.csv(filtered_SE, file = "significant_events_SE.csv")

