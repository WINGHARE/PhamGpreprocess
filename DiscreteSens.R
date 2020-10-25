library(Biobase)
library(SummarizedExperiment)
library(S4Vectors)
library(PharmacoGx)
source("callingWaterFall.R")

df.sens<- read.csv("data/GDSC2_label_9drugs.csv")
reult <- callingWaterfall(df.sens$Cisplatin,type=c("AUC"),intermediate.fold = c(0),plot=TRUE)