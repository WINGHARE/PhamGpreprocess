library(Biobase)
library(SummarizedExperiment)
library(S4Vectors)
library(PharmacoGx)
source("callingWaterFall.R")

df.sens<- read.csv("data/GDSC2_label_9drugs.csv")

for(cd in colnames(df.sens)[2:9]){
  
  result <- callingWaterfall(df.sens[,c(cd)],type=c("AUC"),intermediate.fold = c(0),plot=TRUE)
  
  df.sens[,paste(cd,"binary",sep=".")] <-result
}

df.sens[,9:18]

write.csv(df.sens[,9:18],"data/GDSC2_label_9drugs_binary.csv")