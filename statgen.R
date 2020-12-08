##################################################
## Statistical Genetics / Final project
##################################################
## Johan de Aguas, Leo Pham, Yufan Yang
##################################################

##################################################
## Load libraries
##################################################

library(genetics)
library(tidyverse)

library(randomForest)

##################################################
## Load data
##################################################

# Full data
fullData = read.delim("http://www.biostat.umn.edu/~cavanr/FMS_data.txt",
                     header=T, sep="\t")

# Genotype data
genotypeData = fullData[,4:226] %>% select (-c(b2b, pai1_4g5g))

# Cofounders


##################################################
## Descriptive statistics
##################################################

# Number of NA is SNP
sapply(genoData, function(y) sum(length(which(is.na(y)))))





matData = data.matrix(genoData)
matData[is.na(genoData)] <- 4

PCA = prcomp(matData)
plot(PCA$"x"[,1],PCA$"x"[,2],xlab="PC1",ylab="PC2")

xxx = rfImpute(genoData,fullData$DRM.CH)
