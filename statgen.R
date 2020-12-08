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
library(haplo.stats)

library(qvalue)
library(randomForest)

##################################################
## Load data
##################################################

# Full data
fullData = read.delim("http://www.biostat.umn.edu/~cavanr/FMS_data.txt",
                     header=T, sep="\t")

# Genotype data
genotypeData = fullData[,4:226] %>% 
  dplyr::select (-c(b2b, pai1_4g5g, rs849409))

HWE.exact(genotype(genotypeData$actn3_r577x, sep=''))

# Traits and cofounders
traits = fullData %>% 
  dplyr::select(Gender, Age, Race, pre.BMI, NDRM.CH)

##################################################
## Association Analysis WITHOUT genotype imputation
##################################################

## Exercise 1: Genotypic association analysis, no covariates
## SNPs: actn3_r577x and resistin_c180g
## NDRM.CH is the change of strength in non-dominant arm after exercise

# Significant association is found
summary(aov(NDRM.CH ~ actn3_r577x, data=fullData))
# No significant association is found
summary(aov(NDRM.CH ~ resistin_c180g, data=fullData))



## Exercise 2: Additive association analysis with covariates
## SNPs: actn3_r577x and resistin_c180g
## NDRM.CH is the change of strength in non-dominant arm after exercise

# No significant association is found
summary(lm(NDRM.CH ~ allele.count(genotype(actn3_r577x, sep=""))[,1] 
           + Age + Gender + pre.BMI, data=fullData))

# No significant association is found
summary(lm(NDRM.CH ~ allele.count(genotype(resistin_c180g, sep=""))[,1] 
           + Age + Gender + pre.BMI, data=fullData))



## Example 3: dditive association analysis with covariates
## All SNP
## NDRM.CH is the change of strength in non-dominant arm after exercise

# Extract p value for every SNP
k = ncol(genotypeData)
pvalues = rep(NA, k)
for(i in 1:k){
  mod = lm(NDRM.CH ~ allele.count(genotype(genotypeData[,i], sep=""))[,1] 
           + Age + Gender + pre.BMI, data=fullData)
  pvalues[i] = summary(mod)$coefficients[2,4]  
}

# Manhattan Plot
plot(1:k, -log10(pvalues))
abline(h = -log10(0.05))
abline(h = -log10(1-(1-0.05)^(1/k)))

# Which SNP are associated
snpNaiv = which(pvalues <= 0.05, arr.ind=TRUE, useNames = FALSE)
colnames(genotypeData)[snpNaiv]


##################################################
## Association Analysis WITH unobservable phase
##################################################

# Genotype for selected SNP
geno = cbind.data.frame(
  substr(genotypeData$actn3_r577x,1,1), substr(genotypeData$actn3_r577x,1,1),
  substr(genotypeData$gnb3_rs5443,1,1), substr(genotypeData$gnb3_rs5443,1,1), 
  substr(genotypeData$vdr_taq1,1,1), substr(genotypeData$vdr_taq1,1,1),
  substr(genotypeData$akt2_rs892118,1,1), substr(genotypeData$akt2_rs892118,1,1)
)

# Haplotype estimation
HaploEM = haplo.em(geno, locus.label=c("actn3_r577x","gnb3_rs5443","vdr_taq1",'akt2_rs892118'), 
                   control=haplo.em.control(min.posterior=1e-4))
HaploEM

# Regression on haplotypes
geno = setupGeno(geno)

haplomodel = haplo.glm(NDRM.CH~ geno + Age + Gender + pre.BMI, data=fullData,
                       allele.lev=attributes(geno)$unique.alleles,
                       control=haplo.glm.control(haplo.effect="dominant"))
summary(haplomodel)








# Other side analysis

matData = data.matrix(genoData)
matData[is.na(genoData)] <- 4

PCA = prcomp(matData)
plot(PCA$"x"[,1],PCA$"x"[,2],xlab="PC1",ylab="PC2")