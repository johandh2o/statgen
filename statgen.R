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
library(FactoMineR)
library(factoextra)
library(lattice)
library(qvalue)
library(pegas)
library(haplo.stats)

##################################################
## Auxiliar functions
##################################################

# define a convenience wrapper
updated <- function(object, formula., ..., evaluate = TRUE){
  update(object = object, formula. = formula.,
         data = object$model, ..., evaluate = evaluate)
}

##################################################
## Load data
##################################################

# Full data
fullData = read.delim("http://www.biostat.umn.edu/~cavanr/FMS_data.txt",
                     header=T, sep="\t")

# Binary race variable: Caucasian, non-Caucasian
fullData$Caucasian = ifelse(fullData$Race=='Caucasian', 1, 0)

# Genotype data
genotypeData = fullData[,4:226] %>% 
  dplyr::select (-c(b2b, pai1_4g5g, rs849409,
                    carp_exon3_snp1, carp_exon3_snp2,
                    carp.exon3_snp3, carp_exon5_snp1,
                    gdf8_1225t, gdf8_p198a, gdf8_e164k))

# Traits and cofounders
fullData$respLog = log(1+fullData$NDRM.CH/100)

traits = fullData %>% 
  dplyr::select(Gender, Age, Race, Caucasian, pre.BMI,
                NDRM.CH, respLog)

##################################################
## Subpopulation detection
##################################################

matData = data.matrix(genotypeData)
missingness = apply(matData, 1, function(x) sum(is.na(x)))/ ncol(matData)
matData[is.na(genotypeData)] = 0

PCA = prcomp(matData)

data.frame(PC1 = PCA$"x"[,1],
           PC2 = PCA$"x"[,2],
           Missingness = missingness) %>%
  ggplot(aes(x=PC1, y=PC2, color=missingness)) +
  geom_point() + theme_bw() + 
  theme(legend.position="bottom", legend.box = "horizontal")

data.frame(PC1 = PCA$"x"[,1],
           PC2 = PCA$"x"[,2],
           Race = traits$Race) %>%
  ggplot(aes(x=PC1, y=PC2, color=Race)) +
  geom_point() + theme_bw() + 
  theme(legend.position="bottom", legend.box = "horizontal")

##################################################
## Genotypic association
##################################################

# Extract p value for every SNP
k = ncol(genotypeData)
pvalues = rep(NA, k)
for(i in 1:k){
  # Model 1
  tryCatch( { m1 = lm(respLog ~ genotypeData[,i]+ Caucasian + Age + Gender + pre.BMI,
                      data=fullData) },
            error = function(e) {m1 = lm(respLog ~ genotypeData[,i] + Age + Gender + pre.BMI,
                                         data=fullData)})
  # Model 0
  m0 = updated(m1, .~.-genotypeData[,i])
  # p-values from ANOVA test
  pvalues[i] =  (anova(m1,m0))$`Pr(>F)`[2]
}

# Manhattan Plot
plot(1:k, -log10(pvalues))
abline(h = -log10(0.05))
abline(h = -log10(1-(1-0.05)^(1/k))) # Pvalue correction 

qqmath(~-log10(pvalues),
       distribution=function(x){-log10(qunif(1-x))},
       xlab = '=-log10(expected pvalues)')

# Which SNP are associated using pvalue 0.05
snpNaiv = which(pvalues <= 0.05, arr.ind=TRUE, useNames = FALSE)
colnames(genotypeData)[snpNaiv]

# Q values estimation
qobj = qvalue(pvalues)
min(qobj$pi0.lambda)
min(qobj$pi0.smooth)

# False discovery rate estimation
max(qobj$qvalues[qobj$pvalues <= 0.007])
snpFDR = which(pvalues <= 0.007, arr.ind=TRUE, useNames = FALSE)
colnames(genotypeData)[snpFDR]

##################################################
## Haplotypic association 1
##################################################

# Genotype in gene AKT2
geno1 = genotypeData[,35:38]
# Test for HW equilibrium
hw.test(df2genind(geno1, sep=''))
# Test for Linkage Disequilibrium
LDtable(genetics::LD(makeGenotypes(geno1, sep='')),
        which=c("D", "D'", "r", "X^2","P-value"))

# Haplotype regression
geno1 = cbind.data.frame(
  substr(genotypeData$akt2_7254617,1,1), substr(genotypeData$akt2_7254617,2,2),
  substr(genotypeData$akt2_rs892118,1,1), substr(genotypeData$akt2_rs892118,2,2), 
  substr(genotypeData$akt2_2304186,1,1), substr(genotypeData$akt2_2304186,2,2), 
  substr(genotypeData$akt2_969531,1,1), substr(genotypeData$akt2_969531,2,2)
)
geno1 = setupGeno(geno1)

haplomodel1a = haplo.glm(respLog~ geno1 * Caucasian + Age + Gender + pre.BMI, data=fullData,
                         allele.lev=attributes(geno1)$unique.alleles, 
                         control=haplo.glm.control(haplo.effect="additive"))
aic1a = haplomodel1a$aic
haplomodel1r = haplo.glm(respLog~ geno1 * Caucasian + Age + Gender + pre.BMI, data=fullData,
                         allele.lev=attributes(geno1)$unique.alleles, 
                         control=haplo.glm.control(haplo.effect="recessive"))
aic1r = haplomodel1r$aic
haplomodel1d = haplo.glm(respLog~ geno1 * Caucasian + Age + Gender + pre.BMI, data=fullData,
                         allele.lev=attributes(geno1)$unique.alleles, 
                         control=haplo.glm.control(haplo.effect="dominant"))
aic1d = haplomodel1d$aic

# The recessive model is better according to AIC
c(aic1a, aic1d, aic1r)
summary(haplomodel1r)
qqnorm(haplomodel1r$residuals)
