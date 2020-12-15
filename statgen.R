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
## Functions
##################################################

# define a convenience wrapper
updated <- function(object, formula., ..., evaluate = TRUE){
  update(object = object, formula. = formula., data = object$model, ..., evaluate = evaluate)
}

##################################################
## Load data
##################################################

# Full data
fullData = read.delim("http://www.biostat.umn.edu/~cavanr/FMS_data.txt",
                     header=T, sep="\t")

# Genotype data
genotypeData = fullData[,4:226] %>% 
  dplyr::select (-c(b2b, pai1_4g5g, rs849409,
                    carp_exon3_snp1, carp_exon3_snp2,
                    carp.exon3_snp3, carp_exon5_snp1,
                    gdf8_1225t, gdf8_p198a, gdf8_e164k))


# Check for HWE in one SNP
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
## One SNP
## NDRM.CH is the change of strength in non-dominant arm after exercise

# No significant association is found 
summary(lm(NDRM.CH ~ allele.count(genotype(resistin_c180g, sep=""))[,1] 
           + Age + Gender + pre.BMI, data=fullData))



## Exercise 3: Complete factor association analysis with covariates
## One SNP
## NDRM.CH is the change of strength in non-dominant arm after exercise

# Complete factor model
m1 = lm(NDRM.CH ~ actn3_r577x + Race + Age + Gender + pre.BMI, data=fullData)
m0 = updated(m1, .~.-actn3_r577x)
anova(m1, m0)


## Example 4:  Association analysis with covariates
## All SNP
## NDRM.CH is the change of strength in non-dominant arm after exercise

# Extract p value for every SNP
k = ncol(genotypeData)
pvalues = rep(NA, k)
for(i in 1:k){
  tryCatch( { m1 = lm(NDRM.CH ~ genotypeData[,i]+ Race + Age + Gender + pre.BMI, data=fullData) },
            error = function(e) {m1 = lm(NDRM.CH ~ genotypeData[,i] + Age + Gender + pre.BMI, data=fullData)})
  m0 = updated(m1, .~.-genotypeData[,i])
  pvalues[i] =  (anova(m1,m0))$`Pr(>F)`[2]
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
  substr(genotypeData$akt2_7254617,1,1), substr(genotypeData$akt2_7254617,2,2),
  substr(genotypeData$akt2_2304186,1,1), substr(genotypeData$akt2_2304186,2,2), 
  substr(genotypeData$ampd1_c34t,1,1), substr(genotypeData$ampd1_c34t,2,2)
)

# Haplotype estimation
HaploEM = haplo.em(geno, locus.label=c("akt2_7254617","akt2_2304186",'ampd1_c34t'), 
                   control=haplo.em.control(min.posterior=1e-4))
HaploEM

# Regression on haplotypes
geno = setupGeno(geno)

haplomodel = haplo.glm(NDRM.CH~ geno + Age + Gender + pre.BMI, data=fullData,
                       allele.lev=attributes(geno)$unique.alleles,
                       control=haplo.glm.control(haplo.effect="recessive"))
summary(haplomodel)





texreg(l = list(mod[[1]], mod[[139]], mod[[102]], mod[[210]]))


# Other side analysis

matData = data.matrix(genoData)
matData[is.na(genoData)] <- 4

PCA = prcomp(matData)
plot(PCA$"x"[,1],PCA$"x"[,2],xlab="PC1",ylab="PC2")