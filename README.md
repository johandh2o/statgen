README
================

## Preprocessing

### Packages loading

``` r
library(genetics)
library(tidyverse)
library(FactoMineR)
library(factoextra)
library(lattice)
library(qvalue)
library(pegas)
library(haplo.stats)
```

### Auxiliar functions

``` r
# Update model with same data by changing used covariates
updated = function(object, formula., ..., evaluate = TRUE){
  update(object = object, formula. = formula.,
         data = object$model, ..., evaluate = evaluate)
}
```

### Data loading

``` r
# Complete data
fullData = read.delim("http://www.biostat.umn.edu/~cavanr/FMS_data.txt",
                     header=T, sep="\t")

# Add binary race variable: Caucasian, non-Caucasian
fullData$Caucasian = ifelse(fullData$Race=='Caucasian', 1, 0)
fullData$resp = sqrt(fullData$NDRM.CH)

# Genotype data
genotypeData = fullData[,4:226] %>% 
  dplyr::select (-c(b2b, pai1_4g5g, rs849409,
                    carp_exon3_snp1, carp_exon3_snp2,
                    carp.exon3_snp3, carp_exon5_snp1,
                    gdf8_1225t, gdf8_p198a, gdf8_e164k))

# Traits and cofounders
traits = fullData %>% 
  dplyr::select(Gender, Age, Race, Caucasian, pre.BMI, NDRM.CH)
```

## Subpopulation structure

``` r
# Genotype data in integer form
matData = data.matrix(genotypeData)
missingness = apply(matData, 1, function(x) sum(is.na(x)))/ ncol(matData)
# Missing data is treated as 0 for PCA
matData[is.na(genotypeData)] = 0

# PCA
PCA = prcomp(matData)
fviz_eig(PCA)
```

![](README_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r

# Plot PC1 vs PC2 vs data missingness
data.frame(PC1 = PCA$"x"[,1],
           PC2 = PCA$"x"[,2],
           Missingness = missingness) %>%
  ggplot(aes(x=PC1, y=PC2, color=missingness)) +
  geom_point() + theme_bw() + 
  theme(legend.position="bottom", legend.box = "horizontal")
```

![](README_files/figure-gfm/unnamed-chunk-4-2.png)<!-- -->

``` r
# Plot PC1 vs PC2 vs race
data.frame(PC1 = PCA$"x"[,1],
           PC2 = PCA$"x"[,2],
           Race = traits$Race) %>%
  ggplot(aes(x=PC1, y=PC2, color=Race)) +
  geom_point() + theme_bw() + 
  theme(legend.position="bottom", legend.box = "horizontal")
```

![](README_files/figure-gfm/unnamed-chunk-4-3.png)<!-- -->
