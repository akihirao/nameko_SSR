# Genetic diversity assessmment

## Loading packages

``` r
# Loading packages
library(tidyverse)
library(hierfstat)
library(adegenet)
library(poppr)
```

## Loading data set

``` r
#initializing
rm(list = ls())

# Loading data set
nameko.raw <- read.csv("MLG_Pmicro_123samples.csv",header=T)

# define locus name
locus.names <- unique(str_sub(colnames(nameko.raw)[-c(1:4)],end=-2))
no.locus <- length(locus.names)

# Convert to 6-digit-numeric coded-genotype
nameko.6_digit.genotype.raw = data.frame()
nameko.6_digit.genotype.raw = data.frame(nameko.raw[,c(1:4)])
for (i in 1:no.locus){
  allele.A.position = 3 + i*2
  allele.B.position = 4 + i*2
  target.loci = locus.names[i]
  nameko.6_digit.genotype.raw = data.frame(nameko.6_digit.genotype.raw, target.loci = paste(formatC(nameko.raw[,allele.A.position],width=3, flag="0"),formatC(nameko.raw[,allele.B.position],width=3, flag="0"),sep=""))
}
colnames(nameko.6_digit.genotype.raw)[-c(1:4)] <- locus.names

#filtering out the sample "K23" because of missing alleles as expressed "NA"
nameko.6_digit.genotype  <- na.omit(nameko.6_digit.genotype.raw)

# Convert to genind
nameko.SSR.genind <- df2genind(nameko.6_digit.genotype[,-c(1:4)],ploidy=2,ncode=3,ind.name=nameko.6_digit.genotype$ID,pop=nameko.6_digit.genotype$Pop)

# Convert to genind
nameko.SSR.subpop.genind <- df2genind(nameko.6_digit.genotype[,-c(1:4)],ploidy=2,ncode=3,ind.name=nameko.6_digit.genotype$ID,pop=nameko.6_digit.genotype$Subpop)

strata(nameko.SSR.genind) <- data.frame(nameko.6_digit.genotype[,c(2:3)])

# Convert to genclone
nameko.SSR.genclone <- as.genclone(nameko.SSR.genind)
nameko.SSR.subpop.genclone <- as.genclone(nameko.SSR.subpop.genind)

#Repeat motif of each of the 14 SSR loci
pinfreps <- c(3,2,3,2,2,2,2,3,3,3,2,2,3,2)

#Clone correction
nameko.SSR.MLG.genind <- clonecorrect(nameko.SSR.genind)
nameko.SSR.MLG.subpop.genind <- clonecorrect(nameko.SSR.subpop.genind)
no.MLG.samples <- length(indNames(nameko.SSR.MLG.genind))

nameko.SSR.MLG.three.categories.genind <- nameko.SSR.MLG.genind
subpop.MLG.list <- nameko.SSR.MLG.subpop.genind@pop
MLG.ID.cultivar.indoor <- which(subpop.MLG.list %in% "Cultivar.indoor")
three.categories.MLG.list <- as.character(nameko.SSR.MLG.genind@pop)
three.categories.MLG.list[MLG.ID.cultivar.indoor] <- "Cultivar.indoor"
three.categories.MLG.list[three.categories.MLG.list=="Cultivar"] <- "Cultivar.others"
three.categories.factor.list <- sort(unique(three.categories.MLG.list))
three.categories.MLG.list <- factor(three.categories.MLG.list,levels=three.categories.factor.list)
nameko.SSR.MLG.three.categories.genind@pop <- three.categories.MLG.list

nameko.SSR.MLG.overall.genind <- nameko.SSR.MLG.genind
overall.MLG.list <- rep("overall",no.MLG.samples)
overall.MLG.list <- factor(overall.MLG.list)
nameko.SSR.MLG.overall.genind@pop <- overall.MLG.list
```

## Summary statistics

``` r
# Defining SE function
std_mean <- function(x) sd(x)/sqrt(length(x))

# Summary
summary.out <- summary(nameko.SSR.MLG.three.categories.genind)

# Number of alleles
A.out <- allele.count(nameko.SSR.MLG.three.categories.genind)
A.overall <- nAll(nameko.SSR.MLG.overall.genind,onlyObserved = TRUE)
A.cultivar.indoor <- nAll(nameko.SSR.MLG.three.categories.genind[pop="Cultivar.indoor"],onlyObserved = TRUE)
A.cultivar.others <- nAll(nameko.SSR.MLG.three.categories.genind[pop="Cultivar.others"],onlyObserved = TRUE)
A.wild <- nAll(nameko.SSR.MLG.three.categories.genind[pop="Wild"],onlyObserved = TRUE)

# Percentage of polymophic loci
P.overall <- sum(A.overall > 1)/length(A.overall)
P.cultivar.indoor <- sum(A.cultivar.indoor > 1)/length(A.cultivar.indoor)
P.cultivar.others <- sum(A.cultivar.others > 1)/length(A.cultivar.others)
P.wild <- sum(A.wild > 1)/length(A.wild)

cat("P: cultivar.indoor\n")
```

    ## P: cultivar.indoor

``` r
print(P.cultivar.indoor)
```

    ## [1] 0.3571429

``` r
cat("\n")
```

``` r
cat("P: cultivar.others\n")
```

    ## P: cultivar.others

``` r
print(P.cultivar.others)
```

    ## [1] 0.8571429

``` r
cat("\n")
```

``` r
cat("P: wild\n")
```

    ## P: wild

``` r
print(P.wild)
```

    ## [1] 1

``` r
cat("\n")
```

``` r
cat("P: overall\n")
```

    ## P: overall

``` r
print(P.overall)
```

    ## [1] 1

``` r
cat("\n")
```

``` r
cat("Number of alleles: cultivar.indoor\n")
```

    ## Number of alleles: cultivar.indoor

``` r
print(mean(A.cultivar.indoor))
```

    ## [1] 1.357143

``` r
print(std_mean(A.cultivar.indoor))
```

    ## [1] 0.1328944

``` r
cat("\n")
```

``` r
cat("Number of alleles: cultivar.others\n")
```

    ## Number of alleles: cultivar.others

``` r
print(mean(A.cultivar.others))
```

    ## [1] 2.785714

``` r
print(std_mean(A.cultivar.others))
```

    ## [1] 0.2997906

``` r
cat("\n")
```

``` r
cat("Number of alleles: wild\n")
```

    ## Number of alleles: wild

``` r
print(mean(A.wild))
```

    ## [1] 6.571429

``` r
print(std_mean(A.wild))
```

    ## [1] 0.8815214

``` r
cat("\n")
```

``` r
cat("Number of alleles: overall\n")
```

    ## Number of alleles: overall

``` r
print(mean(A.overall))
```

    ## [1] 6.642857

``` r
print(std_mean(A.overall))
```

    ## [1] 0.9054605

``` r
cat("\n")
```

``` r
# Allelic richness
AR.out <- hierfstat::allelic.richness(nameko.SSR.MLG.three.categories.genind, min.n = 10)
AR.overall.out <- hierfstat::allelic.richness(nameko.SSR.MLG.overall.genind, min.n = 10)

cat("Allelic richness: mean\n")
```

    ## Allelic richness: mean

``` r
print(apply(AR.out$Ar,2,mean))
```

    ## Cultivar.indoor Cultivar.others            Wild 
    ##        1.353985        2.785714        3.322446

``` r
cat("Allelic richness: SE\n")
```

    ## Allelic richness: SE

``` r
print(apply(AR.out$Ar,2,std_mean))
```

    ## Cultivar.indoor Cultivar.others            Wild 
    ##       0.1317393       0.2997906       0.2687665

``` r
cat("\n")
```

``` r
cat("Allelic richness across overall: mean\n")
```

    ## Allelic richness across overall: mean

``` r
print(apply(AR.overall.out$Ar,2,mean)[[1]])
```

    ## [1] 3.237975

``` r
cat("Allelic richness across overall: SE\n")
```

    ## Allelic richness across overall: SE

``` r
print(apply(AR.overall.out$Ar,2,std_mean)[[1]])
```

    ## [1] 0.2768401

``` r
cat("\n")
```

``` r
# Observed heterozygosity
Ho.out <- Ho(nameko.SSR.MLG.three.categories.genind)
cat("Observed heterozygosity\n")
```

    ## Observed heterozygosity

``` r
print(Ho.out)
```

    ## Cultivar.indoor Cultivar.others            Wild 
    ##       0.1785714       0.4285714       0.3583357

``` r
cat("Observed heterozygosity across overall\n")
```

    ## Observed heterozygosity across overall

``` r
print(Ho(nameko.SSR.MLG.overall.genind)[[1]])
```

    ## [1] 0.34345

``` r
cat("\n")
```

``` r
# Expected heterozygosity
Hs.out <- Hs(nameko.SSR.MLG.three.categories.genind)
cat("Expected heterozygosity\n")
```

    ## Expected heterozygosity

``` r
print(Hs.out)
```

    ## Cultivar.indoor Cultivar.others            Wild 
    ##       0.1506696       0.4714286       0.5611409

``` r
cat("Expected heterozygosity across overall\n")
```

    ## Expected heterozygosity across overall

``` r
print(Hs(nameko.SSR.MLG.overall.genind)[[1]])
```

    ## [1] 0.5366123

``` r
cat("\n")
```

``` r
#Baisit diversity
basic.summary <- basic.stats(nameko.SSR.MLG.three.categories.genind)
#cat("\n")
#cat("Basic statistic across overall\n")
#print(basic.summary$overall)
```
