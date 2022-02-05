# Genotype accumulation curve

## Loading packages

``` r
# Loading packages
library(tidyverse)
library(polysat)
library(hierfstat)
library(adegenet)
library(poppr)
library(RColorBrewer)
```

## SE function

``` r
# Defining SE function
std_mean <- function(x) sd(x)/sqrt(length(x))
```

## Loading data set

``` r
# Loading data set
nameko.raw <- read_csv("MLG_Pmicro_123samples.csv")

# Comvert to genotype as expressed in 6-digit-numeric code
nameko.6_digit.genotype.raw = tibble::tibble(nameko.raw[,c(1:3)],
  Phmi01 = str_c(formatC(nameko.raw$Phmi01A,width=3, flag="0"),formatC(nameko.raw$Phmi01B,width=3, flag="0")),
  Phmi02 = str_c(formatC(nameko.raw$Phmi02A,width=3, flag="0"),formatC(nameko.raw$Phmi02B,width=3, flag="0")),
  Phmi03 = str_c(formatC(nameko.raw$Phmi03A,width=3, flag="0"),formatC(nameko.raw$Phmi03B,width=3, flag="0")),
  Phmi05 = str_c(formatC(nameko.raw$Phmi05A,width=3, flag="0"),formatC(nameko.raw$Phmi05B,width=3, flag="0")),
  Phmi07 = str_c(formatC(nameko.raw$Phmi07A,width=3, flag="0"),formatC(nameko.raw$Phmi07B,width=3, flag="0")),
  Phmi08 = str_c(formatC(nameko.raw$Phmi08A,width=3, flag="0"),formatC(nameko.raw$Phmi08B,width=3, flag="0")),
  Phmi09 = str_c(formatC(nameko.raw$Phmi09A,width=3, flag="0"),formatC(nameko.raw$Phmi09B,width=3, flag="0")),
  Phmi10 = str_c(formatC(nameko.raw$Phmi10A,width=3, flag="0"),formatC(nameko.raw$Phmi10B,width=3, flag="0")),
  Phmi13 = str_c(formatC(nameko.raw$Phmi13A,width=3, flag="0"),formatC(nameko.raw$Phmi13B,width=3, flag="0")),
  Phmi14 = str_c(formatC(nameko.raw$Phmi14A,width=3, flag="0"),formatC(nameko.raw$Phmi14B,width=3, flag="0")),
  Phmi17 = str_c(formatC(nameko.raw$Phmi17A,width=3, flag="0"),formatC(nameko.raw$Phmi17B,width=3, flag="0")),
  Phmi20 = str_c(formatC(nameko.raw$Phmi20A,width=3, flag="0"),formatC(nameko.raw$Phmi20B,width=3, flag="0")),
  Phmi23 = str_c(formatC(nameko.raw$Phmi23A,width=3, flag="0"),formatC(nameko.raw$Phmi23B,width=3, flag="0")),
  Phmi24 = str_c(formatC(nameko.raw$Phmi24A,width=3, flag="0"),formatC(nameko.raw$Phmi24B,width=3, flag="0"))
)

#filtering out the sample "K23" because of missing alleles
nameko.6_digit.genotype  <- nameko.6_digit.genotype.raw[-122,]

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
```

## Plotting Genotype accumulation curve

``` r
# genotype accumulate curve
gac <- genotype_curve(nameko.SSR.genclone, sample = 1000, quiet = TRUE)
```

![](GenotypeAccum_files/figure-markdown_github/unnamed-chunk-4-1.png)

## Genetic diversity

``` r
poppr(nameko.SSR.genclone)
```

    ##        Pop   N MLG eMLG   SE    H    G lambda   E.5  Hexp    Ia  rbarD
    ## 1     Wild  72  60 43.3 1.58 3.99 46.3  0.978 0.851 0.583 0.408 0.0345
    ## 2 Cultivar  50  13 13.0 0.00 2.08  5.9  0.830 0.700 0.252 4.381 0.4684
    ## 3    Total 122  73 35.3 2.53 3.89 27.8  0.964 0.561 0.485 1.790 0.1509
    ##                  File
    ## 1 nameko.SSR.genclone
    ## 2 nameko.SSR.genclone
    ## 3 nameko.SSR.genclone

``` r
poppr(nameko.SSR.subpop.genclone)
```

    ##                          Pop   N MLG  eMLG    SE     H     G lambda   E.5  Hexp
    ## 1                   Hokkaido   3   3  3.00 0.000 1.099  3.00  0.667 1.000 0.472
    ## 2                     Aomori   6   6  6.00 0.000 1.792  6.00  0.833 1.000 0.592
    ## 3                      Iwate   2   2  2.00 0.000 0.693  2.00  0.500 1.000 0.256
    ## 4                      Akita   2   2  2.00 0.000 0.693  2.00  0.500 1.000 0.474
    ## 5                     Miyagi   5   3  3.00 0.000 0.950  2.27  0.560 0.802 0.357
    ## 6                   Yamagata  14  14 10.00 0.000 2.639 14.00  0.929 1.000 0.561
    ## 7                    Niigata   5   5  5.00 0.000 1.609  5.00  0.800 1.000 0.597
    ## 8                  Fukushima  13  11  8.69 0.657 2.311  8.89  0.888 0.869 0.558
    ## 9                     Toyama   1   1  1.00 0.000 0.000  1.00  0.000   NaN 0.538
    ## 10                     Fukui   5   2  2.00 0.000 0.500  1.47  0.320 0.725 0.369
    ## 11                    Nagano   9   7  7.00 0.000 1.889  6.23  0.840 0.932 0.595
    ## 12                   Tottori   2   2  2.00 0.000 0.693  2.00  0.500 1.000 0.462
    ## 13                     Kochi   5   2  2.00 0.000 0.500  1.47  0.320 0.725 0.315
    ## 14           Cultivar.indoor  39   8  4.71 1.007 1.654  4.06  0.753 0.723 0.134
    ## 15        CultOutdoorivar.in   1   1  1.00 0.000 0.000  1.00  0.000   NaN 0.385
    ## 16          Cultivar.outdoor   8   3  3.00 0.000 0.736  1.68  0.406 0.630 0.337
    ## 17 Cultivar.indoor.x.outdoor   1   1  1.00 0.000 0.000  1.00  0.000   NaN 0.462
    ## 18                     China   1   1  1.00 0.000 0.000  1.00  0.000   NaN 0.462
    ## 19                     Total 122  73  8.96 0.970 3.885 27.77  0.964 0.561 0.485
    ##        Ia   rbarD                       File
    ## 1  -0.667 -0.0909 nameko.SSR.subpop.genclone
    ## 2  -0.058 -0.0050 nameko.SSR.subpop.genclone
    ## 3      NA      NA nameko.SSR.subpop.genclone
    ## 4      NA      NA nameko.SSR.subpop.genclone
    ## 5   5.863  0.8111 nameko.SSR.subpop.genclone
    ## 6  -0.210 -0.0179 nameko.SSR.subpop.genclone
    ## 7   0.333  0.0327 nameko.SSR.subpop.genclone
    ## 8   1.067  0.0916 nameko.SSR.subpop.genclone
    ## 9      NA      NA nameko.SSR.subpop.genclone
    ## 10  7.333  1.0000 nameko.SSR.subpop.genclone
    ## 11  1.690  0.1447 nameko.SSR.subpop.genclone
    ## 12     NA      NA nameko.SSR.subpop.genclone
    ## 13  6.364  1.0000 nameko.SSR.subpop.genclone
    ## 14  0.716  0.1848 nameko.SSR.subpop.genclone
    ## 15     NA      NA nameko.SSR.subpop.genclone
    ## 16  7.819  0.8349 nameko.SSR.subpop.genclone
    ## 17     NA      NA nameko.SSR.subpop.genclone
    ## 18     NA      NA nameko.SSR.subpop.genclone
    ## 19  1.790  0.1509 nameko.SSR.subpop.genclone
