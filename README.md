## ARTnet: Model Parameterization with the ARTnet Study Data

ARTnet is an anonymous cross-sectional web-based survey conducted from 2017 to 2019 of HIV-related risk behaviors, testing, and use of prevention services among men who have sex with men (MSM) in the United States. It recruited MSM who have completed the American Men’s Internet Survey (AMIS) study, and therefore, the dataset contains variables merged from that study as well. Full access to the dataset from ARTnet will allow the researchers to conduct analyses and disseminate results using the data. 

For further details on the ARTnet Study, you can read the descriptive paper ["Egocentric Sexual Networks of Men Who Have Sex with Men in the United States: Results from the ARTnet Study"](https://www.sciencedirect.com/science/article/pii/S1755436519301409?via%3Dihub) by Weiss et al. in _Epidemics._

Access to the data requires a Memorandum of Understanding (MOU) that outlines the personnel analyzing the data and purposes of the data analyses. This dataset may not be shared without the consent of the ARTnet Study PI (Samuel Jenness, Emory University) as outlined in an MOU. 

### ARTnetData Dependency

The ARTnet package depends on the ARTnetData package, which contains the limited use public dataset. Because of the restrictions of the dataset, the ARTnetData package must be installed separately, before installing the ARTnet package, using the following directions.

#### Installation
The suggested method for accessing the dataset is to directly install the `ARTnetData` package in R, using the `remotes` package as follows:
```r
remotes::install_github("EpiModel/ARTnetData")
```
Because this repository is private, installing this package may require setting up a 
[Github Personal Access Token](https://help.github.com/articles/creating-a-personal-access-token-for-the-command-line/).

This package comes with two included datasets, a wide-form dataset (rows = study participants) and a long-form dataset (rows = partnerships, with multiple rows per unique study participant).

#### Dataset Loading
The suggested method for accessing the dataset is to directly install the `ARTnetData` package in R, using the `remotes` package as follows:
```r
library("ARTnetData")
d <- ARTnet.wide
l <- ARTnet.long
```

The built-in dataset names for the two datasets are `ARTnet.wide` and `ARTnet.long`, and they are "lazy loaded" into global memory when the `ARTnetData` package is loaded. To use or modify the datasets, you might start by assigning those datasets short-hand names. Then any R operations may be performed. 

```r
str(d)
str(l)
```

### The ARTnet Package
The ARTnet package contains standardized scripts to analyze the ARTnet data for the purposes of parameterizing the epidemic modeling with EpiModel and EpiModelHIV. There are three primary functions, detailed below, that conduct statistical analysis of the data for a specific target population of MSM defined by geography, 
age, and race/ethnicity. Users may also conduct analyses of the ARTnet dataset without the ARTnet package, but this package automates several standard analyses needed for many epidemic modeling projects.

#### Installation
The ARTnet package may be installed with the `remotes` package:
```r
remotes::install_github("EpiModel/ARTnet", build_vignettes = TRUE)
```

#### Example Uses
Some of the example uses are then as follows:

```r
# 1. Epistats: Specify geographic features, as well as race stratification 
#              and total age range
epistats <- build_epistats(geog.lvl = "city", 
                           geog.cat = "Atlanta", 
                           race = TRUE, age.limits = c(30, 50),
                           age.breaks = c(35, 45))

# 2. Netparams: Specify age categories if needed, or let ARTnet determine 
                age categories by number of categories desired
netparams <- build_netparams(epistats = epistats, smooth.main.dur = TRUE)

# 3. Netstats: Finalize network setup 
netstats <- build_netstats(epistats, netparams, expect.mort = 0.0005, 
                           network.size = 1000, edges.avg = TRUE)
```

More details of which may be found in the package vignette:
```r
vignette(package = "ARTnet")
```
