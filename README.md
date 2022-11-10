## ARTnet: Model Parameterization with the ARTnet Study Data

ARTnet is an anonymous cross-sectional web-based survey conducted from 2017 to 2019 of HIV-related risk behaviors, testing, and use of prevention services among men who have sex with men (MSM) in the United States. It recruited MSM who have completed the American Men’s Internet Survey (AMIS) study, and therefore, the dataset contains variables merged from that study as well. Full access to the dataset from ARTnet will allow the researchers to conduct analyses and disseminate results using the data. 

For further details on the ARTnet Study, you can read the descriptive paper ["Egocentric Sexual Networks of Men Who Have Sex with Men in the United States: Results from the ARTnet Study"](https://www.sciencedirect.com/science/article/pii/S1755436519301409?via%3Dihub) by Weiss et al. in _Epidemics._ See the **ARTnet Scientific Publications** section below for further details.

### Data Use Agreement
Access to the data requires a Memorandum of Understanding (MOU) that outlines the personnel analyzing the data and purposes of the data analyses. This dataset may not be shared without the consent of the ARTnet Study PI (Samuel Jenness, Emory University) as outlined in an MOU. Please contact the PI by email (mailto:samuel.m.jenness@emory.edu) to request access. A template MOU will be sent; after review access to the ARTnet dataset will be provided via Github.

### ARTnetData Dependency

The ARTnet package depends on the [ARTnetData package](https://github.com/EpiModel/ARTnetData), which contains the limited use public dataset. Because of the restrictions of the dataset, the ARTnetData package must be installed separately, before installing the ARTnet package, using the following directions.

#### Installation
The suggested method for accessing the dataset is to directly install the `ARTnetData` package in R, using the `remotes` package. First, because this repository is private, installing this package requires a [Github Personal Access Token](https://help.github.com/articles/creating-a-personal-access-token-for-the-command-line/).

You can use the `usethis` package to create a new PAT for use in R, like this:
```r
usethis::create_github_token()
```

Copy the PAT, and then put it in your `.Renviron` file as a system variable called `GITHUB_PAT`. You can open your `.Renviron` like this:
```r
usethis::edit_r_environ()
```

Your `.Renviron` file should contain a line like this, but with your own unique PAT on the right-hand side.
```r
GITHUB_PAT=XXXXXX
```

Note that there are other potentially more sophisticated and secure ways to manage your Github PAT, as detailed in the [usethis vignette](https://usethis.r-lib.org/articles/articles/git-credentials.html).

After creating a PAT, restart R/Rstudio and use the `remotes` package to install `ARTnetData`.
```r
remotes::install_github("EpiModel/ARTnetData")
```

#### Dataset Loading
This R package comes with two included datasets, a wide-form dataset (rows = study participants) and a long-form dataset (rows = partnerships, with multiple rows per unique study participant). The suggested method for accessing the dataset is to directly install the `ARTnetData` package in R, using the `remotes` package as follows:
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
The ARTnet package contains standardized scripts to analyze the ARTnet data for the purposes of parameterizing the epidemic modeling with EpiModel and EpiModelHIV. There are three primary functions, detailed below, that conduct statistical analysis of the data for a specific target population of MSM defined by geography, age, and race/ethnicity. Users may also conduct analyses of the ARTnet dataset without the ARTnet package, but this package automates several standard analyses needed for many epidemic modeling projects.

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

### ARTnet Scientific Publications

#### Empirical Analyses

1. Jenness SM, Weiss KM, Prasad P, Zlotorzynska M, Sanchez T. Bacterial STI Screening Rates by Symptomatic Status among Men Who Have Sex with Men in the United States: A Hierarchical Bayesian Analysis. _Sexually Transmitted Diseases._ 2019; 46(1): 25–30. [[PubMed]](https://pubmed.ncbi.nlm.nih.gov/30044334/)

2. Weiss KM, Goodreau SM, Morris M, Prasad P, Ramaraju R, Sanchez T, Jenness SM. Egocentric Sexual Networks of Men Who Have Sex with Men in the United States: Results from the ARTnet Study. _Epidemics._ 2020; 30: 100386. [[PubMed]](https://pubmed.ncbi.nlm.nih.gov/32004795/)

3. Weiss KM, Prasad P, Ramaraju R, Zlotorzynska M, Jenness SM. Estimated Number of Men who have Sex with Men with Indications for HIV Pre-Exposure Prophylaxis in a National Sexual Network Study. _Journal of Acquired Immune Deficiency Syndrome._ 2020; 84(1): 10-17. [[PubMed]](https://pubmed.ncbi.nlm.nih.gov/31939869/)

4. Chandra CL, Weiss KM, Kelley CF, Marcus JL, Jenness SM. Gaps in Screening of Sexually Transmitted Infections among Men Who Have Sex with Men during PrEP Care in the United States. _Clinical Infectious Diseases._ 2021; 73(7): e2261–69. [[PubMed]](https://pubmed.ncbi.nlm.nih.gov/32702116/)

5. Goodreau SM, Maloney KM, Sanchez TH, Morris M, Janulis P, Jenness SM. A Behavioral Cascade of HIV Seroadapation among US Men Who Have Sex with Men in the Era of PrEP and U=U. _AIDS & Behavior._ 2021; 25(12): 3933-3943. [[PubMed]](https://pubmed.ncbi.nlm.nih.gov/33884510/) 

6. Anderson EJ, Weiss KM, Morris M, Sanchez TH, Prasad P, Jenness SM. HIV and Sexually Transmitted Infection Epidemic Potential of Networks of Men Who Have Sex with Men in Two Cities. _Epidemiology._ 2021; 32(5): 681-689. [[PubMed]](https://pubmed.ncbi.nlm.nih.gov/34172692/) 

7. Weiss KM, Prasad P, Sanchez T, Goodreau SM, Jenness SM. Association Between HIV PrEP Indications and Use in a National Sexual Network Study of Men Who Have Sex with Men. _Journal of the International AIDS Society._ 2021; 24(10): e25826. [[PubMed]](https://pubmed.ncbi.nlm.nih.gov/34605174/)

8. 	Maloney KM, Benkeser D, Sullivan PS, Kelley C, Sanchez T, Jenness SM. Sexual Mixing by HIV Status and Pre-exposure Prophylaxis Use Among Men Who Have Sex With Men: Addressing Information Bias. _Epidemiology._ 2022; 33(6): 808-16. [[PubMed]](https://pubmed.ncbi.nlm.nih.gov/35895578/)

9. Mann LM, Le Guillou A, Goodreau SM, Marcus JL, Sanchez T, Weiss KM, Jenness SM. Correlations Between Community-Level HIV Preexposure Prophylaxis Coverage and Individual-Level Sexual Behaviors among US Men Who Have Sex with Men. _AIDS._ 2022. Epub ahead of print. DOI: 10.1097/QAD.0000000000003343. [[PubMed]](https://pubmed.ncbi.nlm.nih.gov/35876641/)

10. Chandra C, Morris M, Van Meter C, Goodreau SM, Sanchez T, Janulis P, Birkett M, Jenness SM. Comparing Sexual Network Mean Active Degree Measurement Metrics among Men Who Have Sex with Men. _Sexually Transmitted Diseases_ 2022. Epub ahead of print. DOI: 10.1097/OLQ.0000000000001708. [[PubMed]](https://pubmed.ncbi.nlm.nih.gov/36112005/)


#### ARTnet Used in HIV/STI Transmission Models

1. Jenness SM, Johnson JA, Hoover KW, Smith DK, Delaney K. Modeling an Integrated HIV Prevention and Care Continuum to Achieve the Ending the HIV Epidemic Goals. _AIDS._ 2020. 34(14): 2103–2113. [[PubMed]](https://pubmed.ncbi.nlm.nih.gov/32910062/)

2. Maloney KM, Driggers R, Sarkar S, Anderson EA, Malik AA, Jenness SM. Projected Impact of Concurrently Available Long-Acting Injectable and Daily-Oral HIV Pre-Exposure Prophylaxis. _Journal of Infectious Diseases._ 2021; 223(1): 72–82. [[PubMed]](https://pubmed.ncbi.nlm.nih.gov/32882043/)

3. Jenness SM, Le Guillou A, Chandra C, Mann L, Sanchez T, Westreich D, Marcus JL. Projected HIV and Bacterial STI Incidence Following COVID-Related Sexual Distancing and Clinical Service Interruption. _Journal of Infectious Diseases._ 2021; 223(6): 1019–28. [[PubMed]](https://pubmed.ncbi.nlm.nih.gov/33507308/)

4. Jenness SM, Knowlton G, Smith DK, Marcus JL, Anderson EJ, Siegler AJ, Jones J, Sullivan PS, Enns E. A Decision Analytics Model to Optimize Investment in Interventions Targeting the HIV PrEP Cascade of Care. _AIDS._ 2021; 35(9): 1479–89. [[PubMed]](https://pubmed.ncbi.nlm.nih.gov/33831910/)

5. Le Guillou A, Buchbinder S, Scott H, Liu A, Havlir D, Scheer S, Jenness SM. Population Impact and Efficiency of Improvements to HIV PrEP Under Conditions of High ART Coverage among San Francisco Men Who Have Sex with Men. _Journal of Acquired Immune Deficiency Syndrome._ 2021; 88(4): 340–347. [[PubMed]](https://pubmed.ncbi.nlm.nih.gov/34354011/)

6. Wheatley MM, Knowlton G, Kao SY, Jenness SM, Enns E. Cost-Effectiveness of Interventions to Improve HIV Pre-Exposure Prophylaxis Initiation, Adherence, and Persistence among Men Who Have Sex with Men. _Journal of Acquired Immune Deficiency Syndrome._ 2022; 90(1): 41-49. [[PubMed]](https://pubmed.ncbi.nlm.nih.gov/35090155/)

7. Jones J, Le Guillou A, Gift TL, Chesson H, Bernstein K, Delaney K, Lyles C, Berruit A, Sullivan PS, Jenness SM. Effect of Screening and Treatment for Gonorrhea and Chlamydia on HIV Incidence among Men who Have Sex with Men in the United States: A Modeling Analysis. _Sexually Transmitted Diseases._ 2022; 49(10):669-676. [[PubMed]](https://pubmed.ncbi.nlm.nih.gov/35921635/)

8. Jenness SM, Le Guillou A, Lyles C, Bernstein KT, Krupinsky K, Enns EA, Sullivan PS, Delaney KP. The Role of HIV Partner Services in the Modern Biomedical HIV Prevention Era: A Network Modeling Study. _Sexually Transmitted Diseases._ 2022. Epub ahead of print. DOI: 10.1097/OLQ.0000000000001711. [[PubMed]](https://pubmed.ncbi.nlm.nih.gov/36194831/)

