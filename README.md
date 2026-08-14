## ARTnet: Model Parameterization with the ARTnet Study Data

<!-- badges: start -->
[![R-CMD-check](https://github.com/EpiModel/ARTnet/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/EpiModel/ARTnet/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

ARTnet is an anonymous cross-sectional web-based survey conducted from 2017 to 2019 of HIV-related risk behaviors, testing, and use of prevention services among men who have sex with men (MSM) in the United States. It recruited MSM who have completed the American Men’s Internet Survey (AMIS) study, and therefore, the dataset contains variables merged from that study as well. Full access to the dataset from ARTnet will allow the researchers to conduct analyses and disseminate results using the data.

For further details on the ARTnet Study, you can read the descriptive paper ["Egocentric Sexual Networks of Men Who Have Sex with Men in the United States: Results from the ARTnet Study"](https://www.sciencedirect.com/science/article/pii/S1755436519301409?via%3Dihub) by Weiss et al. in _Epidemics._ See the **ARTnet Scientific Publications** section below for further details.

### Data Use Agreement

Access to the data requires a Memorandum of Understanding (MOU) that outlines the personnel analyzing the data and purposes of the data analyses. This dataset may not be shared without the consent of the ARTnet Study PI (Samuel Jenness, Emory University) as outlined in an MOU. Please contact the PI by email (mailto:samuel.m.jenness@emory.edu) to request access. A template MOU will be sent; after review access to the ARTnet dataset will be provided via Github.

### ARTnetData Dependency

The ARTnet package depends on the [ARTnetData package](https://github.com/EpiModel/ARTnetData), which contains the limited use public dataset. Because of the restrictions of the dataset, the ARTnetData package must be installed separately, to fully utilize the ARTnet package, using the following directions.

A set of built-in analyses are available within the package as well to experiment before gaining access to the full data-set. See below.

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

### Using the Built-In Models

For `EpiModelHIV` models, the `epistats` and `netstats` object are required. To
experiment with it without access to `ARTnetData`, `ARTnet` provides an example
for each. They are the result of the following calls:

```r
epistats <- build_epistats(
  geog.lvl = "city",
  geog.cat = "Atlanta",
  init.hiv.prev = c(0.33, 0.137, 0.084),
  race = TRUE,
  time.unit = 7
)

netparams <- build_netparams(
  epistats = epistats,
  smooth.main.dur = TRUE
)

netstats <- build_netstats(
  epistats,
  netparams,
  expect.mort = 0.000478213,
  network.size = 1e3
)
```

They can be accessed copied locally with:

```r
file.copy(
  system.file("netstats-example.rds", package = "ARTnet"),
  "netstats.rds"
)

file.copy(
  system.file("epistats-example.rds", package = "ARTnet"),
  "epistats.rds"
)
```

### ARTnet Scientific Publications

Last reconciled 14 August 2026 against PubMed, Europe PMC, OpenAlex, Semantic Scholar, Crossref, and the preprint servers. **48 distinct papers and preprints, 30 distinct first authors, and lead or co-lead institutions in 7 countries, spanning 2019 to 2026. Twenty-one of the 48 have no author from the study team.**

Two notes on citing ARTnet. The canonical descriptive paper's DOI is `10.1016/j.epidem.2020.100386`. Weiss et al. 2020 reports 4,904 participants and 16,198 partnerships; the figure of 15,809 that appears in some downstream work is a post-exclusion analytic count and is not the study total.

If you have published work using ARTnet that is not listed here, please contact the PI so it can be added.

#### Empirical Analyses of ARTnet Data

1. Anderson EJ, Gandhi NR, Sanchez TH, Enns EA, Jenness SM. Predictive Models to Identify Factors Associated with a High One-Time Sexual Partnership Rate Among U.S. Men Who Have Sex with Men. _AIDS and Behavior._ 2026. [[PubMed]](https://pubmed.ncbi.nlm.nih.gov/41934601/)
2. Chandra C, Morris M, Van Meter C, Goodreau SM, Sanchez T, Janulis P, Birkett M, Jenness SM. Comparing Sexual Network Mean Active Degree Measurement Metrics Among Men Who Have Sex With Men. _Sexually Transmitted Diseases._ 2022; 49(12): 808-814. [[PubMed]](https://pubmed.ncbi.nlm.nih.gov/36112005/)
3. Maloney KM, Benkeser D, Sullivan PS, Kelley C, Sanchez T, Jenness SM. Sexual Mixing by HIV Status and Pre-exposure Prophylaxis Use Among Men Who Have Sex With Men: Addressing Information Bias. _Epidemiology._ 2022; 33(6): 808-816. [[PubMed]](https://pubmed.ncbi.nlm.nih.gov/35895578/)
4. Mann LM, Le Guillou A, Goodreau SM, Marcus JL, Sanchez T, Weiss KM, Jenness SM. Correlations Between Community-Level HIV Preexposure Prophylaxis Coverage and Individual-Level Sexual Behaviors Among United States MSM. _AIDS._ 2022; 36(14): 2015-2023. [[PubMed]](https://pubmed.ncbi.nlm.nih.gov/35876641/)
5. Goodreau SM, Maloney KM, Sanchez TH, Morris M, Janulis P, Jenness SM. A Behavioral Cascade of HIV Seroadaptation Among US Men Who Have Sex with Men in the Era of PrEP and U = U. _AIDS and Behavior._ 2021; 25(12): 3933-3943. [[PubMed]](https://pubmed.ncbi.nlm.nih.gov/33884510/)
6. Weiss KM, Prasad P, Sanchez T, Goodreau SM, Jenness SM. Association Between HIV PrEP Indications and Use in a National Sexual Network Study of US Men Who Have Sex with Men. _Journal of the International AIDS Society._ 2021; 24(10): e25826. [[PubMed]](https://pubmed.ncbi.nlm.nih.gov/34605174/)
7. Chandra C, Weiss KM, Kelley CF, Marcus JL, Jenness SM. Gaps in Sexually Transmitted Infection Screening Among Men Who Have Sex with Men in Pre-exposure Prophylaxis (PrEP) Care in the United States. _Clinical Infectious Diseases._ 2021; 73(7): e2261-e2269. [[PubMed]](https://pubmed.ncbi.nlm.nih.gov/32702116/)
8. Anderson EJ, Weiss KM, Morris MM, Sanchez TH, Prasad P, Jenness SM. HIV and Sexually Transmitted Infection Epidemic Potential of Networks of Men Who Have Sex With Men in Two Cities. _Epidemiology._ 2021; 32(5): 681-689. [[PubMed]](https://pubmed.ncbi.nlm.nih.gov/34172692/) (Also fits temporal exponential random graph models to ARTnet and simulates two city networks; listed once, here.)
9. Weiss KM, Prasad P, Ramaraju R, Zlotorzynska M, Jenness SM. Estimated Number of Men Who Have Sex With Men With Indications for HIV Pre-exposure Prophylaxis in a National Sexual Network Study. _Journal of Acquired Immune Deficiency Syndromes._ 2020; 84(1): 10-17. [[PubMed]](https://pubmed.ncbi.nlm.nih.gov/31939869/)
10. Weiss KM, Goodreau SM, Morris M, Prasad P, Ramaraju R, Sanchez T, Jenness SM. Egocentric Sexual Networks of Men Who Have Sex with Men in the United States: Results from the ARTnet Study. _Epidemics._ 2020; 30: 100386. [[PubMed]](https://pubmed.ncbi.nlm.nih.gov/32004795/) (Canonical descriptive paper.)
11. Jenness SM, Weiss KM, Prasad P, Zlotorzynska M, Sanchez T. Bacterial Sexually Transmitted Infection Screening Rates by Symptomatic Status Among Men Who Have Sex With Men in the United States: A Hierarchical Bayesian Analysis. _Sexually Transmitted Diseases._ 2019; 46(1): 25-30. [[PubMed]](https://pubmed.ncbi.nlm.nih.gov/30044334/) (Uses the interim 2017-2018 wave, n = 2,572, and the hyphenated "ART-Net" spelling.)

---

#### Transmission Models Parameterized from ARTnet

Papers with no author from the study team carry the leading group in parentheses. Two papers with a Jenness coauthorship but a senior author outside the team are marked the same way.

1. Clay PA, Jackson DA, Bachmann LH, Spicknall IH. Doxycycline Post-Exposure Prophylaxis May Dramatically Reduce Syphilis Among GBMSM: A Modeling Study. _Sexually Transmitted Diseases._ 2026. [[PubMed]](https://pubmed.ncbi.nlm.nih.gov/42330441/) (CDC. Verified from the parameter table: the annual partner rates defining the three sexual activity groups, 5, 34, and 127 per year, are attributed to `artNET26`, where reference 26 is Weiss 2020.)
2. Prakhova S. Cost-Effectiveness of the Surveillance Strategy for Antimicrobial-Resistant Gonorrhea in the United States: A Modelling Study. _Venereology._ 2026; 5(1): 7. [[DOI]](https://doi.org/10.3390/venereology5010007) (Yale School of Public Health. Structural rather than parametric use: ARTnet justifies excluding main partnership formation from the model.)
3. Gutowska SJ, Hoffman KA, Gurski KF. Improving Adherence to a Daily PrEP Regimen Is Key When Considering Long-Time Partnerships. _Journal of Biological Dynamics._ 2024; 18(1): 2390843. [[PubMed]](https://pubmed.ncbi.nlm.nih.gov/39162356/) (Howard University / University of Maryland, Baltimore County. Verified from the parameter table: mean duration of long-term partnership 3.57 years and average concurrency probability 0.264, both attributed to Weiss 2020.)
4. Pollock ED, Clay PA, Keen A, Currie DW, Carter RJ, Quilter LAS, Gundlapalli AV, Mermin J, Spicknall IH. Potential for Recurrent Mpox Outbreaks Among Gay, Bisexual, and Other Men Who Have Sex with Men - United States, 2023. _MMWR Morbidity and Mortality Weekly Report._ 2023; 72(21): 568-573. [[PubMed]](https://pubmed.ncbi.nlm.nih.gov/37227964/) (CDC. Adapts the network layer of Spicknall 2022 and the Clay preprint, both ARTnet-parameterized; included on PI adjudication.)
5. Steinegger B, Iacopini I, Teixeira AS, Bracci A, Casanova-Ferrer P, Antonioni A, Valdano E. Non-Selective Distribution of Infectious Disease Prevention May Outperform Risk-Based Targeting. _Nature Communications._ 2022; 13(1): 3028. [[PubMed]](https://pubmed.ncbi.nlm.nih.gov/35641538/) (Universitat Rovira i Virgili, Central European University, Aix-Marseille, Universidade de Lisboa, City University of London, Sorbonne/INSERM. Single-statistic use: the ARTnet one-time partnership acquisition rate of 0.16 per week justifies the annealed-network approximation in Supplementary Note 2.)
6. Hendrickx DM, Delva W, Hens N. Influence of Sexual Risk Behaviour and STI Co-Infection Dynamics on the Evolution of HIV Set Point Viral Load in MSM. _Epidemics._ 2021; 36: 100474. [[PubMed]](https://pubmed.ncbi.nlm.nih.gov/34153622/) (Hasselt University, Belgium. Single-statistic use: the ARTnet concurrency estimate of 26.3% is one of two values setting the model's 10% to 30% calibration constraint.)
7. Chandra C, Maloney KM, Le Guillou A, Hoover KW, Jenness SM. Modeling the Potential Impact of Scaling Up Event-Driven PrEP Among Gay, Bisexual, and Other Men Who Have Sex with Men. _Journal of Acquired Immune Deficiency Syndromes._ 2026. (In press; no DOI or PMID assigned. ARTnet use documented in the `EpiModel/EventDrivenPrEP` repository.)
8. Le Guillou AV, Marcus JL, Krakower DS, Violette LR, Jenness SM. Potential Risks and Benefits of Low-Barrier Access and Monitoring for HIV Preexposure Prophylaxis: A Modeling Study. _Journal of Infectious Diseases._ 2026. [[PubMed]](https://pubmed.ncbi.nlm.nih.gov/42467629/)
9. Onwubiko UN, Le Guillou A, Chamberlain AT, Benkeser D, Holland DP, Baral SD, Jenness SM. Evaluating Mental Health Integration in PrEP Programs to Reduce HIV Incidence Among US Gay and Bisexual Men: A Network Modeling Study. _BMC Health Services Research._ 2026. [[PubMed]](https://pubmed.ncbi.nlm.nih.gov/42458454/)
10. Lewnard JA, Paredes MI, Yechezkel M, Davis GS, Hong V, Skela J, Pandey U, Parker NT, Granskog LC, Pomichowski ME, Reyes IAC, Rodriguez-Barraquer I, Müller NF, Tartof SY. Extensive Cryptic Circulation Sustains Mpox Among Men Who Have Sex with Men. _Nature Communications._ 2026; 17(1): 4198. [[PubMed]](https://pubmed.ncbi.nlm.nih.gov/42129165/) (University of California, Berkeley / Kaiser Permanente Southern California)
11. Clay PA, Pollock ED, Saldarriaga EM, Pathela P, Macaraig M, Zucker JR, Crouch B, Kracalik I, Spicknall IH. Modeling the Impact of Vaccine Dose Prioritization Strategies During the 2022 Mpox Outbreak. _American Journal of Epidemiology._ 2026; 195(4): 948-955. [[PubMed]](https://pubmed.ncbi.nlm.nih.gov/40069944/) (CDC / New York City Department of Health and Mental Hygiene)
12. Onwubiko UN, Le Guillou A, Lyles C, Maloney KM, Delaney KP, Jenness SM. Optimizing HIV Partner Services for Gay, Bisexual, and Other Men Who Have Sex With Men Previously Diagnosed With HIV: A Modeling Study. _Sexually Transmitted Diseases._ 2025; 52(8): 495-502. [[PubMed]](https://pubmed.ncbi.nlm.nih.gov/40489800/)
13. Crenshaw EG, Onnela JP. Modeling the 2022 Mpox Outbreak with a Mechanistic Network Model. _BMC Infectious Diseases._ 2025; 25(1): 1507. [[PubMed]](https://pubmed.ncbi.nlm.nih.gov/41194012/) (Harvard T.H. Chan School of Public Health)
14. Clay PA, Markowitz LE, Gopalani SV, Baxter A, Gargano JW, DeSisto CL, Senkomago V, Ekwueme DU, Islam MH, Chesson HW. HPV Vaccination Among Gay, Bisexual, and Other Men Who Have Sex with Men Aged 27-45 Years in the United States Is Potentially Cost-Saving. _Vaccine._ 2025; 65: 127798. [[PubMed]](https://pubmed.ncbi.nlm.nih.gov/41056755/) (CDC)
15. Parker AM, Chang JJ, Chen L, King LM, McCoy SI, Lewnard JA, Bruxvoort KJ. Bacterial Sexually Transmitted Infections and Related Antibiotic Use Among Individuals Eligible for Doxycycline Post-Exposure Prophylaxis in the United States. _Nature Communications._ 2025; 16(1): 9206. [[PubMed]](https://pubmed.ncbi.nlm.nih.gov/41102210/) (University of California, Berkeley / Kaiser Permanente Southern California)
16. Rönn MM, Chesson HW, Grad YH, Reitsma M, Zhu L, Hsu K, Gift TL, Salomon JA. The Influence of Epidemiologic Context on the Success of Partner Notification Programs: Analysis of Gonorrhea Transmission Dynamics. _Journal of Infectious Diseases._ 2025; 232(2): e266-e274. [[PubMed]](https://pubmed.ncbi.nlm.nih.gov/40251014/) (Harvard T.H. Chan School of Public Health / Stanford University)
17. Hamilton DT, Hoover KW, Delaney KP. Reducing HIV Incidence in the Southern United States Through Routine Opt-Out HIV Screening. _AIDS._ 2025; 39(11): 1632-1640. [[PubMed]](https://pubmed.ncbi.nlm.nih.gov/40372031/) (University of Washington / CDC)
18. Fraysse J, Anderson SJ, Smith JC, Matthews DD, Sarkar S, de Aragao F, Blissett R. Achieving the State of Georgia 25% HIV Incidence Reduction Target Among Men Who Have Sex with Men in Atlanta Through Expanded Use of Multimodal Pre-exposure Prophylaxis: A Mathematical Model. _PLoS ONE._ 2025; 20(1): e0312369. [[PubMed]](https://pubmed.ncbi.nlm.nih.gov/39787101/) (ViiV Healthcare / Maple Health Group; inherits the ARTnet network layer through Maloney 2021 and does not cite Weiss 2020)
19. Clay PA, Asher JM, Carnes N, Copen CE, Delaney KP, Payne DC, Pollock ED, Mermin J, Nakazawa Y, Still W, Mangla AT, Spicknall IH. Modelling the Impact of Vaccination and Sexual Behaviour Adaptations on Mpox Cases in the USA During the 2022 Outbreak. _Sexually Transmitted Infections._ 2024; 100(2): 70-76. [[PubMed]](https://pubmed.ncbi.nlm.nih.gov/38050171/) (CDC / DC Health)
20. Kipshidze N, Klein E, Yang W. Understanding the Drivers of Continued Mpox Transmission in the United States: A Modeling Study. _Research Square preprint._ 2024. [[DOI]](https://doi.org/10.21203/rs.3.rs-3817998/v1) (One Health Trust / Columbia University; not indexed as published as of August 2026)
21. Jones J, Jenness SM, Le Guillou A, Sullivan PS, Gift TL, Delaney KP, Chesson H. Estimated Number of Incident HIV Infections in Men Who Have Sex With Men Attributable to Gonorrhea and Chlamydia, Per Gonococcal or Chlamydial Infection, in the United States. _Sexually Transmitted Diseases._ 2023; 50(2): 83-85. [[PubMed]](https://pubmed.ncbi.nlm.nih.gov/36630415/) (Emory / CDC. ARTnet is documented only in the Technical Appendix at `links.lww.com/OLQ/A873`, which carries a section headed "The ARTnet Study"; the research letter itself does not name ARTnet.)
22. Hamilton DT, Wang LY, Hoover KW, Smith DK, Delaney KP, Li J, Hoyte T, Jenness SM, Goodreau SM. Potential Contribution of PrEP Uptake by Adolescents 15-17 Years Old to Achieving the "Ending the HIV Epidemic" Incidence Reduction Goals in the US South. _PLoS ONE._ 2023; 18(11): e0288588. [[PubMed]](https://pubmed.ncbi.nlm.nih.gov/37943869/) (University of Washington)
23. Hamilton DT, Hoover KW, Smith DK, Delaney KP, Wang LY, Li J, Hoyte T, Jenness SM, Goodreau SM. Achieving the "Ending the HIV Epidemic in the U.S." Incidence Reduction Goals Among At-Risk Populations in the South. _BMC Public Health._ 2023; 23(1): 716. [[PubMed]](https://pubmed.ncbi.nlm.nih.gov/37081482/) (University of Washington)
24. Gurski K, Hoffman K. Staged HIV Transmission and Treatment in a Dynamic Model with Long-Term Partnerships. _Journal of Mathematical Biology._ 2023; 86(5): 74. [[PubMed]](https://pubmed.ncbi.nlm.nih.gov/37052718/) (Howard University / University of Maryland, Baltimore County)
25. Jenness SM, Le Guillou A, Lyles C, Bernstein KT, Krupinsky K, Enns EA, Sullivan PS, Delaney KP. The Role of HIV Partner Services in the Modern Biomedical HIV Prevention Era: A Network Modeling Study. _Sexually Transmitted Diseases._ 2022; 49(12): 801-807. [[PubMed]](https://pubmed.ncbi.nlm.nih.gov/36194831/)
26. Jones J, Le Guillou A, Gift TL, Chesson H, Bernstein KT, Delaney KP, Lyles C, Berruti A, Sullivan PS, Jenness SM. Effect of Screening and Treatment for Gonorrhea and Chlamydia on HIV Incidence Among Men Who Have Sex With Men in the United States: A Modeling Analysis. _Sexually Transmitted Diseases._ 2022; 49(10): 669-676. [[PubMed]](https://pubmed.ncbi.nlm.nih.gov/35921635/)
27. Spicknall IH, Pollock ED, Clay PA, Oster AM, Charniga K, Masters N, Nakazawa YJ, Rainisch G, Gundlapalli AV, Gift TL. Modeling the Impact of Sexual Networks in the Transmission of Monkeypox Virus Among Gay, Bisexual, and Other Men Who Have Sex with Men - United States, 2022. _MMWR Morbidity and Mortality Weekly Report._ 2022; 71(35): 1131-1135. [[PubMed]](https://pubmed.ncbi.nlm.nih.gov/36048619/) (CDC. The body text mislabels the source as Atlanta data collected 2016-2019; ARTnet was national and fielded 2017-2019. Do not quote that sentence as a description of ARTnet.)
28. Wheatley MM, Knowlton G, Kao SY, Jenness SM, Enns EA. Cost-Effectiveness of Interventions to Improve HIV Pre-exposure Prophylaxis Initiation, Adherence, and Persistence Among Men Who Have Sex With Men. _Journal of Acquired Immune Deficiency Syndromes._ 2022; 90(1): 41-49. [[PubMed]](https://pubmed.ncbi.nlm.nih.gov/35090155/) (University of Minnesota. ARTnet is documented only in Supplemental Appendix S1-10 at `links.lww.com/QAI/B800`; the main text does not name ARTnet and does not cite Weiss 2020.)
29. Le Guillou A, Buchbinder S, Scott H, Liu A, Havlir D, Scheer S, Jenness SM. Population Impact and Efficiency of Improvements to HIV PrEP Under Conditions of High ART Coverage Among San Francisco Men Who Have Sex With Men. _Journal of Acquired Immune Deficiency Syndromes._ 2021; 88(4): 340-347. [[PubMed]](https://pubmed.ncbi.nlm.nih.gov/34354011/)
30. Jenness SM, Knowlton G, Smith DK, Marcus JL, Anderson EJ, Siegler AJ, Jones J, Sullivan PS, Enns E. A Decision Analytics Model to Optimize Investment in Interventions Targeting the HIV Preexposure Prophylaxis Cascade of Care. _AIDS._ 2021; 35(9): 1479-1489. [[PubMed]](https://pubmed.ncbi.nlm.nih.gov/33831910/)
31. Jenness SM, Le Guillou A, Chandra C, Mann LM, Sanchez T, Westreich D, Marcus JL. Projected HIV and Bacterial Sexually Transmitted Infection Incidence Following COVID-19-Related Sexual Distancing and Clinical Service Interruption. _Journal of Infectious Diseases._ 2021; 223(6): 1019-1028. [[PubMed]](https://pubmed.ncbi.nlm.nih.gov/33507308/)
32. Maloney KM, Le Guillou A, Driggers RA, Sarkar S, Anderson EJ, Malik AA, Jenness SM. Projected Impact of Concurrently Available Long-Acting Injectable and Daily-Oral Human Immunodeficiency Virus Preexposure Prophylaxis: A Mathematical Model. _Journal of Infectious Diseases._ 2021; 223(1): 72-82. [[PubMed]](https://pubmed.ncbi.nlm.nih.gov/32882043/)
33. Jenness SM, Johnson JA, Hoover KW, Smith DK, Delaney KP. Modeling an Integrated HIV Prevention and Care Continuum to Achieve the Ending the HIV Epidemic Goals. _AIDS._ 2020; 34(14): 2103-2113. [[PubMed]](https://pubmed.ncbi.nlm.nih.gov/32910062/)

---

#### Statistical and Software Methods Papers Using ARTnet

Three of the four fit exponential random graph models to the ARTnet limited-use dataset obtained under Memorandum of Understanding, and none of those three has an author with an epidemiologic connection to Emory.

1. Klumb C, Morris M, Goodreau SM, Jenness SM. Improving and Extending STERGM Approximations Based on Cross-Sectional Data and Tie Durations. _Journal of Computational and Graphical Statistics._ 2024; 33(1): 166-180. [[PubMed]](https://pubmed.ncbi.nlm.nih.gov/38455738/) (University of Washington / Emory. Included on PI adjudication: ARTnet-derived degree and duration values set the simulation study's target ranges.)
8. Hod S, Nayak D, Gantenberg JR, Kalemaj I, Trikalinos TA, Smith A. Differentially Private Modeling of Disease Transmission Within Human Contact Networks. _arXiv preprint 2604.07493._ 2026. [[DOI]](https://doi.org/10.48550/arXiv.2604.07493) (Boston University / Brown University. Acknowledgements: "We thank Dr. Samuel M. Jenness for the use of ARTNet data, the egocentric network dataset upon which our generative network models rely.")
9. Brázia J, Kiss IZ, Francisco AP, Teixeira AS. Reconstructing MSM Sexual Networks to Guide PrEP Distribution Strategies for HIV Prevention. _arXiv preprint 2601.04434._ 2026. [[DOI]](https://doi.org/10.48550/arXiv.2601.04434) (Universidade de Lisboa / Northeastern University London. Individual-level reanalysis of 4,667 of the 4,904 ARTnet egos; the data availability statement reproduces the MOU process and the PI's contact address.)
10. Liu Y, Rahimian MA, Yu FY. Seeding with Differentially Private Network Information. _arXiv preprint 2305.16590._ 2023, revised 2026. [[DOI]](https://doi.org/10.48550/arXiv.2305.16590) (University of Pittsburgh / George Mason University. A preliminary version appeared in the Proceedings of AAMAS 2023. The paper reports 4,909 ARTnet participants; the correct figure is 4,904.)

---

#### Use by Groups Outside the Study Team

Twenty-one of the 48 papers have no author from the study team. They come from **14 distinct research groups at 26 institutions in 7 countries** (United States, Portugal, United Kingdom, Belgium, Spain, France, Austria).

| Group or institution | Papers | Subject area |
|---|---|---|
| CDC (Divisions of STD Prevention and HIV Prevention; Center for Forecasting and Outbreak Analytics), with NYC DOHMH and DC Health | 6 | Mpox transmission, vaccination, dose prioritization, and outbreak recurrence (4); HPV vaccination cost-effectiveness (1); doxycycline PEP and syphilis (1) |
| University of California, Berkeley / Kaiser Permanente Southern California | 2 | Mpox cryptic circulation; doxycycline PEP and antibiotic use |
| Harvard T.H. Chan School of Public Health (Onnela group) | 1 | Mpox mechanistic network model |
| Harvard T.H. Chan School of Public Health / Stanford University | 1 | Gonorrhea partner notification |
| University of Washington, with CDC | 1 | Routine opt-out HIV screening in the US South |
| Howard University / University of Maryland, Baltimore County | 2 | HIV pair-formation ordinary differential equation models; PrEP adherence in long-term partnerships |
| Yale School of Public Health | 1 | Cost-effectiveness of antimicrobial-resistant gonorrhea surveillance |
| Hasselt University, Belgium | 1 | HIV set point viral load under STI co-infection |
| Universitat Rovira i Virgili, Central European University, Aix-Marseille, Universidade de Lisboa, City University of London, Sorbonne/INSERM | 1 | Non-selective versus risk-based prevention distribution |
| ViiV Healthcare / Maple Health Group (industry) | 1 | Multimodal PrEP scale-up against a state incidence target |
| One Health Trust / Columbia University | 1 | Drivers of continued mpox transmission (preprint) |
| Boston University / Brown University | 1 | Differential privacy for network epidemic models |
| University of Pittsburgh / George Mason University | 1 | Differentially private seeding on sexual networks |
| Universidade de Lisboa / Northeastern University London | 1 | Sexual network reconstruction and PrEP allocation |

ARTnet is used well beyond the setting it was collected for. It parameterizes published models of pathogens other than HIV, including mpox (six papers), human papillomavirus, gonorrhea and antimicrobial-resistant gonorrhea, syphilis, and doxycycline post-exposure prophylaxis. It is also used outside the EpiModel software stack, in a compartmental susceptible-infectious-susceptible model of HPV (Clay 2025), a pair-formation ordinary differential equation model (Gurski 2023), and a partner notification model (Ronn 2025). Three computer science groups with no epidemiologic connection to Emory have obtained the limited-use dataset under MOU for work on differential privacy and network reconstruction.

This list is assembled from citation and full-text searching rather than from acknowledgements, because the acknowledgement wording requested in the MOU is rarely used in practice. Some downstream uses cite an intermediate modeling paper rather than Weiss et al. 2020 and are therefore missed by citation searching alone, so the list should be read as a lower bound.
