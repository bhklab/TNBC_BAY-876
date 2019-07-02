# TNBC_BAY876

## Introduction
This is a repo that contains the code used to do the bioinfomatics analyses in this paper:

Inhibition of glucose transporter I as a therapy for RB1-overexpressing triple negative breast cancer, 2019, Wu et. al..

The script for the analyses is written in R.


----

## Dependencies:
These R packages need to be installed in order to run the analysis script:
- Biobase
- PharmacoGx
- ggplot2
- GSA
- piano
- fgsea
- dplyr
- EnrichmentBrowser
- ggrepel


----
## Reproducibility of the Analysis:
- Once the project is downloaded to the user computer, the user needs to navigate to the main directory of the project "TNBC_BAY-876-master".
- Inside the main directory, there is an R script file named "run_analysis.R". Running this script will regenerate the different bioinformatics plots used in the paper

**Important Note 1:** the user needs to set the working directory inside the script file before running it, i.e. setting the working directory to "TNBC_BAY-876-master"
**Important Note 2:** the user needs to download all the data needed for this code to work from "figshare url" and put them inside the "data" folder of this repo.


