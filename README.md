# TNBC_BAY876

## Introduction
This is a repo that contains the code used to do the bioinfomatics analyses in this paper:

GLUT1 inhibition blocks growth of RB1-positive Triple Negative Breast Cancer, 2019, Wu et. al..

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
- EnrichmentBrowser (devtools::install_github("lgeistlinger/EnrichmentBrowser")
- ggrepel


----
## Reproducibility of the Analysis:
- Once the project is downloaded to the user computer, the user needs to navigate to the main directory of the project "TNBC_BAY-876-master".
- Inside the main directory, there is an R script file named "run_analysis.R". Running this script will regenerate the different bioinformatics plots used in the paper

<br>
*Important Note 1:* the user needs to set the working directory inside the script file before running it, i.e. setting the working directory to "TNBC_BAY-876-master" or double-click on TNBC_GLUT1.Rproj and it will set it automatically in RStudio
<br>
<br>
*Important Note 2:* the user needs to download all the data needed for this code to work from "https://figshare.com/projects/TNBC_BAY876/65828" and put them inside the "data" folder of this repo. Sensitivity data can be downloaded from "https://figshare.com/s/6e1782e040059d5a34e0" and will be made public in the figshare project once the paper is published.



