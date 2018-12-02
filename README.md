# PAM50

PAM50 package is a trained classifier that predicts PAM50 subtypes using patient's own gene expression data (NanoString, microarray, or RNA-seq). 
It can predict 5 subtypes for breast cancer patients: **Luminal A**, **Luminal B**, **HER2**, **Basal**, and **Normal**. 
Patients with Luminal A subtype usually have the best prognosis. Click [here](https://ww5.komen.org/BreastCancer/SubtypesofBreastCancer.html) for more detailed description of each subtype.

gPAM50 uses 50 genes to make the prediction. As long as a patient has those 50 genes' expression data and Gene ID (Entrez ID), 
you can then predict PAM50 subtype for that patient. **Yes, it works for a single patient's data!**
Interestingly, even the patient is not breast cancer patient, gPAM50 can still produce a prediction which may or may not be meaningful. 
But why don't you give it a try? Maybe you will discovery novel pathway or genes. 


Project page: https://ccchang0111.github.io/PAM50/

Source code: https://github.com/ccchang0111/PAM50

## Installation

You can install gPAM50 from github with:

``` r
# install.packages("devtools")
devtools::install_github("ccchang0111/gPAM50")
```

## Tutorial

Magics can be found here: https://ccchang0111.github.io/PAM50/articles/PAM50_tutorial.html

