# bayesMetaIntegrator
By: Laurynas Kalesinskas

## Bayesian Meta-Analysis for Gene Expression Data
The method is described in [LINK](). It was originally developed for gene expression data, but can be adapted for any type of high dimensional meta-analysis.

## Installation
To install bayesMetaIntegrator, first you must install the JAGs executable for your system. This can be found [here](https://mcmc-jags.sourceforge.io/).  

Following JAGS installation, the bayesMetaIntegrator package can be installed and loaded using devtools:
```
install.packages("devtools")
```
If devtools are already installed, or after devtools installation please run:
```
devtools::install_github('khatrilab/bayesMetaIntegrator', build = TRUE)
```
To get full use of bayesMetaIntegrator, it is recommended you install MetaIntegrator as well. This package can be found on [CRAN](https://cran.r-project.org/web/packages/MetaIntegrator/index.html).

## Vignettes
A vignette with a demo data is included in the package. To open the vignette from R please use:
```
vignette(package = "bayesMetaintegrator", topic = "Example_MetaAnalysis")
```
