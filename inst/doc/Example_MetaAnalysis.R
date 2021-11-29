## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
knitr::opts_chunk$set(fig.width=10, fig.height=6) 



## ----setup--------------------------------------------------------------------
library(bayesMetaintegrator)
#load practice data 
data("test_discovery")

## -----------------------------------------------------------------------------
library(MetaIntegrator)
MetaIntegrator::checkDataObject(test_discovery, "Meta", "Pre-Analysis")

## -----------------------------------------------------------------------------
geneSubsetList <- c("BRCA1", "TBC1D2", "POSTN", "SERPINB2", "PRR4", "CEACAM5", "KRT6A", "SCGB3A1", "C3", "MUC5B", "VEGFA", "SLC6A4", "FANCI")
test_discovery <- probeToGene(test_discovery, geneSubsetList)

## -----------------------------------------------------------------------------
test_discovery <- getDatasetEffectSizes(discovery = test_discovery, cores = 2, steps = 2000)
test_discovery <- filterRhat(test_discovery)
head(test_discovery$bayesianMeta$datasetResults)

## -----------------------------------------------------------------------------
test_discovery <- combineDatasets(discovery = test_discovery, steps = 5000, burnIn = 1000)
head(test_discovery$bayesianMeta$finalResults)

## -----------------------------------------------------------------------------
#this function returns the JAGs model for the hierarchical model
model <- getCombiningModel(discovery = test_discovery, gene = "BRCA1", steps = 5000, burnIn = 1000)
#then can plot the trace (or use other BayesPlot functions)
getCombineTrace(model, parameters = c("MU", "tau"), trace = TRUE)

## -----------------------------------------------------------------------------
test_discovery <- bayesMetaintegrator::filterGenes(test_discovery, probabilityThresh = .05, effectSizeThresh = .5, numberStudiesThresh = 2)

#and now run the MetaIntegrator functions with the resulting filter objects.
MetaIntegrator::multipleROCPlot(test_discovery, test_discovery$filterResults$bayesian_PR_0.05ES_0.5nStudies_2)

## -----------------------------------------------------------------------------
testPlot <- makeBayesianForestPlot(test_discovery, "BRCA1", plotType = "frequentist")

