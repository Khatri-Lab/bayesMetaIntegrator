library(tidyverse)
library(MetaIntegrator)
library(bayesMetaintegrator)
##insert initial MetaIntegrator Object Here

inputFile = ""
outputFile = ""

discovery <- readRDS(inputFile)
#check MetaIntegrator object
print(MetaIntegrator::checkDataObject(discovery, "Meta", "Pre-Analysis"))
discovery <- runMetaAnalysis(discovery)
numStudiesThresh = ceiling(length(discovery$originalData)/2)
discovery$metaAnalysis$pooledResults <- discovery$metaAnalysis$pooledResults %>% filter(numStudies >= numStudiesThresh)

geneSubsetList <- rownames(discovery$metaAnalysis$pooledResults)[discovery$metaAnalysis$pooledResults$numStudies >= numStudiesThresh]
discovery <- probeToGene(discovery, geneSubsetList)
ptm <- proc.time()
invisible(discovery <- getDatasetEffectSizes(discovery = discovery, cores = 30, steps = 2000, burnInSteps = 400, seedNum = 42))
proc.time() - ptm
discovery <- combineDatasets(discovery, seedNum = 42)
saveRDS(discovery, outputFile)
