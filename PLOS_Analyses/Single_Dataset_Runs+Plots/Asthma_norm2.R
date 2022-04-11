library(tidyverse)
library(MetaIntegrator)
library(bayesMetaintegrator)
inputFile = "./AsthmaFullRun.RDS"
outputFile = "./Asthma_norm2.RDS"


discovery <- readRDS(inputFile)
discovery$metaAnalysis <- NULL
discovery$leaveOneOutAnalysis <- NULL
discovery$bayesianMeta <- NULL

#check MetaIntegrator object
print(MetaIntegrator::checkDataObject(discovery, "Meta", "Pre-Analysis"))
discovery <- runMetaAnalysis(discovery)
numStudiesThresh = ceiling(length(discovery$originalData)/2)
discovery$metaAnalysis$pooledResults <- discovery$metaAnalysis$pooledResults %>% filter(numStudies >= numStudiesThresh)

geneSubsetList <-rownames(discovery$metaAnalysis$pooledResults)[discovery$metaAnalysis$pooledResults$numStudies >= numStudiesThresh]
discovery <- probeToGene(discovery, geneSubsetList)
ptm <- proc.time()
invisible(discovery <- getDatasetEffectSizes(discovery = discovery, cores = 20, steps = 2000, burnInSteps = 400, seedNum = 42))
proc.time() - ptm

model <- function() {
  for (i in 1:I) {
    y[i, 1] ~ dnorm(mu[i], pow(se[i, 1], -2))
    mu[i] ~ dnorm(MU, prec.mu)
  }
  MU ~ dnorm(0, pow(2, -2))
  prec.mu <- pow(tau, -2)
  tau ~ dunif(0, 2)
}
jags.data <- list("y", "se", "I")
jags.params <- c("tau", "mu", "MU")
jags.inits <- function() {
  list(tau = runif(1), mu = rnorm(I))
}

discovery <- combineDatasets(discovery, seedNum = 42, definedModel = model)
saveRDS(discovery, outputFile)
