library(tidyverse)
library(CombMSC)
library(bayesMetaintegrator)

#import dataset
discovery <- readRDS("./CardioFullRun.RDS")
discovery$bayesianMeta$finalResults <- NULL
#output path
outPath = "./Results/Cardio/cardio_"
#set combinations
comboNumber = 100
#limit permutations
sizeLimit = ceiling(length(discovery$originalData)/2)

for(simulationSize in seq(1, 12)) {
  #first, generate all combinations - simulationSize
  datasetList <- unique(discovery$bayesianMeta$datasetResults$name)
  allCombinations <- as.data.frame(CombMSC::subsets(length(datasetList), simulationSize, datasetList))
  if (nrow(allCombinations) >= comboNumber) {
    subsetToRun <- allCombinations[sample(nrow(allCombinations), comboNumber), ]
  } else {
    subsetToRun <- allCombinations
  }
  rownames(subsetToRun) <- paste0("Sample_", seq(1:nrow(subsetToRun)))
  originalDiscovery <- discovery
  originalDiscovery$monteCarlo$combinationsToTest <- subsetToRun

  for(i in seq(1, nrow(originalDiscovery$monteCarlo$combinationsToTest))) {
    combinationToTest <- as.character(unname(unlist(originalDiscovery$monteCarlo$combinationsToTest[i,])))
    test <- discovery
    test$bayesianMeta$datasetResults <- test$bayesianMeta$datasetResults %>% filter(name %in% combinationToTest)
    test$bayesianMeta$datasetEffectSizes <- test$bayesianMeta$datasetEffectSizes[,c(combinationToTest, "Gene", "nStudies")]
    test$bayesianMeta$datasetEffectSizes$nStudies <- apply(test$bayesianMeta$datasetEffectSizes, 1, function(x) sum(!is.na(x))-2)
    test$bayesianMeta$datasetStdDev <-  test$bayesianMeta$datasetStdDev[,c(combinationToTest, "Gene", "nStudies")]
    test$bayesianMeta$datasetStdDev$nStudies <- apply(test$bayesianMeta$datasetStdDev, 1, function(x) sum(!is.na(x))-2)

    ptm <- proc.time()
    test$bayesianMeta$finalResults <- NULL
    test$bayesianMeta$pooledResultsDetailed <- NULL
    test <- combineDatasets(test)
    print(proc.time() - ptm)

    samplingName <- paste0("Sample_", i)
    originalDiscovery$monteCarlo$samplings[[samplingName]] <- test$bayesianMeta$finalResults

    #now do regular meta-analysis
    testMeta <- list()
    for(j in seq(1, length(discovery$originalData))) {
      if(discovery$originalData[[j]]$formattedName %in% combinationToTest) {
        testMeta$originalData[[discovery$originalData[[j]]$formattedName]] <- discovery$originalData[[j]]
      }
    }
    testMeta <- MetaIntegrator::runMetaAnalysis(testMeta)
    originalDiscovery$monteCarlo$regularMeta[[samplingName]] <- testMeta$metaAnalysis
    sinkFile <- paste0(outPath, "log_file.txt")
    log = paste0(i, " ", paste(combinationToTest,collapse=" "), "\n")
    cat(log,file=sinkFile,append=TRUE)
    print(i)
  }
  outputFilename = paste0(outPath,simulationSize, ".RDS")
  saveRDS(originalDiscovery, outputFilename)
}
