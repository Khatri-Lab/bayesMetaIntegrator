#' Function that calculates effect size estimate using all cases and controls in a dataset
#'
#' @param dataset MetaIntegrator dataset with gene_expr populated
#' @param steps The number of MCMC steps to run for each gene
#' @param burnInSteps Number of steps for burn-in - this should be less than steps
#' @param seedNum Number for random seed
#' @return Dataset with datasetResults populated
#'
getSingleDatasetES <- function(dataset, steps = 1000, burnInSteps = 200, seedNum = 42) {
  #get expr_set
  print(dataset$formattedName)
  expr_set <- dataset$gene_expr
  classes <- stack(dataset$class)$values
  geneNames <- expr_set$gene
  expr_set$gene <- NULL
  dataset$datasetResults = data.frame()
  if(burnInSteps > steps) {
    stop("MCMC Burn-in steps > steps - please adjust so that burn-in steps < steps.")
  }

  #multi-threaded MCMC - get estimate of effect size per each dataset
  dataset$datasetResults <- foreach(i=1:nrow(expr_set),.combine=rbind, .inorder = FALSE, .packages=c('BEST')) %dopar% {
    cases <- stack(expr_set[i, classes == 1])$values
    controls <- stack(expr_set[i, classes == 0])$values
    a <- tryCatch({
      BESTout <- BESTmcmc(cases, controls, priors=NULL, parallel=FALSE, numSavedSteps = steps, burnInSteps = burnInSteps, rnd.seed = seedNum)
      bestOutput <- summary(BESTout)
      data.frame(Gene = geneNames[i], effSize = bestOutput["effSz", 1], HDIlo = bestOutput["effSz", "HDIlo"],
                 rhatMeanCases = attr(BESTout, "Rhat")[1], rhatMeanControl = attr(BESTout, "Rhat")[2],
                 rhatSigmaCases = attr(BESTout, "Rhat")[4], rhatSigmaControl = attr(BESTout, "Rhat")[5])
    }, error = function(err) {
      data.frame(Gene = geneNames[i], effSize = -100, HDIlo = -100, rhatMeanCases = 0, rhatMeanControl = 0,
                 rhatSigmaCases = 0, rhatSigmaControl = 0)
    })
    return(a)
  }
  #find problematic genes (ones that failed analysis and flag)
  problemGenes <- dataset$datasetResults %>% filter(effSize == -100)
  if (nrow(problemGenes) > 0) {
    problemGenes <- dataset$datasetResults %>% filter(effSize == -100) %>% dplyr::select(Gene)
    problems <- paste(c("Problematic genes found and removed in", dataset$formattedName," (can be found in dataset$flaggedGenes): ", as.vector(problemGenes$Gene)), collapse=" ")
    warning(problems)
    dataset$flaggedGenes <- as.vector(problemGenes$Gene)
  }

  #combine nicely and stick in $bayesianMeta$datasetResults
  dataset$datasetResults <- dataset$datasetResults %>% filter(!effSize == -100)
  colnames(dataset$datasetResults) <- c('Gene', 'effSize', 'HDIlo', 'rhatMeanCases', 'rhatMeanControls', 'rhatSigmaCases', 'rhatSigmaControl')
  dataset$datasetResults$stdDev <- (dataset$datasetResults$effSize - dataset$datasetResults$HDIlo)/1.96
  dataset$datasetResults$HDIlo <- NULL
  dataset$datasetResults$name <- dataset$formattedName
  rownames(dataset$datasetResults) <- dataset$datasetResults$Gene
  dataset$datasetResults <- dataset$datasetResults[, c(1, 2, 7, 3, 4, 5, 6, 8)]
  return(dataset)
}

#' Wrapper function that passes whole collection of datasets to calculate effect size distribution per dataset
#'
#' @param discovery MetaIntegrator object with gene_expr populated
#' @param steps The number of MCMC steps to run for each gene
#' @param burnInSteps Number of steps for burn-in - this should be less than steps
#' @param cores The number of cores allowed to be run
#' @param seedNum Number for random seed
#' @return MetaIntegrator object with $bayesianMeta$datasetEffectSizes and $bayesianData$originalData$dataset$datasetResults populated
#'
getDatasetEffectSizes <- function(discovery, cores = 2, steps = 1000, burnInSteps = 200, seedNum = 42) {
  #set up cluster
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  getDoParWorkers()
  discovery$bayesianMeta$originalData <- lapply(discovery$bayesianMeta$originalData, getSingleDatasetES, steps, burnInSteps, seedNum)
  stopCluster(cl)
  tempFrame <- data.frame()
  for (dataset in discovery$bayesianMeta$originalData) {
    tempFrame <- rbind(tempFrame, dataset$datasetResults)
  }

  #put things in a nice format and create datasetEffectSizes and std deviation dataframes
  discovery$bayesianMeta$datasetResults <- tempFrame
  discovery$bayesianMeta$datasetEffectSizes <-
    tempFrame %>%
    dplyr::select(Gene, effSize, name) %>%
    pivot_wider(names_from = name, values_from = effSize)
  discovery$bayesianMeta$datasetEffectSizes$nStudies <- apply(discovery$bayesianMeta$datasetEffectSizes, 1, function(x) ncol(discovery$bayesianMeta$datasetEffectSizes) - 1 - sum(is.na(x)))

  discovery$bayesianMeta$datasetStdDev <-
    tempFrame %>%
    dplyr::select(Gene, stdDev, name) %>%
    pivot_wider(names_from = name, values_from = stdDev)
  discovery$bayesianMeta$datasetStdDev$nStudies <- apply(discovery$bayesianMeta$datasetStdDev, 1, function(x) ncol(discovery$bayesianMeta$datasetStdDev) - 1 - sum(is.na(x)))
  rownames(discovery$bayesianMeta$datasetResults) <- paste0(discovery$bayesianMeta$datasetResults$Gene, "_", discovery$bayesianMeta$datasetResults$name)
  #add warning for rHats > 1.1
  if(any(discovery$bayesianMeta$datasetResults$rhatMeanCases > 1.1) || any(discovery$bayesianMeta$datasetResults$rhatMeanControls > 1.1)) {
    problemCount <- sum(discovery$bayesianMeta$datasetResults$rhatMeanCases > 1.1) + sum(discovery$bayesianMeta$datasetResults$rhatMeanControls > 1.1)
    percentage <- (problemCount/nrow(discovery$bayesianMeta$datasetResults)) * 100
    warningString <- paste0("Warning: Some mean rhats > 1.1 (", percentage, "%). Please check $bayesianMeta$datasetResults and consider running with more steps.")
    warning(warningString)
  }
  if(any(discovery$bayesianMeta$datasetResults$rhatSigmaCases > 1.1) || any(discovery$bayesianMeta$datasetResults$rhatSigmaControl > 1.1)) {
    problemCount <- sum(discovery$bayesianMeta$datasetResults$rhatSigmaCases > 1.1) + sum(discovery$bayesianMeta$datasetResults$rhatSigmaControl > 1.1)
    percentage <- (problemCount/nrow(discovery$bayesianMeta$datasetResults)) * 100
    warningString <- paste0("Warning: Some sigma rhats > 1.1 (", percentage, "%). Please check $bayesianMeta$datasetResults and consider running with more steps.")
    warning(warningString)
  }
  return(discovery)
}

#' Function that filters out all genes with inidividal dataset Rhats over a certain level
#'
#' @param discovery MetaIntegrator object after getSingleDatasetEffectSizes
#' @param rhatThresh Your rhat threshold
#' @return MetaIntegrator object with $bayesianMeta$filteredSingleDatasetGenes - describing which genes are removed
#'
filterRhat <- function(discovery, rhatThresh = 1.1){
  #get gene names that don't pass
  problemRows <- discovery$bayesianMeta$datasetResults %>% filter(rhatMeanCases > rhatThresh | rhatMeanControls > rhatThresh |
                                                                    rhatSigmaCases > rhatThresh | rhatSigmaControl > rhatThresh)
  discovery$bayesianMeta$datasetResults <- discovery$bayesianMeta$datasetResults %>% filter(!(rhatMeanCases > rhatThresh | rhatMeanControls > rhatThresh |
                                                                                              rhatSigmaCases > rhatThresh | rhatSigmaControl > rhatThresh))
  #now go through problem row pairs and remove them from std. dev and ES dfs
  keysToRemove <- problemRows %>% select("Gene", "name")
  for(rowNum in seq(1, nrow(keysToRemove))) {
    gene <- as.character(keysToRemove[rowNum,1])
    dataset <- as.character(keysToRemove[rowNum,2])
    discovery$bayesianMeta$datasetStdDev[discovery$bayesianMeta$datasetStdDev$Gene == gene, dataset] <- NA
    discovery$bayesianMeta$datasetEffectSizes[discovery$bayesianMeta$datasetStdDev$Gene == gene, dataset] <- NA
  }
  discovery$bayesianMeta$filteredSingleDatasetGenes <- keysToRemove
  return(discovery)
}
