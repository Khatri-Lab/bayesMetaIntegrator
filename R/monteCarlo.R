#' Function that runs Bayesian Monte Carlo Sampling
#'
#' @param discovery Bayesian Meta-Analysis object with basic bayesian meta-analysis run
#' @param numSamplings Number of samplings for monte carlo analysis
#' @param sampleSize Size of the samplings for monte carlo
#' @param steps The number of steps to combine for each sampling
#' @param burnInSteps Number of steps for burn-in for each sampling - this should be less than steps
#' @param seedNum Number for random seed
#' @param model Model for Bayesian combination - default, basic model
#' @return MetaIntegrator object with monteCarlo populated with samplings
#'
bayesianMonteCarloSimulation <- function(discovery, numSamplings, sampleSize, steps = 5000, burnInSteps = 1000, seed = 42, model = NULL) {
  #get dataset names
  datasetList <- unique(discovery$bayesianMeta$datasetResults$name) #note - we are using formattedNames
  allCombinations <- as.data.frame(CombMSC::subsets(length(datasetList), sampleSize, datasetList))
  if(numSamplings > nrow(allCombinations)) {
    stop("There are not enough unique samplings with the parameters selected.")
  }
  subsetToRun <- allCombinations[sample(nrow(allCombinations), numSamplings), ]
  rownames(subsetToRun) <- paste0("Sample_", seq(1:nrow(subsetToRun)))

  originalDiscovery <- discovery
  originalDiscovery$monteCarlo$combinationsToTest <- subsetToRun
  for(i in seq(1, nrow(originalDiscovery$monteCarlo$combinationsToTest))) {
    combinationToTest <- as.character(unname(unlist(subsetToRun[i,])))
    test <- discovery
    test$bayesianMeta$datasetResults <- test$bayesianMeta$datasetResults %>% filter(name %in% combinationToTest)
    test$bayesianMeta$datasetEffectSizes <- test$bayesianMeta$datasetEffectSizes[,c(combinationToTest, "Gene", "nStudies")]
    test$bayesianMeta$datasetEffectSizes$nStudies <- apply(test$bayesianMeta$datasetEffectSizes, 1, function(x) sum(!is.na(x))-2)

    test$bayesianMeta$datasetStdDev <-  test$bayesianMeta$datasetStdDev[,c(combinationToTest, "Gene", "nStudies")]
    test$bayesianMeta$datasetStdDev$nStudies <- apply(test$bayesianMeta$datasetStdDev, 1, function(x) sum(!is.na(x))-2)

    test$bayesianMeta$pooledResultsDetailed <- NULL
    test$bayesianMeta$finalResults <- NULL
    ptm <- proc.time()
    if(!is.null(model)){
      test <- bayesMetaintegrator::combineDatasets(test, steps = steps, burnIn = burnInSteps, seed = seed, definedModel = model)
    } else{
      test <- bayesMetaintegrator::combineDatasets(test, steps = steps, burnIn = burnInSteps, seed = seed)
    }
    print(proc.time() - ptm)
    samplingName <- paste0("Sample_", i)
    originalDiscovery$monteCarlo$samplings[[samplingName]] <- test$bayesianMeta$finalResults
    if(any(test$bayesianMeta$finalResults$ESRhat > 1.1)) {
      problemCount <- sum(test$bayesianMeta$finalResults$ESRhat > 1.1)
      percentage <- (problemCount/nrow(test$bayesianMeta$finalResults)) * 100
      warningString <- paste0("Warning: Some ES Rhats > 1.1 (", percentage, "%). Please check $monteCarlo$samplings$finalResults and consider running with more steps.")
      warning(warningString)
    }
    if(any(test$bayesianMeta$finalResults$TauRhat > 1.1)) {
      problemCount <- sum(test$bayesianMeta$finalResults$TauRhat > 1.1)
      percentage <- (problemCount/nrow(test$bayesianMeta$finalResults)) * 100
      warningString <- paste0("Warning: Some tau Rhats > 1.1 (", percentage, "%). Please check $monteCarlo$samplings$finalResults and consider running with more steps.")
      warning(warningString)
    }

  }
  return(originalDiscovery)
}

#' Function that runs combines across Monte Carlo Samples
#'
#' @param discovery Bayesian Meta-Analysis object with Monte Carlo sampling run
#' @param steps The number of steps to combine for each sampling
#' @param burnInSteps Number of steps for burn-in for each sampling - this should be less than steps
#' @param seedNum Number for random seed
#' @param definedMoel Model for Bayesian combination - default, basic model
#' @return MetaIntegrator object with monteCarlo results in final results
#'
bayesianMonteCarloCombine <- function(discovery, steps = 5000, burnInSteps = 1000, seed = 42, definedModel = NULL) {
  #specify model
  model <- function(){
    for(i in 1:I) { #for each study
      y[i, 1] ~ dnorm(mu[i], pow(se[i,1], -2))   #assume study mean comes from normal with mean mu[i] and given std.dev
      mu[i] ~ dnorm(MU, prec.mu)                 #assume mu[i] comes from overall MU and tau (in precision form)
    }
    MU ~ dnorm(0.0, pow(3, -2))                 #MU is distributed between -3,3 with mean 0
    prec.mu <- pow(tau, -2)                     #Who knows where the precision is
    tau ~ dunif(0, 2)
  }
  #define parameters
  jags.data <- list("y", "se", "I")
  # Define parameters of interest
  jags.params <- c("tau", "mu", "MU")
  # Set initials
  jags.inits <- function() {
    list("tau"=runif(1),
         "mu"=rnorm(I)
    )
  }

  #if model given, replace
  if(!is.null(definedModel)) {
    model <- definedModel
  }

  geneList <- c()
  for(sampling in discovery$monteCarlo$samplings) {
    geneList <- union(geneList, sampling$Gene)
  }
  validGenes = unique(geneList)
  #need to get geneList

  metaAnalysisList <- vector("list", length(validGenes))
  n = 1
  for (gene in validGenes){
    y = c()
    for(i in seq(1, length(discovery$monteCarlo$samplings))){
      y = append(y,(discovery$monteCarlo$samplings[[i]] %>% filter(Gene == gene))[,"ES"])
    }
    se <- c()
    for(i in seq(1, length(discovery$monteCarlo$samplings))){
      se = append(se,(discovery$monteCarlo$samplings[[i]] %>% filter(Gene == gene))[,"Tau"])
    }

    y <- as.matrix(y)
    se <- as.matrix(se)
    I <- length(y)
    print(gene)
    jagsfit <- jags(data=jags.data, inits = jags.inits, jags.params, model.file = model, n.iter = 5000, progress.bar = "none", n.burnin=500, n.thin=10,  DIC=TRUE)
    metaAnalysisList[[n]] <- jagsfit$BUGSoutput$summary
    n = n+1
  }

  #go through meta-analysis list
  for(i in (1:length(metaAnalysisList))) {
    gene = as.character(validGenes[i])
    discovery$monteCarlo$pooledResultsDetailed[[gene]] <- metaAnalysisList[[i]]
  }
  discovery$monteCarlo$finalResults <- bind_rows(lapply(discovery$monteCarlo$pooledResultsDetailed, getProbabilityAndFormat), .id = "Gene")
  discovery$monteCarlo$finalResults <- discovery$monteCarlo$finalResults %>% arrange(Pr0)

  if(any(discovery$monteCarlo$finalResults$ESRhat > 1.1)) {
    problemCount <- sum(discovery$monteCarlo$finalResults$ESRhat > 1.1)
    percentage <- (problemCount/nrow(discovery$monteCarlo$finalResults)) * 100
    warningString <- paste0("Warning: Some ES Rhats > 1.1 (", percentage, "%). Please check $monteCarlo$finalResults and consider running with more steps.")
    warning(warningString)
  }
  if(any(discovery$monteCarlo$finalResults$TauRhat > 1.1)) {
    problemCount <- sum(discovery$monteCarlo$finalResults$TauRhat > 1.1)
    percentage <- (problemCount/nrow(discovery$monteCarlo$finalResults)) * 100
    warningString <- paste0("Warning: Some tau Rhats > 1.1 (", percentage, "%). Please check $monteCarlo$finalResults and consider running with more steps.")
    warning(warningString)
  }
  return(discovery)
}

