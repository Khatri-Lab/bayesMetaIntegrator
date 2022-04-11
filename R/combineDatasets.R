#' Function that formats the combined data
#'
#' @param pooledResultsDetailedList list of results from combineDatafunction
#' @return Nicely formatted final results
#' @keywords internal
#'
getProbabilityAndFormat <- function(pooledResultsDetailed) {
  lowerTail <- pnorm(0, mean=pooledResultsDetailed["MU", "mean"], sd=pooledResultsDetailed["MU", "sd"], lower.tail=TRUE)
  upperTail <- pnorm(0, mean=pooledResultsDetailed["MU", "mean"], sd=pooledResultsDetailed["MU", "sd"], lower.tail=FALSE)
  returnFrame <- data.frame(pooledResultsDetailed["MU", "mean"], pooledResultsDetailed["MU", "sd"],
                            pooledResultsDetailed["tau", "mean"], pooledResultsDetailed["tau", "sd"],
                            min(lowerTail, upperTail), (nrow(pooledResultsDetailed) - 3),
                            pooledResultsDetailed["MU", "Rhat"], pooledResultsDetailed["tau", "Rhat"])
  colnames(returnFrame) <- c("ES", "ES_SD", "Tau", "TauSD", "Pr0", "nStudies", "ESRhat", "TauRhat")
  return(returnFrame)
}


#' Function that takes individual distributions from getDatasetEffectSizes and combines them
#'
#' @param discovery MetaIntegrator object with getDatasetEffectSizes already run.
#' @param steps The number of MCMC steps to run for each gene
#' @param burnInSteps Number of steps for burn-in - this should be less than steps.
#' @param seedNum Number for random seed
#' @param definedModel Optional JAGs model to use for the combine step, if not running standard one
#' @return discovery$bayesianMeta$finalResults - which contains effect size estimate, tau and Pr(ES = 0) and Rhats
#'
combineDatasets <- function(discovery, steps = 5000, burnIn = 1000, seedNum = 42, definedModel = NULL) {
  validGenes <- discovery$bayesianMeta$datasetEffectSizes$Gene[discovery$bayesianMeta$datasetEffectSizes$nStudies > 1]
  if(length(validGenes) < nrow(discovery$bayesianMeta$datasetEffectSizes)) {
    warning('There are genes with only 1 study in the Meta-Analysis. Note: these have been removed.')
  }

  #use passed in model
  if(!is.null(definedModel)) {
    model <- definedModel
  }
  else {
  #define JAGS model
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
  }

  metaAnalysisList <- vector("list", length(validGenes))
  n = 1
  for (gene in validGenes){
    y <- as.matrix((discovery$bayesianMeta$datasetResults %>% filter(Gene == gene))$effSize)
    se <- as.matrix((discovery$bayesianMeta$datasetResults %>% filter(Gene == gene))$stdDev)
    I <- length(y)
    print(gene)
    set.seed(seedNum)
    jagsfit <- jags(data=jags.data, inits = jags.inits, jags.params, model.file = model, n.iter = steps, progress.bar = "none", n.burnin=burnIn, n.thin=10,  DIC=TRUE)
    metaAnalysisList[[n]] <- jagsfit$BUGSoutput$summary
    n = n+1
  }

  for(i in (1:length(metaAnalysisList))) {
    gene = as.character(validGenes[i])
    names <- as.vector(discovery$bayesianMeta$datasetResults %>% filter(Gene == gene) %>% dplyr::select(name))[,1]
    discovery$bayesianMeta$pooledResultsDetailed[[gene]] <- metaAnalysisList[[i]]
    rownames(discovery$bayesianMeta$pooledResultsDetailed[[gene]])[3:(3+length(names)-1)] <- names
  }
  discovery$bayesianMeta$finalResults <- bind_rows(lapply(discovery$bayesianMeta$pooledResultsDetailed, getProbabilityAndFormat), .id = "Gene")
  discovery$bayesianMeta$finalResults <- discovery$bayesianMeta$finalResults %>% arrange(Pr0)
  #check Rhats
  if(any(discovery$bayesianMeta$finalResults$ESRhat > 1.1)) {
    problemCount <- sum(discovery$bayesianMeta$finalResults$ESRhat > 1.1)
    percentage <- (problemCount/nrow(discovery$bayesianMeta$finalResults)) * 100
    warningString <- paste0("Warning: Some ES Rhats > 1.1 (", percentage, "%). Please check $bayesianMeta$finalResults and consider running with more steps.")
    warning(warningString)
  }
  if(any(discovery$bayesianMeta$finalResults$TauRhat > 1.1)) {
    problemCount <- sum(discovery$bayesianMeta$finalResults$TauRhat > 1.1)
    percentage <- (problemCount/nrow(discovery$bayesianMeta$finalResults)) * 100
    warningString <- paste0("Warning: Some tau Rhats > 1.1 (", percentage, "%). Please check $bayesianMeta$finalResults and consider running with more steps.")
    warning(warningString)
  }

  return(discovery)
}

#' Function that returns the MCMC model that was used to combine across datasets
#'
#' @param discovery MetaIntegrator object with getDatasetEffectSizes already run.
#' @param gene The gene we are interested in
#' @param steps The number of MCMC steps to run for each gene
#' @param burnInSteps Number of steps for burn-in - this should be less than steps
#' @param seedNum Number for random seed
#' @param definedModel Optional JAGs model to use for the combine step, if not running standard one
#' @return MCMC model for the gene combination
#'
getCombiningModel <- function(discovery, gene, steps = 5000, burnIn = 1000, seedNum = 42, definedModel = NULL) {
  if(!is.null(definedModel)) {
    model <- definedModel
  }
  else {
  #define JAGS model
    model <- function(){
      for(i in 1:I) { #for each study
        y[i, 1] ~ dnorm(mu[i], pow(se[i,1], -2))   #assume study mean comes from normal with mean mu[i] and given std.dev
        mu[i] ~ dnorm(MU, prec.mu)                 #assume mu[i] comes from overall MU and tau (in precision form)
      }
      MU ~ dnorm(0.0, pow(3, -2))                 #MU is distributed between -3,3 with mean 0
      prec.mu <- pow(tau, -2)                     #Who knows where the precision is
      tau ~ dunif(0, 2)
      #tau ~ dgamma(0, )
    }

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
  y <- as.matrix((discovery$bayesianMeta$datasetResults %>% filter(Gene == gene))$effSize)
  se <- as.matrix((discovery$bayesianMeta$datasetResults %>% filter(Gene == gene))$stdDev)
  I <- length(y)
  set.seed(seedNum)
  jagsfit <- jags(data=jags.data, inits = jags.inits, jags.params, model.file = model, n.iter = steps, progress.bar = "none", n.burnin=burnIn, n.thin=10,  DIC=TRUE)
  return(as.mcmc(jagsfit))
}

#' Function that plots the JAGs traces
#'
#' @param model MCMC model from getCombiningModel
#' @param parameters Parameters we are interested in - likely MU and tau
#' @param trace If true, plots the time series trace. If false, plots the posterior histogram.
#' @return Plot specified by function
#'

getCombineTrace <- function(model, parameters = c("MU", "tau"), trace = TRUE) {
  color_scheme_set("mix-blue-red")
  if (trace == TRUE) {
  return(mcmc_trace(model, pars = parameters,
             facet_args = list(ncol = 1, strip.position = "left")))
  } else {
    return(mcmc_hist_by_chain(model, pars = parameters))
  }
}
