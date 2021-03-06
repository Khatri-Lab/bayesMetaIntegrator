% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/singleDataset.R
\name{getDatasetEffectSizes}
\alias{getDatasetEffectSizes}
\title{Wrapper function that passes whole collection of datasets to calculate effect size distribution per dataset}
\usage{
getDatasetEffectSizes(
  discovery,
  cores = 2,
  steps = 1000,
  burnInSteps = 200,
  seedNum = 42
)
}
\arguments{
\item{discovery}{MetaIntegrator object with gene_expr populated}

\item{cores}{The number of cores allowed to be run}

\item{steps}{The number of MCMC steps to run for each gene}

\item{burnInSteps}{Number of steps for burn-in - this should be less than steps}

\item{seedNum}{Number for random seed}
}
\value{
MetaIntegrator object with $bayesianMeta$datasetEffectSizes and $bayesianData$originalData$dataset$datasetResults populated
}
\description{
Wrapper function that passes whole collection of datasets to calculate effect size distribution per dataset
}
