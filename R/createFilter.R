#' Function that returns a filter in discovery based upon the conditions given
#'
#' @param discovery MetaIntegrator object with combineDatasets already run.
#' @param effectSizeThresh The |ES| threshold level
#' @param numberStudiesThresh Threshold for the number of studies the gene was observed in
#' @param probabilityThresh The probability ES = 0 - lower is more stringent
#' @return discovery object with new filter added
#'
filterGenes <- function(discovery, effectSizeThresh = 0, numberStudiesThresh = 2, probabilityThresh = .05) {
  filter <- list()
  filter$effectSizeThresh <- effectSizeThresh
  filter$numberStudiesThresh <- numberStudiesThresh
  filter$FDRThresh <- probabilityThresh
  filter$filterDescription <- "Bayesian Filtration"
  filter$timestamp <- Sys.time()
  filter$heterogeneityPvalThresh <- 0
  filter$isLeaveOneOut <- FALSE

  filteredFrame <-
    discovery$bayesianMeta$finalResults %>%
    filter(nStudies >= numberStudiesThresh,
           abs(ES) > effectSizeThresh,
           Pr0 < probabilityThresh)
  filter$posGeneNames = filteredFrame[filteredFrame$ES > 0, "Gene"]
  filter$negGeneNames = filteredFrame[filteredFrame$ES < 0, "Gene"]
  filterName <- paste0("bayesian_", "PR_", probabilityThresh, "ES_", effectSizeThresh, "nStudies_", numberStudiesThresh)
  discovery$filterResults[[filterName]] <- filter
  return(discovery)
}
