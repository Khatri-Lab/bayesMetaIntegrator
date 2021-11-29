#' Function that creates a forest plot dataset and pooled effect size estimates
#'
#' @param discovery MetaIntegrator collection of datasets with meta-analysis completed
#' @param geneToPlot The gene we are interested in plotting the forest plot for
#' @param plotType 'bayesian' denotes density forest plot, anything else denotes traditional forest plot
#' @return Forest plot of gene across datasets
#'


makeBayesianForestPlot <- function(discovery, geneToPlot, plotType = "bayesian") {
  if (plotType == "bayesian") {
    listOfPlots <- c()
    geneRowES <- discovery$bayesianMeta$datasetEffectSizes %>% filter(Gene == geneToPlot) %>% select(-Gene, -nStudies)
    geneRowSTD <- discovery$bayesianMeta$datasetStdDev %>% filter(Gene == geneToPlot) %>% select(-Gene, -nStudies)

    for(dataset in colnames(geneRowES)) {
      mean <- as.numeric(geneRowES[,dataset])
      stdDev <- as.numeric(geneRowSTD[,dataset])
      p1 <- ggplot(data = data.frame(x = c(-3, 3)), aes(x)) +
        stat_function(fun = dnorm, n = 101, args = list(mean = mean, sd = stdDev), geom = 'area', fill = 'blue') + ylab(dataset) +
        scale_y_continuous() +
        xlab("Effect Size") +
        geom_vline(xintercept=mean, linetype='dotted', col = 'white') +
        geom_vline(xintercept = 0, size=1) +
        theme_cowplot() +
        theme(aspect.ratio = .2) +
        theme(axis.text.y=element_blank(),
              axis.ticks.y=element_blank(),
              axis.title=element_text(size=10))
      listOfPlots[[dataset]] <- p1
    }

    p1 <- ggplot(data = data.frame(x = c(-3, 3)), aes(x)) +
      stat_function(fun = dnorm, n = 101, args = list(mean = discovery$bayesianMeta$finalResults[discovery$bayesianMeta$finalResults$Gene == geneToPlot, "ES"], sd = discovery$bayesianMeta$finalResults[discovery$bayesianMeta$finalResults$Gene == geneToPlot, "ES_SD"]), geom = 'area', fill = 'red') + ylab("Pooled") +
      scale_y_continuous() + xlab("Pooled Effect Size") +
      geom_vline(xintercept = 0, size=1) +
      theme_cowplot() +
      geom_vline(xintercept=discovery$bayesianMeta$finalResults[discovery$bayesianMeta$finalResults$Gene == geneToPlot, "ES"], linetype='dotted', col = 'white') +
      theme(aspect.ratio = .2) +
      theme(axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            axis.title=element_text(size=10))

    listOfPlots[["pooled"]] <- p1

    return(wrap_plots(listOfPlots, ncol = 1) + plot_annotation(
      title = geneToPlot, theme = theme(plot.title = element_text(hjust = 0.5))))
  }
  else {
    studyData <- data.frame(means = unlist(discovery$bayesianMeta$datasetEffectSizes[discovery$bayesianMeta$datasetEffectSizes$Gene == geneToPlot, !names(discovery$bayesianMeta$datasetEffectSizes) %in% c("nStudies", "Gene")]))
    studyData$SEs <- data.frame(SEs = unlist(discovery$bayesianMeta$datasetStdDev[discovery$bayesianMeta$datasetStdDev$Gene == geneToPlot, !names(discovery$bayesianMeta$datasetStdDev) %in% c("nStudies", "Gene")]))
    getFormattedName <- function(uglyName) {
      if (uglyName %in% names(discovery$originalData)) {
        return(discovery$originalData[uglyName][[1]]$formattedName)
      }
      return("")
    }
    studyData$names <- sapply(rownames(studyData), getFormattedName)
    print(class(studyData$names))
    studyData <- studyData[order(studyData$names), ]
    studyMeans <- studyData$means
    studySEs <- unlist(studyData$SEs)
    studyNames <- studyData$names
    names(studyMeans) <- studyNames
    names(studySEs) <- studyNames


    pooledMean <- discovery$bayesianMeta$finalResults[discovery$bayesianMeta$finalResults == geneToPlot,
                                                        "ES"]
    pooledSE <- discovery$bayesianMeta$finalResults[discovery$bayesianMeta$finalResults == geneToPlot,
                                                      "Tau"]

    rmetaPlots <- rmeta::metaplot(studyMeans, studySEs, labels = studyNames,
                    xlab = "Standardized Mean Difference (log2 scale)",
                    ylab = "", colors = rmeta::meta.colors(box = "royalblue",
                                                           lines = "lightblue", zero = "black", summary = "orange",
                                                           text = "black"), summn = pooledMean, sumse = pooledSE,
                    sumnn = 1/pooledSE^2, main = geneToPlot)
    return(rmetaPlots)

  }
}
