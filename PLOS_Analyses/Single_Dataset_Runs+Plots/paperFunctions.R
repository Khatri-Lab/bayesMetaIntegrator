library(tidyverse)
library(cowplot)
library(ggpubr)

getHealthyDiseaseGeneData <- function(dataset, gene = "FANCI"){
  geneData <- dataset$gene_expr[dataset$gene_expr$gene == gene,] %>% dplyr::select(-gene)
  classesH <- unname(dataset$class == 0)
  classesD <- unname(dataset$class == 1)

  healthy = unlist(unname(geneData[classesH]))
  disease = unlist(unname(geneData[classesD]))

  dataframeH = data.frame(value = healthy, status = "Healthy")
  dataframeD = data.frame(value = disease, status = "Asthma")

  plottingFrame = rbind(dataframeH, dataframeD)
  return(plottingFrame)
}

getJaccardSimPoint <- function(geneListA, geneListB){
  union <- length(union(geneListA, geneListB))
  intersect <- length(intersect(geneListA, geneListB))
  return(intersect/union)
}

#Inputs:
## metaobject: MetaIntegrator Metaobject
## FDR: FDR threshold for regular metaintegrator
## studies: number of studies needed for both bayesian and regular MetaIntegrator
## n: number of genes to calculate Jaccard for [VECTOR]
getJaccardPlot_ESOrder <- function(metaobject, FDR, studies, n){
  #filter down the regular meta-integrator and sort by effect size
  regular <- metaobject$metaAnalysis$pooledResults %>%
    filter(numStudies >= studies) %>%
    filter(effectSizeFDR < FDR) %>%
    arrange(desc(abs(effectSize)))

  #filter down bayesian meta-integrator - sort by Pr0
  bayesian <- metaobject$bayesianMeta$finalResults %>%
    arrange(abs(Pr0))

  jaccardResults <- data.frame(n = numeric(), jaccard = numeric())
  for(nGenes in n){
    metaIntegratorGenes <- rownames(regular)[1:nGenes]
    bayesianGenes <- bayesian$Gene[1:nGenes]
    jaccardResults <- rbind(jaccardResults, data.frame(n = nGenes, jaccard = getJaccardSimPoint(metaIntegratorGenes, bayesianGenes)))
  }
  xlabel <- paste0("Number of Genes (MetaIntegrator FDR = ", FDR, ") - ES Order")
  jac_plot <- ggplot(data=jaccardResults, aes(x=n, y=jaccard)) +
    geom_line()+
    geom_point() +
    theme_cowplot() +
    xlab(xlabel) +
    ylab("Jaccard Similarity") +
    ggtitle("Jaccard Comparison - Bayesian vs. MetaIntegrator")
  return(jac_plot)
}


getJaccardPlot_FDROrder <- function(metaobject, FDR, studies, n){
  #filter down the regular meta-integrator and sort by effect size
  regular <- metaobject$metaAnalysis$pooledResults %>%
    filter(numStudies >= studies) %>%
    arrange((abs(effectSizeFDR)))

  #filter down bayesian meta-integrator - sort by Pr0
  bayesian <- metaobject$bayesianMeta$finalResults %>%
    arrange(abs(Pr0))

  jaccardResults <- data.frame(n = numeric(), jaccard = numeric())
  for(nGenes in n){
    metaIntegratorGenes <- rownames(regular)[1:nGenes]
    bayesianGenes <- bayesian$Gene[1:nGenes]
    jaccardResults <- rbind(jaccardResults, data.frame(n = nGenes, jaccard = getJaccardSimPoint(metaIntegratorGenes, bayesianGenes)))
  }
  xlabel <- paste0("Number of Genes - FDR Order")
  jac_plot <- ggplot(data=jaccardResults, aes(x=n, y=jaccard)) +
    geom_line()+
    geom_point() +
    theme_cowplot() +
    xlab(xlabel) +
    ylab("Jaccard Similarity") +
    ggtitle("Jaccard Similarity")
  return(jaccardResults)
}

#Inputs:
## metaobject: MetaIntegrator Metaobject
## FDR: FDR threshold for regular metaintegrator
## studies: number of studies needed for both bayesian and regular MetaIntegrator
## n: number of genes ES to correlate
getESCorrelation <- function(metaobject, FDR, studies, n = 0){
  #filter down the regular meta-integrator and sort by effect size
  regular <- metaobject$metaAnalysis$pooledResults %>%
    filter(numStudies >= studies) %>%
    filter(effectSizeFDR < FDR) %>%
    arrange(desc(abs(effectSize)))
  regular$Gene <- rownames(regular)
  if (n!=0) {
    regular <- regular[1:n,]
  }

  #filter down bayesian meta-integrator - sort by Pr0
  bayesian <- metaobject$bayesianMeta$finalResults %>%
    arrange(abs(Pr0))

  #join the two together - left_join
  combined <- left_join(regular, bayesian, by = "Gene")
  combined <- combined[complete.cases(combined),]

  cor_plot <- ggscatter(combined, x = "effectSize", y = "ES",
                  add = "reg.line",  # Add regression line
                  add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                  conf.int = TRUE, # Add confidence interval
                  xlab = "MetaIntegrator Effect Size",
                  ylab = "Bayesian Effect Size",
                  title = "MetaIntegrator vs. Bayesian Meta-Analysis Effect Size"
  ) + stat_cor(method = "pearson", label.x = -1, label.y = 1.0)
  return(cor_plot)
}

getTauCorrelation <- function(metaobject, studies, n = 0){
  #filter down the regular meta-integrator and sort by effect size
  regular <- metaobject$metaAnalysis$pooledResults %>%
    filter(numStudies >= studies) %>%
    arrange(desc(abs(effectSize)))
  regular$Gene <- rownames(regular)
  if (n!=0) {
    regular <- regular[1:n,]
  }

  #filter down bayesian meta-integrator - sort by Pr0
  bayesian <- metaobject$bayesianMeta$finalResults %>%
    arrange(abs(Pr0))

  #join the two together - left_join
  combined <- left_join(regular, bayesian, by = "Gene")
  combined <- combined[complete.cases(combined),]
  combined$Tau <- combined$Tau^2

  cor_plot <- ggscatter(combined, x = "tauSquared", y = "Tau",
                        add = "reg.line",  # Add regression line
                        add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                        conf.int = TRUE, # Add confidence interval
                        xlab = "MetaIntegrator Tau Squared",
                        ylab = "Bayesian Tau Squared",
                        title = "MetaIntegrator vs. Bayesian Meta-Analysis Tau"
  ) + stat_cor(method = "pearson", label.x = 1, label.y = .5) + ylim(0, 1.7) + xlim(0, 1.7)
  return(cor_plot)
}

getPCorrelation <- function(metaobject, FDR = .05, studies =2, n = 0){
  #filter down the regular meta-integrator and sort by effect size
  regular <- metaobject$metaAnalysis$pooledResults %>%
    filter(numStudies >= studies) %>%
    filter(effectSizeFDR < FDR) %>%
    filter(fisherFDRDown < FDR | fisherFDRUp < FDR) %>%
    arrange(desc(abs(effectSizeFDR)))
  regular$Gene <- rownames(regular)
  if (n!=0) {
    regular <- regular[1:n,]
  }

  #filter down bayesian meta-integrator - sort by Pr0
  bayesian <- metaobject$bayesianMeta$finalResults %>%
    arrange(abs(Pr0))

  #join the two together - left_join
  combined <- left_join(regular, bayesian, by = "Gene")
  combined <- combined[complete.cases(combined),]

  # cor_plot <- ggscatter(combined, x = "effectSizeFDR", y = "Pr0",
  #                       add = "reg.line",
  #                       add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
  #                       conf.int = TRUE, # Add confidence interval
  #                       xlab = "MetaIntegrator Effect Size FDR",
  #                       ylab = "Bayesian Pr = 0",
  #                       title = "MetaIntegrator vs. Bayesian Meta-Analysis P-Value"
  # ) + stat_cor(method = "pearson")
  return(combined)
}

ROCtuneES <- function(metaobject, nVector, FDR = .05, studies = 2) {
  regular <- metaobject$metaAnalysis$pooledResults %>%
    filter(numStudies >= studies) %>%
    filter(effectSizeFDR < FDR) %>%
    arrange(desc(abs(effectSize)))
  regular$Gene <- rownames(regular)

  #Get names of datasets in metaintegrator object
  rocDataframe <- data.frame(n = numeric(), roc = numeric(), dataset = character())
  for(n in nVector) {
    #we will take top n genes, create a filter object, then create a score
    filter <- list()
    filter$effectSizeThresh <- 0
    filter$numberStudiesThresh <- studies
    filter$FDRThresh <- FDR
    filter$filterDescription <- "testFilter"
    filter$timestamp <- Sys.time()
    filter$heterogeneityPvalThresh <- 0
    filter$isLeaveOneOut <- FALSE
    filteredFrame <-regular %>% slice(1:n)
    filter$posGeneNames = filteredFrame[filteredFrame$effectSize > 0, "Gene"]
    filter$negGeneNames = filteredFrame[filteredFrame$effectSize < 0, "Gene"]


    for(dataset in 1:length(metaobject$originalData)){
      scoreResults <- calculateScore(filter, metaobject$originalData[[dataset]])
      rocRes <- calculateROC(predictions=scoreResults, labels=metaobject$originalData[[dataset]]$class)
      rocDataframe <- rbind(rocDataframe, data.frame(n = n, roc = rocRes$auc, dataset = metaobject$originalData[[dataset]]$formattedName))
    }
  }
  rocDataframeRegular <- rocDataframe
  rocDataframeRegular$analysis <- "MetaIntegrator"
  return(rocDataframeRegular)
}



#Need a function that takes top n genes and gets AUC for those genes
ROCtuneFDR <- function(metaobject, nVector, FDR = .05, studies = 2) {
  regular <- metaobject$metaAnalysis$pooledResults %>%
    filter(numStudies >= studies) %>%
    filter(effectSizeFDR < FDR) %>%
    arrange(abs(effectSizeFDR))
  regular$Gene <- rownames(regular)

  #Get names of datasets in metaintegrator object
  rocDataframe <- data.frame(n = numeric(), roc = numeric(), dataset = character())
  for(n in nVector) {
  #we will take top n genes, create a filter object, then create a score
    filter <- list()
    filter$effectSizeThresh <- 0
    filter$numberStudiesThresh <- studies
    filter$FDRThresh <- FDR
    filter$filterDescription <- "testFilter"
    filter$timestamp <- Sys.time()
    filter$heterogeneityPvalThresh <- 0
    filter$isLeaveOneOut <- FALSE
    filteredFrame <-regular %>% dplyr::slice(1:n)
    filter$posGeneNames = filteredFrame[filteredFrame$effectSize > 0, "Gene"]
    filter$negGeneNames = filteredFrame[filteredFrame$effectSize < 0, "Gene"]


    for(dataset in 1:length(metaobject$originalData)){
      scoreResults <- calculateScore(filter, metaobject$originalData[[dataset]])
      rocRes <- calculateROC(predictions=scoreResults, labels=metaobject$originalData[[dataset]]$class)
      rocDataframe <- rbind(rocDataframe, data.frame(n = n, roc = rocRes$auc, dataset = metaobject$originalData[[dataset]]$formattedName))
    }
  }
  rocDataframeRegular <- rocDataframe
  rocDataframeRegular$analysis <- "MetaIntegrator"
  return(rocDataframeRegular)
}

ROCtuneBayesianES <- function(metaobject, nVector, probability, studies = 2) {
  regular <- metaobject$bayesianMeta$finalResults %>%
    filter(nStudies >= studies) %>%
    filter(Pr0 < probability) %>%
    arrange(desc(abs(ES)))

  #Get names of datasets in metaintegrator object
  rocDataframe <- data.frame(n = numeric(), roc = numeric(), dataset = character())
  for(n in nVector) {
    #we will take top n genes, create a filter object, then create a score
    filter <- list()
    filter$effectSizeThresh <- 0
    filter$numberStudiesThresh <- as.numeric(studies)
    filter$FDRThresh <- 0
    filter$filterDescription <- "testFilter"
    filter$timestamp <- Sys.time()
    filter$heterogeneityPvalThresh <- 0
    filter$isLeaveOneOut <- FALSE
    filteredFrame <- regular %>% slice(1:n)
    filter$posGeneNames = filteredFrame[filteredFrame$ES > 0, "Gene"]
    filter$negGeneNames = filteredFrame[filteredFrame$ES < 0, "Gene"]


    for(dataset in 1:length(metaobject$originalData)){
      scoreResults <- calculateScore(filter, metaobject$originalData[[dataset]])
      rocRes <- calculateROC(predictions=scoreResults, labels=metaobject$originalData[[dataset]]$class)
      rocDataframe <- rbind(rocDataframe, data.frame(n = n, roc = rocRes$auc, dataset = metaobject$originalData[[dataset]]$formattedName))
    }
  }
  rocDataframeBayesian <- rocDataframe
  rocDataframeBayesian$analysis <- "Bayesian"
  return(rocDataframeBayesian)
}

ROCtuneBayesianFilter <- function(metaobject, nVector, studies = 2, ESFilter) {
  regular <- metaobject$bayesianMeta$finalResults %>%
    filter(nStudies >= studies) %>%
    filter(abs(ES) > ESFilter) %>%
    arrange(abs(Pr0))

  #Get names of datasets in metaintegrator object
  rocDataframe <- data.frame(n = numeric(), roc = numeric(), dataset = character())
  for(n in nVector) {
    #we will take top n genes, create a filter object, then create a score
    filter <- list()
    filter$effectSizeThresh <- 0
    filter$numberStudiesThresh <- as.numeric(studies)
    filter$FDRThresh <- 0
    filter$filterDescription <- "testFilter"
    filter$timestamp <- Sys.time()
    filter$heterogeneityPvalThresh <- 0
    filter$isLeaveOneOut <- FALSE
    filteredFrame <- regular %>% dplyr::slice(1:n)
    filter$posGeneNames = filteredFrame[filteredFrame$ES > 0, "Gene"]
    filter$negGeneNames = filteredFrame[filteredFrame$ES < 0, "Gene"]


    for(dataset in 1:length(metaobject$originalData)){
      scoreResults <- calculateScore(filter, metaobject$originalData[[dataset]])
      rocRes <- calculateROC(predictions=scoreResults, labels=metaobject$originalData[[dataset]]$class)
      rocDataframe <- rbind(rocDataframe, data.frame(n = n, roc = rocRes$auc, dataset = metaobject$originalData[[dataset]]$formattedName))
    }
  }
  rocDataframeBayesian <- rocDataframe
  rocDataframeBayesian$analysis <- "Bayesian"
  return(rocDataframeBayesian)
}


ROCtuneBayesian <- function(metaobject, nVector, studies = 2) {
  regular <- metaobject$bayesianMeta$finalResults %>%
    filter(nStudies >= studies) %>%
    arrange(abs(Pr0))

  #Get names of datasets in metaintegrator object
  rocDataframe <- data.frame(n = numeric(), roc = numeric(), dataset = character())
  for(n in nVector) {
    #we will take top n genes, create a filter object, then create a score
    filter <- list()
    filter$effectSizeThresh <- 0
    filter$numberStudiesThresh <- as.numeric(studies)
    filter$FDRThresh <- 0
    filter$filterDescription <- "testFilter"
    filter$timestamp <- Sys.time()
    filter$heterogeneityPvalThresh <- 0
    filter$isLeaveOneOut <- FALSE
    filteredFrame <- regular %>% dplyr::slice(1:n)
    filter$posGeneNames = filteredFrame[filteredFrame$ES > 0, "Gene"]
    filter$negGeneNames = filteredFrame[filteredFrame$ES < 0, "Gene"]


    for(dataset in 1:length(metaobject$originalData)){
      scoreResults <- calculateScore(filter, metaobject$originalData[[dataset]])
      rocRes <- calculateROC(predictions=scoreResults, labels=metaobject$originalData[[dataset]]$class)
      rocDataframe <- rbind(rocDataframe, data.frame(n = n, roc = rocRes$auc, dataset = metaobject$originalData[[dataset]]$formattedName))
    }
  }
  rocDataframeBayesian <- rocDataframe
  rocDataframeBayesian$analysis <- "Bayesian"
  return(rocDataframeBayesian)
}

weightedAUC <- function(discovery, rocDataframe){
  datasetSizes <- data.frame(dataset = character(), size = numeric())
  #first get dataset sizes
  for(datasetNum in seq(1:length(discovery$originalData))){
    size = length(discovery$originalData[[datasetNum]]$class)
    name = discovery$originalData[[datasetNum]]$formattedName
    datasetSizes <- rbind(datasetSizes, data.frame(dataset = name, size = size))
  }
  combinedFrame <- left_join(rocDataframe, datasetSizes)
  combinedFrame$multipliedROC <- combinedFrame$roc * combinedFrame$size
  combinedFrame <- combinedFrame %>% group_by(n, analysis) %>% summarise(finalAUC = sum(multipliedROC)/sum(size))
  return(combinedFrame)
}

# metaobject <- discovery
# nVector <- seq(from = 5, to = 200, by = 5)
# bayesian <- ROCtuneBayesian(metaobject, nVector)
# regular <- ROCtuneRegular(metaobject, nVector, .05)
# rocDataframe <- rbind(bayesian, regular)
# ggplot(data = rocDataframe, aes(x=n, y=roc, group = interaction(dataset, analysis), linetype=analysis, color=dataset)) +
#  geom_line() +
#  theme_cowplot() +
#  xlab("Number Genes") +
#  ylab("AUC per Dataset") +
#  ggtitle("Tuning Curve - AUC - MetaIntegrator/Bayesian - Effect Size FDR") +
#  ylim(.5, 1)
