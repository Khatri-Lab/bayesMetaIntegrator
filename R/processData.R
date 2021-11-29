#' Function that converts probes to gene expression matrix
#'
#' @param dataset Dataset in MetaIntegrator format
#' @param geneSubsetList A vector with genes that we are interested in -
#' these will be converted from probe expression values to gene expression values
#' @return MetaIntegrator datasetObject with gene_expr matrix added
#'

getGeneTable <- function(dataset, geneSubsetList) {
  exprMatrix <- as.data.frame(dataset$expr, check.rows = FALSE)
  #make complete cases
  exprMatrix <- exprMatrix[complete.cases(exprMatrix),]
  keys <- as.data.frame(stack(dataset$keys))
  keys <- tidyr::separate_rows(keys, 1, sep = ",")
  colnames(keys) <- c("gene", "probe")
  exprMatrix$probe <- rownames(exprMatrix)
  #Filter gather probes, genesand filter
  longGenes <- exprMatrix %>% gather(sample, value, -probe) %>% dplyr::left_join(keys) %>% dplyr::filter(gene %in% geneSubsetList)
  #take average across probes
  dm <- longGenes %>% group_by(gene, sample) %>% summarize(mean = mean(value))
  geneSampleFrame <- dm %>% spread(sample, mean) %>% as.data.frame
  geneSampleFrame <- geneSampleFrame[!grepl(",", geneSampleFrame$gene),]
  #gene sample frame check
  geneSampleFrame <- geneSampleFrame[complete.cases(geneSampleFrame), ]
  rownames(geneSampleFrame) <- geneSampleFrame$gene
  geneSampleFrame$gene <- NULL
  #check overall variance
  if(sum(apply(geneSampleFrame, 1, sd) > 0) < nrow(geneSampleFrame)){
    warnMessage <- paste0("Note: Genes with zero variance found in dataset ", dataset$formattedName, ". These have been removed.")
    warning(warnMessage)
    keep <- apply(geneSampleFrame, 1, sd) > 0
    geneSampleFrame <- geneSampleFrame[keep,]
  }
  #within group variance
  if(sum(apply(geneSampleFrame[,unname(dataset$class==1)], 1, sd) > 0) < nrow(geneSampleFrame)){
    warnMessage <- paste0("Note: Genes with zero variance found in dataset, ", dataset$formattedName, " in class 1. These have been removed.")
    warning(warnMessage)
    keep <- apply(geneSampleFrame[,unname(dataset$class==1)], 1, sd) > 0
    geneSampleFrame <- geneSampleFrame[keep,]
  }
  if(sum(apply(geneSampleFrame[,unname(dataset$class==0)], 1, sd) > 0) < nrow(geneSampleFrame)){
    warnMessage <- paste0("Note: Genes with zero variance found in dataset, ", dataset$formattedName, " in class 0. These have been removed.")
    warning(warnMessage)
    keep <- apply(geneSampleFrame[,unname(dataset$class==0)], 1, sd) > 0
    geneSampleFrame <- geneSampleFrame[keep,]
  }
  #knock out genes that have < 3 unique values
  if(sum(apply(geneSampleFrame, 1, function(x)length(unique(x))) < 3) > 0){
    numberFiltered <- sum(apply(geneSampleFrame, 1, function(x)length(unique(x))) < ncol(exprMatrix)/2)
    warnMessage <- paste0("Note:", numberFiltered, " Genes with less than 3 unique values found in dataset, ", dataset$formattedName, ". These have been removed.")
    warning(warnMessage)
    keep <- apply(geneSampleFrame, 1, function(x)length(unique(x))) >= 3
    geneSampleFrame <- geneSampleFrame[keep,]
  }
  #set var and return
  geneSampleFrame$gene <- rownames(geneSampleFrame)
  #set var and return
  geneSampleFrame$gene <- rownames(geneSampleFrame)
  dataset$gene_expr <- geneSampleFrame
  dataset$expr <- NULL
  return(dataset)
}

#' Wrapper function that will loop over a MetaIntegrator discovery object and convert probes to genes for all datasets.
#'
#' @param discovery MetaIntegrator collection of datasets
#' @param geneSubsetList A vector with genes that we are interested in -
#' these will be converted from probe expression values to gene expression values
#' @return MetaIntegrator collection of datasets with gene_expr added to each one.
#' @examples probesToGene(discovery, c('BRCA1', 'APOE'))

probeToGene <- function(discovery, geneSubsetList) {
  #loop through discovery$originalData and convert genes to probes
  discovery$bayesianMeta$originalData <-
    lapply(discovery$originalData, getGeneTable, geneSubsetList = geneSubsetList)
  return(discovery)
}
