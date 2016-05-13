
#' Filters a data frame in wide format that contains metabolites present
#' in every sample
#'
#' @param cMetDataWide a count matrix of Metabolon data in Wide format
#' @return a Data Frame which has metabolites missing in *any* sample removed
#' @examples
#' allMetsPresentWide <- AllMetsPresent(metDataWide)
#' @export

AllMetsPresent <- function(cMetDataWide) {
    allMetsPresent <- (colSums(is.na(cMetDataWide)) == 0)
    return(cMetDataWide[, allMetsPresent])
}


#' Filters a data frame in long format that contains metabolites present
#' in every sample
#'
#' @param cMetDataLong a count matrix of Metabolon data in Wide format
#' @return a Data Frame which has metabolites missing in *any* sample removed
#' @examples
#' allMetsPresentLong <- AllMetsPresentLong(cMetDataLong)
#' @export

AllMetsPresentLong <- function(cMetDataLong) {
    # function returns the metabolites present in every sample in the supplied long frame
    completeMets <- cMetDataLong %>% plyr::ddply(.(METABOLITE_NAME), plyr::summarise, NO_MISSING = sum(is.na(METABOLITE_COUNT)))
    completeMets <- completeMets %>% dplyr::filter(NO_MISSING == 0) %>% dplyr::select(METABOLITE_NAME) %>% unlist %>% as.vector
    allPresentMetsInDataset <- dplyr::filter(cMetDataLong, METABOLITE_NAME %in% completeMets)
}

#' Filters a metabolite count data frame in long format.  Removes metabolites missing in samples at a greater
#' frequency than a threshold set by the user
#'
#' @param cMetDataLong a count matrix of Metabolon data in Long format
#' @param cMetDataLong proportion of metabolite frequency to filter out
#' @return a Data Frame which has metabolites present at a lower frequency than 'proportion' removed
#' @examples
#' allMetsPresentLong <- FilterMissingMets(cMetDataLong, 0.2)
#' removes all metabolites missing in more than 20% of the samples
#' @export

FilterMissingMets <- function(cMetDataLong, proportion) {
    # function that filters mets by a missing percentage i.e > 20% of mets are missing then remove this metabolite
    filteredMets <- cMetDataLong %>% plyr::ddply(.(METABOLITE_NAME), plyr::summarise, PROP_MISSING = sum(is.na(METABOLITE_COUNT))/length(METABOLITE_COUNT))
    filteredMets <- filteredMets %>% dplyr::filter(PROP_MISSING <= proportion) %>% dplyr::select(METABOLITE_NAME) %>% unlist %>% as.vector
    filteredMetsInDataset <- dplyr::filter(cMetDataLong, METABOLITE_NAME %in% filteredMets)
}

#' Transforms the metabolite count in long format by taking log10 of each value
#'
#' @param cMetDataLong a count matrix of Metabolon data in Long format
#' @return a Data Frame where the metabolite counts are now log10 values
#' @examples
#' logMetCounts <- TransformLog10Mets(cMetDataLong)
#' @export

TransformLog10Mets <- function(cMetDataLong) {
    return(dplyr::mutate(cMetDataLong, METABOLITE_COUNT = log10(METABOLITE_COUNT)))
}

#' Transforms the metabolite count in long format by exponentiating each value to base 10
#'
#' @param cMetDataLong a log10 count matrix of Metabolon data in Long format
#' @return a Data Frame where the metabolite counts are original values
#' @examples
#' expMetCounts <- TransformExp10Mets(logMetCounts)
#' @export

TransformExp10Mets <- function(cMetDataLong) {
    # transforms the metabolite data by taking 10^ of the count
    return(dplyr::mutate(cMetDataLong, METABOLITE_COUNT = 10^METABOLITE_COUNT))
}

#' Calculates the z-score of each metabolite count within each metabolite
#' Typically requires a log-transformed count matrix for normality
#'
#' @param cMetDataLong a count matrix of Metabolon data in Long format
#' @return a Data Frame where the metabolite counts within each metabolite are z-scores
#' @examples
#' zMetCounts <- TransformZscoreMets(logLongMetCounts)
#' @export

TransformZscoreMets <- function(cMetDataLong) {
    # transforms the metabolite data calculating the z-score for each sample within each metabolite
    return(plyr::ddply(cMetDataLong, "METABOLITE_NAME", transform, METABOLITE_COUNT = (METABOLITE_COUNT - mean(METABOLITE_COUNT, na.rm = TRUE))/sd(METABOLITE_COUNT, na.rm = TRUE)))
}

#' Convenience function for transforming a Long count Data Frame to a Wide count Data Frame
#'
#' @param cMetDataLong a count matrix of Metabolon data in Long format
#' @return a count matrix of Metabolon data in Wide format
#' @examples
#' wideMetCounts <- MetDataLongToWide(longMetCounts)
#' @export

MetDataLongToWide <- function(cMetDataLong) {
    # presents the metabolite data in wide format
    return(reshape2::dcast(cMetDataLong, formula = SAMPLE_NAME ~ METABOLITE_NAME, value.var = "METABOLITE_COUNT"))
}

#' Convenience function for transforming a Wide count Data Frame to a Long count Data Frame
#'
#' @param cMetDataWide a count matrix of Metabolon data in Long format
#' @return a count dataframe of Metabolon data in Long format
#' @examples
#' wideMetCounts <- MetDataLongToWide(longMetCounts)
#' @export

MetDataWideToLong <- function(cMetDataWide) {
    # melts the wide data format to the long data format
    cMetDataLong <- reshape2:::melt.data.frame(cMetDataWide, id.vars = "SAMPLE_NAME", variable.name = "METABOLITE_NAME", value.name = "METABOLITE_COUNT")
    return(cMetDataLong)
}

#' Calculates a scaled and centred PCA analysis of a Metabolon dataset
#'
#' @param sampleInfo a dataframe of sample metadata
#' @param cMetDataLong a count dataframe of Metabolon data in Long format
#' @return an S3 pca object
#' @examples
#' wideMetCounts <- CalculateMetPCA(longMetCounts)
#' @export

CalculateMetPCA <- function(cMetDataLong) {

  nrowsOrig <- dim(cMetDataLong)[1]

  cMetDataCompleteLong <- AllMetsPresentLong(cMetDataLong)

  nrowsNew <- dim(cMetDataCompleteLong)[1]

  if(nrowsOrig != nrowsNew){
    message(paste("[WARN]: pca_plot - datapoints removed from",nrowsOrig,"to",nrowsNew,"due to missing data"))
  }

  cMetData <- MetDataLongToWide(cMetDataCompleteLong)
  cMetDataNumbers <- dplyr::select(cMetData, -SAMPLE_NAME)

  pca <- prcomp(cMetDataNumbers, center=T, scale. = T)
  pca$SAMPLE_NAME <- cMetData$SAMPLE_NAME
  return(pca)

}
