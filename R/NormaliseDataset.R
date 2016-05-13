
#' Convenience function for the median normalisation.
#'
#' The function assigns a specific value to members of a specific group. Used in order
#' to group Metabolon run days together, e.g. 1, 2, 3 are actually the same run day.
#'
#' @param groupList A list of arrays of groups
#' @param valInGroup An value to find in one of the arrays in the groupList
#' @return the index in the list that valInGroup belongs to
#' @examples
#' AssignLabelsToGroups(list(c(1, 2, 3), c(4, 5, 6)), 6)
#' #returns 2 as item 6 is in the second group in the list


AssignLabelsToGroups <- function(groupList, valInGroup) {
    # returns the group in which the value resides
    for (i in 1:length(groupList)) {
        groupItems <- groupList[[i]]
        if (valInGroup %in% groupItems) {
            return(i)
        }
    }
    return(-1)
}


#' Performs Metabolon median normalisation, typically by 'Run Day'
#'
#' @param cMetInfo A DataFrame of metabolite information
#' @param cSampleInfo A DataFrame of sample information
#' @param cMetLog10DataLong A DataFrame of log 10 Metabolite Counts in Long Format
#' @param normaliseField A field in cSampleInfo over which the Normalisation should be performed
#' @param normaliseGroups A list of groupings of the levels of normaliseField to treat as normalisation blocks
#' @param scaleToMedianRawCount A boolean that asks whether to adjust to the global median metabolite level
#' or to adjust to a median of 1
#' @return A dataframe in Long format normalised by ''normaliseField''.
#' @importFrom magrittr "%>%"
#' @importFrom plyr "." "here"
#'
#' @export



MedianNormaliseDataset <- function(cSampleInfo, cMetLog10DataLong, normaliseField , normaliseGroups, scaleToMedianRawCount = TRUE) {

    # check normaliseField is in cSampleInfo
    if (!(normaliseField %in% colnames(cSampleInfo))) {
        stop(paste("median_within_rundays:Normalisation field [", normaliseField, "]", " is not included in this sample information dataset"))
    }

    # pull the relevant samples out
    sampleNames <- cSampleInfo %>% dplyr::select(SAMPLE_NAME) %>% unlist %>% as.vector
    metData <- dplyr::filter(cMetLog10DataLong, SAMPLE_NAME %in% sampleNames)

    # check normaliseGroups are consistent with cSampleInfo and all are accounted for
    normLevels <- as.character(sort(unname(unlist(unique(dplyr::select_(cSampleInfo, normaliseField))))))
    normaliseGroupsTest <- as.character(sort(as.numeric(unlist((normaliseGroups)))))

    if (!isTRUE(all.equal(normaliseGroupsTest, normLevels))) {
        print("all levels in Normalisation Groups not accounted for:")
        print("Proposed Groups: ")
        print(normaliseGroupsTest)
        print("Groups in data:")
        print(normLevels)
        stop()
    }

    # we have checked the inputs - lets perform the normalisation by subtracting the median of each metabolite in the group from each metabolite

    # create the groups over which data is to be normalised
    cSampleInfo$NORMALISATION_GROUP <- sapply(cSampleInfo$RUN_DAY, FUN = function(x) {AssignLabelsToGroups(normaliseGroups, x)})
    normalisationData <- dplyr::left_join(metData, dplyr::select(cSampleInfo, SAMPLE_NAME, NORMALISATION_GROUP), by = "SAMPLE_NAME")

    # work out the function (e,g median) on each metabolite
    metData <- plyr::ddply(normalisationData, .(METABOLITE_NAME), dplyr::mutate, GLOBAL_METABOLITE_MEDIAN = median(METABOLITE_COUNT, na.rm = TRUE))

    # convert to logs for stablility
    metData$GLOBAL_METABOLITE_MEDIAN <- log2(metData$GLOBAL_METABOLITE_MEDIAN)
    metData$METABOLITE_COUNT <- log2(metData$METABOLITE_COUNT)

    if (scaleToMedianRawCount) {
        metData <- plyr::ddply(metData, .(METABOLITE_NAME, NORMALISATION_GROUP), dplyr::mutate, NORMALISED_BY_MEDIAN_GROUP = 2^((METABOLITE_COUNT - median(METABOLITE_COUNT, na.rm = TRUE)) +
            GLOBAL_METABOLITE_MEDIAN))
    } else {
        metData <- plyr::ddply(metData, .(METABOLITE_NAME, NORMALISATION_GROUP), dplyr::mutate, NORMALISED_BY_MEDIAN_GROUP = 2^(METABOLITE_COUNT - median(METABOLITE_COUNT, na.rm = TRUE)))
    }
    return(dplyr::select(metData, SAMPLE_NAME, METABOLITE_NAME, METABOLITE_COUNT = NORMALISED_BY_MEDIAN_GROUP))
}


#' Adjust the metabolite counts for differing volumes in supplied biological material
#'
#' @param cSampleInfo A Data Frame containing cleaned information about the samples
#' @param cMetInfo A Data Frame  containing cleaned information about the detected metabolites
#' @param cMetDataLong A Data Frame  containing cleaned data raw ion-counts in Long format
#' @param normalVolume The volume (in microL) to adjust the ion-counts to
#' @param sampleType The labelling of the sampletype to apply the adjustment to - i.e. not MTRX
#' @return A Data Frame containing volume adjusted count information
#' @examples
#' volumeAdjustedData <- VolumeAdjustMetabolonData(sampleInfo, metInfo, metData, normalVolume = 0.95)
#'
#' @export

VolumeNormaliseDataset <- function(cSampleInfo, cMetInfo, cMetDataLong, normalVolume, sampleType) {

  #choose just the experimental samples
  samplesToBeAdjusted <- cSampleInfo %>% dplyr::filter(SAMPLE_TYPE %in% sampleType) %>% dplyr::select(SAMPLE_NAME) %>% unlist %>% as.vector

  #ignore any internal or recovery standards
  metsToBeAdjusted <- cMetInfo %>% dplyr::filter(METABOLITE_TYPE %in% c("KNOWN", "UNKNOWN")) %>% dplyr::select(METABOLITE_NAME) %>% unlist %>% as.vector
  dataToBeAdjusted <- cMetDataLong %>% dplyr::filter(SAMPLE_NAME %in% samplesToBeAdjusted & METABOLITE_NAME %in% metsToBeAdjusted)

  # add the volume field back to the frame
  dataToBeAdjusted <- dplyr::left_join(dataToBeAdjusted, dplyr::select(cSampleInfo, SAMPLE_NAME, VOLUME_EXTRACTED_UL), by = "SAMPLE_NAME")
  volumeAdjustedData <- dplyr::mutate(dataToBeAdjusted, VOL_ADJ_COUNT = (normalVolume/VOLUME_EXTRACTED_UL) * METABOLITE_COUNT)
  volumeAdjustedData <- volumeAdjustedData %>% dplyr::select(SAMPLE_NAME, METABOLITE_NAME, METABOLITE_COUNT = VOL_ADJ_COUNT)
  unadjustedData <- dplyr::filter(cMetDataLong, !(SAMPLE_NAME %in% samplesToBeAdjusted), !(METABOLITE_NAME %in% metsToBeAdjusted))

  # need to add back any internal standards data in the adjusted samples
  metsToBeAddedBack <- cMetInfo %>% dplyr::filter(!(METABOLITE_NAME %in% metsToBeAdjusted)) %>% dplyr::select(METABOLITE_NAME) %>% unlist %>% as.vector
  dataToBeAddedBack <- cMetDataLong %>% dplyr::filter(SAMPLE_NAME %in% samplesToBeAdjusted & METABOLITE_NAME %in% metsToBeAddedBack)

  # internal standards and volume adjusted data of the adjusted samples
  cMetDataVolAdjusted <- rbind(dataToBeAddedBack, volumeAdjustedData)

  # we want to add the missing unadjusted rows unadjusted samples back onto the long data
  cMetDataVolAdjusted <- rbind(unadjustedData, cMetDataVolAdjusted)
}

#' Normalises a Metabolon dataset using the Cross-contribution Compensating
#' Multiple Internal standard normalistaion" (CCMN) method, Redestig et al. 2009
#'
#' @param metInfo A dataframe of metabolite metadata (must include Internal Standards)
#' @param sampleInfo A dataframe containing the sample metadata
#' @param metDataLong A dataframe of metabolite data in Long Format
#' @param biologicalFactor A chararacter factor in the samples metadata to calulcate cross-contamination with
#' @param matchOrder flag for re-ordering
#'
#' @importFrom plyr "."
#' @importFrom dplyr "%>%"
#'
#' @export

CCMNNormaliseDataset <- function(metInfo, sampleInfo, metDataLong, biologicalFactor, matchOrder = NULL) {

  #check normaliseField is in cSampleInfo
  if (!(biologicalFactor %in% colnames(sampleInfo))) {
    stop(paste("normalise_by_ccmn:Normalisation field [", normaliseField, "]", " is not included in this sample information dataset"))
  }

  #check normaliseGroups are consistent with cSampleInfo and all are accounted for
  normLevels <- as.character(unname(unlist(unique(dplyr::select_(sampleInfo, normaliseField)))))
  normaliseGroupsTest <- unlist(normaliseGroups)

  if ( ! isTRUE(all.equal(normaliseGroupsTest, normLevels) )) {
    print("all levels in Normalisation Groups not accounted for:")
    print("Proposed Groups: ")
    print(normaliseGroupsTest)
    print("Groups in data:")
    print(normLevels)
    stop()
  }



  #create the groups over which data is to be normalised
  cSampleInfo <- plyr::ddply(cSampleInfo, .(SAMPLE_NAME), plyr::mutate, NORMALISATION_GROUP = AssignLabelsToGroups(normaliseGroups, RUN_DAY))
  normalisationData <- plyr::join(metData, dplyr::select(cSampleInfo,SAMPLE_NAME,NORMALISATION_GROUP), by="SAMPLE_NAME")

  # get the metabolites which are labelled as internal standards
  ISMetabolites <- dplyr::metInfo %>% dplyr::filter(METABOLITE_TYPE == "IS")

  #check the list inst empty
  if (nrow(ISMetabolites) == 0){
    stop("normalise_by_NOMIS : Cannot normalise by nomis - no internal standards detected in MetInfo")
  }

  #pull out the metabolite data pertaining to these standards
  ISMetaboliteData <- metDataLong %>% dplyr::filter(METABOLITE_NAME %in% ISMetabolites$METABOLITE_NAME)

  #check we have internal standard data
  if (nrow(ISMetaboliteData) == 0){
    stop("normalise_by_NOMIS : Cannot normalise by nomis - no internal standards detected in MetData")
  }

  #pick out the metabolites and their platform by which they are to be normalised (known and unknown)
  expMetabolites <- metInfo %>% dplyr::filter(METABOLITE_TYPE == "KNOWN" | METABOLITE_TYPE == "UNKNOWN") %>% dplyr::select(METABOLITE_NAME)
  expMetaboliteData <- cMetDataLong %>% dplyr::filter(METABOLITE_NAME %in% expMetabolites$METABOLITE_NAME)

  #need to create a design matrix which specifies which samples are in which biological factor of interest
  design <- dplyr::select_(cSampleInfo,"SAMPLE_NAME", biologicalFactor)
  design <- reshape2::dcast(design,paste("SAMPLE_NAME", "~", biologicalFactor) )
  design <- design %>% dplyr::mutate_each(funs(ifelse(is.na(.), 0, 1)), -SAMPLE_NAME)

  #check design matrix - each row should sum to 1 ie. there should be one biological factor for each sample
  if(!all(rowSums(select(design, -SAMPLE_NAME))==1)){
    stop("[ERROR]: CCMN - invalid design matrix")
  }

  #step 1:scale and mean centre the internal standards data
  ISMetaboliteDataScaled <- plyr::ddply(ISMetaboliteData, "METABOLITE_NAME" , plyr::mutate, METABOLITE_COUNT = (METABOLITE_COUNT - mean(METABOLITE_COUNT)) / sd(METABOLITE_COUNT))

  #regress the internal standards on the design matrix

  ISMetaboliteDataScaled <- MetabolonR::met_data_long_to_wide(ISMetaboliteDataScaled)

  if(!is.null(matchOrder)){
    ISMetaboliteDataScaled <- ISMetaboliteDataScaled[match( matchOrder , ISMetaboliteDataScaled$SAMPLE_NAME),]
  }

  #save sample names for later
  sampleNames <- ISMetaboliteDataScaled[1]
  ISMetaboliteDataScaled <- ISMetaboliteDataScaled[-1]

  lmFit <- lm(as.matrix(ISMetaboliteDataScaled) ~ -1 + as.matrix(design[-1]))
  residualDesign <- resid(lmFit)

  #decompose the residuals into their principal components
  pc <- prcomp(residualDesign, center = F, scale = F)

  #predict the scores of the first n prinicpal components
  pcScores <- predict(pc, residualDesign)
  pcScores <- pcScores[, 1:2]

  expMetaboliteDataScaled <- MetabolonR::met_data_long_to_wide(expMetaboliteData)

  if(!is.null(matchOrder)){
    expMetaboliteDataScaled <- expMetaboliteDataScaled[match(matchOrder , expMetaboliteDataScaled$SAMPLE_NAME),]
  }

  expMetaboliteDataScaled <- expMetaboliteDataScaled[-1]
  expMetaboliteDataScaled <- scale(expMetaboliteDataScaled)
  means <- attr(expMetaboliteDataScaled, "scaled:center")
  sd <- attr(expMetaboliteDataScaled, "scaled:scale")

  lmFit2 <- lm(as.matrix(expMetaboliteDataScaled) ~ -1 + as.matrix(pcScores))
  scaledNorm <- expMetaboliteDataScaled - predict(lmFit2, data.frame(pcScores))

  #reapply the mean and scaling
  scaledNorm <- sweep(scaledNorm, 2, sd, "*")
  scaledNorm <- sweep(scaledNorm, 2, means, "+")

  #stick the SAMPLE_NAME back on.
  normalisedData <- data.frame(sampleNames, scaledNorm)
  MetabolonR::met_data_wide_to_long(normalisedData)
}


#' Normalises a Metabolon dataset using the "Normalization using Optimal
#' selection of Multiple Internal Standards" method of Sysi-Aho et al 2007
#'
#' @param metInfo A dataframe of metabolite metadata (must include the Internal Standards)
#' @param metDataLong A dataframe of metabolite data in Long Format
#'
#' @return A dataframe of NOMIS normalised data in Long Format
#'
#' @export


NOMISNormaliseDataset <- function(metInfo, metDataLong) {
  ISMetabolites <- metInfo %>% dplyr::filter(METABOLITE_TYPE == "IS")

  #check the list isn't empty
  if (nrow(ISMetabolites) == 0){
    stop("normalise_by_NOMIS : Cannot normalise by nomis - no internal standards detected in MetInfo")
  }

  #pull out the metabolite data pertaining to these standards
  ISMetaboliteData <- metDataLong %>% dplyr::filter(METABOLITE_NAME %in% ISMetabolites$METABOLITE_NAME)

  #check we have internal standard data
  if (nrow(ISMetaboliteData) == 0){
    stop("normalise_by_NOMIS : Cannot normalise by nomis - no internal standards detected in MetData")
  }

  expMetabolites <- metInfo %>% dplyr::filter(METABOLITE_TYPE == "KNOWN" | METABOLITE_TYPE == "UNKNOWN") %>% dplyr::select(METABOLITE_NAME, PLATFORM)
  expMetaboliteData <- metDataLong %>% dplyr::filter(METABOLITE_NAME %in% expMetabolites$METABOLITE_NAME)

  #create the response and covariate data frames
  response <- as.matrix(MetabolonR::MetDataLongToWide(expMetaboliteData)[-1]) #remove sample name

  #we need to construct the covariates based on the individual platforms of each IS
  platforms <- unique(ISMetabolites$PLATFORM)

  normalised <- matrix(NA,dim(response)[1],dim(response)[2])

  for (ii in seq(1,dim(response)[2])) {
    metaboliteName <- colnames(response)[ii]
    platformForMetabolite <- metInfo %>% dplyr::filter(METABOLITE_NAME == metaboliteName) %>% dplyr::select(PLATFORM) %>% unlist %>% as.vector
    internalStandardsForMetabolite <- ISMetabolites %>% dplyr::filter(PLATFORM == platformForMetabolite) %>% dplyr::select(METABOLITE_NAME) %>% unlist %>% as.vector
    covariates <- as.matrix(MetabolonR::MetDataLongToWide(dplyr::filter(ISMetaboliteData,METABOLITE_NAME %in% internalStandardsForMetabolite))[-1])
    regress <- lm(response[,ii] ~ covariates, na.action = na.exclude)

    normalised[,ii] <- residuals(regress) #want to keep the NAs in the output
    normalised[,ii] <- normalised[,ii] + mean(response[,ii], na.rm = TRUE)
  }

  normalised <- as.data.frame(normalised)
  colnames(normalised) <- colnames(response)
  normalisedData <- cbind(MetabolonR::MetDataLongToWide(expMetaboliteData)[1] , normalised)
  normalisedData <- MetabolonR::MetDataWideToLong(normalisedData)
}

