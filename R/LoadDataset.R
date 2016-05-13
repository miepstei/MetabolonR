
#' Load Metabolon datasets from preconstructed sample file, metabolite info file and a data file
#'
#' @param sampleFile A filename containing information about the samples
#' @param metaboliteFile A filename containing information about the detected metabolites
#' @param dataFile A filename containing the data, typically raw ion-counts, for the samples and metabolites
#' @return A list of DataFrames containing the sample, metabolite and count information
#' @examples
#' metabolonDataSet <- ReadMetabolonData('sampleInfo.csv', 'metaboliteInfo.csv', 'data.csv')
#' @export

LoadMetabolonData <- function(sampleFile, metaboliteFile, dataFile) {

    # returns: list containing dataframes of the samples, metabolites and metabolite counts

    # read in the metabolites metadata and preprocess. We want the sortorder to match the sample metabolite sortorder
    metabolites <- read.csv(paste(metaboliteFile, sep = ""), header = TRUE, sep = ",", stringsAsFactors = FALSE)
    metabolites$PATHWAY_SORTORDER <- gsub("^([0-9]+)", "M\\1", metabolites$PATHWAY_SORTORDER, perl = TRUE)

    metabolites$PATHWAY_SORTORDER <- gsub("^X - ([0-9]+)", "X\\1", metabolites$PATHWAY_SORTORDER, perl = TRUE)

    # note - there are no unique keys (i.e. PATHWAY SORTORDERS) with the raw data for metabolites containing IS and RS standards so assign a dummy
    metabolites$METABOLITE_NAME <- paste("M", 1:dim(metabolites)[1], sep = "")

    metabolites$METABOLITE_TYPE <- "KNOWN"
    metabolites$METABOLITE_TYPE[grepl('^RS', metabolites$BIOCHEMICAL, perl=TRUE) ] <- "RS"
    metabolites$METABOLITE_TYPE[grepl('^IS', metabolites$BIOCHEMICAL, perl=TRUE)] <- "IS"
    metabolites$METABOLITE_TYPE[grepl('^X', metabolites$BIOCHEMICAL, perl=TRUE)] <- "UNKNOWN"

    # we want to classify all blanks in the SUPER_PATHWAY and SUB_PATHWAYS so change the LEVEL
    levels(metabolites$SUPER_PATHWAY)[1] <- "Unknown"
    levels(metabolites$SUB_PATHWAY)[1] <- "Unknown"

    # read in sample IDs
    sampleIds <- read.csv(paste(sampleFile, sep = ""), header = TRUE, sep = ",", stringsAsFactors = FALSE)
    sampleIds$SAMPLE_NAME <- gsub("\\.|\\-", "", sampleIds$SAMPLE_NAME)


    # read in the volume adjusted ion-counts and preprocess
    rawCounts <- read.csv(paste(dataFile, sep = ","), stringsAsFactors = FALSE, header = TRUE, sep = ",")
    rawCounts$SAMPLE_NAME <- gsub("\\.|\\-", "", rawCounts$SAMPLE_NAME)
    colnames(rawCounts)[2:ncol(rawCounts)] <- metabolites$METABOLITE_NAME

    # convert to numeric all metabolite fields -> causes NAs in blanks or character in cells
    rawCounts[2:ncol(rawCounts)] <- suppressWarnings(sapply(rawCounts[2:ncol(rawCounts)], as.numeric))

    # check for consistency between fields
    if (!setequal(rawCounts$SAMPLE_NAME, sampleIds$SAMPLE_NAME) | !setequal(colnames(rawCounts)[2:ncol(rawCounts)], metabolites$METABOLITE_NAME)) {
        stop("Data is inconsistent in ReadMetabolonData()")
    }

    return(list(sampleIds, metabolites, rawCounts))
}

