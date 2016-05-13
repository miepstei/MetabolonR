#' creates or cleans a given directory (non recursively). If the dir doesn't exist
#' the function will create it, if it does exist, the function will remove all items
#' from the directory. Used to create ''Output'' directories whose items are expendible
#'
#' @param directory - the dir to create or to remove all items from non recursively
#' @examples
#' CreateOrClean('path/to/my/directory')
#'
#' @export

CreateOrClean <- function(directory) {
    # creates or cleans a given directory (non recursively)
    if (dir.exists(directory)) {
        unlink(file.path(directory, "*"), recursive = FALSE, force = TRUE)
    } else {
        dir.create(directory, showWarnings = FALSE)
    }
}

#' Calculates the number of NA entries in a vector x
#'
#' @param x - a vector of numeric values, potentially with some NAs
#' @return sum of NA entries in x
#'
#' @examples
#' HasNA(c(3,4,NA,8))
#' @export

HasNA <- function(x) sum(is.na(x))


#' Converts a factor with numeric character levels into a factor with numeric levels
#'
#' @param x - a factor with numeric character levels
#' @return a factor with levels that are numeric, not character
#'
#' @examples
#' TBD
#' @export
#'
AsNumericFactor <- function(x) {
    as.numeric(levels(x))[x]
}


#' Removes the mean from a numeric vector, removes NAs
#'
#' @param x - numeric vector
#' @return a mean centred numeric vector
#'
#' @examples
#' MeanCentered(c(2, 1, 3)) returns c(0, -1, 1)
#'
#' @export

MeanCentred <- function(x) x - mean(x, na.rm = T)


#' Calculates the percentage Coefficent of Variation from a numeric vector
#'
#' @param x - numeric vector
#' @return % Coefficient of Variation of the vector x
#'
#' @examples
#' CoVar(c(2, 1, 3)) returns 50


CoVar <- function(x) {
    100 * sd(x)/mean(x)
}

#'Renames column names in a dataframe, so each column name matching
#'columnsBy[[i]] becomes columnsTo[[i]]
#'
#'@param df A data frame with which to rename the colnames
#'@param columnsBy A vector of colnames to be replaced
#'@param columnsTo A vector of names to replace with
#'@param fill Values to insert into the now missing columns
#'
#'@export

RenameColumns <- function(df,columnsBy, columnsTo, fill=NA){

  #select the columns and data of the dataset you want to use
  dataToMove <- dplyr::select_(df,.dots = columnsBy)

  #relabel the column names to the names you want to replace
  colnames(dataToMove) <- columnsTo

  remain <-  dplyr::select_(df,.dots = setdiff(colnames(df),columnsTo))
  remain <- cbind(remain, dataToMove)

  setToFill <- setdiff(columnsBy,columnsTo)
  remain[[setToFill]] <- fill

  return(remain)
}
