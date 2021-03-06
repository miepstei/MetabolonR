% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ManuipulateDataset.R
\name{MetDataWideToLong}
\alias{MetDataWideToLong}
\title{Convenience function for transforming a Wide count Data Frame to a Long count Data Frame}
\usage{
MetDataWideToLong(cMetDataWide)
}
\arguments{
\item{cMetDataWide}{a count matrix of Metabolon data in Long format}
}
\value{
a count dataframe of Metabolon data in Long format
}
\description{
Convenience function for transforming a Wide count Data Frame to a Long count Data Frame
}
\examples{
wideMetCounts <- MetDataLongToWide(longMetCounts)
}

