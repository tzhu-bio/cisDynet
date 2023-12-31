% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TimeSeries.R
\name{getTimeRNA}
\alias{getTimeRNA}
\title{Plot the RNA dynamic along the pseudo time.}
\usage{
getTimeRNA(
  peak_time,
  peak2gene,
  rna_matrix,
  corr_cutoff = 0.4,
  Palette = NA,
  return_matrix = FALSE
)
}
\arguments{
\item{peak_time}{The result of **getTimeATAC** with the parameter return_matrix = TRUE.}

\item{peak2gene}{The rds file obtained by **getPeak2Gene**.}

\item{rna_matrix}{The RNA-seq quantification matrix.}

\item{corr_cutoff}{The cutoff of correlation to get the reliable Peak2gene links. Default: 0.4.}

\item{return_matrix}{Whether to return the result as data frame. Default: FALSE.}
}
\description{
Plot the RNA dynamic along the pseudo time.
}
\examples{
  getTimeRNA(peak_time = peak_time, peak2gene = peak2gene, rna_matrix = rna_matrix, corr_cutoff = 0.4, return_matrix = T)

}
