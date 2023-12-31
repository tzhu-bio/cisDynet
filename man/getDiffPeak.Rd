% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DifferentialAnalysis.R
\name{getDiffPeak}
\alias{getDiffPeak}
\title{A function to identify the differential peaks between two groups.}
\usage{
getDiffPeak(
  norm_data,
  condition,
  rep_N,
  control,
  experment,
  log2fc = 1,
  padj = 0.05
)
}
\arguments{
\item{norm_data}{The normalizated data obtained by **quantification**.}

\item{condition}{Your designed conditions, like the control and experiment.}

\item{rep_N}{The replicates number.}

\item{control}{Which condition is the control.}

\item{experment}{Which condition is the experment.}

\item{log2fc}{The log2 fold change to determine the significant differential peaks.}

\item{padj}{The adjusted P values to determine the significant differential peaks.}
}
\description{
A function to identify the differential peaks between two groups.
}
\examples{
  diff_peak <- getDiffPeak(data, condition=c("control","exp1"), rep_N = 2, control="control", experment="exp1")

}
