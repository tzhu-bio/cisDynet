% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FootprintAnalysis.R
\name{plotFootprintScore}
\alias{plotFootprintScore}
\title{A function to plot all samples footprint with scatter / heatmap.}
\usage{
plotFootprintScore(
  score_result,
  cluster_N = 6,
  label_motif,
  limit_Zscore = FALSE
)
}
\arguments{
\item{score_result}{The footprint score result obtained by **getFootprintScore**.}

\item{cluster_N}{The cluster number for heatmap. Default: 6.}

\item{label_motif}{The motif list to label in the right of the heatmap.}

\item{limit_Zscore}{Whether to limit zscore to -2 to 2. Default: FALSE.}
}
\description{
A function to plot all samples footprint with scatter / heatmap.
}
\examples{
  plotFootprintScore(score_result=res, label_motif("motif1","motif2","motif3"...))

}
