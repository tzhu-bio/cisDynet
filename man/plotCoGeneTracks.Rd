% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotTracks.R
\name{plotCoGeneTracks}
\alias{plotCoGeneTracks}
\title{Plot specific gene ATAC signal and co-accessible with tracks.}
\usage{
plotCoGeneTracks(
  samples_path,
  samples_suffix,
  gene_name,
  left = 30000,
  right = 30000,
  peaks,
  colist,
  coaccess_cutoff = 0.4,
  line_size = 0.8,
  curvature = 0.3,
  color = NA,
  back.color = TRUE
)
}
\arguments{
\item{samples_path}{The bigwig (bw) file folder.}

\item{samples_suffix}{The suffix of the samples. e.g. ".cpm.bw".}

\item{gene_name}{The gene name.}

\item{left}{The gene left flanking to plot. Default: 30000.}

\item{right}{The gene right flanking to plot. Default: 30000.}

\item{peaks}{The merged peaks.}

\item{colist}{The co-accessible file obtained by **getCoaccessible**.}

\item{coaccess_cutoff}{The correlation cutoff of co-accessible.  Default: 0.4.}

\item{line_size}{The curve line with.  Default: 0.8.}

\item{curvature}{The degree of curvature. Default: 0.3.}

\item{color}{Color list to plot tracks.}

\item{back.color}{Whether to plot backgroud color.  Default: TRUE.}
}
\description{
Plot specific gene ATAC signal and co-accessible with tracks.
}
\examples{
 plotCoGeneTracks(samples_path = "./signal", samples_suffix = ".cpm.bw", gene_name="ENSG00000103888" , peaks="final_res/meged.peaks.bed",colist="./coaccess/tsv")
}
