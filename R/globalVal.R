
#' A function to add the necessary annotation for **CAT**.
#'
#' @param gene_bed The bed file for genes annotation.
#' @param gtf The GTF file.
#' @param genome_size The genome size file.
#'
#' @return
#' @export
#'
#' @examples   addAnnotation(gene_bed="./gene.bed", gtf="hg19.gtf", genome_size="hg19.chrom.size")
addAnnotation <- function(gene_bed, gtf, genome_size) {
  # gene <- valr::read_bed(gene_bed, n_fields = 6)
  gene <- valr::read_bed(gene_bed)  ## valr 0.7.0 update
  tss <- gene
  tss$start <- ifelse(tss$strand=="+", tss$start, tss$end)
  tss$end <- ifelse(tss$strand=="+", tss$start + 1, tss$end+1)
  genome_size <- valr::read_genome(genome_size)
  CATAnno <- list(gene, tss, genome_size, gtf)
  names(CATAnno) <- c("gene", "tss", "genome","gtf")
  assign("CATAnno", CATAnno, envir = .GlobalEnv)
}

##############################################################

#' A function to add necessary annotation for foorprinting analysis.
#'
#' @param corrected_signal The absolute path of **ATACorrect** folder.
#' @param bindetect_result The absolute path of **BINDetect** folder.
#'
#' @return
#' @export
#'
#' @examples  addFootprint(corrected_signal = "./signal", bindetect_result = "./res/")
addFootprint <- function(corrected_signal, bindetect_result){
  FTAnno <- list(corrected_signal, bindetect_result)
  names(FTAnno) <- c("signal", "bindetect")
  assign("FTAnno", FTAnno, envir = .GlobalEnv)
}
