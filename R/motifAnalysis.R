
#' Get the motifs from JASPAR database (2018).
#'
#' @param species Species.
#' @param collection Collection type.
#' @param ...
#'
#' @return
#' @export
#'
#' @examples  getJasparMotifs2018(species = "Homo sapiens", collection = "CORE")
getJasparMotifs2018 <- function(species = "Homo sapiens", collection = "CORE", ...) {
  opts <- list()
  opts["species"] <- species
  opts["collection"] <- collection
  opts <- c(opts, list(...))
  out <- TFBSTools::getMatrixSet(JASPAR2018::JASPAR2018, opts)
  if (!isTRUE(all.equal(TFBSTools::name(out), names(out))))
    names(out) <- paste(names(out), TFBSTools::name(out), sep = "_")
  return(out)
}


#' Make SummarizedExperiment for chromVAR input.
#'
#' @param count_data The raw count data obtained by **getCounts**.
#'
#' @return
#' @export
#'
#' @examples makeSummarizedExperiment(count_data)
#'
makeSummarizedExperiment <- function(count_data){
  peak <- data.frame(peak = rownames(count_data))
  peak$chr <- sapply(strsplit(peak$peak,":"), `[`, 1 )
  peak$bed <- sapply(strsplit(peak$peak,":"), `[`, 2 )
  mat2 <- peak %>% tidyr::separate(bed, c("start", "end"), "-") %>% dplyr::select(-peak)
  peak_range <- GenomicRanges::makeGRangesFromDataFrame(mat2)
  logfile("Making SummarizedExperiment object...")
  fragment_counts <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = as.matrix(count_data)),
                                                                rowRanges = peak_range)
  return(fragment_counts)
}


#' Rank the motif deviation score.
#'
#' @param deviation_result   The deviation result obtained by **chromVAR**.
#' @param label_topN  The top N to label. Default is 30.
#' @param return_matrix  Wheather to return ranked matrix. Default is FALSE.
#'
#' @return
#' @export
#'
#' @examples   getRankMotif(deviation_result = res)
#'
getRankMotif <- function(deviation_result, label_topN = 30, return_matrix = F){
  logfile("Computing the variability...")
  variability <- chromVAR::computeVariability(deviation_result)
  res <- variability[order(variability$variability, decreasing = T),]
  res$index <- 1:nrow(res)
  res$group <- ifelse(res$index <= 30,"yes","no")
  p <- ggplot2::ggplot(res, aes(x=index, y=variability)) +
       ggplot2::geom_point(aes(size=variability,color= group), shape = 1) +
       ggrepel::geom_text_repel(data= res[res$index <= 30,], aes(label = name), size = 3, min.segment.length = 0,seed = 42,
                 box.padding = 0.5, max.overlaps = Inf, arrow = arrow(length = unit(0.010, "npc")), nudge_x = .15, nudge_y = .5, color = "grey50") +
       ggpubr::theme_pubr() + ggplot2::xlab("Sorted TF") + ggplot2::ylab("Variability")
  if(return_matrix){
    return(res[res$group=="yes", ])
  }else{
    return(p)
  }
}


