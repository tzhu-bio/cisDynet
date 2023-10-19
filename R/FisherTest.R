#' A function to calculate the two group peaks enrichment based on Fisher Test.
#'
#' @param peak_a  The peak sets A you want to enrich with peak sets B.
#' @param peak_b The peak sets B you want to enrich with peak sets A.
#' @param genome_size The genome size file. Two columns: chr chr_length.
#'
#' @return
#' @export
#'
#' @examples calFisherPeaks(peak_a = c("a1.bed", "a2.bed", "a3.bed"),
#'                          peak_b = c("b1.bed", "b2.bed", "b3.bed"), genome_size= "genome.txt")
#'
calFisherPeaks <- function(peak_a, peak_b){
  checkGeAnno()
  genome <- CATAnno$genome
  if (length(peak_a) == 1){
    peaka <- valr::read_bed(peak_a)
    res <- list()
    for (pb in peak_b){
      peakb <- valr::read_bed(pb)
      temp <- valr::bed_fisher(peaka, peakb, genome)
      temp$Peak1 <- basename(pa)
      temp$Peak2 <- basename(pb)
      #temp$group <- sprintf("%s-%s",basename(peak_a),basename(pb))
      res[[pb]] <- temp
    }
    final <- as.data.frame(dplyr::bind_rows(res))
    return(final)
  }
  else{
    logfile(sprintf("Fisher testing for %s x %s.",length(peak_a), length(peak_b)))
    all <- list()
    for (pa in peak_a){
      peaka <- valr::read_bed(pa)
      res <- list()
      for (pb in peak_b){
        peakb <- valr::read_bed(pb)
        temp <- valr::bed_fisher(peaka, peakb, genome)
        temp$Peak1 <- basename(pa)
        temp$Peak2 <- basename(pb)
        #temp$group <- sprintf("%s-%s",basename(pa),basename(pb))
        res[[pb]] <- temp
      }
      all[[pa]] <- as.data.frame(dplyr::bind_rows(res))
    }
    final <- as.data.frame(dplyr::bind_rows(all))
    return(final)
  }
}


#' A function to plot the enrichment rsult based on **calFisherPeaks**.
#'
#' @param fisher_res The enrichment result obtained by **calFisherPeaks**.
#' @param scale_size The relative size between big and small dot size. default: c(2,12)
#'
#' @return
#' @export
#'
#' @examples  plotEnrich(res)
#'
plotEnrich <- function(fisher_res, scale_size=c(2,12)){
  df <- fisher_res
  df[sapply(df, is.infinite)] <- 100
  df$p.value <- ifelse(df$p.value==0,1e-100,df$p.value)
  df$score <- -log10(df$p.value)
  p <- ggplot2::ggplot(df, aes(x=Peak1, y=Peak2)) + ggplot2::geom_point(aes(size= estimate, fill = score),shape=21,alpha=0.9)+
       ggplot2::scale_fill_gradientn(colours =paletteer::paletteer_d("ggsci::pink_material")[1:6]) + ggplot2::theme_bw()+
       ggplot2::scale_size_continuous(range = scale_size) + ggplot2::ylab('') + ggplot2::xlab("") +
       ggplot2::theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggplot2::labs(fill = "-log10(p.value)", size="Ratio")
  return(p)
}

