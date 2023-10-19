#' A function to identify the super enahncers.
#'
#' @param sample_list A vectot contains sample names.
#' @param peak_path The path contains the sample peak files.
#' @param peak_suffix The suffix of peak files. Default: ".narrowPeak.bed".
#' @param N_top The top number of peakset return. Default: 1000.
#' @param blacklist A blacklist would prevent us from identifying super enhancers in those areas. Default: NA.
#'
#' @return
#' @export
#'
#' @examples  getSuperEnhancer(sample_list=c("s1","s2"),
#'                             peak_path="~/peak/",
#'                             blacklist="~/peak/balcklist.bed")
#'
getSuperEnhancer <- function(sample_list, peak_path, peak_suffix=".narrowPeak.bed", N_top=1000, blacklist = NA){
  peak <- list()
  for (i in sample_list){
    peak[[i]] <- read.table(sprintf("%s/%s%s",peak_path,i, peak_suffix))
  }
  all_peaks <- dplyr::bind_rows(peak)[,c(1:3)]
  colnames(all_peaks) <- c("chrom", "start", "end")
  merged_peaks <- valr::bed_merge(all_peaks)
  #res <- as.data.frame(valr::bed_cluster(merged_peaks))
  merged_peaks$peak_len <- merged_peaks$end - merged_peaks$start
  merged_peaks$peak <- sprintf("%s:%s-%s",merged_peaks$chrom, merged_peaks$start, merged_peaks$end)
  final <- merged_peaks[order(merged_peaks$peak_len,decreasing=T),]
  if(is.na(blacklist)){
    return(as.data.frame(head(final,n = N_top)))
  }
  else{
    black <- valr::read_bed(blacklist)
    final_f <- valr::bed_intersect(final, black, invert = T)
    return(as.data.frame(head(final_f,n = N_top)))
  }
}
