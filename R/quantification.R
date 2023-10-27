## Merge the sample peaks, quantification and normalization.
#'
#'
#' @param sample_list Input a vector containing the names of the samples.
#' @param cut_path The cuts directory obtained by running snakmake.
#' @param peak_path The cuts directory obtained by running snakmake.
#' @param save_file_path Storing the final data.
#' @param peak_suffix The suffix of peak file. Default: "_peaks_unique.narrowPeak.bed".
#'
#' @return
#' @export
#'
#' @examples  quantification(sample_list=c("s1","s2","s3","s4"), cut_path = "~/snakemake/cuts/", peak_path = "~/snakemake/peaks/", save_file_path = "~/quant/res")
#'
#'
quantification <- function(sample_list, cut_path, peak_path, peak_suffix="_peaks_unique.narrowPeak.bed", save_file_path=NA){

  peak <- list()
  peak <- lapply(sample_list, function(x){
    read.table(sprintf("%s/%s%s", peak_path, x, peak_suffix))
  })
  all_peaks <- dplyr::bind_rows(peak)[,c(1:3)]
  colnames(all_peaks) <- c("chrom", "start", "end")

  merged_peaks <- valr::bed_merge(all_peaks)
  coverage <- lapply(sample_list, function(x) {
    cut <- valr::read_bed(sprintf("%s/%s_q30_cut_sites.bed",cut_path,x))
    colnames(cut)[4] <- "name"
    res <- valr::bed_map(merged_peaks,cut,sum=sum(name))[4]
    return(res)
  })
  all_cov <- rlist::list.cbind(coverage)
  colnames(all_cov) <- sample_list
  all_cov[is.na(all_cov)] <- 0
  peak_len <- merged_peaks$end - merged_peaks$start
  norm_data <- quantile_normalization(as.matrix(all_cov / peak_len))
  rownames(norm_data) <- sprintf("%s:%s-%s",merged_peaks$chrom, merged_peaks$start,merged_peaks$end)
  cpm_data <- sweep(norm_data,2,colSums(norm_data),`/`) * 1000000
  cpm_data <- cpm_data[rowSums(cpm_data)!=0, ]
  if(!is.na(save_file_path)){
    write.table(cpm_data, file=sprintf("%s/ATAC_CPM_Norm_Data.tsv",save_file_path), sep='\t',quote=F,col.names=T,row.names=T)
  }
  return(cpm_data)
}


#' A function to get the count file for **chromVAR**.
#'
#' @param sample_list Input a vector containing the names of the samples.
#' @param cut_path The cuts directory obtained by running snakmake.
#' @param peak_path The cuts directory obtained by running snakmake.
#' @param peak_suffix The suffix of the peaks. Default: "_peaks_unique.narrowPeak.bed".
#' @param save_file_path The path to save files.
#'
#' @return
#' @export
#'
#' @examples   getCount(sample_list=c("s1","s2","s3","s4"), cut_path = "~/snakemake/cuts/", peak_path = "~/snakemake/peaks/", save_file_path = "~/count")
getCount <- function(sample_list, cut_path, peak_path, peak_suffix="_peaks_unique.narrowPeak.bed", save_file_path=NA){

  peak <- list()
  peak <- lapply(sample_list, function(x){
    read.table(sprintf("%s/%s%s", peak_path, x, peak_suffix))
  })
  all_peaks <- dplyr::bind_rows(peak)[,c(1:3)]
  colnames(all_peaks) <- c("chrom", "start", "end")

  merged_peaks <- valr::bed_merge(all_peaks)
  coverage <- lapply(sample_list, function(x) {
    cut <- valr::read_bed(sprintf("%s/%s_q30_cut_sites.bed",cut_path,x))
    res <- valr::bed_map(merged_peaks,cut,sum=sum(X4))[4]
    return(res)
  })
  all_cov <- rlist::list.cbind(coverage)
  colnames(all_cov) <- sample_list
  all_cov[is.na(all_cov)] <- 0
  rownames(all_cov) <- sprintf("%s:%s-%s",merged_peaks$chrom, merged_peaks$start,merged_peaks$end)
  all_cov <- all_cov[rowSums(all_cov) !=0, ]
  if(!is.na(save_file_path)){
    write.table(all_cov, file=sprintf("%s/ATAC_Counts_Data.tsv",save_file_path), sep='\t',quote=F,col.names=T,row.names=T)
  }
  return(all_cov)
}


#' A function to get the count file for **GWASNorm**.
#'
#' @param sample_list Input a vector containing the names of the samples.
#' @param cut_path The cuts directory obtained by running snakmake.
#' @param peak_path The cuts directory obtained by running snakmake.
#' @param peak_suffix The suffix of the peaks. Default: "_peaks_unique.narrowPeak.bed".
#' @param save_file_path The path to save files.
#'
#' @return
#' @export
#'
#' @examples  getCountSplit(sample_list=c("s1","s2","s3","s4"), cut_path = "~/snakemake/cuts/", peak_path = "~/snakemake/peaks/", save_file_path = "~/count")
getCountSplit <- function(sample_list, cut_path, peak_path, peak_suffix="_peaks_unique.narrowPeak.bed", save_file_path){
  peak <- list()
  peak <- lapply(sample_list, function(x){
    read.table(sprintf("%s/%s%s", peak_path, x, peak_suffix))
  })
  all_peaks <- dplyr::bind_rows(peak)[,c(1:3)]
  colnames(all_peaks) <- c("chrom", "start", "end")

  merged_peaks <- valr::bed_merge(all_peaks)
  coverage <- lapply(sample_list, function(x) {
    cut <- valr::read_bed(sprintf("%s/%s_q30_cut_sites.bed",cut_path,x))
    res <- valr::bed_map(merged_peaks,cut,sum=sum(X4))[4]
    return(res)
  })
  all_cov <- rlist::list.cbind(coverage)
  colnames(all_cov) <- sample_list
  all_cov[is.na(all_cov)] <- 0
  count_res <- cbind(merged_peaks, all_cov)
  lapply(sample_list, function(x) {
    temp <- count_res[,c("chrom","start","end",x)]
    write.table(temp, sprintf("%s/%s.bed", save_file_path, x),sep='\t', quote=F, col.names = F, row.names = F)
  })
  return(invisible(NULL))
}
