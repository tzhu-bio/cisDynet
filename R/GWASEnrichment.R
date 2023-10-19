# Read input count files and generate a data frame
read_files <- function(indir) {
  files <- list.files(indir, pattern = "\\.bed$", full.names = TRUE)
  data_list <- lapply(files, function(file) {
    condition_name <- tools::file_path_sans_ext(basename(file))
    data <- read.table(file, sep = "\t", header = FALSE)
    colnames(data) <- c("chr", "start", "end", condition_name)
    return(data)
  })
  count_data <- Reduce(function(x, y) merge(x, y, by = c("chr", "start", "end"), all = TRUE), data_list)
  return(count_data)
}


# Normalize counts to account for sequencing depth
normalize_counts <- function(count_data) {
  count_sums <- colSums(count_data[, !colnames(count_data) %in% c("chr", "start", "end")], na.rm = TRUE)
  norm_factors <- max(count_sums) / count_sums
  normalized_data <- count_data[, !colnames(count_data) %in% c("chr", "start", "end")] * norm_factors
  normalized_data <- cbind(count_data[, c("chr", "start", "end")], normalized_data)
  return(normalized_data)
}


# Drop the bottom 10th percentile of peaks by total count
drop_peaks <- function(norm_data, drop_quantile) {
  peak_sums <- rowSums(norm_data[, !colnames(norm_data) %in% c("chr", "start", "end")], na.rm = TRUE)
  tenth_percentile <- quantile(peak_sums, drop_quantile)
  return(norm_data[peak_sums >= tenth_percentile, ])
}

# Normalize each peak's count by the Euclidean norm across all conditions
euclidean_norm <- function(norm_data) {
  norm_data_squared <- norm_data[, -c(1:3)]^2
  euclidean_norm <- sqrt(rowSums(norm_data_squared, na.rm = TRUE))
  return(norm_data[, -c(1:3)] / euclidean_norm)
}

# Quantile normalize each condition's peak count distributions
quantile_normalisation <- function(df){
  df_rank <- apply(df,2,rank,ties.method="min")
  df_sorted <- data.frame(apply(df, 2, sort))
  df_mean <- apply(df_sorted, 1, mean)

  index_to_mean <- function(my_index, my_mean){
    return(my_mean[my_index])
  }

  df_final <- apply(df_rank, 2, index_to_mean, my_mean=df_mean)
  rownames(df_final) <- rownames(df)
  return(df_final)
}


##################################################
#               Data Normalization               #
##################################################

#' Preprocessing the data for GWAS enrichment.
#'
#' @param indir The directory containing the counts files obtained by the **getCountSplit** function.
#' @param outdir  The directory save the normalized data.
#' @param prefix  The prefix of save file.
#'
#' @return
#' @export
#'
#' @examples    GWASNorm(indir = "./gwas/indir", outdir = "./gwas_out", prefix = "Immune_Cell")
GWASNorm <- function(indir, outdir, prefix, drop_quantile = 0.1){
  logfile("Reading counts data...")
  count_data <- read_files(indir)
  logfile("Normalizing...")
  norm_data <- normalize_counts(count_data)
  norm_data <- drop_peaks(norm_data, drop_quantile = drop_quantile)
  folder_path <- outdir
  if (!dir.exists(folder_path)) {
    dir.create(folder_path)
  }
  norm_data_path <- file.path(folder_path, paste0(prefix, "_counts_normToMax.txt"))
  write.table(norm_data, norm_data_path, sep='\t',quote=F, row.names=F)
  quantile_norm_data <- quantile_normalisation(norm_data[,-c(1,2,3)])
  quantile_norm_data <- quantile_norm_data / colSums(quantile_norm_data) *15000
  quantile_norm_data <- cbind(norm_data[,c(1,2,3)], quantile_norm_data)
  quantile_norm_data_path <- file.path(folder_path, paste0(prefix, "_counts_normToMax_quantileNorm.txt"))
  write.table(quantile_norm_data, quantile_norm_data_path, sep='\t', quote=F, row.names=F)
  euclidean_norm_data_path <- file.path(folder_path, paste0(prefix, "_counts_normToMax_quantileNorm_euclideanNorm.txt"))
  euclidean_norm_data <- euclidean_norm(quantile_norm_data)
  euclidean_norm_data  <- round(euclidean_norm_data, 4)
  euclidean_norm_data <- cbind(norm_data[,c(1,2,3)],euclidean_norm_data)
  write.table(euclidean_norm_data, euclidean_norm_data_path, sep='\t', quote=F,row.names=F)
}


##################################################
#                  GWAS Enrichment               #
##################################################

#' Do enrichment analysis with GWAS variants.
#'
#' @param euclideanNorm_file The euclideanNorm file obtained by **GWASNorm** function.
#' @param snp   The SNPs list. It consists of two columns, the first being the chromosome and the second the position. Tabs between the two columns and no header.
#' @param trait The traits name you want to set.
#' @param output_path The path to save the enrichment results.
#'
#' @return
#' @export
#'
#' @examples    GWASEnrichment(euclideanNorm_file = "Immune_Cell_counts_normToMax_quantileNorm_euclideanNorm.txt",
#'                              snp = "test.snp", trait = "test", output_path = "./")
GWASEnrichment <- function(euclideanNorm_file, snp, trait, output_path){

  # Load data
  load_data <- function(file_path) {
    data <- read.table(file_path, sep = "\t", head=T)
    data <- data %>% dplyr::arrange(chr, start, end)
    return(data)
  }

  # Rank peaks within each condition
  rank_data <- function(norm_data) {
    peak_info <- norm_data[, c("chr", "start", "end"), drop = FALSE]
    rank_matrix <- lapply(norm_data[, -c(1:3), drop = FALSE], rank, ties.method = "min")
    rank_data <- cbind(peak_info, rank_matrix)
    return(rank_data)
  }
  # Find the SNP overlapping Peaks
  overlap_snp_peak <- function(snp, rank_data){
    snp_dt <- data.table::fread(snp)
    colnames(snp_dt) <- c("chromosome","position")
    quant_dt <- ranked_data
    snp_gr <- GenomicRanges::GRanges(seqnames = snp_dt$chromosome,
                                     ranges = IRanges::IRanges(start = snp_dt$position, end = snp_dt$position))
    quant_gr <- GenomicRanges::GRanges(seqnames = quant_dt$chr,
                                       ranges = IRanges::IRanges(start = quant_dt$start, end = quant_dt$end))
    overlap_idx <- GenomicRanges::findOverlaps(snp_gr, quant_gr)
    if (length(overlap_idx) > 0) {
      overlap_snp <- snp_dt[S4Vectors::queryHits(overlap_idx),]
      overlap_quant <- quant_dt[S4Vectors::subjectHits(overlap_idx),]

      # merge the data
      overlap_df <- cbind(overlap_snp, overlap_quant)
      rownames(overlap_df) <- NULL
      unique_df <- overlap_quant[!duplicated(overlap_quant[c(1,2,3)]),]
      return(list(overlap_df, unique_df))
    } else{
      return(NULL)
  }
}
  # Calculate enrichment
  calc_enrichment <- function(unique_peaks, num_peaks) {

    # Define function to calculate standard deviation
    calc_sd <- function(n, num_peaks) {
      mean_sd <- sqrt((num_peaks^2 - 1) / (12 * n))
      return(mean_sd)
    }

    # Calculate observed rank means

    observed_rank_means <- unique_peaks %>%
      select(-c("chr", "start", "end")) %>%
      colMeans()

    # Calculate number of unique peaks and p-values
    n <- nrow(unique_peaks)
    if (n == 0) {
      p_values <- rep(NA, length(observed_rank_means))
      names(p_values) <- names(observed_rank_means)
      return(list(observed_rank_means, p_values, num_peaks, n, NaN, NaN))
    } else {
      mean_mean <- (1 + num_peaks) / 2
      mean_sd <- calc_sd(n, num_peaks)
      p_values <- pnorm(observed_rank_means, mean = mean_mean, sd = mean_sd, lower.tail = FALSE)
      #p_values <- format(p_values, scientific = TRUE)
      names(p_values) <- names(observed_rank_means)
      return(list(observed_rank_means, p_values, num_peaks, n, mean_mean, mean_sd))
    }
  }

  # Write log file
  write_log <- function(enrich_res, out_path, trait){
    fileConn<- file(sprintf("%s/%s_Enrichment.log", out_path, trait))
    writeLines(c(sprintf("Total number of peaks: %s",enrich_res[[3]]),
                 sprintf("Number of overlapping peaks: %s",enrich_res[[4]]),
                 sprintf("Number of SNPs overlapping peaks: %s", nrow(overlap_matrix) -1),
                 sprintf("Distribution mean: %s", enrich_res[[5]]),
                 sprintf("Distribution sd: %s", enrich_res[[6]])), fileConn)
    close(fileConn)
  }
  # Load data
  norm_data <- load_data(euclideanNorm_file)

  # Rank peaks
  ranked_data <- rank_data(norm_data)

  # Identify overlapping peak-SNP pairs and unique peaks
  overlap_results <- overlap_snp_peak(snp, ranked_data)
  if (is.null(overlap_results)){
    print("No 0verlapping SNPs!")
  } else{
  overlap_matrix <- overlap_results[[1]]
  unique_peaks <- overlap_results[[2]]

  # Calculate enrichment
  num_peaks <- nrow(ranked_data)
  enrichment_results <- calc_enrichment(unique_peaks, num_peaks)
  observed_rank_means <- enrichment_results[[1]]
  p_values <- enrichment_results[[2]]
  num_peaks <- enrichment_results[[3]]
  n <- enrichment_results[[4]]
  enrichment_scores <- enrichment_results[[5]]
  expected_rank_means <- enrichment_results[[6]]
  write_log(enrichment_results, output_path, trait)
  write.table(data.frame(p_values), sprintf("%s/%s_Enrichment_Pvalues.txt", output_path, trait),sep='\t',col.names=F, row.names=T,quote=F)
  write.table(data.frame(overlap_matrix), sprintf("%s/%s_SNP_Overlap_Peaks.txt", output_path, trait),sep='\t',col.names=T, row.names=F,quote=F)
  write.table(data.frame(unique_peaks), sprintf("%s/%s_Unique_Peaks.txt", output_path, trait),sep='\t',col.names=T, row.names=F,quote=F)
}
}
##################################################
#     Integrate the Enrichment Results.          #
##################################################

#' The function to integrate the GWAS enrichment results.
#'
#' @param enrichment_result_path The GWAS enrichment result path you set in **GWASEnrichment**
#' @param return_matrix Whether to plot.
#'
#' @return
#' @export
#'
#' @examples  getGER(enrichment_result_path = "./out")
getGER <- function(enrichment_result_path, return_matrix = FALSE){
  lst <- list()
  for (i in Sys.glob(sprintf("%s/*_Enrichment_Pvalues.txt", enrichment_result_path))){
    pval <- read.table(i, row.names = 1)
    colnames(pval) <- gsub("_Enrichment_Pvalues.txt", "", basename(i))
    lst[[i]] <- pval
  }
  all_res <- rlist::list.cbind(lst)
  all_res <- -log10(all_res)
  df <- data.frame(t(all_res))
  p <- ComplexHeatmap::pheatmap(as.matrix(df), scale="row", color=c(rep("white",5), paletteer::paletteer_d("grDevices::blues9")[1:7]),
                                cellwidth=16, cellheight=16,cluster_row=T,border_color = "#767F8B",
                                display_numbers = ifelse(df > 1.3, "*", ""),
                                fontsize_number=20,treeheight_row=30,treeheight_col=30, name = "Z-score")
  if(return_matrix){
    return(all_res)
  }else{
    return(p)
  }
}


##################################################
#          Plot the Enrichment Results.          #
##################################################

#' A function to plot the GWAS Enrichment Results.
#'
#' @param result The GWAS enrichment result obtained by **getGER**.
#' @param plot_type c("radar","lollipop").
#'
#' @return
#' @export
#'
#' @examples
plotGER <- function(result, plot_type = "lollipop"){
  col1 <- paletteer::paletteer_d("ggthemes::Tableau_10")
  col2 <- paletteer::paletteer_d("ggthemes::Classic_Cyclic")
  result$sample <- rownames(result)
  long_df <- tidyr::pivot_longer(data =result,cols = -sample, names_to = "GWAS",values_to = "P") %>%as.data.frame()
  if(plot_type=="lollipop"){
    p <-  ggplot2::ggplot(long_df, aes(x = sample, y = P, color=GWAS)) +
      ggplot2::geom_segment(aes(xend = sample, yend = 0), color = "grey", size = 0.6) +
      ggplot2::geom_point(size = 3) +
      ggplot2::coord_flip() + ggplot2::facet_wrap(~GWAS) +
      ggplot2::scale_color_manual(values = c(col1, col2)) +
      ggpubr::theme_pubr() + ggplot2::geom_hline(yintercept = -log10(0.05), color = "#E7298A", linetype = "dashed")
  }else{
    p <- ggplot2::ggplot(long_df, aes(x = sample, y = P, color=GWAS, group = GWAS)) + ggplot2::geom_point(size=3, alpha=1)+
      ggplot2::geom_polygon(alpha=0.1,linewidth=0.5) + ggplot2::coord_polar() + ggplot2::labs(x="",y="-log10(Pvalue)") + ggplot2::theme_bw()+
      ggplot2::theme(axis.text = element_text(color = "black"),
                     panel.border = element_rect(size = 0),
                     panel.grid = element_line(linewidth = 0.5),
                     legend.text.align = 0) + ggplot2::geom_hline(yintercept = -log10(0.05), color = "#E7298A", linetype = "dashed")+
      ggplot2::scale_color_manual(values = c(col1, col2))
  }
  return(p)
}
