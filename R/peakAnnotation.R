
##############################################################################
#         Annotating the peaks with ChIPSeeker.                              #
##############################################################################
#' Annotating the peaks with ChIPSeeker.
#' @param peak_path The path containing the peak bed files.
#' @param sample Sample names.
#' @param peak_suffix Provide a suffix for the peak file. Default "_peaks_unique.narrowPeak.bed".
#' @param annoDb_name An annotation package for your species.
#' @param plot_percent Plotting proportions rather than specific numbers. Default: TRUE.
#'
#' @return
#' @export
#'
#' @examples annoPeaks("~/H3ac/peaks/",c("S1","S2","S3","S4","S5","S6"),
#'                     "_peaks_unique.narrowPeak.bed","org.Hs.eg.db",plot_percent=F)
#'
annoPeaks <- function(peak_path, sample, peak_suffix="_peaks_unique.narrowPeak.bed", annoDb_name, plot_percent=T){
  checkGeAnno()
  lst <- list()
  for (i in sample){
    lst[[i]] <-  sprintf("%s/%s%s", peak_path, i, peak_suffix)
  }
  peak_files <- unlist(lst)
  txdb <- GenomicFeatures::makeTxDbFromGFF(CATAnno$gtf)
  suppressMessages({
  peak_annot_list <- lapply(lst, function(x) ChIPseeker::annotatePeak(x, TxDb = txdb, annoDb = annoDb_name, verbose=FALSE))
  })
  names(peak_annot_list) <- names(peak_files)
  peak_annot_summary_list <- lapply(peak_annot_list, function(x) {
    anno_df <- as.data.frame(x@anno)
    anno_df$annotation <- ifelse(grepl("Exon",anno_df$annotation),"Exon",anno_df$annotation)
    anno_df$annotation <- ifelse(grepl("Intron",anno_df$annotation),"Intron",anno_df$annotation)
    anno_df <- anno_df %>% dplyr::group_by(annotation) %>% dplyr::summarise(count = n()) %>% as.data.frame()
    return(anno_df)
  })
  peak_annot_summary_df <- dplyr::bind_rows(peak_annot_summary_list,.id = "sample")
  if(plot_percent){
    peak_annot_summary_df <- dplyr::group_by(peak_annot_summary_df, sample) %>% dplyr::mutate(percent = count/sum(count))
    p <- ggplot2::ggplot(peak_annot_summary_df, aes(x = sample, y = percent*100, fill = annotation)) +
         ggplot2::geom_bar(stat = "identity") +
         ggplot2::theme_bw() + ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
         ggplot2::labs(x = "", y = "Percent (%)", fill = "Annotation")+
         ggplot2::scale_fill_manual(values=paletteer::paletteer_d("ggthemes::Tableau_20"))

  }else{
    p <- ggplot2::ggplot(peak_annot_summary_df, aes(x = sample, y = count/1000, fill = annotation)) +
         ggplot2::geom_bar(stat = "identity") +
         ggplot2::theme_bw() + ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
         ggplot2::labs(x = "", y = "Count (x 1000)", fill = "Annotation")+
         ggplot2::scale_fill_manual(values=paletteer::paletteer_d("ggthemes::Tableau_20"))
  }
  return(p)
}

##############################################################################
#              Plot the distance of peak summit to TSS                       #
##############################################################################

#' Plot the distance of peak summit to TSS.
#' @param peak_path  The path containing the peak bed files.
#' @param sample  Sample names.
#' @param suffix Provide a suffix for the peak file.
#' @return
#' @export
#'
#' @examples   plotSummitDis("~/peaks", c("Control1","Control2","Control3","Control4","Control5","IDD1","IDD2","IDD3","IDD4","IDD5"),
#'                           ".narrowPeak.bed")
#'
plotSummitDis <- function(peak_path, sample, suffix){
  checkGeAnno()
  sample_lst <- list()
  for (i in sample){
    sample_lst[[i]] <- sprintf("%s/%s%s", peak_path, i, suffix)
  }
  tss <- CATAnno$tss
  dis <- lapply(sample_lst, function(x) {
    a <- valr::read_narrowpeak(x)
    a$start <- a$start + a$peak
    a$end <- a$start + 1
    a <- a[,c(1:4)]
    res <- valr::bed_closest(a, tss)
    return(res)
  })
  all_res <- as.data.frame(dplyr::bind_rows(dis, .id = "group"))
  all_res$.dist <- abs(all_res$.dist) + 1
  if(length(sample) > 1){
    p <- ggplot2::ggplot(all_res, aes(x = .dist, y = group, fill=group)) +
      ggridges::geom_density_ridges() +
      ggplot2::scale_y_discrete(expand = c(0, 0)) +
      ggplot2::scale_x_continuous(expand = c(0, 0)) +
      ggplot2::coord_cartesian(clip = "off") + ggplot2::xlab("Distance to closest gene (bp)") + ggplot2::ylab("")+
      ggpubr::theme_pubr() +ggplot2::scale_fill_manual(values=c(as.character(paletteer::paletteer_d("ggthemes::Hue_Circle")),"#1F77B4"))+
      ggplot2::scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),labels = scales::trans_format("log10", scales::math_format(10^.x)))
    return(p)
  }else{
    p <- ggplot2::ggplot(all_res, aes(x=.dist,color=group)) +
      ggplot2::geom_density()+ggpubr::theme_pubr(border=T)+
      ggplot2::xlab("Distance to closest gene (bp)") +ggplot2::ylab("Density")+
      ggplot2::scale_color_manual(values=paletteer::paletteer_d("ggthemes::Tableau_10"))
    return(p)
  }
}

##############################################################################
#      Determining the peak category (Proximal/Distal/Intragenic)            #
##############################################################################

#' Determining the peak category (Proximal/Distal/Intragenic).
#' @param peak_path  The path containing the peak bed files.
#' @param sample Sample names.
#' @param suffix Provide a suffix for the peak file. Default: "_peaks_unique.narrowPeak.bed".
#' @param cutoff The distance of peak summit to TSS (bp) to determine the proximal / distal.
#' @param save_path Save file path.
#' @param save_name Save file name.
#' @param tss_flank How many bp downstream of the TSS is considered proximal.
#'
#' @return
#' @export
#'
#' @examples plotPDI("F:/CAT/example/peaks/", c("Bulk_B", "Mem_B","Naive_B"),
#'                   "_peaks_unique.narrowPeak.bed",
#'                   "F:/CAT/example/hg19_gene_standard.bed", cutoff=3000, tss_flank=1000)
#'
#'
plotPDI <- function(peak_path, sample, suffix="_peaks_unique.narrowPeak.bed", tss_flank, cutoff, save_path=NA, save_name=NA){
  checkGeAnno()
  sample_lst <- list()
  for (i in sample){
    sample_lst[[i]] <- sprintf("%s/%s%s", peak_path, i, suffix)
  }
  gene <- CATAnno$gene
  tss <- CATAnno$tss
  gene$start <- ifelse(gene$strand == "+", gene$start + tss_flank, gene$start)
  gene$end <- ifelse(gene$strand == "+", gene$end, gene$end - tss_flank)
  genome <- CATAnno$genome
  promoter <- valr::bed_flank(tss, genome, left = cutoff, right = tss_flank, strand = TRUE)
  dis <- lapply(sample_lst, function(x) {
    a <- valr::read_narrowpeak(x)
    a$start <- a$start + a$peak
    a$end <- a$start + 1
    b <- a[,c(1:4)]
    pro <- valr::bed_intersect(b, promoter)[,c(1:4)]
    pro <-  pro[!duplicated(pro[,c(1,2,3)]),]
    pro$Type <- "Proximal"
    colnames(pro)[2:4] <- c("start","end","name")
    intra <- valr::bed_intersect(b, promoter,invert = TRUE)[,c(1:4)]
    intra <-  intra[!duplicated(intra[,c(1,2,3)]),]
    intra_res <- valr::bed_intersect(intra, gene)[,c(1:4)]
    intra_res <- intra_res[!duplicated(intra_res[,c(1,2,3)]),]
    intra_res$Type <- "Intragenic"
    colnames(intra_res)[2:4] <- c("start","end","name")
    re <- valr::bed_intersect(intra, intra_res,invert = TRUE)[,c(1:4)]
    re <-  re[!duplicated(re[,c(1,2,3)]),]
    re$Type <- "Distal"
    intra_pe <- rbind(re[,c(1,2,3,4,5)], intra_res, pro)
    #num <- as.data.frame(table(intra_pe$Type))
    return(intra_pe)
  })
  all_res <- as.data.frame(dplyr::bind_rows(dis, .id = "group"))
  fre <- as.data.frame(table(all_res$group, all_res$Type))
  if(!is.na(save_path)){
    write.table(all_res, sprintf("%s/%s_Proximal_Distal_Intragenic_Annotation.tsv", save_path, save_name),sep='\t',row.names=F,quote=F)
  }
  p <- ggplot2::ggplot(fre, aes(x = Var1, y = Freq/1000, fill = Var2)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggplot2::labs(x = "", y = "Frequency (x 1000)", fill = "Annotation")+
    ggplot2::scale_fill_manual(values=c("#1F78B4","#33A02C","#FF7F00"))
  return(p)
}


##############################################################################
#    Determining the merged peak category (Proximal/Distal/Intragenic)       #
##############################################################################

#' Determining the merged peak category (Proximal/Distal/Intragenic).
#' Determining the merged peak category (Proximal / Distal / Intragenic).
#' @param quant_data Quantified data obtained by **quantification**.
#' @param tss_flank  How many bp downstream of the TSS is considered proximal.
#' @param cutoff The distance of peak summit to TSS (bp) to determine the proximal / distal.
#' @param save_path Save file path.
#' @param save_name Save file name.
#'
#' @return
#' @export
#'
#' @examples    annoMergedPeaks("F:/CAT/example/ATAC_CPM_Norm_Data.tsv,tss_flank=1000, cutoff=3000)
#'
annoMergedPeaks <- function(quant_data, tss_flank, cutoff, save_path=NA, save_name=NA){
  checkGeAnno()
  peak <- read.table(quant_data, header = T, row.names = 1)
  merged<- data.frame(peak=rownames(peak))
  merged_peaks <- merged %>% tidyr::separate(peak, c("chrom", "tem1"), ":") %>% tidyr::separate(tem1, c("start", "end"), "-")
  merged_peaks$start <- as.integer(merged_peaks$start)
  merged_peaks$end <- as.integer(merged_peaks$end)

  gene <- CATAnno$gene
  tss <- CATAnno$tss
  gene$start <- ifelse(gene$strand == "+", gene$start + tss_flank, gene$start)
  gene$end <- ifelse(gene$strand == "+", gene$end, gene$end - tss_flank)
  genome <- CATAnno$genome
  promoter <- valr::bed_flank(tss, genome, left = cutoff, right = tss_flank, strand = TRUE)

  merged_peaks$Start <- merged_peaks$start
  merged_peaks$End <- merged_peaks$end
  merged_peaks$start <- as.integer((merged_peaks$start + merged_peaks$end) / 2)
  merged_peaks$end <- merged_peaks$start + 1
  pro <- valr::bed_intersect(merged_peaks, promoter)[,c(1:5)]
  pro <-  pro[!duplicated(pro[,c(1,2,3)]),]
  pro$Type <- "Proximal"
  colnames(pro)[2:5] <- c("start","end","Start","End")
  intra <- valr::bed_intersect(merged_peaks, promoter,invert = TRUE)
  intra <-  intra[!duplicated(intra[,c(1,2,3)]),]
  intra_res <- valr::bed_intersect(intra, gene)[,c(1:5)]
  intra_res <- intra_res[!duplicated(intra_res[,c(1,2,3)]),]
  intra_res$Type <- "Intragenic"
  colnames(intra_res)[2:5] <- c("start","end","Start","End")
  re <- valr::bed_intersect(intra, intra_res,invert = TRUE)[,c(1:5)]
  re <-  re[!duplicated(re[,c(1,2,3)]),]
  re$Type <- "Distal"
  intra_pe <- rbind(re, intra_res, pro)
  intra_pe <- as.data.frame(intra_pe[!duplicated(intra_pe[,c(1,2,3)]), ])
  intra_pe <- intra_pe[,c(1,4,5,6,2)]
  colnames(intra_pe) <- c("Chromosome","Start","End","Type","summit")
  rownames(intra_pe) <- sprintf("%s:%s-%s", intra_pe$Chromosome, intra_pe$Start, intra_pe$End)
  if(!is.na(save_path)){
  write.table(sprintf("%s/%s_Merged_Peaks_Annotations.tsv", save_path, save_name), sep='\t', quote=F, col.names=T, row.names=T)
  }
  return(intra_pe)
}
