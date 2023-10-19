#' A function to identify the differential peaks between two groups.
#'
#' @param norm_data The normalizated data obtained by **quantification**.
#' @param condition Your designed conditions, like the control and experiment.
#' @param control Which condition is the control.
#' @param experment Which condition is the experment.
#' @param log2fc The log2 fold change to determine the significant differential peaks.
#' @param padj The adjusted P values to determine the significant differential peaks.
#' @param rep_N The replicates number.
#'
#' @return
#' @export
#'
#' @examples   diff_peak <- getDiffPeak(data, condition=c("control","exp1"), rep_N = 2, control="control", experment="exp1")
#'
getDiffPeak <- function(norm_data, condition, rep_N, control, experment, log2fc=1, padj=0.05){
  coldata <- data.frame(condition = factor(rep(condition, each = rep_N), levels = condition))
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = round(norm_data), colData = coldata, design= ~condition)
  dds1 <- DESeq2::DESeq(dds, fitType = 'mean', minReplicatesForReplace = 7, parallel = TRUE)
  res <-  DESeq2::results(dds1, contrast = c('condition', experment, control))
  res1 <- data.frame(res, stringsAsFactors = FALSE, check.names = FALSE)
  res1[which(res1$log2FoldChange >= log2fc & res1$padj < padj),'sig'] <- 'Up'
  res1[which(res1$log2FoldChange <= -log2fc & res1$padj < padj),'sig'] <- 'Down'
  res1[which(abs(res1$log2FoldChange) <= log2fc | res1$padj >= padj),'sig'] <- 'None'
  res1 <- res1[complete.cases(res1),]
  return(res1)
}


#' A function to plot the differential result.
#'
#' @param diff_data The differential data obtained by **getDiffPeak**.
#' @param save_path The path to save the figure.
#' @param file_prefix The save file prefix.
#' @param figure_height The PDF figure height. Default: 8.12
#' @param figure_width The PDF figure width. Default: 8.12
#'
#' @return
#' @export
#'
#' @examples  plotVolcano(diff_peak)
plotVolcano <- function(diff_data, save_path=NA, file_prefix, figure_height=8.12,figure_width=8.12){
  p <- ggplot2::ggplot(data = diff_data, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
    ggplot2::geom_point(size = 1) + ggplot2::scale_color_manual(values = c('#F02720', 'gray', '#2C69B0'), limits = c('Up', 'None', 'Down')) +
    ggplot2::labs(x = 'log2 Fold Change', y = '-log10 adjust p-value', color = '') +
    ggplot2::theme(plot.title = element_text(hjust = 0.5, size = 14), panel.grid = element_blank(),
          panel.background = element_rect(color = 'black', fill = 'transparent'),
          legend.key = element_rect(fill = 'transparent')) +
    ggplot2::geom_vline(xintercept = c(-1, 1), lty = 3, color = 'black') +
    ggplot2::geom_hline(yintercept = -log10(0.05), lty = 3, color = 'black') +
    ggplot2::xlim(-round(max(abs(diff_data$log2FoldChange))),round(max(abs(diff_data$log2FoldChange)))) + ggplot2::ylim(0, round(max(-log10(diff_data$padj))))
  return(p)
  if(!is.na(save_plot)){
    pdf(sprintf("%s/%s_Diff_Peaks_Volcano.pdf",save_path,file_prefix),height=figure_height, width=figure_width)
    print(p)
    dev.off()
  }
}


#' A function to call the differential peak target genes.
#'
#' @param diff_data The differential data obtained by **getDiffPeak**.
#' @param save_path The path to save the result.
#' @param file_prefix The save file prefix.
#' @return
#' @export
#'
#' @examples   getDiffTargetGenes(diff_data = diff_peak)
#'
getDiffTargetGenes <- function(diff_data, save_path=NA, file_prefix){
  checkGeAnno()
  diff_peak <- diff_data
  diff_peak$chrom <- sapply(strsplit(rownames(diff_peak),":"), `[`, 1)
  diff_peak$start1 <- sapply(strsplit(rownames(diff_peak),":"), `[`, 2)
  diff_peak <- diff_peak %>% tidyr::separate(start1, c("start", "end"), "-")
  diff_peak$peak_center1 <- as.integer((as.integer(diff_peak$start) + as.integer(diff_peak$end))/2)
  diff_peak$peak_center2 <- diff_peak$peak_center1 + 1
  peak_summit <- diff_peak[,c("chrom","peak_center1","peak_center2")]
  colnames(peak_summit) <- c("chrom","start","end")
  tss <- CATAnno$tss
  target <- valr::bed_closest(peak_summit,tss)
  target$dist_strand <- ifelse(target$strand.y == "+", target$.dist, -(target$.dist))
  res <- target[,c(1,2,3,4,7,8,11)]
  final <- merge(diff_peak, res, by.x = c("chrom","peak_center1","peak_center2"),by.y=c("chrom","start.x","end.x"),all.x=T)
  final_res <- final[,c(1,11,12,4,5,6,7,8,9,10,2,13,14,15,16)]
  colnames(final_res)[11] <- "peak_center"
  colnames(final_res)[12] <- "TSS"
  colnames(final_res)[13] <- "gene"
  colnames(final_res)[14] <- "strand"
  colnames(final_res)[15] <- "distance"
  return(final_res)
  if(!is.na(save_path)){
    write.table(final_res, sprintf("%s/%s_Diff_Peak_Target.tsv",save_path,file_prefix),sep='\t',quote=F,col.names=T,row.names=T)
  }
}


#' A fucntion to plot MA plot
#'
#' @param target The target genes obtained by **getDiffTargetGenes**
#'
#' @return
#' @export
#'
#' @examples plotMA(target)
plotMA <- function(target){
  p <- ggplot2::ggplot(data = target, aes(x = log10(abs(distance)), y = log2FoldChange, color = -log10(padj))) +
    ggplot2::geom_point(size = 1) +  ggplot2::labs(x = 'log10(Distance to TSS)', y = 'log2(Fold change)', color = '')+
    ggplot2::theme(plot.title = element_text(hjust = 0.5, size = 14), panel.grid = element_blank(),
                   panel.background = element_rect(color = 'black', fill = 'transparent'),
                   legend.key = element_rect(fill = 'transparent')) + ggplot2::scale_color_distiller(palette = "Spectral",name="-log10(FDR)")+
    ggplot2::geom_hline(yintercept = 0, lty = 5, color = '#A6A6A6')
  return(p)
}



#' A function to compare the chromatin accessibility across thress samples.
#'
#' @param group_name The group names.
#' @param quant_df The normalized quantification matrix obtained by **quantification**.
#' @param return_matrix Whether to return a matrix instead of a plot. Default is FALSE.
#' @param point_size The point size. Default is 0.1.
#' @param color The color list.
#'
#' @return
#' @export
#'
#' @examples   getTriads(group_name=c("A","B","C"), quant_df = quant_df)
getTriads <- function(group_name, quant_df, return_matrix = FALSE, point_size=0.1, color = NA){
  triadpos <- data.frame(A=c(1,0,0,0,0.5,0.5,0.33), B=c(0,1,0,0.5,0,0.5,0.33),C=c(0,0,1,0.5,0.5,0,0.33))
  rownames(triadpos) <- c(sprintf("%s dominant",group_name[1]),
                          sprintf("%s dominant",group_name[2]),
                          sprintf("%s dominant",group_name[3]),
                          sprintf("%s suppressed", group_name[1]),
                          sprintf("%s suppressed",group_name[2]),
                          sprintf("%s suppressed",group_name[3]),"Balanced")

  quant_df1 <- quant_df[rowSums(quant_df) > 1, ]
  quant_df2 <- quant_df1 / rowSums(quant_df1)
  quant_df2$group <- apply(quant_df2, 1, function(x){
    rownames(triadpos)[which.min(rdist::cdist(t(x), triadpos))]
  })
  if(!is.na(color)){
    triadcol <- setNames(color[1:7], rownames(triadpos))
  }else{
    triadcol <- setNames(paletteer::paletteer_d("ggthemes::Classic_10_Medium")[1:7], rownames(triadpos))
  }
  triadcol['Balanced'] <- 'lightgrey'
  options(repr.plot.width = 8, repr.plot.height = 8, repr.plot.res = 100)
  p <- ggtern::ggtern(data = quant_df2, aes_string(x=group_name[1],y=group_name[2],z=group_name[3], color = "group")) +
    ggplot2::geom_point(size=point_size) + ggplot2::scale_color_manual(values = triadcol) + ggtern::theme_rgbw() +
    ggplot2::theme(plot.title = element_text(hjust = 0.5,size = 16)) + ggtern::geom_mask()+
    ggplot2::guides(colour = guide_legend(override.aes = list(size=6))) +
    ggplot2::theme(legend.title=element_blank())+
    ggplot2::theme(panel.background = element_blank(), legend.key = element_rect(fill="transparent"))+
    ggtern::theme_custom(col.T = "#dc6c50", col.L = "#3f7fa7", col.R = "#62b19a")
  if (return_matrix){
    mat <- merge(quant_df1, quant_df2, by = 0, all = T)[,c(1,2,3,4,8)]
    rownames(mat) <- mat$Row.names
    mat <- mat[,c(2:5)]
    colnames(mat)[1:3] <- colnames(quant_df1)
    return(mat)
  }else{
    return(p)
  }
}


#' A function to plot the triads result.
#'
#' @param triads_data The triads result obtained by **getTriads**.
#'
#' @return
#' @export
#'
#' @examples  plotTriads(triads_data = res)
plotTriads <- function(triads_data){
  df <- reshape2::melt(triads_data)
  df$value <- log2(df$value + 1)
  p <- ggpubr::ggboxplot(df, x = "group", y = "value", fill = "variable", palette = c("#00AFBB", "#E7B800", "#FC4E07"))+
       ggplot2::theme(axis.text.x = element_text(angle = 30, vjust = 0.5, hjust=1)) + ggplot2::xlab("") + ggplot2::ylab("log2 (CPM + 1)")
  return(p)
}


#' A function to get the target genes of triads peak.
#'
#' @param triads_data The triads result obtained by **getTriads**.
#' @param save_path The path to save the result.
#' @param file_prefix The prefix of the file.
#'
#' @return
#' @export
#'
#' @examples   getTriadsTargetGenes(triads_data = res)
getTriadsTargetGenes <- function(triads_data, save_path=NA, file_prefix=NA){
  checkGeAnno()
  diff_peak <- triads_data
  diff_peak$chrom <- sapply(strsplit(rownames(diff_peak),":"), `[`, 1)
  diff_peak$start1 <- sapply(strsplit(rownames(diff_peak),":"), `[`, 2)
  diff_peak <- diff_peak %>% tidyr::separate(start1, c("start", "end"), "-")
  diff_peak$peak_center1 <- as.integer((as.integer(diff_peak$start) + as.integer(diff_peak$end))/2)
  diff_peak$peak_center2 <- diff_peak$peak_center1 + 1
  peak_summit <- diff_peak[,c("chrom","peak_center1","peak_center2")]
  colnames(peak_summit) <- c("chrom","start","end")
  tss <- CATAnno$tss
  target <- valr::bed_closest(peak_summit, tss)
  target$dist_strand <- ifelse(target$strand.y == "+", target$.dist, -(target$.dist))
  res <- target[,c(1,2,3,4,6,8,11)]
  final <- merge(diff_peak, res, by.x = c("chrom","peak_center1","peak_center2"),by.y=c("chrom","start.x","end.x"),all.x=T)
  final_res <- final[,-c(3)]
  colnames(final_res)[2] <- "Peak_center"
  colnames(final_res)[7:8] <- c("Peak_start","Peak_end")
  colnames(final_res)[9] <- "TSS"
  colnames(final_res)[10] <- "gene"
  colnames(final_res)[11] <- "strand"
  colnames(final_res)[12] <- "distance"
  if(!is.na(save_path)){
    write.table(final_res, sprintf("%s/%s_Diff_Peak_Target.tsv",save_path,file_prefix),sep='\t',quote=F,col.names=T,row.names=T)
  }
  return(final_res)
}

