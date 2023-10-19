
##############################################################################
#                 ATAC-to-RNA / Peak-to-Gene links.                          #
##############################################################################

#' Calculate the ATAC-to-RNA / Peak-to-Gene links.
#' @param atac_matrix  The normalized ATAC-seq data frame obtained by **quantification**.
#' @param rna_matrix  The normalized RNA-seq data frame. Note: The RNA and ATAC matrix have the same column names and orders.
#' @param peak_annotation  Peak annotation obtained by **plotPDI**. The file suffix: _Proximal_Distal_Intragenic_Annotation.tsv
#' @param max_distance  The max distance (bp) around the gene to calculate the correlations.
#'                       If you set 20000, it will find the peak-to-gene within 20kb upstream and downstream of the gene. Default: 20000.
#' @param N_permutation Simulate how many times pseudo-data to measure the P value of each peak-to-gene. Default: 10000.
#' @param save_path The path to save the result.
#'
#' @return
#' @export
#'
#' @examples   getPeak2Gene(atac_matrix="~/zhutao/CAT/ATAC_norm_quant.tsv", rna_matrix="~/zhutao/CAT/RNA_norm_quant.tsv",
#'                          tss_file="~/nip_tss.bed", peak_annotation="~/zhutao/CAT/peak_anno.tsv", max_distance=20000, N_permutation=10000)
getPeak2Gene <- function(atac_matrix, rna_matrix, peak_annotation,
                          max_distance=20000, N_permutation=10000, save_path=NA){
  checkGeAnno()
  atac_paired_norm <- read.table(atac_matrix, head=T, row.names = 1)
  rna_paired <- read.table(rna_matrix,head=T, row.names = 1)
  logfile("Remove the gene with all expression value is 0.")
  rna_paired <- rna_paired[rowSums(rna_paired) > 0,]
  tss <- CATAnno$tss %>% as.data.frame()
  rownames(tss) <- tss$name
  tss_df <- tss[,c(1,2,3,4,6)]
  colnames(tss_df) <- c("chr","start","end","gene","strand")
  gr.tss <- GenomicRanges::makeGRangesFromDataFrame(tss_df,keep.extra.columns=TRUE)
  peak_bed <- read.table(peak_annotation,head=T)
  peak_bed$Peak <- sprintf("%s:%s-%s",peak_bed$Chromosome, peak_bed$Start, peak_bed$End)
  rownames(peak_bed) <- NULL
  peak_gr <- GenomicRanges::makeGRangesFromDataFrame(peak_bed, keep.extra.columns=TRUE)

  ## Make the pseudo data for permutation
  logfile("Make the pseudo data for permutation...")
  test_random <- atac_paired_norm[sample(nrow(atac_paired_norm), N_permutation), ]


  ## Calculate the correlation coefficient and p-value for a given gene.
  calculate_pvalue <- function(row_data, df, gene, test_random=test_random){
    corr_random <- apply(test_random, 1, function(row){
      corr <- cor(as.numeric(row_data), row)
    })
    correlations <- apply(df, 1, function(row){
      corr <- cor(as.numeric(row_data), row)
      if (is.numeric(corr)){
      p <- BSDA::z.test(corr_random, mu = corr, sigma.x = 15)$p.value
      } else{
        corr = 0
        p <- BSDA::z.test(corr_random, mu = corr, sigma.x = 15)$p.value
      }
      return(list(corr = corr, p = p))
    })
    result <- data.frame(Peak = rownames(df), Gene = gene)
    result$correlations <- sapply(correlations, function(x) x$corr)
    result$p.value <- sapply(correlations, function(x) x$p)
    return(result)
  }

  ## get the gene2peak links function
  getGenePeaks <- function(gene, max_distance){
    gr.bm_gs <- unique(gr.tss[gr.tss$gene %in% gene])
    GenomicRanges::start(gr.bm_gs) <- GenomicRanges::start(gr.bm_gs) - max_distance
    GenomicRanges::end(gr.bm_gs) <- GenomicRanges::end(gr.bm_gs) + max_distance
    ol <- GenomicRanges::findOverlaps(peak_gr, gr.bm_gs)
    peak_in_window <- peak_gr[unique(IRanges::from(ol))]
    if(length(peak_in_window)!=0){
      df1 <- data.frame(chr=GenomicRanges::seqnames(peak_in_window), start=GenomicRanges::start(peak_in_window), end=GenomicRanges::end(peak_in_window))
      df1$peak <- sprintf("%s:%s-%s", df1$chr, df1$start, df1$end)
      ATAC <- atac_paired_norm[df1$peak,]
      RNA <- rna_paired[gene,]
      res <- calculate_pvalue(RNA, ATAC, gene, test_random=test_random)
      return(res)
    }else{
      return(NA)
    }
  }
  ##  calculate the all gene peak2gene links
  gene_list <- rownames(rna_paired)
  pb <- progress::progress_bar$new(total = length(gene_list))
  result_list <- lapply(c(1:length(gene_list)), function(x) {
    result <- getGenePeaks(gene_list[x],max_distance=max_distance)
    pb$tick()
    return(result)
  })
  na.omit.list <- function(y) { return(y[!sapply(y, function(x) all(is.na(x)))]) }
  res_list <- na.omit.list(result_list)
  all_res <- dplyr::bind_rows(res_list)
  final <- merge(all_res,peak_bed,by="Peak",all.x=T)
  final_res <- merge(final,tss_df, by.x="Gene", by.y="gene", all.x=T)
  final_res$Summit2TSS <- final_res$summit - final_res$start
  final_res$orientation <- ifelse((final_res$strand=="+" & final_res$Summit2TSS>=0)|(final_res$strand=="-" & final_res$Summit2TSS<0), "Downstream", "Upstream")
  final_res <- final_res[,c("Peak","Gene","correlations","p.value","Type","summit","start","Summit2TSS","strand","orientation")]
  colnames(final_res)[6:7] <- c("PeakSummit","TSS")
  final_res <- tibble::add_column(final_res, FDR = p.adjust(final_res$p.value, method = "fdr"), .after = 4)
  if(!is.na(save_path)){
    saveRDS(final_res, sprintf("%s/Peak2Gene_All_Links.rds",save_path))
  }
  return(final_res)
}


##############################################################################
#                 Plot all Peak2Gene links with heatmap.                     #
##############################################################################
#' Plot all the Peak2Gene links with heatmap.
#' @param p2g_res The peak2gene result obtained by **getPeak2Gene**
#'
#' @param cor_cutoff  Pearson's correlation coefficient threshold for determining significance Peak2Gene links.
#' @param atac_matrix The normalized ATAC-seq data frame obtained by **quantification**. Same as **getPeak2Gene** parameter.
#' @param rna_matrix The normalized RNA-seq data frame. Same as **getPeak2Gene** parameter.
#' @param cluster_N Clustering rows into N modules. Default: 4.
#' @param palATAC The color of ATAC heatmap.
#' @param palRNA The color of RNA heatmap.
#' @param save_path The path to save result.
#' @param fig_width The PDF width. Default: 12.
#' @param fig_height The PDF height. Default: 12.
#' @param raster Wheather to raster the plot. Default: FALSE.
#'
#' @return
#' @export
#'
#' @examples   plotP2GHeatmap("./P2G_all_links.rds",cor_cutoff=0.4,atac_matrix="./ATAC_norm_quant.tsv",rna_matrix="./RNA_norm_quant.tsv",save_path="./")
#'
plotP2GHeatmap <- function(p2g_res, cor_cutoff, atac_matrix, rna_matrix, cluster_N =4, raster = F,
                           palATAC = NA, palRNA= NA, save_path=NA, fig_width=12,fig_height=12){
  checkGeAnno()
  blueYellow <- c("1"="#352A86","2"="#343DAE","3"="#0262E0","4"="#1389D2","5"="#2DB7A3","6"="#A5BE6A","7"="#F8BA43","8"="#F6DA23","9"="#F8FA0D")
  solarExtra <- c("5"='#3361A5', "7"='#248AF3', "1"='#14B3FF', "8"='#88CEEF', "9"='#C1D5DC', "4"='#EAD397', "3"='#FDB31A',"2"= '#E42A2A', "6"='#A31D1D')
  p2g <- readRDS(p2g_res)
  p2g <- p2g[p2g$correlations >= cor_cutoff, ]
  logfile(sprintf("Depending on the cutoff. The Peak2Gene links number: %s",dim(p2g)[1]))
  #logfile("This step may take a while, please wait...")
  rna <- read.table(rna_matrix, header = T, row.names = 1)
  p2g$idx_atac <- sprintf("ATAC_%s",1:nrow(p2g))
  p2g$idx_rna <- sprintf("RNA_%s",1:nrow(p2g))
  p2g_atac <- p2g[,c("Peak","idx_atac")]
  rownames(p2g_atac) <- p2g_atac$idx_atac
  atac <- read.table(atac_matrix, row.names = 1, header = T)
  atac_mat <- merge(p2g_atac, atac, by.x="Peak",by.y=0,all=F)
  rownames(atac_mat) <- atac_mat$idx_atac
  atac_mat <- atac_mat[,-c(1,2)]
  atac_scale <- rowZscores(as.matrix(atac_mat),limit=T)

  df1 <- rowZscores(as.matrix(atac_scale), limit = TRUE) %>% as.data.frame()
  logfile("Calculating ATAC matrix.")
  set.seed(2023)
  row_dend <- hclust(dist(df1))
  #mat <- df1[row_dend$order, ]
  group <- data.frame(C=cutree(row_dend, k = as.integer(cluster_N)))
  group$Cluster <- paste0("Cluster",group$C)
  group <- group[order(group$Cluster),]
  group$Cluster <- factor(group$Cluster, levels = sort(unique(group$Cluster)))
  mat <- df1[rownames(group),]
  group <- subset(group,select=-(C))
  annotation <- mat[match(rownames(group),rownames(as.data.frame(mat))),]
  mycol <- c(paletteer::paletteer_d("ggthemes::Tableau_10"), paletteer::paletteer_d("ggthemes::Tableau_20"))
  gcols <- setNames(as.character(mycol[1:cluster_N]),unique(group$Cluster))
  gcol <- list(Cluster=gcols)

  N <- cluster_N - 1
  if(!is.na(palATAC)){
    p3 <- ComplexHeatmap::pheatmap(as.matrix(annotation), show_rownames=F,cluster_row=F,cluster_col=F,border_color=NA,use_raster=raster,
                                   annotation_colors =gcol, color = colorRampPalette(palATAC)(256),
                                   annotation_row = group, gaps_row = cumsum(as.numeric(table(group$Cluster)))[1:N],
                                   main="ATAC-seq",name="ATAC Z-score")
  } else{
    p3 <- ComplexHeatmap::pheatmap(as.matrix(annotation), show_rownames=F,cluster_row=F,cluster_col=F,border_color=NA,use_raster=raster,
                                   annotation_colors =gcol, color = colorRampPalette(blueYellow)(256),
                                   annotation_row = group, gaps_row = cumsum(as.numeric(table(group$Cluster)))[1:N],
                                   main="ATAC-seq",name="ATAC Z-score")
  }
  # Set the RNA matrix gene order
  rownames(p2g) <- p2g$idx_atac
  logfile("Calculating RNA matrix.")
  p2gnew <- p2g[rownames(group),]
  p2g_rna <- p2gnew[,c("Gene","idx_rna")]
  rna_order <- merge(p2g_rna, rna, by.x="Gene",by.y=0,all.x =T)
  rownames(rna_order) <- rna_order$idx_rna
  rna_order <- rna_order[p2gnew$idx_rna,-c(1,2)]
  rna_order <- rna_order[,colnames(annotation)]
  rna_scale <- rowZscores(as.matrix(rna_order),limit=T)
  if(!is.na(palRNA)){
    p4 <- ComplexHeatmap::pheatmap(as.matrix(rna_scale), show_rownames=F,cluster_row=F,cluster_col=F,border_color=NA,use_raster=raster,
                                   annotation_colors =gcol, color = colorRampPalette(palRNA)(256),
                                   annotation_row = group, gaps_row = cumsum(as.numeric(table(group$Cluster)))[1:N],
                                   main="RNA-seq", name="RNA Z-score")
  }else{
    p4 <- ComplexHeatmap::pheatmap(as.matrix(rna_scale), show_rownames=F,cluster_row=F,cluster_col=F,border_color=NA,use_raster=raster,
                                   annotation_colors =gcol, color = colorRampPalette(solarExtra)(256),
                                   annotation_row = group, gaps_row = cumsum(as.numeric(table(group$Cluster)))[1:N],
                                   main="RNA-seq", name="RNA Z-score")
  }
  pcor <- ComplexHeatmap::Heatmap(p2gnew$correlations, cluster_rows=F, name = "Correlations", col = paletteer::paletteer_d("RColorBrewer::Oranges")[1:7])
  ppvalue <- ComplexHeatmap::Heatmap(-log10(p2gnew$FDR), cluster_rows=F, name = "-log10(FDR)", col = paletteer::paletteer_d("RColorBrewer::Greens")[1:7])
  gene_type_col = setNames(as.character(paletteer_d("ggthemes::Classic_10")[1:length(unique(p2gnew$Type))]),unique(p2gnew$Type))
  ptype <- ComplexHeatmap::Heatmap(p2gnew$Type, cluster_rows=F, name = "Type", col = gene_type_col)
  pdis <- ComplexHeatmap::Heatmap(log10(abs(p2gnew$Summit2TSS)+1), cluster_rows=F, name = "log10(Summit2TSS + 1)", col = paletteer::paletteer_d("grDevices::blues9")[1:7])
  if(!is.na(save_path)){
    write.table(rna_scale, sprintf("%s/Peak2Gene_Links_RNA_matrix.tsv",save_path),sep='\t',quote=F)
    write.table(annotation, sprintf("%s/Peak2Gene_Links_ATAC_matrix.tsv",save_path),sep='\t',quote=F)
    write.table(p2g, sprintf("%s/Peak2Gene_Links_matrix.tsv",save_path),sep='\t',quote=F)
    pdf(sprintf("%s/ATAC_RNA_All_Links_Heatmap.pdf",save_path),width=fig_width,height=fig_height)
    ComplexHeatmap::draw(p3 + p4 + pdis + ptype + pcor + ppvalue, column_title = sprintf("%s Peak-to-Gene linkages",dim(p2g)[1]))
    dev.off()
  }
  return(ComplexHeatmap::draw(p3 + p4 + pdis + ptype + pcor + ppvalue, column_title = sprintf("%s Peak-to-Gene linkages",dim(p2g)[1])))
}



##############################################################################
#                 Plot one gene / one peak scatter plot.                     #
##############################################################################
#' Plot one gene / one peak scatter plot.
#' @param peak The peak coordinate you want to plot.
#' @param gene The gene ID you want to plot.
#' @param atac_matrix The normalized ATAC-seq data frame obtained by **quantification**. Same as **getPeak2Gene** parameter.
#' @param rna_matrix The normalized RNA-seq dataframe. Same as **getPeak2Gene** parameter.
#' @param color The color to plot.
#' @param legend Wheather to plot legend. default: TRUE.
#'
#' @return
#' @export
#'
#' @examples   plotPeakGene("1:39081-39341","gene","ATAC_norm_quant.tsv","RNA_norm_quant.tsv",legend=F)
plotPeakGene <- function(peak, gene, atac_matrix, rna_matrix, color=NA, legend = TRUE){
  if(!is.na(color)){
    col <- color
  }else{
    col <- c(as.character(paletteer_d("ggthemes::Hue_Circle")),as.character(paletteer_d("ggthemes::Tableau_20")),
             as.character(c("1"="#D51F26","2"="#272E6A","3"="#208A42","4"="#89288F","5"="#F47D2B", "6"="#FEE500",
                            "7"="#8A9FD1","8"="#C06CAB","19"="#E6C2DC","10"="#90D5E4", "11"="#89C75F","12"="#F37B7D",
                            "13"="#9983BD","14"="#D24B27","15"="#3BBCA8", "16"="#6E4B9E","17"="#0C727C", "18"="#7E1416",
                            "9"="#D8A767","20"="#3D3D3D")))
    atac_paired_norm <- read.table(atac_matrix, row.names=1, header=T)
    #rownames(atac_paired_norm) <- gsub("chr","",rownames(atac_paired_norm))
    rna_paired <- read.table(rna_matrix, row.names=1, header=T)
    if(dim(atac_paired_norm)[2]!= dim(rna_paired)[2]){
      print("The ATAC matrix and the RNA have different matrix columns")
    } else{
      samples <- colnames(rna_paired)
      col <- setNames(col[1:length(colnames(rna_paired))], samples)
      df <- data.frame(ATAC=as.numeric(as.matrix(atac_paired_norm)[peak,]),RNA=as.numeric(rna_paired[gene,]))
      df$sample <- colnames(atac_paired_norm)
      options(repr.plot.width = 5, repr.plot.height = 5, repr.plot.res = 100)
      pp <- ggplot2::ggplot(df, aes(x=ATAC, y=RNA)) +
            ggplot2::geom_point(aes(colour = sample))+
            ggplot2::geom_smooth(method=lm,color = "#F28E2B") + ggplot2::scale_color_manual(values = col)+
            ggpubr::theme_pubr(legend = c("right")) + ggplot2::theme(panel.border = element_rect(colour = "black", fill=NA))+
            ggplot2::xlab('ATAC-seq  log2(CPM+1)') + ggplot2::ylab("RNA-seq  log2(TPM+1)") + ggplot2::ggtitle(sprintf("%s\n%s", gene, peak))
      if(!legend){
        pp <- pp + ggplot2::theme(legend.position = "none")
      }
      return(pp)
    }
  }
}


##############################################################################
#         Plot specific gene with tracks and peak2gene linkages.             #
##############################################################################

#' Plot specific gene with tracks and peak2gene linkages.
#'
#' @param samples_path The bigwig (bw) file folder.
#' @param samples_suffix The suffix of the samples. Like the ".cpm.bw"
#' @param gene_name The gene name.
#' @param left The gene left flanking to plot. Default: 10000.
#' @param right The gene right flanking to plot. Default: 10000.
#' @param peaks The merged peaks.
#' @param peak2gene The RDS file obtained by **getPeak2Gene**.
#' @param cor_cutoff The correlation cutoff of peak-to-gene links. Default: 0.4.
#' @param line_size The curve line with.  Default: 0.8.
#' @param curvature The degree of curvature. Default: 0.3.
#' @param color Color list to plot tracks.
#' @param back.colorWhether to plot backgroud color.  Default: TRUE.
#'
#' @return
#' @export
#'
#' @examples  plotP2GTracks(samples_path="./signal", samples_suffix=".cpm.bw", gene_name="ENST00000220244", peaks="./merged_peaks.bed", peak2gene="./All_Peak2Gene_links.rds")
#'
plotP2GTracks <-  function(samples_path, samples_suffix, gene_name, left = 10000, right = 10000,
                          peaks, peak2gene, cor_cutoff = 0.4,line_size = 0.8, curvature = 0.3, color = NA, back.color = TRUE){
  checkGeAnno()
  coords <- getGeneBed(gene_name = gene_name, left = left, right = right)
  chr <- coords[1]
  start <- as.integer(coords[2])
  end <- as.integer(coords[3])
  link <- readRDS(peak2gene)
  link$chr <- sapply(strsplit(link$Peak,":"), `[`, 1)
  link <- link[,c("chr","PeakSummit","TSS","correlations")]
  colnames(link) <- c("chr","start","end","correlations")
  link$correlations <- round(link$correlations,2)
  link_filter <- link[abs(link$correlations)>=cor_cutoff & link$chr==chr & link$start>= start & link$end <= end, ]
  colora <- as.character(paletteer::paletteer_d("ggthemes::Hue_Circle"))
  colorb <- as.character(paletteer::paletteer_d("ggthemes::Tableau_20"))
  colorc <- as.character(paletteer::paletteer_d("ggthemes::Tableau_10"))
  if(nrow(link_filter)==0){
    print("Your provide region without significant peak2gene links!")
  }
  else{
    plink <- linkVis(linkData = link_filter,
                     start = "start",
                     end = "end",
                     link.aescolor = "correlations",
                     link.color = c('grey80','#FF800E'),
                     line.size=line_size,
                     facet = FALSE,
                     curvature = curvature,
                     xAixs.info = FALSE)

    bwlist <- list()
    for (bw in Sys.glob(sprintf("%s/*%s",samples_path, samples_suffix))){
      bwlist[[length(bwlist) + 1]] <- bw
    }
    bwfile <- unlist(bwlist)
    if (!is.na(color)){
      if(length(color) < length(bwfile)){
        print("The color list you provided is less than the number of samples.")
      }
      else{
        mybw <- loadBigWig(bwfile, chr = chr, start = start, end = end)
        mybw$score <- round(mybw$score,1)
        ## plot
        ptrack <- trackVis(bWData = mybw,
                           chr = chr,
                           region.min = start,
                           region.max = end,
                           color = color,
                           theme = "bw",
                           yAxis.info = F,new.yaxis = T,
                           pos.ratio = c(0.02,0.8),
                           back.color = back.color)

        gtf <- rtracklayer::import(CATAnno$gtf,format = "gtf") %>% data.frame()
        trans <- trancriptVis(gtfFile = gtf,
                              Chr = chr,
                              posStart = start - 3000,
                              posEnd = end + 3000,
                              addNormalArrow = FALSE,
                              newStyleArrow = T,
                              absSpecArrowLen = T,
                              speArrowRelLen = 0.2,
                              textLabel = "gene_name",
                              textLabelSize = 3,
                              relTextDist = 0.4,
                              exonWidth = 0.3,
                              collapse = F,
                              selecType = "lt"
        )
        peakvis <- bedVis(bdFile = peaks,
                          chr = chr, region.min = start, region.max = end, fill = "#006BA4", show.legend=F)

        ptrack %>% aplot::insert_bottom(peakvis,height = 0.03)%>% aplot::insert_bottom(plink,height = 0.1) %>% aplot::insert_bottom(trans,height = 0.08)
      }
    }
    else{
      mybw <- loadBigWig(bwfile, chr = chr, start = start, end = end)
      mybw$score <- round(mybw$score,1)
      ## plot
      ptrack <- trackVis(bWData = mybw,
                         chr = chr,
                         region.min = start,
                         region.max = end,
                         color = c(colora, colorb, colorc),
                         theme = "bw",
                         yAxis.info = F,new.yaxis = T,
                         pos.ratio = c(0.02,0.8),
                         back.color = back.color)

      gtf <- rtracklayer::import(CATAnno$gtf,format = "gtf") %>% data.frame()
      trans <- trancriptVis(gtfFile = gtf,
                            Chr = chr,
                            posStart = start - 3000,
                            posEnd = end + 3000,
                            addNormalArrow = FALSE,
                            newStyleArrow = T,
                            absSpecArrowLen = T,
                            speArrowRelLen = 0.2,
                            textLabel = "gene_name",
                            textLabelSize = 3,
                            relTextDist = 0.4,
                            exonWidth = 0.3,
                            collapse = F,
                            selecType = "lt"
      )
      peakvis <- bedVis(bdFile = peaks,
                        chr = chr, region.min = start, region.max = end, fill = "#006BA4", show.legend=F)

      ptrack %>% aplot::insert_bottom(peakvis,height = 0.03)%>% aplot::insert_bottom(plink,height = 0.1) %>% aplot::insert_bottom(trans,height = 0.08)
    }
  }
}
