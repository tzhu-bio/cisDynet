
#' A function to get the footprint score.
#'
#' @return  A matrix with all samples footprint score.
#' @export
#'
#' @examples  getFootprintScore(samples=c("s1","s2","s3"...))
#'
getFootprintScore <- function(){
  checkFTAnno()
#  all_lst <- list()
#  for (sample in samples){
#    df <- read.table(sprintf("%s/%s/bindetect_results.txt", FTAnno$bindetect, sample), head=T)
#    df <- df[,c(1,6)]
#    all_lst[[sample]] <- df
#  }
#  all <- Reduce(function(x, y) merge(x, y, by = 1, all = TRUE), all_lst)
  ft <- read.table(sprintf("%s/bindetect_results.txt",FTAnno$bindetect), head = T)
  ftscore <- ft[, grepl(paste(c("output_prefix", "mean_score"),collapse="|"), colnames(ft))]
  colnames(ftscore) <- gsub("_footprints_mean_score","",colnames(ftscore))
  return(ftscore)
}

######################################################################################

#' A function to plot all samples footprint with scatter / heatmap.
#'
#' @param score_result The footprint score result obtained by **getFootprintScore**.
#' @param cluster_N The cluster number for heatmap. Default: 6.
#' @param label_motif The motif list to label in the right of the heatmap.
#' @param limit_Zscore Whether to limit zscore to -2 to 2. Default: FALSE.
#'
#' @return
#' @export
#'
#' @examples   plotFootprintScore(score_result=res, label_motif("motif1","motif2","motif3"...))
#'
plotFootprintScore <- function(score_result, cluster_N = 6, label_motif, limit_Zscore = FALSE){
  if(ncol(score_result)==3){
    s1 <- colnames(score_result)[2]
    s2 <- colnames(score_result)[3]
    score_result$group <- ifelse(score_result[,2] > score_result[,3], s1, s2)
    score_result$diff <- abs(score_result[,2] - a[,3])
    score_result$label <- ifelse(score_result$diff > quantile(score_result$diff, probs = 0.98), score_result$output_prefix, NA)
    p <- ggscatter(score_result, x = s1, y = s2, color = "group",palette = c("#00AFBB",  "#FC4E07"),size=1, label = "label",repel = TRUE) + geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red")
    return(p)
  } else {
    rownames(score_result) <- score_result$output_prefix
    df <- score_result[,c(2:ncol(score_result))]
    df <- quantile_normalization(df)
    #df <- preprocessCore::normalize.quantiles(as.matrix(df))
    df1 <- rowZscores(as.matrix(df), limit = limit_Zscore)
    df1 <- df1[complete.cases(df1), ]
    row_dend <- hclust(dist(df1))
    mat <- df1[row_dend$order, ]
    group <- data.frame(C=cutree(row_dend, k = cluster_N))
    group$Cluster <- paste0("Cluster",group$C)
    group <- group[order(group$Cluster),]
    group$Cluster <- factor(group$Cluster, levels = sort(unique(group$Cluster)))
    mat <- df1[rownames(group),]
    group <- subset(group,select=-(C))
    annotation <- mat[match(rownames(group),rownames(as.data.frame(mat))),]
    gcols <- setNames(as.character(paletteer::paletteer_d("ggthemes::Tableau_10")[1:cluster_N]),unique(group$Cluster))
    gcol <- list(Cluster=gcols)
    anno <- data.frame(motif = label_motif)
    anno$id <- sapply(strsplit(anno$motif,"_"), `[`, 1)
    ha <- ComplexHeatmap::rowAnnotation(foo = ComplexHeatmap::anno_mark(at=match(anno$motif,rownames(as.data.frame(mat))), labels=anno$id))
    p <- ComplexHeatmap::pheatmap(mat, color = rev(paletteer::paletteer_c("ggthemes::Red-Blue-White Diverging", 30)[5:25]), show_rownames = F, name="Z-score", cluster_rows = F,
                                  annotation_colors =gcol, annotation_row = group, gaps_row = cumsum(as.numeric(table(group$Cluster)))[1:(cluster_N-1)]) + ha
  }
  return(p)
}


#' A function to compare the Tn5 signal around the motif in the whole genome.
#'
#' @param samples Samples used to compare the foorprint score.
#' @param motif Specific motifs ID. Note the inclusion of the suffix (_MAxxx).
#' @param smooth_window The size of the window used to smooth the lines.
#' @param flanking_length The flanking length to calculate and visualisation.
#' @param color The color list.
#' @param save_path  The path to save the plot.
#' @param file_prefix The prefix of the PDF file.
#'
#' @return
#' @export
#'
#' @examples   compareFootprints(samples = c("CD8pos_T","Monocytes"),motif = "Arntl_MA0603.1",bindetect_path = "F:/CAT/example/BINDetect/",signal_path = "F:/CAT/example/ATACcorrect/")
#'
compareFootprints <- function(samples, motif, smooth_window=5, flanking_length=200, N_cores = 1, color = NA, save_path=NA, file_prefix=NA){
  checkFTAnno()
  cal_footprint <- function(gene_bed, coverage_bed, group_name) {
    logfile("Reading footprint...")
    #gene <- genomation::readBed(gene_bed)
    gene <- valr::read_bed(gene_bed, n_fields = 6)
    gene <- GenomicRanges::makeGRangesFromDataFrame(gene, keep.extra.columns=T)
    logfile("Reading bigwig...")
    coverage <- loadFullBigWig(coverage_bed)
    colnames(coverage) <- c("chr","start","end","coverage","sample")
    cov <- GenomicRanges::makeGRangesFromDataFrame(coverage, keep.extra.columns=T)
    logfile("Normalize to matrix...")
    mat1 = EnrichedHeatmap::normalizeToMatrix(cov, gene, value_column = "coverage",
                             extend = flanking_length, mean_mode = "w0", w = smooth_window,smooth = TRUE,background = NA,trim=c(0.01, 0.01))

    res <- data.frame(mat1)
    colnames(res) <- 1:ncol(res)
    aa <- reshape2::melt(res)
    quant_smo <- aa %>% dplyr::group_by(variable) %>% dplyr::summarise(
              mean_value = mean(value, na.rm = TRUE),
              lower_ci = mean(value, na.rm = TRUE) - qt(0.975, df = sum(!is.na(value))) * sd(value, na.rm = TRUE) / sqrt(sum(!is.na(value))),
              upper_ci = mean(value, na.rm = TRUE) + qt(0.975, df = sum(!is.na(value))) * sd(value, na.rm = TRUE) / sqrt(sum(!is.na(value)))
  )
    quant_smo$idx <- 1:nrow(quant_smo)
    quant_smo$group <- group_name
    #quant_smo <- data.frame(x = c(1:length(colnames(mat1))), density=colMeans(mat1,na.rm=T),group=group_name)
    return(quant_smo)
  }
  lista <- list()
  for (a in samples){
    lista[[a]] <- sprintf("%s/%s/beds/%s_all.bed",FTAnno$bindetect, motif, motif)
  }
  listb <- list()
  for (b in samples){
    listb[[b]] <- sprintf("%s/%s_q30_corrected.bw",FTAnno$signal, b)
  }
  c <- as.list(samples)
  res <- list()

  res <- parallel::mcmapply(function(gene_bed, coverage_bed, group_name){
    return(cal_footprint(gene_bed, coverage_bed, group_name))
  }, gene_bed=lista, coverage_bed=listb, group_name=c, mc.cores=N_cores, SIMPLIFY = FALSE)

  # res <- parallel::mclapply(function(gene_bed, coverage_bed, group_name){
  #   return(cal_footprint(gene_bed, coverage_bed, group_name))
  # }, gene_bed=lista, coverage_bed=listb, group_name=c, mc.cores=length(samples),SIMPLIFY = FALSE)

  tt <- data.frame(do.call(rbind, res))
  motif_len <- read.table(lista[[1]])
  motif_len$len <- motif_len$V3 - motif_len$V2
  motif_length <- as.integer(unique(motif_len$len))
  linea <- flanking_length / smooth_window
  lineb <- linea + (motif_length / smooth_window)+2
  scale <- (flanking_length + motif_length)/ smooth_window
  options(repr.plot.width = 7, repr.plot.height = 6, repr.plot.res = 100)
  if(is.na(color)){
    p3 <- ggplot2::ggplot(tt, aes(x=idx, y=mean_value, color=group, fill = group)) +
          ggplot2::geom_line(size=0.5) +
          ggplot2::scale_color_manual(values = paletteer::paletteer_d("ggthemes::Classic_10"))+
          ggplot2::scale_linetype_manual(values=c(rep("solid",5)))+
          ggpubr::theme_pubr(base_size = 12,base_family = "", border = TRUE, margin = TRUE, legend = c("right"), x.text.angle = 0)+
          ggplot2::geom_vline(xintercept = c(linea, lineb), linetype="dashed", color = "#5C6068", size=0.5)+
          ggplot2::coord_cartesian(expand = T, ylim = c(quantile(tt$mean_value, 0.001), 1.2*quantile(tt$mean_value, 0.999)))+
          ggplot2::scale_x_continuous(breaks=c(0,0.5 * scale, scale,1.5 * scale, 2 * scale),
                         labels=c(sprintf("-%s",flanking_length), sprintf("-%s",flanking_length/2), "0",sprintf("%s",flanking_length/2),sprintf("%s",flanking_length)))+
          ggplot2::labs(x = "Distance to motif center (bp)", y = 'Tn5 bias-corrected normalized insertions') + ggplot2::ggtitle(motif)+
          ggplot2::geom_ribbon(aes(ymin=lower_ci, ymax=upper_ci, fill=group,color=NULL), alpha = 0.2)+ggplot2::theme(plot.title = element_text(hjust = 0.5))+
          ggplot2::scale_fill_manual(values = paletteer::paletteer_d("ggthemes::Classic_10"))
  }
  else{
    p3 <- ggplot2::ggplot(tt, aes(x=x, y=density, color=group)) +
          ggplot2::geom_line(size=0.5) +
          ggplot2::scale_color_manual(values = color)+scale_linetype_manual(values=c(rep("solid",5)))+
          ggpubr::theme_pubr(base_size = 12,base_family = "", border = TRUE, margin = TRUE, legend = c("right"), x.text.angle = 0)+
          ggplot2::geom_vline(xintercept = c(linea, lineb), linetype="dashed", color = "#5C6068", size=0.5)+
          ggplot2::coord_cartesian(expand = T, ylim = c(quantile(tt$density, 0.001), 1.2*quantile(tt$density, 0.999)))+
          ggplot2::scale_x_continuous(breaks=c(0,0.5 * scale, scale,1.5 * scale, 2 * scale),
                         labels=c(sprintf("-%s",flanking_length), sprintf("-%s",flanking_length/2), "0",sprintf("%s",flanking_length/2),sprintf("%s",flanking_length)))+
          ggplot2::labs(x = "Distance to motif center (bp)", y = 'Tn5 bias-corrected normalized insertions') + ggplot2::ggtitle(motif)+
          ggplot2::geom_ribbon(aes(ymin=lower_ci, ymax=upper_ci, fill=group,color=NULL), alpha = 0.2)+ ggplot2::theme(plot.title = element_text(hjust = 0.5))+
          ggplot2::scale_fill_manual(values = colors)
  }
  if(!is.na(save_path)){
    pdf(sprintf("%s/%s_%s_Footprint.pdf",save_path, file_prefix, motif),width=7,height=6)
    print(p3)
    dev.off()
  }
  return(p3)
}



#' Plot the differential footprint with Volcano, Bar or Lollipop.
#'
#' @param samples  The two sample names to compare.
#' @param cutoff  The cutoff of fold change. Default is NA.
#' @param plot_type Plot type. Three types are available: "Volcano","Bar","Lollipop".
#' @param save_path The path to save PDF file.
#' @param file_prefix The prefix of PDF file.
#' @param figure_width PDF width. Default is 8.21.
#' @param figure_height PDF height. Default is 8.21.
#' @param scale
#'
#' @return
#' @export
#'
#' @examples   plotDiffFootprint(samples=c("NIP_YP1","NIP_YP2"), bindetect_path="/public/workspace/zhutao/encode/CAT",plot_type="Bar",cutoff=0.1)
#'
plotDiffFootprint <- function(samples,cutoff=NA, plot_type=c("Volcano","Bar","Lollipop"),scale = 5, save_path=NA,
                              file_prefix=NA,figure_width=8.21, figure_height=8.21){
  checkFTAnno()
  plot_type <- match.arg(plot_type)
  bindetect <- read.table(sprintf("%s/bindetect_results.txt",FTAnno$bindetect),head=T)
  check_name <- dim(bindetect[,grepl(sprintf("%s_footprints_%s_footprints",samples[1], samples[2]),colnames(bindetect))])[2]
  if (check_name==0){
    bindetect <- bindetect[,c("name",sprintf("%s_footprints_%s_footprints_change",samples[2], samples[1]),
                              sprintf("%s_footprints_%s_footprints_pvalue",samples[2], samples[1]))]
    colnames(bindetect) <- c("motif","change","p")
    if(is.na(cutoff)){
      cutoff <- as.numeric(quantile(bindetect$change, 0.99))
      bindetect[which(bindetect$change >= cutoff & bindetect$p < 0.05),'sig'] <- sprintf('%s_Up',samples[2])
      bindetect[which(bindetect$change <= - cutoff & bindetect$p < 0.05),'sig'] <- sprintf('%s_Up',samples[1])
      bindetect[which(abs(bindetect$change) <= cutoff | bindetect$p >= 0.05),'sig'] <- 'n.s.'
      res <- bindetect[bindetect$sig!="n.s.",]
      rownames(res) <- NULL
      s1 <- sprintf('%s_Up',samples[2])
      s2 <- sprintf('%s_Up',samples[1])
      #return(bindetect)
    }
    else{
      cutoff <- cutoff
      bindetect[which(bindetect$change >= cutoff & bindetect$p < 0.05),'sig'] <- sprintf('%s_Up',samples[2])
      bindetect[which(bindetect$change <= - cutoff & bindetect$p < 0.05),'sig'] <- sprintf('%s_Up',samples[1])
      bindetect[which(abs(bindetect$change) <= cutoff | bindetect$p >= 0.05),'sig'] <- 'n.s.'
      res <- bindetect[bindetect$sig!="n.s.",]
      rownames(res) <- NULL
      s1 <- sprintf('%s_Up',samples[2])
      s2 <- sprintf('%s_Up',samples[1])
      if(dim(bindetect[bindetect$sig!="n.s.",][1]==0)){
        logfile("Cutoff is too tight, please lower the cutoff!")
      }
    }
  }
  if (check_name == 3){
    bindetect <- bindetect[,c("name",sprintf("%s_footprints_%s_footprints_change",samples[1], samples[2]),
                              sprintf("%s_footprints_%s_footprints_pvalue",samples[1], samples[2]))]
    colnames(bindetect) <- c("motif","change","p")
    if(is.na(cutoff)){
      cutoff <- as.numeric(quantile(bindetect$change, 0.99))
      bindetect[which(bindetect$change >= cutoff & bindetect$p < 0.05),'sig'] <- sprintf('%s_Up',samples[1])
      bindetect[which(bindetect$change <= - cutoff & bindetect$p < 0.05),'sig'] <- sprintf('%s_Up',samples[2])
      bindetect[which(abs(bindetect$change) <= cutoff | bindetect$p >= 0.05),'sig'] <- 'n.s.'
      res <- bindetect[bindetect$sig!="n.s.",]
      rownames(res) <- NULL
      s1 <- sprintf('%s_Up',samples[1])
      s2 <- sprintf('%s_Up',samples[2])
      #return(bindetect)
    }
    else{
      cutoff <- cutoff
      bindetect[which(bindetect$change >= cutoff & bindetect$p < 0.05),'sig'] <- sprintf('%s_Up',samples[1])
      bindetect[which(bindetect$change <= - cutoff & bindetect$p < 0.05),'sig'] <- sprintf('%s_Up',samples[2])
      bindetect[which(abs(bindetect$change) <= cutoff | bindetect$p >= 0.05),'sig'] <- 'n.s.'
      res <- bindetect[bindetect$sig!="n.s.",]
      rownames(res) <- NULL
      s1 <- sprintf('%s_Up',samples[1])
      s2 <- sprintf('%s_Up',samples[2])
      if(dim(res)[1]==0){
        logfle("Cutoff is too tight, please lower the cutoff!")
      }
    }
  }
  if(plot_type=="Volcano"){
    options(repr.plot.width = 8, repr.plot.height = 7, repr.plot.res = 100)
    p <- ggplot2::ggplot(data = bindetect, aes(x = change, y = -log10(p), color = sig)) +
         ggplot2::geom_point(size = abs(bindetect$change)*scale) +
         ggplot2::scale_color_manual(values = c('#C51B7D', '#E0E0E0', '#4D9221'),limits = c(s1, 'n.s.', s2)) +
         ggplot2::labs(x = 'Differential binding score', y = '-log10 adjust p-value',title = sprintf('%s vs. %s',gsub("_Up","",s1),gsub("_Up","",s2)), color = '') +
         ggplot2::theme(plot.title = element_text(hjust = 0.5, size = 14), panel.grid = element_blank(),
            panel.background = element_rect(color = 'black', fill = 'transparent'),
            legend.key = element_rect(fill = 'transparent')) + ggplot2::geom_vline(xintercept = c(-cutoff, cutoff), lty = 3, color = 'black') +
         ggrepel::geom_label_repel(data=subset(bindetect, sig!="n.s."), aes(label=motif,color=sig), size=3,label.size=NA,max.overlaps=30,box.padding = 0.4, show.legend=FALSE)
  }
  if(plot_type=="Bar"){
    res$score <- -log10(res$p)
    p <- ggpubr::ggbarplot(res, x = "motif", y = "change",fill = "score",color="white",
                   sort.by.groups = FALSE,sort.val ="asc",
                   x.text.angle = 90, ylab = "Differential binding score",xlab = '',
                   rotate = TRUE,ggtheme = theme_pubr()) +
         ggplot2::scale_fill_distiller(palette = "Spectral",name="-log10(p)")+ggpubr::theme_pubr(base_size = 12, base_family = "", border = TRUE, margin = TRUE,
         legend = c("right"), x.text.angle = 0)+ theme(plot.title = element_text(hjust = 0.5,size = 16))+
         ggplot2::geom_hline(yintercept = 0, color = "#7A7A7A", lwd = 0.5) + ggplot2::ggtitle(sprintf('%s vs. %s',gsub("_Up","",s1),gsub("_Up","",s2)))+
         ggplot2::geom_hline(yintercept = 0, color = "#7A7A7A", lwd = 0.5)
  }
  if(plot_type=="Lollipop"){
    res$score <- -log10(res$p)
    p <- ggplot2::ggplot(res, aes(x=reorder(motif,change), y=change,color=score)) +
         ggplot2::geom_segment( aes(x=reorder(motif,change), xend=motif, y=0, yend=change),color=ifelse(res$change>0,"#FF7F0E","#1F77B4")) +
         ggplot2::geom_point( size=abs(res$change)*scale, alpha=1)+ ggplot2::theme_light() + ggplot2::coord_flip() +
         ggplot2::scale_color_distiller(palette = "Spectral",name="-log10(p)") + ggplot2::theme(plot.title = element_text(hjust = 0.5,size = 16)) +
         ggplot2::ggtitle(sprintf('%s vs. %s',gsub("_Up","",s1),gsub("_Up","",s2)))+
         ggplot2::xlab('') + ggplot2::ylab("Differential binding score")
  }
  if(!is.na(save_path)){
    pdf(sprintf("%s/%s_Differential_Footprint_%s.pdf",save_path, file_prefix,plot_type),width=figure_width,height=figure_height)
    print(p)
    dev.off()
  }
  else{
    return(p)
  }
}




#' Convert motif ID to gene ID.
#'
#' @param motif_file The PWM matrix file.
#' @param orgdb The orgdb name.
#' @param tf_family The TF family files. Thee columns are included:Symbol	Ensembl	Family.
#'
#' @return
#' @export
#'
#' @examples   getMotif2Gene(motif_file = "jaspar.motif", orgdb = "org.Hs.eg.db", tf_family="tf.family.tsv")
getMotif2Gene <- function(motif_file, orgdb, tf_family){
  file_lines <- readLines(motif_file)
  header_lines <- grep("^>", file_lines, value = TRUE)
  df <- data.frame(header_lines) %>% tidyr::separate(col = header_lines, into = c("Motif_ID", "Motif_Name"), sep = "\t", remove = TRUE)
  df$Motif_ID <- gsub(">","",df$Motif_ID)
  id1 <- AnnotationDbi::mapIds( x = orgdb, keys = df$Motif_Name, column = c("ENSEMBL"), keytype = c('SYMBOL'),fuzzy = TRUE, multiVals = "first", ignore.case = TRUE)
  df1 <- data.frame(ENSEMBL=id1, Gene_Name=names(id1))
  res <- cbind(df, df1)
  res <- res[complete.cases(res), ]
  family <- read.table(tf_family, head=T)
  res1 <- merge(res, family,by.x="ENSEMBL",by.y="Ensembl",all.x=T)
  res1[is.na(res1)] <- "unknow"
  return(res1[,c(1,2,3,4,6)])
}


#' Plot the correlation of footprint score and TF expression.
#'
#' @param motif2gene The motif to gene annotation file obtained by **getMotif2Gene**.
#' @param gene_exp The normalized gene expression matrix.
#' @param footprint_score The footprint score file obtained by **getFootprintScore**.
#' @param return_matrix Whether to return matrix. Default is FALSE.
#' @param cor_cutoff The correlation cutoff of TF footprint score and TF expression. Default is 0.5.
#' @param color Color list.
#'
#' @return
#' @export
#'
#' @examples  plotFootprintScoreExp(motif2gene=m2g, gene_exp="rsem_exp.tsv", footprint_score = ftscore)
plotFootprintScoreExp <- function(motif2gene, gene_exp, footprint_score, return_matrix = FALSE, cor_cutoff = 0.5, color=NA){
  motif <- motif2gene
  #rownames(motif) <- motif$output_prefix
  #motif <- motif[,2:ncol(motif)]
  motif$id <- sprintf("%s_%s", motif$Motif_Name, motif$Motif_ID)
  rownames(footprint_score) <- footprint_score$output_prefix
  footprint_score <- footprint_score[,2:ncol(footprint_score)]
  ftscore <- footprint_score[motif$id,]
  gene_exp <- read.table(gene_exp, row.names=1, head=T)
  rna <- gene_exp[motif$ENSEMBL,colnames(ftscore)]
  res <- lapply(1:nrow(ftscore),function(i){ data.frame(TF = rownames(ftscore)[i], Correlation = cor(as.numeric(rna[i,]), as.numeric(ftscore[i,])), Family = motif$Family[i])})
  corres <- dplyr::bind_rows(res)
  corres[is.na(corres)] <- 0
  corres$Correlation <- round(corres$Correlation,2)
  if(return_matrix){
  return(corres)
  } else if (is.na(color)){
    res <- corres %>% mutate(is_annotate=ifelse(abs(corres$Correlation) >= cor_cutoff, "yes", "no"))
    freq <- table(res$Family)
    res$new <- ifelse(res$Family %in% names(freq[freq > 10]), res$Family, "Others")
    res <- res[order(res$new),]
    res$label <- sapply(strsplit(res$TF,"_"), `[`, 1)
    dfcol <- data.frame(x=unique(res$new), y=0, label=c(1:length(unique(res$new))))
    p <- ggplot2::ggplot(res, aes(x=new, y=Correlation)) + ggplot2::geom_hline(yintercept=-0.5, linetype="dashed", color = "grey")+
      ggplot2::geom_hline(yintercept=0.5, linetype="dashed", color = "grey")+
      ggplot2::geom_jitter( aes(color=as.factor(new)),size=abs(res$Correlation)*4, width=0.4)+ ggplot2::scale_color_manual(values = paletteer::paletteer_d("ggthemes::Classic_Cyclic"))+
      ggpubr::theme_pubr()+ggplot2::labs(x = "", y = 'Correlation of TF footprint score and TF expression')+
      ggrepel::geom_label_repel(data=subset(res, is_annotate=="yes"), aes(label=label), size=3,label.size=NA,max.overlaps=30,box.padding = 0.4)+
      ggplot2::geom_tile(data = dfcol,aes(x=x,y=y),height=0.2,color = "black",fill = paletteer::paletteer_d("ggthemes::Classic_Cyclic"),alpha = 0.8,show.legend = F)+
      ggplot2::geom_text(data=dfcol,aes(x=x,y=y,label=x),size =4,color ="white")+ggplot2::theme(panel.border = element_rect(colour = "black", fill=NA))+
      ggplot2:: theme(legend.position="none")
    return(p)
  }
  else if(!is.na(color)){
    res <- corres %>% mutate(is_annotate=ifelse(abs(corres$Correlation) >= cor_cutoff, "yes", "no"))
    freq <- table(res$Family)
    res$new <- ifelse(res$Family %in% names(freq[freq > 10]), res$Family, "Others")
    res <- res[order(res$new),]
    res$label <- sapply(strsplit(res$TF,"_"), `[`, 1)
    dfcol <- data.frame(x=unique(res$new), y=0, label=c(1:length(unique(res$new))))
    p <- ggplot2::ggplot(res, aes(x=new, y=Correlation)) + ggplot2::geom_hline(yintercept=-0.5, linetype="dashed", color = "grey")+
      ggplot2::geom_hline(yintercept=0.5, linetype="dashed", color = "grey")+
      ggplot2::geom_jitter( aes(color=as.factor(new)),size=abs(res$Correlation)*4, width=0.4)+ ggplot2::scale_color_manual(values = color)+
      ggpubr::theme_pubr()+ggplot2::labs(x = "", y = 'Correlation of TF footprint score and TF expression')+
      ggrepel::geom_label_repel(data=subset(res, is_annotate=="yes"), aes(label=label), size=3,label.size=NA,max.overlaps=30,box.padding = 0.4)+
      ggplot2::geom_tile(data = dfcol,aes(x=x,y=y),height=0.2,color = "black",fill = paletteer::paletteer_d("ggthemes::Classic_Cyclic"),alpha = 0.8,show.legend = F)+
      ggplot2::geom_text(data=dfcol,aes(x=x,y=y,label=x),size =4,color ="white")+ggplot2::theme(panel.border = element_rect(colour = "black", fill=NA))+
      ggplot2:: theme(legend.position="none")
  }
}






# getAllFPD <- function(samples, flanking_length= 50, N_cores = 10){
#   result <- read.table(sprintf("%s/%s/bindetect_results.txt", FTAnno$bindetect, samples[1]),head = T)
#   motif_list <- result$output_prefix
#   all_lst <- lapply(samples, function(sample){
#     logfile(sprintf("Footprint Depth with %s...", sample))
#     coverage <- valr::read_bigwig(sprintf("%s/%s_q30_corrected.bw", FTAnno$signal, sample))
#     all_lst <- list()
#     #pb <- progress::progress_bar$new(total = length(motif_list))
#     result_list <- parallel::mcmapply(function(x){
#       motif_center <- valr::read_bed(sprintf("%s/%s/%s/beds/%s_all.bed", FTAnno$bindetect, sample, x, x), n_fields = 3)
#       motif_flank <- valr::bed_flank(motif_center, CATAnno$genome, both = flanking_length)
#       flank_score <- valr::bed_map(motif_flank, coverage, .sum = sum(score))
#       avg_flank <- sum(flank_score$.sum, na.rm = T) / (flanking_length * 2)
#       center_score <- valr::bed_map(motif_center, coverage, .sum = sum(score))
#       motif_center$len <- motif_center$end - motif_center$start
#       motif_length <- as.integer(unique(motif_center$len))
#       #avg_center <- sum(center_score$.sum, na.rm = T) / motif_length
#       bg <- (flank_score + center_score) / (flanking_length * 2 + motif_length)
#       fpd <- avg_flank - bg
#       #fpd <- (avg_flank - avg_center) / (avg_flank + avg_center)
#       #pb$tick()
#       return(data.frame(Sample = sample, Motif = x, FPD = round(fpd,2)))
#       }, x = motif_list, mc.cores = N_cores, SIMPLIFY = FALSE)
#     return(do.call(rbind, result_list))
#     })
#   res <- do.call(rbind, all_lst)
#   ## Remove the value <0
#   res$FPD <- ifelse(res$FPD>=0, res$FPD, 0)
#   rownames(res) <- NULL
#   return(res)
# }
#
#
# plotAllFPD <- function(FPD_result){
#   mycolor <- c(paletteer_d("ggthemes::Classic_Cyclic"),rev(paletteer_d("ggthemes::Tableau_20")))
#   FPD_result$is_annotate <- ifelse(FPD_result$FPD>= 2, "yes", "no")
#   FPD_result$gene <- sapply(strsplit(FPD_result$Motif,"_"), `[`, 1)
#   dfcol<-data.frame(x=unique(FPD_result$Sample), y = 0, label=c(1:13))
#   p <- ggplot2::ggplot(FPD_result, aes(x=Sample, y=FPD)) + #ggplot2::geom_hline(yintercept = -1, linetype="dashed", color = "grey") +
#        ggplot2::geom_hline(yintercept = 2, linetype="dashed", color = "grey") +
#        ggplot2::geom_jitter(aes(color=as.factor(Sample)), size = abs(FPD_result$FPD) * 2, width=0.4)+ ggplot2::scale_color_manual(values = paletteer_d("ggthemes::Classic_Cyclic")) +
#        ggpubr::theme_pubr() + ggplot2::labs(x = "", y = 'Footprint Depth Score') +
#        ggplot2::geom_label_repel(data=subset(FPD_result, is_annotate=="yes"), aes(label=gene), size=3,label.size=NA, max.overlaps=30,box.padding = 0.4) +
#        ggplot2::geom_tile(data = dfcol,aes(x=x,y=y), height=0.2,color = "black",fill = mycolor,alpha = 0.8,show.legend = F) +
#        ggplot2::geom_text(data=dfcol,aes(x=x,y=y,label=x),size =4,color ="white") + ggplot2::theme(panel.border = element_rect(colour = "black", fill=NA)) + guides(fill="none")
#   return(p)
# }
#
# compareFPD <- function(FPD_result, samples){
#   fpd1<- FPD_result[FPD_result$Sample %in% c(samples[1], samples[2]),]
#   wide_matrix <- dcast(fpd1, Motif ~ Sample, value.var = "FPD")
#   p<- ggplot(wide_matrix, aes(x = "Bulk_B", y = "CD8pos_T")) +
#   geom_point(color = "steelblue", size = 3) +
#   geom_abline(intercept = 30, slope = -5, linetype = "dashed", color = "darkred") +
#   labs(x = "Weight (1000 lbs)", y = "Miles per Gallon") +
#   theme_minimal()
#   ggscatter(wide_matrix, x = "Bulk_B", y = "CD8pos_T")
# }
