
##############################################################################
#                 Calculating the SPM matrix.                                #
##############################################################################

#' A function to quantify tissue-specificity as the Specificity Measure (SPM).
#'
#' @param quant_data A normalized data matrix obtained by **quantification**.
#' @param plot Wheather to plot figure. Default is TRUE.
#'
#' @return
#' @export
#'
#' @examples calSPM(quant_data)
#'
calSPM <- function(quant_data, plot = TRUE){
  color <- c(as.character(paletteer::paletteer_d("ggthemes::Hue_Circle")),"#D51F26","#272E6A","#208A42","#89288F","#F47D2B", "#FEE500","#8A9FD1","#C06CAB","#E6C2DC",
             "#90D5E4", "#89C75F","#F37B7D","#9983BD","#D24B27","#3BBCA8", "#6E4B9E","#0C727C", "#7E1416","#D8A767","#3D3D3D")
  spm <- function(vector){
    norm_value <- sqrt(sum(vector**2))
    vector[vector!=0] <- vector[vector!=0]**2 / (norm_value*vector[vector!=0])
    return(vector)
  }
  res <- as.data.frame(t(apply(quant_data, 1, function(x) spm(x))))
  df <- reshape2::melt(res)
  if (plot){
    p <- ggplot2::ggplot(df, aes(x=variable, y=value, fill = variable)) + ggplot2::geom_violin(trim=TRUE)+
         ggplot2::geom_boxplot(width=0.1) + ggplot2::theme_minimal() +
         ggplot2::scale_fill_manual(values=color) +
         ggplot2::guides(fill="none")+ggplot2::theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
         ggplot2::xlab("Samples") + ggplot2::ylab("SPM score")
    return(p)
  }
  else{
    return(res)
  }
}


##############################################################################
#                 Get sample top N specific peaks.                           #
##############################################################################


#' Get the top N specific peaks in each sample.
#'
#' @param spm_data The Specificity Measure(SPM) result obtained by **calSPM**.
#' @param norm_data A normalized data matrix obtained by **quantification**.
#' @param top_N The top N specific peaks to return. Default is 1000.
#' @param save_path The path to save file.
#' @param file_prefix The save file prefix.
#'
#' @return
#' @export
#'
#' @examples    getTopSpecifcPeaks(spm, quant_data, top_N = 1000)
getTopSpecifcPeaks <- function(spm_data, norm_data, top_N = 1000, save_path=NA, file_prefix=NA){
  sorted_row_names <- unique(unlist(lapply(spm_data, function(x) {
    sorted <- sort(x)
    row_names <- row.names(spm_data)[order(x)]
    row_names[1:top_N]
  })))
  res <- norm_data[rownames(norm_data) %in% sorted_row_names, ]
  if(!is.na(save_path)){
    write.table(res, sprintf("%s/%s_top%s", save_path, file_prefix, top_N), sep='\t', quote=F)
  }
  return(res)
}

##############################################################################
#                 Calculating the JSD score.                                 #
##############################################################################


#' A function to calculate the JSD score.
#'
#' @param norm_data A normalized data matrix obtained by **quantification**.
#' @param save_path The path to save the JSD score result.
#' @param file_prefix The save file prefix.
#'
#' @return
#' @export
#'
#' @examples  calJSDscore(norm_data = quant_data)
#'
calJSD <- function(norm_data, save_path=NA, file_prefix){
  logfile("Calculating JSD...")
  entropy <- apply(norm_data, 1, function(x){philentropy::H(x/sum(x))})
  H_specificity <- 1 - entropy / log2(dim(norm_data)[2])
  tissue_sp <- as.data.frame(H_specificity)
  return(tissue_sp)
  if(!is.na(save_path)){
    write.table(tissue_sp, sprintf("%s/%s_All_Peaks_JSD.tsv",save_path, file_prefix),sep='\t',quote=F,col.names=T)
  }
}

##############################################################################
#                 Get sample specific peaks.                                 #
##############################################################################


#' A function to identify sample specific peaks.
#'
#' @param norm_data A normalized data matrix obtained by **quantification**.
#' @param jsd_file The JSD result obtained by **calJSD**.
#' @param cutoff The JSD cutoff to determine the sample specific peaks. Default is 0.2.
#' @param save_path The path to save the sample specific peak result.
#' @param file_prefix The save file prefix.
#'
#' @return
#' @export
#'
#' @examples  getSpecificPeak(norm_data, jsd_file, cutoff=0.2)
getSpecificPeak <- function(norm_data, jsd_file, cutoff = 0.2, save_path=NA, file_prefix=NA){
  jsd_res1 <- jsd_file
  jsd_res1$peak <- rownames(jsd_res1)
  logfile(sprintf("Get the TOP %.2f%% specifci peaks.",cutoff*100))
  sp_peak <- head(jsd_res1[order(jsd_res1$H_specificity, decreasing = TRUE),], n = as.integer(cutoff * (dim(jsd_res1)[1])))
  sp_data <- norm_data[rownames(sp_peak),]
  if(!is.na(save_path)){
    write.table(sp_data,sprintf("%s/%s_Specific_Peaks_Quant.tsv",save_path,file_prefix),sep='\t', quote=F)
  }
  return(sp_data)
}

##############################################################################
#                 Plot sample specific peaks.                                #
##############################################################################


#' Plot the sample specific peaks.
#'
#' @param specific_peak A sample-specific peak file obtained by **getSpecificPeak**.
#' @param zscore_min The maximum value of the zscore. Default: 2.
#' @param zscore_max The minium value of the zscore. Default: -2.
#' @param cluster_col Wheather to clustering columns. Default: FALSE.
#' @param color Heatmap color.
#' @param save_path The path to save the sample specific peak pdf file.
#' @param file_prefix The save file prefix.
#' @param figure_height The PDF height. Default: 14.
#' @param figure_width The PDF width. Default: 7.
#'
#' @return
#' @export
#'
#' @examples  plotSpecificPeak(sp)
plotSpecificPeak <- function(specific_peak, zscore_min = -2, zscore_max = 2, cluster_col=FALSE, color=NA, save_path=NA, file_prefix=NA,figure_height=14,figure_width=7){
  scale_mat <- rowZscores(as.matrix(specific_peak), min= -2, max=2, limit=T)
  options(repr.plot.width = 6, repr.plot.height = 14, repr.plot.res = 100)
  #scale_mat <- orderMatrix(mat = as.matrix(specific_peak))
  p <- ComplexHeatmap::pheatmap(scale_mat, cluster_row = T, cluster_cols = T,
                          show_rownames = F, show_colnames = T,
                          color = colorRampPalette(blueYellow)(256), main = "Specific Peaks", border_color = NA,name="Z-score")
  return(p)
  if(!is.na(save_path)){
    pdf(sprintf("%s/%s_Specific_Peaks.pdf", save_path, file_prefix), width = figure_width, height = figure_height)
    print(p)
    dev.off()
  }
}

##############################################################################
#                 Get the optimal cluster number.                            #
##############################################################################


#' A function to help you determine the cluster number before you use **plotClusterSpecificPeak**.
#'
#' @param specific_peak A sample specific peak file obtained by **getSpecificPeak**.
#'
#' @return
#' @export
#'
#' @examples  getClusterNum(sp_data)
getClusterNum <- function(specific_peak){
  scale_mat <- rowZscores(as.matrix(specific_peak), min= -2, max=2, limit=T)
  p <- factoextra::fviz_nbclust(scale_mat, kmeans, method = "wss") +
       labs(subtitle = "Elbow method")
  return(p)
}


##############################################################################
#             Plot sample specific peaks with cluster.                      #
##############################################################################



#' A function to help you plot the sample specific peaks with heatmap.
#'
#' @param specific_peak A sample specific peak file obtained by **getSpecificPeak**.
#' @param zscore_min The maximum value of the zscore. Default: -2.
#' @param zscore_max The minium value of the zscore. Default: 2.
#' @param cluster_num The cluster number obtained by **getClusterNum**.
#' @param cluster_col Wheather to clustering columns. Default: FALSE.
#' @param color Providing your color.
#' @param save_path The path to save the sample specific peak pdf file.
#' @param file_prefix The save file prefix.
#' @param figure_height The PDF height. Default: 14.
#' @param figure_width The PDF width. Default: 7.
#'
#' @return
#' @export
#'
#' @examples  plotClusterSpecificPeak(sp_peaks, cluster_num = 7)
plotClusterSpecificPeak <- function(specific_peak, zscore_min = -2, zscore_max = 2, cluster_N = 10,
                                    cluster_col=FALSE, color=NA, save_path=NA,
                                    file_prefix=NA,figure_height=14,figure_width=7){
  df1 <- specific_peak[complete.cases(specific_peak), ]
  df1 <- rowZscores(as.matrix(df1), limit = TRUE, min = zscore_min, max = zscore_max) %>% as.data.frame()
  row_dend <- hclust(dist(df1))
  #mat <- df1[row_dend$order, ]
  group <- data.frame(C=cutree(row_dend, k = as.integer(cluster_N)))
  group$Cluster <- paste0("Cluster",group$C)
  group$Cluster <- factor(group$Cluster, levels = gtools::mixedsort(unique(group$Cluster)))
  group <- group[order(group$Cluster),]
  #group$Cluster <- factor(group$Cluster, levels = gtools::mixedsort(unique(group$Cluster)))
  mat <- df1[rownames(group),]
  group <- subset(group,select=-(C))
  annotation <- mat[match(rownames(group),rownames(as.data.frame(mat))),]
  mycol <- c(paletteer::paletteer_d("ggthemes::Tableau_10"), paletteer::paletteer_d("ggthemes::Tableau_20"))
  gcols <- setNames(as.character(mycol[1:cluster_N]),unique(group$Cluster))
  gcol <- list(Cluster=gcols)
  # matr <- rowZscores(as.matrix(sp_data1), min = -2, max = 2, limit = T)
  # p <- ComplexHeatmap::pheatmap(matr, show_rownames = F, row_km = cluster_N, border_color = NA)
  # p2 = ComplexHeatmap::draw(p)
  # r.dend <- ComplexHeatmap::row_dend(p2) # If needed, extract row dendrogram
  # names(r.dend) <- c(1:cluster_N)
  # rcl.list <- ComplexHeatmap::row_order(p2)
  # clu_df <- lapply(c(1:cluster_N), function(i){
  #   out <- data.frame(GeneID = rownames(matr[rcl.list[[i]],]),
  #                     Cluster = paste0("Cluster", i),
  #                     stringsAsFactors = FALSE)
  #   return(out)
  # }) %>% do.call(rbind, .)
  # group <- clu_df
  # group$Cluster <- factor(group$Cluster, levels = unique(group$Cluster))
  # rownames(group) <- group$GeneID
  # group <- subset(group, select = -(GeneID))
  # annotation <- matr[match(rownames(group),rownames(matr)), ]
  # annotation <- annotation[complete.cases(annotation), ]
  # gcols <- setNames(as.character(paletteer::paletteer_d("ggthemes::Tableau_20"))[1:cluster_N], unique(group$Cluster))
  # gcol <- list(Cluster = gcols)
  p3 <- ComplexHeatmap::pheatmap(as.matrix(mat), show_rownames = F, cluster_row = F, cluster_col = T, border_color = NA,
                                 use_raster = F, annotation_colors = gcol,
                                 annotation_row = group, gaps_row = cumsum(as.numeric(table(group$Cluster)))[1:(cluster_N-1)],
                                 color =colorRampPalette(blueYellow)(256), annotation_names_row=F,name="ATAC Z-score")
  if(!is.na(save_path)){
      write.table(group,sprintf("%s/%s_Cluster_Group.txt",save_path, file_prefix),sep='\t',col.names=F, row.names=T, quote=F)
      write.table(mat,sprintf("%s/%s_Clustered_Scaled_Matr.tsv",save_path, file_prefix),sep='\t',col.names=T, row.names=T, quote=F)
      pdf(sprintf("%s/%s_Clsuter_Specific_Peaks.pdf",save_path, file_prefix),width=figure_width, height=figure_height)
      print(p3)
      dev.off()
  }
  return(p3)
}

##############################################################################
#                 Get the target genes in each cluster.                      #
##############################################################################


#' A function to find the target genes based on specific peaks.
#'
#' @param cluster_group A group-cluster file obtained by **plotClusterSpecificPeak**.
#' @param save_path  The path to save the sample specific peak pdf file.
#' @param file_prefix  The save file prefix.
#'
#' @return
#' @export
#'
#' @examples  getClusterTargetGenes("F:/CAT/data/Specific_Peaks_Cluster_Group.txt",tss_file = "F:/CAT/data/nip_tss2.bed")
#'
getClusterTargetGenes <- function(cluster_group, save_path=NA, file_prefix){
  checkGeAnno()
  diff_peak <- read.table(cluster_group,row.names=1)
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
  res <- target[, c(1,2,3,4,6,8,11)]
  final <- merge(diff_peak, res, by.x = c("chrom","peak_center1","peak_center2"), by.y=c("chrom","start.x","end.x"),all.x=T)
  final_res <- final[, c(1,5,6,2,4,7,8,9,10)]
  colnames(final_res)[4:9] <- c("peak_center","cluster","TSS","gene","strand","distance")
  if(!is.na(save_path)){
    write.table(final_res, sprintf("%s/%s_Cluster_Peak_Target.tsv",save_path,file_prefix),sep='\t',quote=F,col.names=T,row.names=T)
  }
  return(final_res)
}

##############################################################################
#                           Plot the GO results.                             #
##############################################################################

#' A function to plot the every cluster target genes GO enrichment.
#'
#' @param target The target genes result obtained by **getClusterTargetGenes**.
#' @param orgdb The orgdb annotation package name.
#' @param N_top The top N GO terms showing in the plot. Default: 8.
#' @param scale_size The max and min value of dot size. Default: c(2,12).
#' @param ont The GO type to enrich. Three type are available: "BP", "CC", "MP". Default: BP.
#' @param clustering Whether to clutsering the GO terms. Default: TRUE.
#'
#' @return
#' @export
#'
#' @examples   plotGO(target = res, orgdb= org.Hs.eg.db, N_top = 8, scale_size=c(2,12), ont = "BP", clustering = TRUE)
#'
plotGO <- function(target, orgdb, N_top = 8, scale_size=c(2,12), ont = "BP", clustering = TRUE){
golist <- list()
for (i in unique(target$cluster)){
     logfile(sprintf("GO with %s", i))
     flush.console()
    go <- clusterProfiler::enrichGO(OrgDb=orgdb,
                 gene = unique(target[target$cluster==i,]$gene),
                 pvalueCutoff = 1,
                 qvalueCutoff = 1,
                 keyType = 'ENSEMBL',
                 pAdjustMethod = 'fdr',
                 ont = ont)
    go1 <- data.frame(go)
    go1$cluster <- i
    golist[[i]] <- go1

}
res <- dplyr::bind_rows(golist)
res <- res[res$p.adjust < 0.05,]
res$score <- -log10(res$p.adjust)
res <- res %>% dplyr::group_by(cluster) %>% dplyr::slice_max(order_by = score, n = N_top, with_ties=FALSE)
res$cluster <- factor(res$cluster, levels = gtools::mixedsort(unique(res$cluster)))
if(clustering == FALSE){
p <- ggplot2::ggplot(res, aes(x=cluster, y=Description)) +
     ggplot2::geom_point(aes(size= Count, fill = score),shape=21,alpha=0.9)+
     ggplot2::scale_fill_gradientn(colours =paletteer::paletteer_d("ggsci::pink_material")[1:6])+
     ggplot2::theme_bw()+ggplot2::theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1))+
     ggplot2::scale_size_continuous(range = scale_size)+ ggplot2::ylab('') + ggplot2::xlab("") +
     ggplot2::labs(fill = "-log10(p.adjust)", size="Count")
} else if(clustering){
  mat <- res %>%
  dplyr::select(-ID, -GeneRatio, -BgRatio,  -pvalue, -p.adjust, -qvalue, -geneID, -score) %>%
  tidyr::pivot_wider(names_from = cluster, values_from = Count) %>% data.frame()
  row.names(mat) <- mat$Description  # put gene in `row`
  mat <- mat[,-1] #drop gene column as now in rows
  mat[is.na(mat)] <- 0
  clust <- hclust(dist(mat %>% as.matrix())) # hclust with distance matrix
  res$Description <- factor(res$Description, levels = clust$labels[clust$order])
  v_clust <- hclust(dist(mat %>% as.matrix() %>% t())) # hclust with distance matrix
  res$cluster <- factor(res$cluster, levels = v_clust$labels[v_clust$order])
  p <- ggplot2::ggplot(res, aes(x=cluster, y=Description)) +
     ggplot2::geom_point(aes(size= Count, fill = score),shape=21,alpha=0.9)+
     ggplot2::scale_fill_gradientn(colours =paletteer::paletteer_d("ggsci::pink_material")[1:6])+
     ggplot2::theme_bw()+ggplot2::theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1))+
     ggplot2::scale_size_continuous(range = scale_size)+ ggplot2::ylab('') + ggplot2::xlab("") +
     ggplot2::labs(fill = "-log10(p.adjust)", size="Count")
}
return(p)
}


##############################################################################
#              Plot the GO results with network.                             #
##############################################################################


#' A function to compare the GO enrichment result with Network.
#'
#' @param target  The target genes result obtained by **getClusterTargetGenes**.
#' @param orgdb    The orgdb annotation package name.
#' @param ont  The GO type to enrich. Three type are available: "BP", "CC", "MP". Default: BP.
#' @param showCategory  The max category number to show. Default: 50.
#' @param category_node  The size of category nodes. Default: 0.5.
#' @param line  The width of the line. Default: 0.5.
#' @param category_labelThe size of category label. Default: 1.2.
#'
#' @return
#' @export
#'
#' @examples   plotGONetwork(target = res, orgdb = org.Hs.eg.db)
#'
plotGONetwork <- function(target, orgdb, ont = "BP", showCategory = 50, category_node = 0.5, line=0.5, category_label =1.2){
     numbers <- as.integer(sub("Cluster", "", unique(target$cluster)))
     sorted_data <- unique(target$cluster)[order(numbers)]
     target$cluster <- factor(target$cluster, levels = sorted_data)
     logfile("Compare the GO results with clusters...")
     xx <- clusterProfiler::compareCluster(gene~cluster, data=target,
                     pvalueCutoff = 0.05,
                     keyType = 'ENSEMBL',
                     pAdjustMethod = 'fdr',
                     ont = ont,
                     OrgDb = orgdb)
     xx <- enrichplot::pairwise_termsim(xx)
     p <- enrichplot::emapplot(xx, showCategory = showCategory, cex.params = list(category_node =category_node,
          category_label=category_label, line=line))+
        scale_fill_manual(values = paletteer::paletteer_d("ggthemes::Tableau_10"))
        return(p)
}

