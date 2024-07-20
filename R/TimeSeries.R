
#' Order the rows along with pseudo time.
#'
#' @param mat
#'
#' @return
#' @export
#'
#' @examples
#'
orderMatrix <- function(mat){
  getQuantiles <- function(v = NULL, len = length(v)){
    if(length(v) < len){
      v2 <- rep(0, len)
      v2[seq_along(v)] <- v
    }else{
      v2 <- v
    }
    p <- trunc(rank(v2))/length(v2)
    if(length(v) < len){
      p <- p[seq_along(v)]
    }
    return(p)
  }
  scaled_rows <- mat #t(scale(t(as.matrix(mat))))
  scaled_rows <- scaled_rows[,colSums(scaled_rows)!= 0]
  varQ <- getQuantiles(matrixStats::rowVars(scaled_rows))
  mat <- scaled_rows[order(varQ, decreasing = TRUE), ]
  idx <- order(apply(mat, 1, which.max))
  res <- mat[idx, ]
  return(rowZscores(res,limit = T))
}


#' Fitting the ATAC signal along the time.
#'
#' @param norm_data Input the normalized dynamic matrix.
#' @param Palette Providing the custom Palette from the https://r-charts.com/color-palettes/.
#' @param return_matrix Whether return matrix. Default: FALSE.
#'
#' @return
#' @export
#'
#' @examples
#'
getTimeATAC <- function(norm_data, Palette = NA, return_matrix = FALSE){

  time_points <- seq_along(1:ncol(norm_data))

  interpolate_expression <- function(gene_expression) {
    loess_fit <- loess(gene_expression ~ time_points)
    interpolated_values <- predict(loess_fit, newdata = seq(0, ncol(norm_data), length.out = 2000))
    interpolated_values[!is.finite(interpolated_values)] <- 0
    return(interpolated_values)
  }

  ## Interpolation for each peak
  interpolated_expression <- apply(norm_data, 1, interpolate_expression)
  res <- t(interpolated_expression)

  ## plot
  df <- orderMatrix(res)
  if(return_matrix){
    return(as.data.frame(df))
  }else{
  if(!is.na(Palette[1])){
    p <- ComplexHeatmap::pheatmap(df, col = Palette,
                                  cluster_cols=F,cluster_row=F,show_rownames=F, border_color=NA,
                                  name = "Z-score", top_annotation = ComplexHeatmap::HeatmapAnnotation(Time = 1:ncol(df)),
                                  show_legend=F)+
         ComplexHeatmap::colAnnotation(foo = 1:ncol(df))
  }else{
    p <- ComplexHeatmap::pheatmap(df, col = colorRampPalette(c("#1F77B4", "white", "#FF7F0E"))(100),
                                  cluster_cols=F,cluster_row=F,show_rownames=F, border_color=NA,
                                  name = "Z-score", top_annotation = ComplexHeatmap::HeatmapAnnotation(Time = 1:ncol(df),
                                                                                                       show_legend=F))
  }
  return(p)
}
}

#' Plot the RNA dynamic along the pseudo time.
#'
#' @param peak_time The result of **getTimeATAC** with the parameter return_matrix = TRUE.
#' @param peak2gene The rds file obtained by **getPeak2Gene**.
#' @param rna_matrix The RNA-seq quantification matrix.
#' @param corr_cutoff The cutoff of correlation to get the reliable Peak2gene links. Default: 0.4.
#' @param return_matrix Whether to return the result as data frame. Default: FALSE.
#'
#' @return
#' @export
#'
#' @examples   getTimeRNA(peak_time = peak_time, peak2gene = peak2gene, rna_matrix = rna_matrix, corr_cutoff = 0.4, return_matrix = T)
#'
getTimeRNA <- function(peak_time, peak2gene, rna_matrix, corr_cutoff = 0.4, Palette = NA, return_matrix = FALSE){
  peak2gene <- readRDS(peak2gene)
  p2g1 <- peak2gene[peak2gene$correlations >= corr_cutoff, ]
  p2g1 <- p2g1 %>% dplyr::group_by(Peak) %>% dplyr::slice_max(correlations, n=1) %>% as.data.frame()
  p2g1 <- p2g1[!duplicated(p2g1$Peak),]
  rownames(p2g1) <- p2g1$Peak
  peak <- p2g1[rownames(peak_time), ]
  peak <- na.omit(peak)
  gene_exp <- read.table(rna_matrix, head = T, row.names = 1)
  rna_meta <- peak[,c("Peak","Gene")]
  rna <- merge(rna_meta, gene_exp, by.x="Gene", by.y=0, all.x = T)# %>% dplyr::select(-Peak, -Gene)
  rownames(rna) <- rna$Peak
  res <- rna[rna_meta$Peak,]
  res1 <- res %>% dplyr::select(-Peak, -Gene)
  res1 <- res1[,colSums(res1)!= 0]
  time_points <- seq_along(1:ncol(res1))

  interpolate_expression <- function(gene_expression) {
    loess_fit <- loess(gene_expression ~ time_points)
    interpolated_values <- predict(loess_fit, newdata = seq(0, ncol(res1), length.out = 2000))
    interpolated_values[!is.finite(interpolated_values)] <- 0
    return(interpolated_values)
  }

  # Interpolation for each peak
  interpolated_expression_raw <- as.data.frame(apply(res1, 1, interpolate_expression))
  interpolated_expression <- interpolated_expression_raw[rowSums(interpolated_expression_raw) > 0,] %>% as.matrix()
  res2 <- t(interpolated_expression)
  df1 <- rowZscores(res2,limit = T)
  #column_idx <- which(df[1, ] > 0)[1]
  #df1 <- df[,c(column_idx:ncol(df))]
  if(!is.na(Palette[1])){
    p <- ComplexHeatmap::pheatmap(df1, col = Palette,
                                  cluster_cols=F,cluster_row=F,show_rownames=F, border_color=NA,show_colnames = F,
                                  name = "Z-score", top_annotation = ComplexHeatmap::HeatmapAnnotation(Time = 1:ncol(df1),
                                                                                                       show_legend=F))
  } else{
    p <- ComplexHeatmap::pheatmap(df1, col = colorRampPalette(c("#4575B4", "white", "#D73027"))(100),
                                cluster_cols=F,cluster_row=F,show_rownames=F, border_color=NA, show_colnames = F,
                                name = "Z-score", top_annotation = ComplexHeatmap::HeatmapAnnotation(Time = 1:ncol(df1), show_legend=F))
  }
  if(return_matrix){
    return(as.data.frame(res2))#[,c("Gene","Peak")]))
  }else{
    return(p)
  }
}


#' Get the fitting time genes.
#'
#' @param peak_time The result of **getTimeATAC** with the parameter return_matrix = TRUE.
#' @param peak2gene The rds file obtained by **getPeak2Gene**.
#' @param rna_matrix The RNA-seq quantification matrix.
#' @param corr_cutoff The cutoff of correlation to get the reliable Peak2gene links. Default: 0.4.
#'
#' @return
#' @export
#'
#' @examples   getTimeGene(peak_time = peak_time, peak2gene = peak2gene,rna_matrix = rna_matrix,corr_cutoff = 0.4)
getTimeGene <- function(peak_time, peak2gene, rna_matrix, corr_cutoff = 0.4){
  peak2gene <- readRDS(peak2gene)
  p2g1 <- peak2gene[peak2gene$correlations >= corr_cutoff, ]
  p2g1 <- p2g1 %>% dplyr::group_by(Peak) %>% dplyr::slice_max(correlations, n=1) %>% as.data.frame()
  p2g1 <- p2g1[!duplicated(p2g1$Peak),]
  rownames(p2g1) <- p2g1$Peak
  peak <- p2g1[rownames(peak_time), ]
  peak <- na.omit(peak)
  gene_exp <- read.table(rna_matrix, head = T, row.names = 1)
  rna_meta <- peak[,c("Peak","Gene")]
  rna <- merge(rna_meta, gene_exp, by.x="Gene", by.y=0, all.x = T)# %>% dplyr::select(-Peak, -Gene)
  rownames(rna) <- rna$Peak
  res <- rna[rna_meta$Peak,]
  rownames(res) <- 1:nrow(res)
  #res1 <- res %>% dplyr::select(-Peak, -Gene)
  return(res[,c("Gene","Peak")])
}


#' Plot the fitting time ATAC and RNA.
#'
#' @param norm_data Input the normalized dynamic matrix.
#' @param peak2gene The peak2gene RDS file obtained by **getPeak2Gene**.
#' @param rna_matrix The RNA expression matrix.
#' @param corr_cutoff The cutoff of correlation to get the reliable Peak2gene links. Default: 0.4.
#' @param label The label annotation file.
#' @param ATAC_Palette Color for ATAC heatmap.
#' @param RNA_Palette Color for RNA heatmap.
#'
#' @return
#' @export
#'
#' @examples  plotTimeAll(norm_data = "dynamic_ocrs.tsv", peak2gene = "F:/CAT/example/All_Peak2Gene_links.rds" , rna_matrix = "RSEM_matrix.tsv")
plotTimeAll <- function(norm_data, peak2gene, rna_matrix, corr_cutoff = 0.4, label = NA, ATAC_Palette = NA, RNA_Palette = NA){
  logfile("Calculating ATAC time...")
  atac <- getTimeATAC(norm_data, return_matrix = TRUE)
  rna <- getTimeRNA(peak_time = atac, peak2gene = peak2gene, rna_matrix = rna_matrix, corr_cutoff = corr_cutoff, return_matrix = T)
  pg <- getTimeGene(peak_time = atac, peak2gene = peak2gene, rna_matrix = rna_matrix,corr_cutoff = corr_cutoff)
  atac_mat <- atac[pg$Peak, ]
  rna_matrix <- read.table(rna_matrix, head =T, row.names = 1)
  rna_mat <- rna_matrix[pg$Gene, ]
  time_points <- seq_along(1:ncol(rna_mat))

  interpolate_expression <- function(gene_expression) {
    loess_fit <- loess(gene_expression ~ time_points)
    interpolated_values <- predict(loess_fit, newdata = seq(0, ncol(rna_mat), length.out = 2000))
    interpolated_values[!is.finite(interpolated_values)] <- 0
    return(interpolated_values)
  }
    logfile("Calculating RNA time...")
  # Interpolation for each peak
  interpolated_expression_raw <- as.data.frame(apply(rna_mat, 1, interpolate_expression))
  interpolated_expression <- interpolated_expression_raw[rowSums(interpolated_expression_raw) > 0,] %>% as.matrix()
  res2 <- t(interpolated_expression)
  df1 <- rowZscores(res2,limit = T)
  #column_idx <- which(df[1, ] > 1.9)[1]
  #df1 <- df[,c(column_idx:ncol(df))]
  if(!is.na(ATAC_Palette)){
    p1 <- ComplexHeatmap::pheatmap(rowZscores(as.matrix(atac_mat)), col = ATAC_Palette,
                                   cluster_cols=F, cluster_row=F,show_rownames=F, border_color=NA,show_colnames=F,
                                   name = "ATAC Z-score", top_annotation = ComplexHeatmap::HeatmapAnnotation(Time = 1:ncol(atac_mat), show_legend=F))
    }else{
  p1 <- ComplexHeatmap::pheatmap(rowZscores(as.matrix(atac_mat)), col = colorRampPalette(c("#1F77B4", "white", "#FF7F0E"))(100),
                                 cluster_cols=F, cluster_row=F,show_rownames=F, border_color=NA,show_colnames=F,
                                 name = "ATAC Z-score", top_annotation = ComplexHeatmap::HeatmapAnnotation(Time = 1:ncol(atac_mat),show_legend=F))
    }
  if(!is.na(label)){
    anno <- read.table(label, head=T)
    ha <- ComplexHeatmap::rowAnnotation(foo=ComplexHeatmap::anno_mark(at=match(anno$gene,rownames(as.data.frame(df1))),
                                                                      labels=anno$name))
    if(!is.na(RNA_Palette)){
      p2 <- ComplexHeatmap::pheatmap(df1, col = RNA,
                                     cluster_cols=F, cluster_row=F,show_rownames=F, border_color=NA,show_colnames=F,
                                     name = "RNA Z-score", top_annotation = ComplexHeatmap::HeatmapAnnotation(Time = 1:ncol(df1),show_legend=F))
    }else{
      p2 <- ComplexHeatmap::pheatmap(df1, col = colorRampPalette(c("#4575B4", "white", "#D73027"))(100),
                                 cluster_cols=F, cluster_row=F,show_rownames=F, border_color=NA,show_colnames=F,
                                 name = "RNA Z-score", top_annotation = ComplexHeatmap::HeatmapAnnotation(Time = 1:ncol(df1),show_legend=F))
      }
  }
  p3 <- p1 + p2 + ha
  return(p3)
}



#' Plot the specific gene expression along the time.
#'
#' @param peak_time The result of **getTimeATAC** with the parameter return_matrix = TRUE.
#' @param peak2gene The rds file obtained by **getPeak2Gene**.
#' @param rna_matrix The RNA-seq quantification matrix.
#' @param gene The gene to plot.
#' @param corr_cutoff The cutoff of correlation to get the reliable Peak2gene links. Default: 0.4.
#'
#' @return
#' @export
#'
#' @examples
plotTimeGene <- function(peak_time, peak2gene, rna_matrix, gene, corr_cutoff = 0.4){
    peak2gene <- readRDS(peak2gene)
    p2g1 <- peak2gene[peak2gene$correlations >= corr_cutoff, ]
    p2g1 <- p2g1 %>% dplyr::group_by(Peak) %>% dplyr::slice_max(correlations, n=1) %>% as.data.frame()
    p2g1 <- p2g1[!duplicated(p2g1$Peak),]
    rownames(p2g1) <- p2g1$Peak
    peak <- p2g1[rownames(peak_time), ]
    peak <- na.omit(peak)
    gene_exp <- read.table(rna_matrix, head = T, row.names = 1)
    rna_meta <- peak[,c("Peak","Gene")]
    rna <- merge(rna_meta, gene_exp, by.x="Gene", by.y=0, all.x = T)# %>% dplyr::select(-Peak, -Gene)
    rownames(rna) <- rna$Peak
    res <- rna[rna_meta$Peak,]
    res1 <- res %>% dplyr::select(-Peak, -Gene)
    res1 <- res1[,colSums(res1)!= 0]
    time_points <- seq_along(1:ncol(res1))

    interpolate_expression <- function(gene_expression) {
      loess_fit <- loess(gene_expression ~ time_points)
      interpolated_values <- predict(loess_fit, newdata = seq(0, ncol(res1), length.out = 2000))
      interpolated_values[!is.finite(interpolated_values)] <- 0
      return(interpolated_values)
    }

    # Interpolation for each peak
    interpolated_expression <- apply(res1, 1, interpolate_expression)
    #interpolated_expression <- interpolated_expression_raw[rowSums(interpolated_expression_raw) > 0,] %>% as.matrix()
    res2 <- t(interpolated_expression)
    column_idx <- which(res2[1, ] > 0)[1]
    df <- res2[,c(column_idx:ncol(res2))]
    peakid <- res[res$Gene == gene, "Peak"][1]

    # Plotting
    dfG <- data.frame(exp=as.numeric(df[peakid,]),time=c(1:ncol(df)))
    p <- ggplot2::ggplot(dfG, aes(x=time/10, y=exp)) + ggplot2::geom_jitter(data=dfG, aes(x=time/10,y=exp,color=time/10),size=1,width=7) + ggplot2::geom_smooth(color="#FF7F00") +
         ggplot2::scale_color_viridis_c() + ggplot2::theme_bw()+ ggplot2::xlab("Time") + ggplot2::ylab("Expression (TPM)") + ggplot2::ggtitle(gene) + ggplot2::labs(color = "Time")
    return(p)
}
