
#' A function to plot the insertion fragment size.
#'
#' @param frag_path Provide a directory of frags from the snakemake flow.
#' @param sample Input a vector containing the names of the samples.
#' @param plot_type Option to draw ridges / density plot. Default: "ridges".
#'
#' @return
#' @export
#'
#' @examples  plotFragments("~/snakemake/fragments_size",sample=c("s1","s2"),plot_type="ridges")
#'
plotFragments <- function(frag_path, sample, plot_type = "ridges"){
  if (length(sample) == 1){
    combined_data <- data.table::fread(sprintf("%s/%s_q30_frag_size.txt",frag_path,sample))
    combined_data$group <- sample
    p <- ggplot2::ggplot(combined_data, ggplot2::aes(x=V3,fill=group)) +
      ggplot2::geom_density(color="white",alpha=0.8)+ggpubr::theme_pubr(border=T)+
      ggplot2::xlab("Fragment size (bp)") +ggplot2::ylab("Density")+
      ggplot2::scale_fill_manual(values=paletteer::paletteer_d("ggthemes::Tableau_10"))
    return(p)
  }
  if (length(sample) > 1){
    cuts_lst <- lapply(sample, function(x){ a <- data.table::fread(sprintf("%s/%s_q30_frag_size.txt",frag_path,x))
    a$group <- x
    return(a)})
    combined_data <- do.call(rbind, cuts_lst)
    options(repr.plot.width = 10, repr.plot.height =5, repr.plot.res = 100)
    if (plot_type=="ridges"){
      p <- ggplot2::ggplot(combined_data, ggplot2::aes(x = V3, y = group, fill=group)) +
        ggridges::geom_density_ridges() +
        ggplot2::scale_y_discrete(expand = c(0, 0)) +
        ggplot2::scale_x_continuous(expand = c(0, 0)) +
        ggplot2::coord_cartesian(clip = "off") + ggplot2::xlab("Fragment size (bp)") + ggplot2::ylab("")+
        ggpubr::theme_pubr() +ggplot2::scale_fill_manual(values=paletteer::paletteer_d("ggthemes::Tableau_10"))
    }else if(plot_type=="density"){
      p <- ggplot2::ggplot(combined_data, ggplot2::aes(x=V3,fill = group)) +
        ggplot2::geom_density() + ggpubr::theme_pubr(border=T)+
        ggplot2::xlab("Fragment size (bp)") +ggplot2::ylab("Density")+
        ggplot2::scale_fill_manual(values=paletteer::paletteer_d("ggthemes::Classic_Cyclic"))+
        ggplot2::facet_wrap(~group)
    }
    return(p)
  }
}


#' A function to plot the Tn5 signal around the TSSs.
#'
#' @param tss_path Provide a directory of TSS from the snakemake flow.
#' @param sample Input a vector containing the names of the samples.
#' @param split_group Wheather to split groups. Dafault: FALSE.
#'
#' @return
#' @export
#'
#' @examples   plotTSS("~/snakemake/tss/",c("s1,"s2"))
#'
plotTSS <- function(tss_path, sample, split_group = FALSE){
  if (length(sample) == 1){
    tss <- data.table::fread(sprintf("%s/%s_tss.csv", tss_path, sample), header=T)
    tss$group <- sample
    p <- ggplot2::ggplot(data=tss, aes(x=pos, y=cut_sum,color=group)) +
      ggplot2::geom_line() + ggplot2::theme_bw() + ggplot2::xlab("Distance to TSS (bp)") +
      ggplot2::ylab("Density") + ggplot2::scale_color_manual(values=paletteer_d("ggthemes::Tableau_10"))
    return(p)
  }
  if (length(sample) > 1){
    tss_lst <- lapply(sample, function(x){ a <- data.table::fread(sprintf("%s/%s_tss.csv", tss_path, x))
    a$group <- x
    return(a)})
    tss <- do.call(rbind, tss_lst)
    p <- ggplot2::ggplot(data = tss, aes(x = pos, y = cut_sum, color = group)) +
      ggplot2::geom_line() + ggplot2::theme_bw() + ggplot2::xlab("Distance to TSS (bp)") +
      ggplot2::ylab("TSS enrichment score") + ggplot2::scale_color_manual(values = paletteer::paletteer_d("ggthemes::Tableau_10"))
    if (split_group){
      p <- ggplot2::ggplot(data = tss, aes(x = pos, y = cut_sum)) +
        ggplot2::geom_line(aes(color = group)) + ggplot2::theme_bw() + ggplot2::xlab("Distance to TSS (bp)") +
        ggplot2::ylab("TSS enrichment score") + ggplot2::scale_color_manual(values = paletteer::paletteer_d("ggthemes::Tableau_10")) + ggplot2::facet_wrap(~group)
    }
    return(p)
  }
}


#' Plot the PCA plot based on quantification.
#'
#' @param norm_data The normalized data obtained by **quantification**.
#' @param save_path The path to save the result.
#' @param figure_height PDF height. Default: 8.27.
#' @param figure_width  PDF width. Default: 8.27.
#'
#' @return
#' @export
#'
#' @examples   plotPCA(norm_data = quant, save_path="~/CAT/")
#'
plotPCA <- function(norm_data, save_path=NA,figure_height=8.27,figure_width=8.27){
  df <- tibble::rownames_to_column(as.data.frame(t(norm_data)), "group")

  pca_res <- stats::prcomp(t(norm_data), scale. = TRUE)

  options(repr.plot.width = 10, repr.plot.height =10, repr.plot.res = 100)
  p <- ggplot2::autoplot(pca_res, frame = TRUE, frame.type = 'norm', label = F, colour = 'group',data=df) +
       ggpubr::theme_pubr() + ggplot2::scale_color_manual(values = pcol30)+
       ggrepel::geom_text_repel(aes(label = group, color = group))+
       ggplot2::theme(legend.position="none")
  if (!is.na(save_path)){
    pdf(sprintf("%s/Samples_PCA.pdf",save_path),height=figure_height, width=figure_width)
    print(p)
    dev.off()
  }
  return(p)
}


#' A function used to plot tree clustering result.
#'
#' @param norm_data The normalized data obtained by **quantification**.
#' @param k   How many clusters do you want to get.
#' @param save_path The path to save the result.
#' @param figure_height PDF height. Default: 8.27.
#' @param figure_width PDF width. Default: 8.27.
#'
#' @return
#' @export
#'
#' @examples   plotDendrogram(norm_data = quant, k=10, save_path="~/CAT/")
#'
plotDendrogram <- function(norm_data, k, save_path=NA, figure_height=8.12, figure_width=8.12){
  #options(repr.plot.width = 10, repr.plot.height =10, repr.plot.res = 100)
  t(norm_data) %>% dist() %>%
    stats::hclust() %>% as.dendrogram() -> dend
  if (!is.na(save_path)){
    pdf(sprintf("%s/Samples_Dendrogram.pdf",save_path),height=figure_height, width=figure_width)
    if (k <= 10){
      dend %>% dendextend::set("labels_col", value = paletteer::paletteer_d("ggthemes::Tableau_10"), k = k) %>%
        dendextend::set("branches_k_color", value = paletteer::paletteer_d("ggthemes::Tableau_10"), k = k)%>% plot(horiz=TRUE, axes=FALSE)
      }
    if (k>10 & k <=20){
      p <- dend %>% dendextend::set("labels_col", value = paletteer::paletteer_d("ggthemes::Tableau_20"), k = k) %>%
        dendextend::set("branches_k_color", value = paletteer::paletteer_d("ggthemes::Tableau_20"), k = k)%>% plot(horiz=TRUE, axes=FALSE)
      }
    dev.off()
  }else{
    par(mar = c(5, 7, 4, 15))
    if (k <= 10){
      dend %>% dendextend::set("labels_col", value = paletteer::paletteer_d("ggthemes::Tableau_10"), k = k) %>%
        dendextend::set("branches_k_color", value = paletteer::paletteer_d("ggthemes::Tableau_10"), k = k)%>% plot(horiz=TRUE, axes=FALSE) %>% return()
    }
    if (k>10 & k <=20){
      p <- dend %>% dendextend::set("labels_col", value = paletteer::paletteer_d("ggthemes::Tableau_20"), k = k) %>%
        dendextend::set("branches_k_color", value = paletteer::paletteer_d("ggthemes::Tableau_20"), k = k)%>% plot(horiz=TRUE, axes=FALSE) %>% return()
    }
  }
  if (k > 20){
    print("Please set the k <= 20 !")
  }
}

