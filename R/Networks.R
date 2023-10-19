
##############################################################################
#                 Get the regulatory network.                                #
##############################################################################

#'  Get the regulatory network.
#' @param sample The sample name. Only one sample support.
#' @param tf The motif ID. Only one motif support.
#' @param peak2gene  The RDS file obtained by **getPeak2Gene**.
#' @param p2g_cutoff  The absolute cutoff of peak2gene correlation. Default: 0.5.
#' @param motif2gene The annotation file obtained by **getMotif2Gene**.
#' @param targets The target only contain TF(1) or contain genes and TFs (2). Default: 1.
#' @param level Consider only the target genes of your TF of interest (1). Consider not only the target genes of your TF of interest, but also more downstream target genes (2). Default: 2.
#'
#' @return
#' @export
#'
#' @examples  getNetwork(sample="Bulk_B", tf="Atf1_MA0604.1", peak2gene="All_Peak2Gene_links.rds", p2g_cutoff = 0.5, motif2gene = m2g)
getNetwork <- function(sample, tf, peak2gene, p2g_cutoff = 0.5, motif2gene, targets = 1, level = 2){
  checkFTAnno()
  getTarget <- function(tf, level){
    ft <- read.table(sprintf("%s/%s/beds/%s_%s_footprints_bound.bed",FTAnno$bindetect, tf, tf, sample))
    ft$peak <- sprintf("%s:%s-%s", ft$V7, ft$V8, ft$V9)
    ft <- ft[,c(4,11)]
    res <- merge(ft,p2g, by.x="peak", by.y="Peak", all = F)
    res <- res[,c(2,3,4)]
    if(targets==1){   ## targets==1, only create the TF regulatory networks.
      net <- merge(res, motif2gene, by.x = "Gene", by.y = "ENSEMBL", all = F)[,c(2,5,3,4)]
      colnames(net) <- c("Source", "Target", "P2G_Correlation", "Motif_ID")
      net$Motif_ID <- sprintf("%s_%s",net$Target, net$Motif_ID)
      net <- net %>% dplyr::group_by(Target) %>% dplyr::mutate(P2G_Correlation = max(P2G_Correlation)) %>% dplyr::distinct() %>% dplyr::mutate(Level = level)
    } else if(targets==2){
      net <- merge(res, motif2gene, by.x = "Gene", by.y = "ENSEMBL", all.x = T)
      net$Motif_Name <- ifelse(is.na(net$Motif_Name), net$Gene, net$Motif_Name)
      net <- net[,c(2,5,3,4)]
      colnames(net) <- c("Source", "Target", "P2G_Correlation", "Motif_ID")
      net$Motif_ID <- sprintf("%s_%s",net$Target, net$Motif_ID)
      net <- net %>% dplyr::group_by(Target) %>% dplyr::mutate(P2G_Correlation = max(P2G_Correlation)) %>% dplyr::distinct() %>% dplyr::mutate(Level = level)
    }
    return(net)
  }

  p2g <- readRDS(peak2gene)
  p2g <- p2g[abs(p2g$correlations) >= p2g_cutoff,]
  network <- getTarget(tf, level=1)
  if(level==1){
      network <- network
    }else if(level==2){
      tarlist <- lapply(grep("_NA", unique(network$Motif_ID), invert=TRUE, value = TRUE), function(x){
        return(getTarget(tf=x, level = 2))
      })
      res1 <- dplyr::bind_rows(tarlist)
      network <- rbind(network, res1) %>% as.data.frame()
    }
  network <- dplyr::distinct(network, Source, Target, .keep_all = TRUE)
  network$Source <- sapply(strsplit(network$Source,"_"), `[`, 1)
  network$Target <- sapply(strsplit(network$Target,"_"), `[`, 1)
  network$Motif_ID <- gsub("_NA", "", network$Motif_ID)
  return(network)
  }



##############################################################################
#                 Plot the regulatory network.                              #
##############################################################################
#' Plot the regulatory network.
#' @param network   The network result obtained by **getNetwork**.
#' @param seed Set seed number to keep the same result. Default: 1.
#'
#' @return
#' @export
#'
#' @examples
plotNetwork <- function(network, seed=1){
  g <- igraph::graph_from_data_frame(d=network, directed=T)
  igraph::E(g)$arrow.size <- 0.5
  igraph::E(g)$arrow.width <- 1
  igraph::V(g)$label.cex <- 0.8
  igraph::V(g)$label.font <- 0.6
  igraph::V(g)$label.dist <-1
  igraph::V(g)$size <- 5

  #igraph::V(g)$size <- ifelse(igraph::degree(g) > 5, 7,4)
  #igraph::V(g)$size <- ifelse(network$Level==1, 7,4)
  set.seed(seed)
  igraph::E(g)$color[igraph::E(g)$P2G_Correlation > 0] <- "#1BA3C6"
  igraph::E(g)$color[igraph::E(g)$P2G_Correlation < 0] <- "#FC719E"
  igraph::E(g)$width <- abs(igraph::E(g)$P2G_Correlation)*2.5
  coords <- igraph::layout_with_fr(g, niter=9999)

  plot(g, layout=coords, vertex.frame.color="#7E756D",
       vertex.label.color = "black",vertex.label.family = "Arial",vertex.color=ifelse(net$Level==1, "#6FB899","#F7C480"), vertex.frame.width=1)
  ## Set Legend.
  legend("bottomright", legend = c("Level 1", "Level 2"),
         col = c("#6FB899","#F7C480"), pch = 20, cex = 0.8, bty = "n")

  legend("bottomleft", legend = c("Positive","Negative"), lwd = 1.5, lty = 1, col = c("#1BA3C6","#FC719E"), cex = 0.8, bty = "n")

}
