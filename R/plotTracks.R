
##############################################################################
#         Plot the ATAC signal with tracks.                                  #
##############################################################################
#' Plot the genome signal with tracks.
#' @param samples_path  The bigwig (bw) file folder.
#' @param samples_suffix  The suffix of the samples. e.g. ".cpm.bw".
#' @param chr  The chromosome ID.
#' @param start   The start position.
#' @param end   The end position.
#' @param peaks  The merged peaks file.
#' @param color  color list to plot tracks.
#' @param back.color  Whether to plot backgroud color.  Default: TRUE.
#'
#' @return
#' @export
#'
#' @examples  plotTracks(samples_path = "./signal", samples_suffix = ".cpm.bw", chr="chr1", start=10000, end=20000, peaks="final_res/meged.peaks.bed")
#'
plotTracks <-  function(samples_path, samples_suffix, chr, start, end, peaks, color=NA, back.color=TRUE){
  checkGeAnno()
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
      mybw$score <- round(mybw$score, 1)
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
      ptrack %>% aplot::insert_bottom(peakvis,height = 0.03) %>% aplot::insert_bottom(trans,height = 0.08)
    }
  }
  else{
    mybw <- loadBigWig(bwfile, chr = chr, start = start, end = end)
    mybw$score <- round(mybw$score, 1)
    ## plot
    colora <- as.character(paletteer::paletteer_d("ggthemes::Hue_Circle"))
    colorb <- as.character(paletteer::paletteer_d("ggthemes::Tableau_20"))
    colorc <- as.character(paletteer::paletteer_d("ggthemes::Tableau_10"))
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

    ptrack %>% aplot::insert_bottom(peakvis,height = 0.03)%>% aplot::insert_bottom(trans,height = 0.08)
  }
}

##############################################################################
#         Plot genome tracks with co-accessible information                  #
##############################################################################

#' Plot genome tracks with co-accessible information.
#' @param samples_path  The bigwig (bw) file folder.
#' @param samples_suffix  The suffix of the samples. e.g. ".cpm.bw".
#' @param chr  The chromosome ID.
#' @param start  The start position.
#' @param end   The end position.
#' @param peaks  The merged peaks file.
#' @param colist  The co-accessibile result obtained by **getCoaccessible**.
#' @param coaccess_cutoff The correlation cutoff of co-accessible.  Default: 0.4.
#' @param line_size The curve line with.  Default: 0.8.
#' @param curvature  The degree of curvature. Default: 0.3.
#' @param color Color list to plot tracks.
#' @param back.color Whether to plot backgroud color.  Default: TRUE.
#'
#' @return
#' @export
#'
#' @examples  plotCoTracks(samples_path = "./signal", samples_suffix = ".cpm.bw", chr="chr1", start=10000, end=20000, peaks="final_res/meged.peaks.bed",colist="./coaccess/tsv").
#'
plotCoTracks <-  function(samples_path, samples_suffix, chr, start, end,
                          peaks, colist, coaccess_cutoff=0.4,line_size=0.8, curvature = 0.3, color=NA, back.color=TRUE){
  checkGeAnno()
  link <- readRDS(colist)
  link$chr <- sapply(strsplit(link$Peak1,":"), `[`, 1)
  link <- link[,c("chr","Center1","Center2","coaccess")]
  colnames(link) <- c("chr","start","end","coaccess")
  link$coaccess <- round(link$coaccess,2)
  link_filter <- link[abs(link$coaccess)>=coaccess_cutoff & link$chr==chr & link$start>= start & link$end <= end, ]
  if(nrow(link_filter)==0){
    print("Your provide region without significant co-accessible links.")
  }
    plink <- linkVis(linkData = link_filter,
                     start = "start",
                     end = "end",
                     link.aescolor = "coaccess",
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
                              selecType = "lt",
        )
        peakvis <- bedVis(bdFile = peaks,
                          chr = chr, region.min = start, region.max = end, fill = "#006BA4", show.legend=F)

        ptrack %>% aplot::insert_bottom(peakvis,height = 0.03)%>% aplot::insert_bottom(plink,height = 0.06*max(abs(link$coaccess))) %>% aplot::insert_bottom(trans,height = 0.08)
      }
    }
    else{
      mybw <- loadBigWig(bwfile, chr = chr, start = start, end = end)
      mybw$score <- round(mybw$score,1)
      ## plot
      colora <- as.character(paletteer::paletteer_d("ggthemes::Hue_Circle"))
      colorb <- as.character(paletteer::paletteer_d("ggthemes::Tableau_20"))
      colorc <- as.character(paletteer::paletteer_d("ggthemes::Tableau_10"))
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

      ptrack %>% aplot::insert_bottom(peakvis,height = 0.03)%>% aplot::insert_bottom(plink,height = 0.06*max(abs(link$coaccess))) %>% aplot::insert_bottom(trans,height = 0.08)
    }
  }


##############################################################################
#         Get the specific gene coordinates and flanks.                      #
##############################################################################

#' Get the specific gene coordinates and flanks.
#' @param gene_name The gene name.
#' @param left  Left flanking. Default: TRUE.
#' @param right Right flanking. Default: TRUE.
#'
#' @return
#' @export
#'
#' @examples
getGeneBed <- function(gene_name, left = 10000, right = 10000){
  checkGeAnno()
  geneBed <- CATAnno$gene
  geneBed <- geneBed[geneBed$name== gene_name, ]
  genome <- CATAnno$genome
  geneFlank <- valr::bed_slop(geneBed, genome, left = left, right = right)
  chr <- geneFlank$chrom
  start <- geneFlank$start
  end <- geneFlank$end
  gene_coord <- c(chr, start, end)
  return(gene_coord)
}

##############################################################################
#         Plot specific gene ATAC signal with tracks.                        #
##############################################################################

#' Plot specific gene ATAC signal with tracks.
#' @param samples_path The bigwig (bw) file folder.
#' @param samples_suffix  The suffix of the samples.
#' @param gene_name The gene name.
#' @param left The gene left flanking to plot. Default is 10000.
#' @param right The gene right flanking to plot. Default is 10000.
#' @param peaks The merged peaks.
#' @param color Color list to plot tracks.
#' @param back.color Whether to plot backgroud color.  Default if TRUE.
#'
#' @return
#' @export
#'
#' @examples  plotGeneTracks(samples_path = "./signal", samples_suffix = ".cpm.bw", gene_name="ENSG00000103888" , peaks="final_res/meged.peaks.bed")
plotGeneTracks <- function(samples_path, samples_suffix, gene_name, left = 10000, right = 10000, peaks, color=NA, back.color=TRUE){
  coords <- getGeneBed(gene_name = gene_name, left = left, right = right)
  p <- plotTracks(samples_path = samples_path, chr = coords[1], start = as.integer(coords[2]), end = as.integer(coords[3]), samples_suffix = samples_suffix, peaks = peaks, color = color,back.color = back.color)
  return(p)
}


##############################################################################
#         Plot specific gene ATAC signal and co-accessible with tracks.      #
##############################################################################

#' Plot specific gene ATAC signal and co-accessible with tracks.
#' @param samples_path The bigwig (bw) file folder.
#' @param samples_suffix The suffix of the samples. e.g. ".cpm.bw".
#' @param gene_name The gene name.
#' @param left The gene left flanking to plot. Default: 30000.
#' @param right The gene right flanking to plot. Default: 30000.
#' @param peaks The merged peaks.
#' @param colist The co-accessible file obtained by **getCoaccessible**.
#' @param coaccess_cutoff The correlation cutoff of co-accessible.  Default: 0.4.
#' @param line_size The curve line with.  Default: 0.8.
#' @param curvature The degree of curvature. Default: 0.3.
#' @param color  Color list to plot tracks.
#' @param back.color Whether to plot backgroud color.  Default: TRUE.
#'
#' @return
#' @export
#'
#' @examples  plotCoGeneTracks(samples_path = "./signal", samples_suffix = ".cpm.bw", gene_name="ENSG00000103888" , peaks="final_res/meged.peaks.bed",colist="./coaccess/tsv")
plotCoGeneTracks <- function(samples_path, samples_suffix, gene_name, left = 30000, right = 30000,
                             peaks, colist, coaccess_cutoff=0.4,line_size=0.8, curvature = 0.3, color=NA, back.color=TRUE){
  coords <- getGeneBed(gene_name = gene_name, left = left, right = right)
  p <- plotCoTracks(samples_path = samples_path, chr = coords[1], start = as.integer(coords[2]), end = as.integer(coords[3]),
                  samples_suffix = samples_suffix, peaks = peaks, color = color,back.color = back.color,colist=colist, coaccess_cutoff=coaccess_cutoff,line_size=line_size, curvature = curvature)
  return(p)
}

