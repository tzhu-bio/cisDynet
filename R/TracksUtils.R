#' The following functions were adapted from https://github.com/junjunlab/transPlotR
#' We would like to thank the authors for providing such a great code resource.
#' @title bedVis
#' @name bedVis
#' @author JunZhang
#' @description visualize peaks(bed files).
#'
#' @param bdFile the bed file path, default(NULL).
#' @param chr the chromesome of peak, default(NULL).
#' @param region.min the peak start coordinate, default(NULL).
#' @param region.max the peak end coordinate, default(NULL).
#' @param track.width track width, default(0.1).
#' @param collapse whether collapse the track, default(FALSE).
#' @param fill track fill colors, default(NULL).
#' @param show.legend whether show fill color legend, default(TRUE).
#' @param add.label whether add peak name, default(FALSE).
#' @param label.column the peak name column name, default(NULL).
#' @param label.vjsut the peak label vjust, default(0.1).
#'
#' @return a ggplot object.
#'
#' @export

globalVariables(c("sn","ymin"))

bedVis <- function(bdFile = NULL,
                   chr = NULL,
                   region.min = NULL,
                   region.max = NULL,
                   track.width = 0.1,
                   collapse = FALSE,
                   fill = NULL,
                   show.legend = TRUE,
                   add.label = FALSE,
                   label.column = NULL,
                   label.vjsut = 0.1){
  # loop read bed
  purrr::map_df(1:length(bdFile),function(x){
    tmp <- rtracklayer::import.bed(bdFile[x]) %>%
      data.frame()

    # add name
    tmp$fileName <- strsplit(bdFile[x],split = ".bed") %>% unlist()

    # add sn
    tmp$sn <- x
    return(tmp)
  }) -> bdData

  # filter region data
  if(!is.null(region.min) & !is.null(region.max)){
    regeion.bd <- bdData %>%
      dplyr::filter(seqnames == chr) %>%
      dplyr::filter(start >= region.min & end <= region.max)
  }else{
    regeion.bd <- bdData %>%
      dplyr::filter(seqnames == chr)
  }

  # whether collapse track
  if(collapse == TRUE){
    regeion.bd$sn <- 1
  }

  # add y region
  regeion.bd <- regeion.bd %>%
    dplyr::mutate(ymin = sn - track.width,
                  ymax = sn + track.width)

  # plot
  bed <-
    ggplot2::ggplot(regeion.bd) +
    ggplot2::geom_rect(ggplot2::aes_string(xmin = 'start',xmax = 'end',
                                           ymin = 'ymin',ymax = 'ymax',
                                           fill = 'fileName'),
                       show.legend = show.legend) +
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   axis.text = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(),
                   axis.title = ggplot2::element_blank()) +
    ggplot2::guides(fill = ggplot2::guide_legend(title = ''))

  if(!is.null(fill)){
    bed.col <- bed +
      ggplot2::scale_fill_manual(values = fill)
  }else{
    bed.col <- bed
  }

  # ========================================
  # add peak label
  if(add.label == TRUE){
    bed.label <- bed.col +
      ggplot2::geom_text(ggplot2::aes(x = (start + end)/2,y = ymin - label.vjsut,
                                      label = get(label.column)),
                         check_overlap = TRUE)
  }else{
    bed.label <- bed.col
  }
  return(bed.label)
}

# title calcuRotatedCoord
#' @name calcuRotatedCoord
#' @author JunZhang
#' @description calculate the rotated rectangle coordinate with specified degree.
#'
#' @param data data.frame
#' @param theta rotate degree, default(45).
#' @param workers how many workers for parallel calculation, default(1).
#' @param rx x variable name, default(NULL).
#' @param ry y variable name, default(NULL).
#'
#' @return a data.frame
#' @export

globalVariables(c('.data','multisession','xr','yr'))

calcuRotatedCoord <- function(data = NULL,
                              rx = NULL,
                              ry = NULL,
                              theta = 45,
                              workers = 1){
  # Set a "plan" for how the code should run.
  future::plan(future::multisession, workers = workers)

  # get coord
  furrr::future_map_dfr(1:nrow(data),function(i){
    tmp <- data[i,]

    x <- rx
    y <- ry

    tmp <- tmp %>%
      dplyr::mutate(xr = .data[[x]]*cos(pi*(theta/180)) + .data[[y]]*sin(pi*(theta/180)),
                    yr = .data[[y]]*cos(pi*(theta/180)) - .data[[x]]*sin(pi*(theta/180)))

    tmp <- tmp %>%
      dplyr::mutate(xr = xr*cos(pi*(theta/180)),
                    yr = yr*sin(pi*(theta/180)))
    return(tmp)
  }) -> data

  return(data)
}


#' @title linkVis
#' @name linkVis
#' @author JunZhang
#' @description visualize the coordinate relation like chromtin accessbility or peak sites correlation.
#'
#' @param linkData the link data with data.frame format, default(NULL).
#' @param start link start position, default(NULL).
#' @param end link end position, default(NULL).
#' @param group facet group variable name, default(NULL).
#' @param link.aescolor link line color or mapping variable name, default(NULL).
#' @param link.color colors to change the link line colors when "link.aescolor" is a mapping variable, default(NULL).
#' @param line.size link line size, default(0.5).
#' @param curvature the link line curvature, default(0.5).
#' @param yshift the space upper the link line, default(0.1).
#' @param legend.title the legend title, default("").
#' @param facet whether show the plot with facet plot, default(TRUE).
#' @param facet.placement the facet label placement, default("outside").
#' @param facet.fill facet rectangle fill color, default("grey90").
#' @param facet.color facet rectangle border color, default("black").
#' @param facet.text.angle facet text angle, default(90).
#' @param facet.text.size facet text size, default(14).
#' @param xAixs.info whether remove X axis info, default(FASLE).
#'
#' @param base_size theme base_size, default(14).
#'
#' @return a ggplot object.
#' @export

linkVis <- function(linkData = NULL,
                    start = NULL,
                    end = NULL,
                    group = NULL,
                    base_size = 14,
                    link.aescolor = NULL,
                    link.color = NULL,
                    line.size = 0.5,
                    curvature = 0.5,
                    yshift = 0.1,
                    legend.title = "",
                    facet = TRUE,
                    facet.placement = "outside",
                    facet.fill = "grey90",
                    facet.color = "black",
                    facet.text.angle = 90,
                    facet.text.size = 14,
                    xAixs.info = TRUE){
  # get input data
  data <- linkData
  colname <- colnames(data)

  # whether facet by groups
  if(facet == TRUE){
    if(link.aescolor %in% colname){
      p1 <-
        ggplot2::ggplot(linkData) +
        ggplot2::geom_curve(ggplot2::aes_string(x = start,xend = end,
                                                y = group,yend = group,
                                                color = link.aescolor),
                            curvature = curvature,
                            lwd = line.size)
    }else{
      p1 <-
        ggplot2::ggplot(linkData) +
        ggplot2::geom_curve(ggplot2::aes_string(x = start,xend = end,
                                                y = group,yend = group),
                            color = link.aescolor,
                            curvature = curvature,
                            lwd = line.size)
    }
  }else{
    linkData$yc <- as.character(1)

    if(link.aescolor %in% colname){
      p1 <-
        ggplot2::ggplot(linkData) +
        ggplot2::geom_curve(ggplot2::aes_string(x = start,xend = end,
                                                y = "yc",yend = "yc",
                                                color = link.aescolor),
                            curvature = curvature,
                            lwd = line.size)
    }else{
      p1 <-
        ggplot2::ggplot(linkData) +
        ggplot2::geom_curve(ggplot2::aes_string(x = start,xend = end,
                                                y = "yc",yend = "yc"),
                            color = link.aescolor,
                            curvature = curvature,
                            lwd = line.size)
    }
  }

  # ajust y
  p1 <- p1 +
    ggplot2::scale_y_discrete(expand = ggplot2::expansion(mult = c(0,-yshift))) +
    ggplot2::theme_bw(base_size) +
    ggplot2::xlab('') + ggplot2::ylab('')

  # whether mapping color
  if(link.aescolor %in% colname){
    if(is.null(link.color)){
      if(is.numeric(linkData[,link.aescolor])){
        p2 <- p1 +
          ggplot2::scale_color_gradient(low = "#339933",
                                        high = "#FF9900",
                                        name = legend.title)
      }else{
        p2 <- p1 +
          ggplot2::scale_color_discrete(name = legend.title)
      }
    }else{
      if(is.numeric(linkData[,link.aescolor])){
        p2 <- p1 +
          ggplot2::scale_color_gradient(low = link.color[1],
                                        high = link.color[2],
                                        name = legend.title)
      }else{
        p2 <- p1 +
          ggplot2::scale_color_manual(values = link.color,
                                      name = legend.title)
      }
    }
  }else{
    p2 <- p1
  }

  # facet
  if(facet == TRUE){
    p3 <- p2 +
      ggplot2::theme(axis.ticks.y = ggplot2::element_blank(),
                     axis.text.y = ggplot2::element_blank(),
                     panel.grid = ggplot2::element_blank(),
                     strip.placement = facet.placement,
                     strip.background = ggplot2::element_rect(fill = facet.fill,
                                                              color = facet.color),
                     strip.text.y.left = ggplot2::element_text(angle = facet.text.angle,
                                                               size = facet.text.size)) +
      ggplot2::facet_wrap(facets = group,
                          ncol = 1,scales = "free_y",
                          strip.position = "left")
  }else{
    p3 <- p2 +
      ggplot2::theme(axis.ticks.y = ggplot2::element_blank(),
                     axis.text.y = ggplot2::element_blank(),
                     panel.grid = ggplot2::element_blank())
  }

  # whether remove X axis info
  if(xAixs.info == FALSE){
    p4 <- p3 +
      ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                     axis.ticks.x = ggplot2::element_blank()) + ylab("Co-accessibility")
  }else{
    p4 <- p3 + ylab("Co-accessibility")
  }
  return(p4)
}



#' @title loadBigWig
#' @name loadBigWig
#' @author JunZhang
#' @description read bigwig files.
#'
#' @param bwFile the path of bigwig files, default(NULL). bigwig files should end with ".bw" or ".bigwig" and directory should not be named "bw" or"bigwig".
#'
#' @return a data.frame
#' @export

loadBigWig <- function(bwFile = NULL, chr = NULL, start, end){
  # loop read bed
  which <- GenomicRanges::GRanges(c(chr), IRanges::IRanges(c(start), c(end)))
  purrr::map_df(1:length(bwFile),function(x){
    tmp <- rtracklayer::import.bw(bwFile[x],which=which) %>%
      data.frame() %>%
      dplyr::select(-width,-strand)

    # sampe name
    spt <- strsplit(bwFile[x],split = "/|.bw|.bigwig") %>% unlist()
    sname <- spt[length(spt)]
    # add name
    tmp$fileName <- sname

    return(tmp)
  }) -> bWData
  return(bWData)
}



#' @title loadBigWig
#' @name loadBigWig
#' @author JunZhang
#' @description read bigwig files.
#'
#' @param bwFile the path of bigwig files, default(NULL). bigwig files should end with ".bw" or ".bigwig" and directory should not be named "bw" or"bigwig".
#'
#' @return a data.frame
#' @export

loadFullBigWig <- function(bwFile = NULL){
  # loop read bed
  purrr::map_df(1:length(bwFile),function(x){
    tmp <- rtracklayer::import.bw(bwFile[x]) %>%
      data.frame() %>%
      dplyr::select(-width,-strand)

    # sampe name
    spt <- strsplit(bwFile[x],split = "/|.bw|.bigwig") %>% unlist()
    sname <- spt[length(spt)]
    # add name
    tmp$fileName <- sname

    return(tmp)
  }) -> bWData
  return(bWData)
}




#' @title trackVis
#' @name trackVis
#' @author JunZhang
#' @description visualize bigwig files.
#'
#' @param bWData the data.frame bigwig data, default(NULL).
#' @param gtf.file whether supply gtf annotation file, default(NULL).
#' @param gene.name the gene name to be choosed for visualization, default(NULL).
#' @param chr chr the chromesome of peak, default(NULL).
#' @param region.min the start coordinate, default(NULL).
#' @param region.max the end coordinate, default(NULL).
#' @param show.legend whether show color legend, default(FALSE).
#' @param legend.position the legend position, default("right").
#' @param color the track color, default(NULL).
#' @param extend.up extend for upstream of start site, default(3000).
#' @param extend.dn extend for downstream of start site, default(3000).
#' @param base_size theme base size, default(14).
#' @param label.angle the facet label angle, default(0).
#' @param label.face the facet label face, default("bold").
#' @param space.y the facet panel space, default(0.5).
#' @param sample.order the sample order to be plotted in graph, default(NULL).
#' @param sampleName.dist the facet label distance from Y axis, default(0.4).
#' @param sampleName.hjust the facet label hjust, default(1).
#' @param sampleName.vjust the facet label vjust, default(0.5).
#' @param xAxis.info whether retain X axis info, default(TRUE).
#' @param yAxis.info whether retain Y axis info, default(TRUE).
#' @param ticks.y.len the y axis ticks length, default(0.3).
#' @param theme plot theme, "bw" or "classic", default("classic").
#' @param scales the facet scales settings, default("fixed").
#' @param ncol the columns to be arranged, default(1).
#' @param mark.region whether highlight regions in plot, default(FALSE).
#' @param mark.col the colors of marked regions, default(NULL).
#' @param mark.alpha the color alpha of marked regions, default(0.5).
#' @param new.yaxis whether add new style Y axis, default(FALSE).
#' @param pos.ratio the new style Y axis relative position, default(c(0.01,0.8)).
#' @param yinfo.text.size the new style Y axis text size, default(5).
#'
#' @param back.color whether add panel background color, default(FALSE).
#' @param back.color.alpha panel background color alpha, default(0.15).
#' @param y.max the ylim, default(NULL).
#' @param new.label whether add label in plot, default(FALSE).
#' @param label.color the label color, default(NULL).
#' @param pos.label.ratio the new label relative position, default(c(0.99,0.8)).
#' @param label.text.size the new label text size, default(5).
#' @param label.hjust the new label text hjust, default(1).
#' @param yinfo.hjust the new style Y axis text hjust, default(0).
#' @param facetGroup the annotation for samples, default(NULL).
#' @param annoLine.size the annotation line size, default(1).
#' @param line.arrow the annotation line arrow, default(NULL).
#'
#' @return a ggplot object.
#' @export

globalVariables(c("fileName","group","label","score","x","y"))

trackVis <- function(bWData = NULL,
                     gtf.file = NULL,
                     gene.name = NULL,
                     chr = NULL,
                     region.min = NULL,
                     region.max = NULL,
                     show.legend = FALSE,
                     legend.position = "right",
                     color = NULL,
                     extend.up = 3000,
                     extend.dn = 3000,
                     base_size = 14,
                     label.angle = 0,
                     label.face = "bold",
                     space.y = 0.5,
                     sample.order = NULL,
                     sampleName.dist = 0,
                     sampleName.hjust = 1,
                     sampleName.vjust = 0.5,
                     facetGroup = NULL,
                     annoLine.size = 1,
                     line.arrow = NULL,
                     y.max = NULL,
                     xAxis.info = TRUE,
                     yAxis.info = TRUE,
                     ticks.y.len = 0.3,
                     theme = "classic",
                     scales = "fixed",
                     ncol = 1,
                     mark.region = NULL,
                     mark.col = NULL,
                     mark.alpha = 0.5,
                     new.yaxis = FALSE,
                     pos.ratio = c(0.01,0.8),
                     yinfo.text.size = 5,
                     yinfo.hjust = 0,
                     new.label = FALSE,
                     label.color = NULL,
                     pos.label.ratio = c(0.99,0.8),
                     label.text.size = 5,
                     label.hjust = 1,
                     back.color = FALSE,
                     back.color.alpha = 0.15){
  # whether supply gene name
  if(!is.null(gtf.file) & !is.null(gene.name)){
    gene <- gtf.file %>% dplyr::filter(gene_name == gene.name)
    chr <- unique(gene$seqnames) %>% as.character()
    region.min <- min(gene$start)
    region.max <- max(gene$start)
  }

  # filter specified region
  regeion.bw <- bWData %>%
    dplyr::filter(seqnames == chr) %>%
    dplyr::filter(start >= (region.min - extend.up) & end <= (region.max + extend.dn))

  # whether change order
  if(!is.null(sample.order)){
    regeion.bw$fileName <- factor(regeion.bw$fileName,levels = sample.order)
  }else{
    regeion.bw$fileName <- factor(regeion.bw$fileName,levels = unique(regeion.bw$fileName))
  }

  # panel background data
  dback <- data.frame(fileName = factor(unique(regeion.bw$fileName),
                                        levels = levels(regeion.bw$fileName)))
  # dback <- regeion.bw[,7:ncol(regeion.bw)] %>% unique() %>% data.frame()

  # whether add background
  if(back.color == TRUE){
    p0 <-
      ggplot2::ggplot(regeion.bw) +
      ggplot2::geom_rect(data = dback,
                         ggplot2::aes(xmin = -Inf,xmax = Inf,ymin = 0,ymax = Inf,
                                      fill = fileName),
                         show.legend = FALSE,
                         alpha = back.color.alpha)
  }else{
    p0 <-
      ggplot2::ggplot(regeion.bw)
  }

  # plot
  p1 <- p0 +
    ggplot2::geom_rect(ggplot2::aes(xmin = start,xmax = end,
                                    ymin = 0,ymax = score,
                                    fill = fileName,color = fileName),
                       show.legend = show.legend) +
    ggplot2::geom_segment(ggplot2::aes(x = min(start),xend = max(end),
                                       y = 0,yend = 0),
                          size = 1) +
    ggplot2::xlab('') + ggplot2::ylab('') +
    ggplot2::coord_cartesian(expand = 0) +
    ggplot2::guides(fill = ggplot2::guide_legend(title = '')) +
    ggplot2::guides(color = ggplot2::guide_legend(title = ''))

  if(is.null(facetGroup)){
    p1 <- p1 +
      ggplot2::facet_wrap(~fileName,ncol = ncol,
                          strip.position = 'left',scales = scales)
  }else{
    p1 <- p1 +
      ggh4x::facet_nested_wrap(facets = c(facetGroup,"fileName"),
                               nest_line = ggplot2::element_line(colour = 'black',
                                                                 size = annoLine.size,
                                                                 arrow = line.arrow),
                               strip.position = 'left',
                               scales = scales,ncol = ncol)
  }

  # y labels
  if(is.null(y.max)){
    if(scales == "fixed"){
      ylimit <- range(regeion.bw$score)[2]
    }else{
      # free_y
      if(is.null(facetGroup)){
        ylimit <- regeion.bw %>%
          dplyr::group_by(fileName) %>%
          dplyr::summarise(maxScore = max(score))
      }else{
        ylimit <- regeion.bw %>%
          dplyr::group_by(.data[[facetGroup]],fileName) %>%
          dplyr::summarise(maxScore = max(score))
      }
    }
  }else{
    if(scales == "fixed"){
      ylimit <- y.max
    }else{
      # free_y
      if(is.null(facetGroup)){
        ylimit <- regeion.bw %>%
          dplyr::group_by(fileName) %>%
          dplyr::summarise(maxScore = max(score))

        # add new range
        ylimit$minScore <- unlist(y.max[[1]])
        ylimit$maxScore <- unlist(y.max[[2]])
      }else{
        ylimit <- regeion.bw %>%
          dplyr::group_by(.data[[facetGroup]],fileName) %>%
          dplyr::summarise(maxScore = max(score))

        # add new range
        ylimit$minScore <- unlist(y.max[[1]])
        ylimit$maxScore <- unlist(y.max[[2]])
      }
    }
  }

  # change y range
  if(scales == "fixed"){
    p1 <- p1 +
      ggplot2::scale_y_continuous(breaks = c(0,ylimit),
                                  limits = c(0,ylimit))
  }else{
    if(!is.null(y.max)){
      low <- unlist(y.max[[1]])
      high <- unlist(y.max[[2]])

      # facet new scales
      lapply(1:length(unique(regeion.bw$fileName)), function(x){
        ggplot2::scale_y_continuous(limits = c(low[x],high[x]),
                                    breaks = c(low[x],high[x]))
      }) -> new.scales

      # add new scales
      p1 <- p1 +
        ggh4x::facetted_pos_scales(y = new.scales)
    }else{
      p1 <- p1
    }
  }

  # choose theme
  if(theme == "bw"){
    p1.1 <- p1 +
      ggplot2::theme_bw(base_size = base_size) +
      ggplot2::theme(strip.text.y.left = ggplot2::element_text(angle = label.angle,
                                                               face = label.face,
                                                               size = base_size,
                                                               hjust = sampleName.hjust,
                                                               vjust = sampleName.vjust),
                     panel.grid = ggplot2::element_blank(),
                     legend.position = legend.position,
                     strip.placement = 'outside',
                     strip.background = ggplot2::element_rect(fill = NA,colour = NA),
                     strip.switch.pad.wrap = ggplot2::unit(sampleName.dist,'cm'),
                     panel.spacing = ggplot2::unit(space.y,'cm'))

  }else{
    p1.1 <- p1 +
      ggplot2::theme_classic(base_size = base_size) +
      ggplot2::theme(strip.text.y.left = ggplot2::element_text(angle = label.angle,
                                                               face = label.face,
                                                               size = base_size,
                                                               hjust = sampleName.hjust,
                                                               vjust = sampleName.vjust),
                     panel.grid = ggplot2::element_blank(),
                     legend.position = legend.position,
                     strip.placement = 'outside',
                     strip.background = ggplot2::element_rect(fill = NA,colour = NA),
                     strip.switch.pad.wrap = ggplot2::unit(sampleName.dist,'cm'),
                     panel.spacing = ggplot2::unit(space.y,'cm'),
                     axis.ticks.length.y = ggplot2::unit(ticks.y.len,"cm"))
  }

  # whether supply own colors
  if(!is.null(color)){
    p2 <- p1.1 +
      ggplot2::scale_fill_manual(values = color) +
      ggplot2::scale_color_manual(values = color)
  }else{
    p2 <- p1.1
  }

  # whether retain X axis info
  if(xAxis.info == FALSE){
    p3 <- p2 +
      ggplot2::theme(axis.ticks.x = ggplot2::element_blank(),
                     axis.text.x = ggplot2::element_blank())
  }else{
    p3 <- p2
  }

  # whether retain Y axis info
  if(yAxis.info == FALSE){
    p4 <- p3 +
      ggplot2::theme(axis.ticks.y = ggplot2::element_blank(),
                     axis.text.y = ggplot2::element_blank())
  }else{
    p4 <- p3
  }

  # whether mark some regions
  if(!is.null(mark.region)){
    if(is.list(mark.region)){
      mark.df.tmp <- data.frame(start = unlist(mark.region[[1]]),
                                end = unlist(mark.region[[2]]),
                                group = as.character(1:length(unlist(mark.region[[1]]))))

      lapply(1:length(unique(regeion.bw$fileName)), function(x){
        tmp <- mark.df.tmp
        tmp$fileName <- unique(regeion.bw$fileName)[x]
        return(tmp)
      }) %>% do.call("rbind",.) %>% data.frame() -> mark.df

      mark.df$fileName <- factor(mark.df$fileName,levels = levels(regeion.bw$fileName))

      # add mark
      p5 <- p4 +
        ggnewscale::new_scale_fill() +
        ggplot2::geom_rect(data = mark.df,
                           ggplot2::aes(xmin = start,xmax = end,
                                        ymin = 0,ymax = ylimit,
                                        fill = group),alpha = mark.alpha,
                           show.legend = FALSE)

      # change mark colors
      if(!is.null(mark.col)){
        p5 <- p5 +
          ggplot2::scale_fill_manual(values = mark.col)
      }else{
        p5 <- p5
      }
    }else{
      print("Please supply list object!")
    }
  }else{
    p5 <- p4
  }

  # whether add new yaxis
  if(new.yaxis == TRUE){
    if(scales == "fixed"){
      yinfo <- data.frame(label = paste("[0-",ylimit,"]",sep = ''),
                          x = range(regeion.bw$start)[1] + pos.ratio[1]*(range(regeion.bw$start)[2] - range(regeion.bw$start)[1]),
                          y = pos.ratio[2]*ylimit)
    }else{
      if(is.null(y.max)){
        yinfo <- data.frame(label = paste("[0-",ylimit$maxScore,"]",sep = ''),
                            fileName = ylimit$fileName,
                            group = if("group" %in% colnames(ylimit)){ylimit$group}else{1},
                            x = range(regeion.bw$start)[1] + pos.ratio[1]*(range(regeion.bw$start)[2] - range(regeion.bw$start)[1]),
                            y = pos.ratio[2]*ylimit$maxScore)
      }else{
        yinfo <- data.frame(label = paste("[",ylimit$minScore,"-",ylimit$maxScore,"]",sep = ''),
                            fileName = ylimit$fileName,
                            group = if("group" %in% colnames(ylimit)){ylimit$group}else{1},
                            x = range(regeion.bw$start)[1] + pos.ratio[1]*(range(regeion.bw$start)[2] - range(regeion.bw$start)[1]),
                            y = pos.ratio[2]*ylimit$maxScore)
      }
    }

    # add text label
    p6 <- p5 +
      ggplot2::geom_text(data = yinfo,
                         ggplot2::aes(x = x,y = y,label = label),
                         size = yinfo.text.size,
                         hjust = yinfo.hjust)
  }else{
    p6 <- p5
  }

  # whether add new sample label
  if(new.label == TRUE){
    if(is.null(facetGroup)){
      group <- NA
      fileName <- unique(regeion.bw$fileName)
    }else{
      tmp <- regeion.bw %>%
        dplyr::select(fileName,group) %>%
        unique()
      fileName <- tmp$fileName
      group <- tmp$group
    }

    # label info
    labelinfo <- data.frame(fileName = fileName,
                            group = group,
                            x = range(regeion.bw$start)[1] + pos.label.ratio[1]*(range(regeion.bw$start)[2] - range(regeion.bw$start)[1]),
                            y = pos.label.ratio[2]*ylimit)

    labelinfo$fileName <- factor(labelinfo$fileName,levels = levels(regeion.bw$fileName))

    # add text label
    if(is.null(label.color)){
      label.color <- rep('black',nrow(labelinfo))
    }else{
      label.color <- label.color
    }

    # plot
    p7 <- p6 +
      ggnewscale::new_scale_color() +
      ggplot2::geom_text(data = labelinfo,
                         ggplot2::aes(x = x,y = y,label = fileName,color = fileName),
                         show.legend = FALSE,
                         size = label.text.size,
                         hjust = label.hjust,
                         fontface = label.face) +
      ggplot2::scale_color_manual(values = label.color)
    # ggplot2::theme(strip.text.y.left = ggplot2::element_blank())
    if(is.null(facetGroup)){
      p7 <- p7 +
        ggplot2::theme(strip.text.y.left = ggplot2::element_blank())
    }else{
      p7 <- p7
    }
  }else{
    p7 <- p6
  }

  return(p7)
}




#' @title trancriptVis
#' @name trancriptVis
#' @author JunZhang
#' @description This package is to visualize gene diffrent isoforms.
#' @param gtfFile GTF file.
#' @param gene Target gene to plot.
#' @param myTranscript Specify which transcripts to plot use transcipt id.
#' @param Chr Chromosome number.
#' @param posStart Region start position on genome.
#' @param posEnd Region end position on genome.
#' @param collapse Whether to collapse multiple transcripts into one, default(FALSE).
#' @param exonWidth Exon width to plot, default(0.3).
#' @param relTextDist Transcripts name or gene name relative to exon, default(0.3).
#' @param intronSize Intron line size, default(0.5).
#' @param arrowBreak How many gap distance to draw arrows, the smaller the more arrows, default(0.15).
#' @param exonColorBy Whether color group by "transcript_id" or "gene_name", default(NULL).
#' @param exonFill Exon fill color, default('#333399').
#' @param circle Whether make plot into a circle plot, default(FALSE).
#' @param cicStart Circle plot start position, default(pi).
#' @param circSegCol Circle sgement color, default('#333399').
#' @param text_only When circle plot labeled by gene name, whether remove the line connected with gene name, default(FALSE).
#' @param ylimLow The Y axis lower limitation of Circle plot, default(-10).
#' @param openAngle The gap of the circle plot, default(0.5).
#' @param arrowCol Normal arrow color, default('#333399').
#' @param arrowAngle Normal arrow angle, default(30).
#' @param arrowLength Normal arrow length, default(0.1).
#' @param arrowType Normal arrow type, default('open').
#' @param addNormalArrow Whether add normal arrow on plot, default(TRUE).
#' @param newStyleArrow Whether add new style arrow on plot, default(FALSE).
#' @param absSpecArrowLen Whether make new style arrow length to be relative to each transcript length or absolute length to the longest transcript, default(FALSE).
#' @param speArrowRelPos The relative position to the transcript on horizontal direction of new style arrow, default(0).
#' @param speArrowRelLen The relative length to the transcript length of new style arrow, default(0.05).
#' @param speArrowStart The new style arrow start position on the vertical direction, default(-0.15).
#' @param speArrowRelHigh The relative height of new style arrow to the vertical length, default(2).
#' @param speArrowLineSize The new style arrow line size, default(0.5).
#' @param speArrowCol The new style arrow line color, default('black').
#' @param speArrowAngle The new style arrow angle, default(30).
#' @param speArrowLen The new style arrow length, default(0.1).
#' @param speArrowType The new style arrow type, default('closed').
#' @param textLabel The text label aesthetic mappings, default('transcript_id').
#' @param textLabelSize The text label size, default(5).
#' @param textLabelColor The text label color, default('black').
#' @param base_size Theme basesize, default(14).
#' @param marginX Plot left and right margins, default(0.2).
#' @param marginY Plot top and bottomn margins, default(0.2).
#' @param aspect.ratio Plot ratio, default(NULL).
#' @param facetByGene Whether facet by gene to plot, this useful for your genes which are far away from each other or not located on the same chromosome, default(FALSE).
#' @param ncolGene The column numbers to plot, default(NULL).
#' @param scales Facet plot scales, same as "facet_wrap" function, default('free').
#' @param strip.position Facet plot strip.position, same as "facet_wrap" function, default('top').
#' @param forcePosRel Whether force the genome coordinate to relative position to transcript start/end position, default('FALSE').
#' @param panel.spacing Facet plot panel space, default(0.3).
#' @param revNegStrand Whether reverse the negtive strand when set "forcePosRel=TRUE", default('FALSE').
#'
#' @param xAxis.info Whether retain X axis ticks and text, default(TRUE).
#' @param reverse.y whether reverse the Y axis, default(FALSE).
#' @param text.pos the label position(left/right), default(middle).
#' @param selecType choose the representative transcript to show("lt(longest transcript)" or "lcds(longest CDS)"), default(NULL).
#' @param topN the top number representative transcript to be shown, default(1).
#' @param show.legend whether show color legend, default(FALSE).
#'
#' @import tidyverse
#' @import cowplot
#' @import stats
#'
#' @return A ggplot object.
#'
#' @export
#' @examples
#' ##############################################################
#' # test function
#'
#' ########################################################
#' # load data
#' data(gtf)
#'
#' # non-coding gene
#' trancriptVis(gtfFile = gtf,
#'              gene = 'Xist')
#'
#' # coding gene
#' trancriptVis(gtfFile = gtf,
#'              gene = 'Nanog')
#'
#' # change fill color
#' trancriptVis(gtfFile = gtf,
#'              gene = 'Nanog',
#'              exonFill = '#CCFF00')
#'
#' # change inrton line size
#' trancriptVis(gtfFile = gtf,
#'              gene = 'Nanog',
#'              intronSize = 1)
#'
#' # change label size,color and position
#' trancriptVis(gtfFile = gtf,
#'              gene = 'Nanog',
#'              textLabelSize = 4,
#'              textLabelColor = 'red',
#'              relTextDist = 0)
#'
#' # aes by gene name
#' trancriptVis(gtfFile = gtf,
#'              gene = 'Nanog',
#'              textLabel = 'gene_name')
#'
#' # color aes by transcript
#' trancriptVis(gtfFile = gtf,
#'              gene = 'Tpx2',
#'              exonColorBy = 'transcript_id')
#'
#' # change arrow color and type
#' trancriptVis(gtfFile = gtf,
#'              gene = 'Nanog',
#'              arrowCol = 'orange',
#'              arrowType = 'closed')
#'
#' # no intron gene and add arrow color
#' # change arrow color and type
#' trancriptVis(gtfFile = gtf,
#'              gene = 'Jun',
#'              textLabel = 'gene_name',
#'              arrowCol = 'white',
#'              arrowType = 'closed') +
#'   theme_void()
#'
#' # add arrow breaks
#' trancriptVis(gtfFile = gtf,
#'              gene = 'Nanog',
#'              arrowCol = 'orange',
#'              arrowType = 'closed',
#'              arrowBreak = 0.1)
#'
#' # draw specific transcript
#' p1 <- trancriptVis(gtfFile = gtf,
#'                    gene = 'Commd7')
#'
#' p2 <- trancriptVis(gtfFile = gtf,
#'                    gene = 'Commd7',
#'                    myTranscript = c('ENSMUST00000071852','ENSMUST00000109782'))
#'
#' # combine
#' cowplot::plot_grid(p1,p2,ncol = 2,align = 'hv')

# global variables
globalVariables(c('end', 'gene_id', 'gene_name','seqnames',".",
                  'start', 'strand','transcript_id','transcript_name',
                  'type', 'vl_x1' ,'width', 'yPos','.env','cdslen','tlen'))

# use_package("tidyverse", type = "depends")

# define function
trancriptVis <- function(gtfFile = NULL,
                         gene = NULL,
                         selecType = NULL,
                         topN = 1,
                         myTranscript = NULL,
                         Chr = NULL,
                         posStart = NULL,
                         posEnd = NULL,
                         collapse = FALSE,
                         exonWidth = 0.3,
                         relTextDist = 0.3,
                         intronSize = 0.5,
                         arrowBreak = 0.15,
                         exonColorBy = NULL,
                         show.legend = FALSE,
                         exonFill = "#333399",
                         circle = FALSE,
                         cicStart = pi,
                         circSegCol = '#333399',
                         text_only = FALSE,
                         ylimLow = -10,
                         openAngle = 0.5,
                         arrowCol = '#333399',
                         arrowAngle = 30,
                         arrowLength = 0.1,
                         arrowType = 'open',
                         addNormalArrow = TRUE,
                         newStyleArrow = FALSE,
                         absSpecArrowLen = FALSE,
                         speArrowRelPos = 0,
                         speArrowRelLen = 0.05,
                         speArrowStart = -0.15,
                         speArrowRelHigh = 2,
                         speArrowLineSize = 0.5,
                         speArrowCol = 'black',
                         speArrowAngle = 30,
                         speArrowLen = 0.1,
                         speArrowType = "closed",
                         textLabel = 'transcript_id',
                         text.pos = "middle",
                         textLabelSize = 5,
                         textLabelColor = 'black',
                         base_size = 14,
                         marginX = 10,
                         marginY = 10,
                         aspect.ratio = NULL,
                         facetByGene = FALSE,
                         ncolGene = NULL,
                         scales = 'free',
                         strip.position = 'top',
                         forcePosRel = FALSE,
                         panel.spacing = 0.3,
                         revNegStrand = FALSE,
                         xAxis.info = TRUE,
                         reverse.y = FALSE){
  ##############################################################################
  # test whether with a given specific gene or region

  # select columns
  # if("transcript_name" %in% colnames(gtfFile)){
  #   sln <- c('seqnames','start','end','width','strand','type','gene_id','gene_name','transcript_id','transcript_name')
  # }else{
  #   sln <- c('seqnames','start','end','width','strand','type','gene_id','gene_name','transcript_id')
  # }

  # load GTF file
  if(is.character(gtfFile)){
    gtfFile <- rtracklayer::import(gtfFile,format = "gtf") %>%
      data.frame()
  }else{
    gtfFile <- gtfFile
  }

  # check colnames
  if("transcript_name" %in% colnames(gtfFile)){
    if("gene_name" %in% colnames(gtfFile)){
      sln <- c('seqnames','start','end','width','strand','type','gene_id','gene_name','transcript_id','transcript_name')
    }else{
      sln <- c('seqnames','start','end','width','strand','type','gene_id','transcript_id','transcript_name')
    }
  }else{
    if("gene_name" %in% colnames(gtfFile)){
      sln <- c('seqnames','start','end','width','strand','type','gene_id','gene_name','transcript_id')
    }else{
      sln <- c('seqnames','start','end','width','strand','type','gene_id','transcript_id')
    }
  }

  # select data
  if(is.null(gene)){
    # filter gene by region
    myGene <- gtfFile %>%
      dplyr::filter(seqnames == Chr & start >= posStart & end <= posEnd) %>%
      dplyr::filter(type != 'gene') %>%
      dplyr::select(sln)
  }else{
    # use gene_id stands for gene_name
    if("gene_name" %in% colnames(gtfFile)){
      # filter gene by gene name
      myGene <- gtfFile %>%
        dplyr::filter(gene_name %in% .env$gene) %>%
        dplyr::filter(type != 'gene') %>%
        dplyr::select(sln)
    }else{
      # filter gene by gene id
      myGene <- gtfFile %>%
        dplyr::select(sln) %>%
        dplyr::mutate(gene_name = gene_id) %>%
        dplyr::filter(gene_name %in% .env$gene) %>%
        dplyr::filter(type != 'gene')
    }
  }

  ##############################################################################
  # whether plot specific transcript
  if(is.null(myTranscript)){
    myData <- myGene
  }else{
    myData <- myGene %>%
      dplyr::filter(transcript_id %in% myTranscript)
  }

  # whether type column contains "transcript" info
  typeInfo <- unique(myData$type)
  if("transcript" %in% typeInfo){
    myData <- myData
  }else{
    purrr::map_df(unique(myData$transcript_id),function(x){
      tmp <- myData %>% dplyr::filter(transcript_id == x)
      tinfo <- tmp[1,]

      # assign start and end
      tinfo <- tinfo %>%
        dplyr::mutate(start = min(tmp$start),
                      end = max(tmp$end),
                      width = abs(max(tmp$end) - min(tmp$start)) + 1,
                      type = "transcript")

      # combine
      tData <- rbind(tinfo,tmp)
      return(tData)
    }) -> myData
  }

  ##############################################################################
  # select representive transcript
  ##############################################################################
  # function
  filterRepTrans <- function(data,selecType,topN = 1){
    tg <- data

    # calculate transcript/CDS length
    purrr::map_df(unique(tg$transcript_id),function(x){
      tmp <- tg %>% dplyr::filter(transcript_id == x & type == "exon")
      tmp2 <- tg %>% dplyr::filter(transcript_id == x & type == "CDS")
      tranLength <- data.frame(tid = x,tlen = sum(tmp$width),cdslen = sum(tmp2$width))
      return(tranLength)
    }) -> lenInfo

    # choose type
    if(selecType == "lt"){
      rankTran <- lenInfo %>%
        dplyr::arrange(dplyr::desc(tlen),dplyr::desc(cdslen)) %>%
        dplyr::slice_head(n = topN)
      return(rankTran$tid)
    }else if(selecType == "lcds"){
      rankTran <- lenInfo %>%
        dplyr::arrange(dplyr::desc(cdslen),dplyr::desc(tlen)) %>%
        dplyr::slice_head(n = topN)
      return(rankTran$tid)
    }else{
      print("Please choose 'lt' or 'lcds'!")
    }
  }

  # get rep trans
  if(!is.null(selecType)){
    purrr::map_df(unique(myData$gene_id),function(x){
      tmp <- myData %>% dplyr::filter(gene_id == x)
      res <- tmp %>% dplyr::filter(transcript_id %in% filterRepTrans(tmp,selecType,topN))
    }) -> myData

  }else{
    myData <- myData
  }

  ##############################################################################

  ##############################################################################
  # get gene id
  gid <- unique(myData$gene_id)

  ##############################################################################
  # add y axis position
  if(collapse == FALSE){
    # expand gene
    purrr::map_df(1:length(gid),function(x){
      tmp <- myData %>%
        dplyr::filter(gene_id == gid[x])
      # loop for tid
      trans_tmp <- tmp %>% dplyr::filter(type == 'transcript') %>%
        dplyr::arrange(width)

      # whether has transcript in gene info
      if(nrow(trans_tmp) > 0){
        tid <- unique(trans_tmp$transcript_id)

        # assign y position
        purrr::map_df(1:length(tid),function(x){
          tmp1 <- myData %>%
            dplyr::filter(transcript_id == tid[x]) %>%
            dplyr::mutate(yPos = x)
          tmp1$ymin <- ifelse(tmp1$type == 'CDS',
                              tmp1$yPos - exonWidth/2,
                              tmp1$yPos - exonWidth/6)
          tmp1$ymax <- ifelse(tmp1$type == 'CDS',
                              tmp1$yPos + exonWidth/2,
                              tmp1$yPos + exonWidth/6)
          return(tmp1)
        }) -> exon_ypos
      }
      # return(exon_ypos)
    }) -> mul_exon_ypos
  }else{
    # collapse gene
    mul_exon_ypos <- myData %>% dplyr::mutate(yPos = 1)
    mul_exon_ypos$ymin <- ifelse(mul_exon_ypos$type == 'CDS',
                                 mul_exon_ypos$yPos - exonWidth/2,
                                 mul_exon_ypos$yPos - exonWidth/6)
    mul_exon_ypos$ymax <- ifelse(mul_exon_ypos$type == 'CDS',
                                 mul_exon_ypos$yPos + exonWidth/2,
                                 mul_exon_ypos$yPos + exonWidth/6)
  }

  ##############################################################################
  # whether transform coordinate
  if(forcePosRel == TRUE){
    purrr::map_df(gene,function(x){
      tmp <- mul_exon_ypos %>%
        dplyr::filter(gene_name == x)
      purrr::map_df(unique(tmp$transcript_id),function(t){
        tmp1 <- tmp %>% dplyr::filter(transcript_id == t) %>%
          dplyr::arrange(start,end)

        # whether reverse negtive strand
        if(revNegStrand == FALSE){
          # start coord
          startPos <- min(tmp1$start)

          # add new pos
          tmp1 <- tmp1 %>% dplyr::mutate(start = start - startPos,
                                         end = end - startPos)
        }else{
          if(unique(tmp1$strand) == '-'){
            # end coord
            endPos <- max(tmp1$end)

            # add new pos
            tmp1 <- tmp1 %>% dplyr::mutate(start = endPos - start,
                                           end = endPos - end)
          }else{
            # start coord
            startPos <- min(tmp1$start)

            # add new pos
            tmp1 <- tmp1 %>% dplyr::mutate(start = start - startPos,
                                           end = end - startPos)
          }
        }
        return(tmp1)
      }) -> relPos_tmp
      return(relPos_tmp)
    }) -> exonNewPos
  }else{
    exonNewPos <- mul_exon_ypos
  }

  ##############################################################################
  # extarct data
  exon <- exonNewPos %>% dplyr::filter(type != 'transcript')
  trans <- exonNewPos %>% dplyr::filter(type == 'transcript')

  # add text x/y pos
  if(revNegStrand == FALSE){
    # define label position
    if(text.pos == "left"){
      text.pos.hjust = 1
      trans$textX <- trans$start
      trans$textY <- trans$yPos + relTextDist
    }else if(text.pos == "right"){
      text.pos.hjust = 0
      trans$textX <- trans$end
      trans$textY <- trans$yPos + relTextDist
    }else{
      text.pos.hjust = 0.5
      trans$textX <- ifelse(trans$strand == '+',
                            trans$start + trans$width/2,
                            trans$end - trans$width/2)
      trans$textY <- trans$yPos + relTextDist
    }
  }else{
    # define label position
    if(text.pos == "left"){
      text.pos.hjust = 1
      trans$textX <- trans$start
      trans$textY <- trans$yPos + relTextDist
    }else if(text.pos == "right"){
      text.pos.hjust = 0
      trans$textX <- trans$end
      trans$textY <- trans$yPos + relTextDist
    }else{
      text.pos.hjust = 0.5
      trans$textX <- (trans$start + trans$end)/2
      trans$textY <- trans$yPos + relTextDist
    }
  }

  ##############################################################################
  # whether add specific arrow
  if(newStyleArrow == FALSE){
    arrow_trans <- trans
  }else{
    # whether supply gene or coordinate
    if(is.null(gene)){
      gene <- unique(trans$gene_name)
    }

    # loop
    purrr::map_df(gene,function(gen){
      genTrans <- trans %>% dplyr::filter(gene_name == gen) %>%
        dplyr::arrange(dplyr::desc(width))

      # define longest transcript length
      longestWidth = genTrans$width[1]

      # add special arrow
      purrr::map_df(genTrans$transcript_id,function(x){
        tmp <- genTrans %>%
          dplyr::filter(transcript_id == x)

        # test strand
        if(absSpecArrowLen == TRUE){
          # add absolute specArrow
          if(tmp$strand == '+'){
            tmp <- tmp %>% dplyr::mutate(vl_x1 = start + speArrowRelPos*width,
                                         vl_x2 = vl_x1 + speArrowRelLen*longestWidth,
                                         vl_y1 = yPos + speArrowStart,
                                         vl_y2 = yPos + speArrowStart*speArrowRelHigh)
          }else if(tmp$strand == '-'){
            tmp <- tmp %>% dplyr::mutate(vl_x1 = end - speArrowRelPos*width,
                                         vl_x2 = vl_x1 - speArrowRelLen*longestWidth,
                                         vl_y1 = yPos + speArrowStart,
                                         vl_y2 = yPos + speArrowStart*speArrowRelHigh)
          }
        }else{
          # add relative specArrow
          if(tmp$strand == '+'){
            tmp <- tmp %>% dplyr::mutate(vl_x1 = start + speArrowRelPos*width,
                                         vl_x2 = vl_x1 + speArrowRelLen*width,
                                         vl_y1 = yPos + speArrowStart,
                                         vl_y2 = yPos + speArrowStart*speArrowRelHigh)
          }else if(tmp$strand == '-'){
            tmp <- tmp %>% dplyr::mutate(vl_x1 = end - speArrowRelPos*width,
                                         vl_x2 = vl_x1 - speArrowRelLen*width,
                                         vl_y1 = yPos + speArrowStart,
                                         vl_y2 = yPos + speArrowStart*speArrowRelHigh)
          }
        }
        return(tmp)
      }) -> arrow_trans1
      return(arrow_trans1)
    }) -> arrow_trans
  }

  # change reversed specArrow direction
  if(revNegStrand == FALSE){
    arrow_trans <- arrow_trans
  }else{
    arrow_trans$vl_x2 <- abs(arrow_trans$vl_x2)
  }

  # strand control arrow direction
  arrow_trans$ad <- ifelse(arrow_trans$strand == '+','last','first')

  # arrow breaks
  arrow_seq = c(0.05,seq(0,0.95,arrowBreak)[-1],1)

  ##############################################################################
  # first layer
  if(is.null(exonColorBy)){
    p1 <- ggplot2::ggplot(exon) +
      ggplot2::geom_rect(ggplot2::aes_(xmin = ~start,xmax = ~end,
                                       ymin = ~ymin,ymax = ~ymax),
                         fill = exonFill)
  }else{
    p1 <- ggplot2::ggplot(exon) +
      ggplot2::geom_rect(ggplot2::aes_(xmin = ~start,xmax = ~end,
                                       ymin = ~ymin,ymax = ~ymax,
                                       fill = ~get(exonColorBy)),
                         show.legend = show.legend)
  }

  ##############################################################################
  if(newStyleArrow == FALSE){
    p2 <- p1
  }else{
    p2 <- p1 +
      # add vertical line
      ggplot2::geom_segment(data = arrow_trans,
                            ggplot2::aes_string(x = "vl_x1",xend = "vl_x1",
                                                y = "vl_y1",yend = "vl_y2"),
                            color = speArrowCol,
                            size = speArrowLineSize) +
      # add horizotal line and arrow
      ggplot2::geom_segment(data = arrow_trans,
                            ggplot2::aes_string(x = "vl_x1",xend = "vl_x2",
                                                y = "vl_y2",yend = "vl_y2"),
                            color = speArrowCol,
                            size = speArrowLineSize,
                            arrow = ggplot2::arrow(angle = speArrowAngle,
                                                   length = ggplot2::unit(speArrowLen, "inches"),
                                                   ends = "last",
                                                   type = speArrowType))
  }


  ##############################################################################
  # whether facet by gene when gene far away each other or not on same chromosome
  if(facetByGene == FALSE){
    # whether draw ploar plot
    if(circle == FALSE){
      # add arrow and segment line with geom_arrowsegment
      if(addNormalArrow == TRUE){
        p3 <- p2 +
          ggarchery::geom_arrowsegment(data = arrow_trans,
                                       ggplot2::aes_(x = ~start,xend = ~end,
                                                     y = ~yPos,yend = ~yPos),
                                       color = arrowCol,
                                       size = intronSize,
                                       arrow_positions = arrow_seq,
                                       arrow_fills = rep(arrowCol,length(arrow_seq)),
                                       arrows = list(ggplot2::arrow(angle = arrowAngle,
                                                                    length = ggplot2::unit(arrowLength, "inches"),
                                                                    ends = arrow_trans$ad,
                                                                    type = arrowType)))
      }else{
        p3 <- p2 +
          ggplot2::geom_segment(data = arrow_trans,
                                ggplot2::aes_(x = ~start,xend = ~end,
                                              y = ~yPos,yend = ~yPos),
                                color = arrowCol,
                                size = intronSize)
      }
    }else{
      # add arrow and segment line geom_segment
      if(addNormalArrow == TRUE){
        p3 <- p2 +
          ggplot2::geom_segment(data = arrow_trans,
                                ggplot2::aes_(x = ~start,xend = ~end,
                                              y = ~yPos,yend = ~yPos),
                                color = circSegCol,
                                size = intronSize,
                                arrow = ggplot2::arrow(angle = arrowAngle,
                                                       length = ggplot2::unit(arrowLength, "inches"),
                                                       ends = arrow_trans$ad,
                                                       type = arrowType))
      }else{
        p3 <- p2 +
          ggplot2::geom_segment(data = arrow_trans,
                                ggplot2::aes_(x = ~start,xend = ~end,
                                              y = ~yPos,yend = ~yPos),
                                color = circSegCol,
                                size = intronSize)
      }
    }
  }else{
    # define columns to facet
    if(is.null(ncolGene)){
      ncol = length(gene)
    }else{
      ncol = ncolGene
    }
    # facet plot

    # whether draw ploar plot
    if(circle == FALSE){
      # add arrow and segment line with geom_arrowsegment
      if(addNormalArrow == TRUE){
        # loop add arrow
        for (g in gene) {
          tmpTrans <- arrow_trans %>% dplyr::filter(gene_name == g)
          p2 <- p2 +
            ggarchery::geom_arrowsegment(data = tmpTrans,
                                         ggplot2::aes_(x = ~start,xend = ~end,
                                                       y = ~yPos,yend = ~yPos),
                                         color = arrowCol,
                                         size = intronSize,
                                         arrow_positions = arrow_seq,
                                         arrow_fills = rep(arrowCol,length(arrow_seq)),
                                         arrows = list(ggplot2::arrow(angle = arrowAngle,
                                                                      length = ggplot2::unit(arrowLength, "inches"),
                                                                      ends = tmpTrans$ad,
                                                                      type = arrowType)))
        }
        p3_tmp <- p2
      }else{
        p3_tmp <- p2 +
          ggplot2::geom_segment(data = arrow_trans,
                                ggplot2::aes_(x = ~start,xend = ~end,
                                              y = ~yPos,yend = ~yPos),
                                color = arrowCol,
                                size = intronSize)
      }
    }else{
      # add arrow and segment line geom_segment
      if(addNormalArrow == TRUE){
        # loop add arrow
        for (g in gene) {
          tmpTrans <- arrow_trans %>% dplyr::filter(gene_name == g)
          p2 <- p2 +
            ggplot2::geom_segment(data = tmpTrans,
                                  ggplot2::aes_(x = ~start,xend = ~end,
                                                y = ~yPos,yend = ~yPos),
                                  color = circSegCol,
                                  size = intronSize,
                                  arrow = ggplot2::arrow(angle = arrowAngle,
                                                         length = ggplot2::unit(arrowLength, "inches"),
                                                         ends = tmpTrans$ad,
                                                         type = arrowType))
        }
        p3_tmp <- p2
      }else{
        p3_tmp <- p2 +
          ggplot2::geom_segment(data = arrow_trans,
                                ggplot2::aes_(x = ~start,xend = ~end,
                                              y = ~yPos,yend = ~yPos),
                                color = circSegCol,
                                size = intronSize)
      }
    }

    # facet
    p3 <- p3_tmp +
      ggplot2::facet_wrap(facets = "gene_name",
                          scales = scales,
                          ncol = ncol,
                          strip.position = strip.position)
  }

  ##############################################################################
  # add text label
  if(circle == FALSE){
    # geom_text
    p4 <- p3 +
      ggplot2::geom_text(data = arrow_trans,
                         ggplot2::aes_string(x = 'textX',y = 'textY',
                                             label = textLabel),
                         size = textLabelSize,
                         color = textLabelColor,
                         hjust = text.pos.hjust,
                         check_overlap = T) +
      ggplot2::theme_bw(base_size = base_size) +
      ggplot2::theme(panel.grid = ggplot2::element_blank(),
                     axis.ticks.y = ggplot2::element_blank(),
                     axis.text.y = ggplot2::element_blank(),
                     axis.title.y = ggplot2::element_blank(),
                     plot.margin = ggplot2::margin(t = marginY,r = marginX,b = marginY ,l = marginX)) +
      ggplot2::xlab('Positions on genome')
  }else{
    # geom_textpath
    p4 <- p3 +
      geomtextpath::geom_textpath(data = arrow_trans,
                                  ggplot2::aes_string(x = "textX",y = "textY",label = textLabel),
                                  size = textLabelSize,
                                  color = textLabelColor,
                                  hjust = text.pos.hjust,
                                  text_only = text_only) +
      ggplot2::coord_polar(theta = 'x',start = cicStart) +
      ggplot2::scale_y_continuous(limits = c(ylimLow,nrow(arrow_trans) + 1)) +
      ggplot2::scale_x_continuous(expand = c(0,openAngle*max(trans$width))) +
      ggplot2::theme_void()
  }

  ##############################################################################
  # add ratio
  if(is.null(aspect.ratio)){
    p5 <- p4
  }else{
    p5 <- p4 +
      ggplot2::theme(aspect.ratio = aspect.ratio)
  }

  ##############################################################################
  # facet background
  if(facetByGene == TRUE){
    p6 <- p5 +
      ggplot2::theme(strip.text.x = ggplot2::element_text(size = base_size + 2),
                     strip.background = ggplot2::element_rect(fill = 'grey90'),
                     panel.spacing = ggplot2::unit(panel.spacing,'cm'))
  }else{
    p6 <- p5
  }

  ##############################################################################
  # xlabel
  if(forcePosRel == TRUE){
    p7 <- p6 +
      ggplot2::xlab('Positions on genome')
  }else{
    p7 <- p6 +
      ggplot2::xlab('Position')
  }

  # X aixs text and ticks
  if(xAxis.info == F){
    p8 <- p7 +
      ggplot2::theme(axis.ticks.x = ggplot2::element_blank(),
                     axis.text.x = ggplot2::element_blank())
  }else{
    p8 <- p7
  }

  # whether reverse y axis
  if(reverse.y == TRUE){
    p9 <- p8 +
      ggplot2::scale_y_reverse()
  }else{
    p9 <- p8
  }
  return(p9)
}
