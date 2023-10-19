#' Draw venn plot of bed intersection.
#'
#' @param input_df A dataframe contains key and value(From get_combain_result())
#' @param bed_list A list contains records of all bed files(Read by function get_beds())
#' @param colors Color list. Default: Tableau_10
#' @param opacity Degree of opacity for the color(s) specified with colors (less opacity, more transparency). Default: 0.8
#' @param plotsize Default: 15
#'
#' @return Ggplot object of venn
#' @importFrom venn venn
#'
#' @description
#' This function facilitates the generation of Venn plot based on the intersection
#' results obtained from the get_beds() function and get_combine_result()
#'
#'
#' @export
#'
#' @examples
#' # Get beds using get_beds
#' beds <- get_beds(c('path1', 'path2', ... ))
#'
#' # Use get_combine_result() to get intersection result
#' # The processing time may vary depending on the number of sets involved
#' inputs <- get_combine_result(beds)
#'
#' draw_venn(inputs, beds)
plotVenn <- function(input_df, bed_list, colors=c('#4E79A7FF','#F28E2BFF','#E15759FF','#76B7B2FF','#59A14FFF','#EDC948FF','#B07AA1FF'), opacity = 0.8, plotsize = 15){
  set_size <- log2(dim(input_df)[1]+1)
  if(set_size>7){
    warning("The number of sets greater than 7, please use upset plot")
  }else{
    venn_input <- data.frame(key=sapply(input_df$key, get_venn_key, key_len=set_size), value=input_df$value)
    plot <- venn::venn(set_size,
                       counts=c(0,venn_input[order(venn_input$key),]$value),
                       zcolor=colors,
                       snames = unlist(lapply(bed_list, function(record) record$name)),
                       ggplot=T, opacity=opacity, plotsize=plotsize)
    return(plot)
  }

}


#' Draw upset plot of bed intersection.
#'
#' @param input_df A dataframe contains key and value(From get_combain_result())
#' @param bed_list A list contains records of all bed files(Read by function get_beds())
#' @param decreasing See UpSetR for details. Boolean value, Default: TRUE.
#' @param angles The angle of main bar text, Default: 45.
#' @param text_scale The scale of the text. Default: 0.65.
#' @param point_size The point size. Default: 2.8.
#' @param line_size The line width. Default: 1.
#'
#' @return Upset plot of intersection result
#' @importFrom UpSetR upset fromExpression
#' @importFrom paletteer paletteer_d
#'
#' @description
#' This function facilitates the generation of Upset plot based on the intersection
#' results obtained from the get_beds() function and get_combine_result()
#'
#'
#' @export
#'
#' @examples
#' # Get beds using get_beds
#' beds <- get_beds(c('path1', 'path2', ... ))
#'
#' # Use get_combine_result() to get intersection result
#' # The processing time may vary depending on the number of sets involved
#' inputs <- get_combine_result(beds)
#'
#' draw_upset(inputs, beds)
plotUpset <- function(input_df, bed_list, decreasing=T, angles=45, text_scale=0.65, point_size=2.8, line_size=1){
  beds_name <- unlist(lapply(bed_list, function(record) record$name))
  upset_input <- get_upset_input(input_df, beds_name)

  plot <- UpSetR::upset(UpSetR::fromExpression(upset_input),
                        nsets=length(beds_name),
                        order.by = "freq",
                        decreasing = decreasing,
                        number.angles = angles,
                        text.scale = text_scale,
                        point.size = point_size,
                        line.size = line_size,
                        main.bar.color= "#1170AA",matrix.color ="#FC7D0B",
                        sets.bar.color=paletteer::paletteer_d("ggthemes::Tableau_10")[c(1:length(beds_name))])
  return(plot)
}


#' Draw flower plot which represents the number of merged intervals.
#'
#' @param bedpaths The vector contains the paths of bed files, can be dir_path.
#' @param names Define the names of input files.
#' @param pattern The suffix of the bed files. Default: '.bed'.
#' @param start Start degree. Default: 90.
#' @param a Width of the petal. Default: 0.5.
#' @param b Length of the petal. Default: 2.
#' @param r Radius of the core. Default: 1.5.
#' @param width Width of the plot. Default: 10.
#' @param height Heigth of the plot. Default: 10.
#' @param ellipse_col Color of the petal.
#' @param circle_col Color of the core.
#' @param circle_text_cex Text size of label around the petal. Default: 1
#'
#' @return Ggplot object of venn
#' @importFrom valr bed_merge
#'
#' @description
#' Uses can use this function to get the number of merged intervals via Flower plot
#'
#'
#' @export
#'
#' @examples
#' draw_flower(c('path1','path2', ...))
plotFlower <- function(bedpaths, names=NULL, pattern='.bed', start=90, a=0.5, b=2, r=1.5, width=10, height=10, circle_text_cex=1,
                        circle_col=rgb(0, 162, 214, max = 255), ellipse_col=rgb(135, 206, 235, 150, max = 255)){
  beds <- get_beds(bedpaths, names, pattern)
  single_bed <- get_singlebed(bedpaths, pattern)

  flower_lables <- unlist(lapply(beds, function(record) record$name))
  flower_values <- unlist(lapply(beds, function(record) dim(record$interval)[1]))

  core_value <- dim(valr::bed_merge(single_bed))[1]

  plot <- flower_plot(flower_lables, flower_values,labels=paste0(core_value), start=start, a=a, b=b, r=r,
                      width=width, height=height, circle_text_cex=circle_text_cex, circle_col = circle_col, ellipse_col = ellipse_col)
  return(plot)
}


#' Calculate pairwised intersection of bed files.
#'
#' @param bed_list A list contains records of all bed files(Read by function get_beds()).
#' @param picture Reture plot object. Default: TRUE.
#' @param method Method of corrplot. Default: 'pie', can be "circle", "square", "ellipse", "number", "shade", "color", "pie".
#' @param color Color of glyphs.
#' @param tl_col Color of textg.
#'
#' @return A matrix/corrplot of pairwised intersection
#' @importFrom corrplot corrplot COL1
#' @importFrom valr bed_intersect
#'
#' @description
#' The purpose of the cal_pairwise() function is to perform a pairwise
#' intersection calculation among multiple BED files and then display with
#' corrplot
#'
#'
#' @export
#'
#' @examples
#' # Get beds using get_beds
#' beds <- get_beds(c('path1', 'path2', ... ))
#'
#' plot <- cal_pairwise(beds)
calPairwise <- function(bed_list, picture=T, method='pie', color="Blues",  tl_col="black"){
  n_beds <- length(bed_list)
  comb_df <- data.frame(combn(c(1:n_beds), 2))
  result_mat <- matrix(rep(1, n_beds**2), ncol=n_beds)

  for(comb in comb_df){
    a_i <- comb[1]
    b_i <- comb[2]

    bed_a <- bed_list[[a_i]]$interval
    bed_b <- bed_list[[b_i]]$interval

    size_a <- dim(bed_a)[1]
    size_b <- dim(bed_b)[1]

    result_mat[a_i, b_i] <- dim(get_unique(valr::bed_intersect(bed_a, bed_b)))[1]/size_a
    result_mat[b_i, a_i] <- dim(get_unique(valr::bed_intersect(bed_b, bed_a)))[1]/size_b
  }

  names <- unlist(lapply(bed_list, function(record) record$name))
  colnames(result_mat) <- names
  rownames(result_mat) <- names

  if(picture){
    corrplot::corrplot(result_mat, method=method, type='full', diag = FALSE, col.lim=c(0,1), col=corrplot::COL1(color),  tl.col= tl_col)
    return(invisible(NULL))
  }else{
    return(result_mat)
  }
}
