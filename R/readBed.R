#' Parse bedpaths as file based path.
#'
#' @param bedpaths The vector contains the paths of bed files, can be dir_path.
#' @param pattern The suffix of the bed files. Default: '.bed'.
#'
#' @return A vector which contains parsed file based paths
#' @export
#'
#' @examples
#' parse_bedpath(c('./Encode_bed', './peak/xxx.bed','./peak/zzz.bed'))
parse_bedpath <- function(bedpaths, pattern='.bed'){
  result_path <- c()
  for(path in bedpaths){
    if(dir.exists(path)){
      beds_path <- list.files(path, pattern = pattern)
      if(length(beds_path)>0)
        result_path <- c(result_path, paste(path, beds_path, sep='/'))
    }else{
      if(file.exists(path)){
        result_path <- c(result_path, path)
      }
    }
  }
  return(result_path)
}


#' Read beds into list, can be arbitary dir_path.
#'
#' @param bedpaths The vector contains the paths of bed files, can be dir_path.
#' @param names Define the names of input files.
#' @param pattern The suffix of the bed files. Default: '.bed'.
#'
#' @return A list which can be used for calculating intersection of each beds
#' @export
#'
#' @importFrom valr read_bed
#'
#' @description
#' Users are required to utilize this function in order to convert BED format
#' files into a list, which is an essential prerequisite for subsequent analyses
#'
#' @examples
#' get_beds(c('./Encode_bed/ccc.bed', './peak/xxx.bed','./peak/zzz.bed'), names=c('A1','A2','A3'))
get_beds <- function(bedpaths, names=NULL, pattern='.bed'){
  bedpaths <- parse_bedpath(bedpaths, pattern=pattern)
  rec <- list()
  if(!is.null(names) && length(names) != length(bedpaths)){
    warning("Assigned names can not match the number of bed files, use file name...")
    names <- NULL
  }
  for(i in c(1:length(bedpaths))){
    current_path <- bedpaths[i]
    name <- ifelse(is.null(names), gsub("_peaks_unique.narrowPeak.bed","",tail(strsplit(current_path,"/")[[1]],1)), names[i])
    rec[[i]] <- list(name=name, interval=valr::read_bed(current_path))
    rec[[i]]$interval <- dplyr::mutate(rec[[i]]$interval, len=end-start, ident=as.character(i))
    rec[[i]]$res <- dplyr::mutate(rec[[i]]$interval[,1:3], len=end-start, ident=as.character(i))
  }
  return(rec)
}


#' Read beds into a single tibble.
#'
#' @param bedpaths The vector contains the paths of bed files, can be dir_path.
#' @param sort A boolean value. Default: TRUE.
#' @param pattern The suffix of the bed files. default: '.bed'
#'
#' @return A tibble contains all the beds interval records
#' @importFrom dplyr arrange
#'
#' @export
#'
#' @examples
#' get_singlebed(c('./Encode_bed/ccc.bed', './peak/xxx.bed','./peak/zzz.bed'))
get_singlebed <- function(bedpaths, sort = TRUE, pattern='.bed'){
  bedpaths <- parse_bedpath(bedpaths, pattern=pattern)
  for(i in c(1:length(bedpaths))){
    if(i==1){
      df <- valr::read_bed(bedpaths[i])[,1:3]
    }else{
      df <- rbind(df, valr::read_bed(bedpaths[i])[,1:3])
    }
  }

  if(sort==TRUE){
    df <- dplyr::arrange(df, chrom, start)
  }
  return(df)
}
