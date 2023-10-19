#' Get the co-accessible peaks in genome wide.
#'
#' @param quant_matrix   The normalized matrix obtained by **quantification**.
#' @param window   The window size to shift.
#' @param peak_anno The merged peaks annotation.
#'
#' @return
#' @export
#'
#' @examples   getCoaccessible(peak_anno="peak_anno.txt", quant_matrix="quant_mat.tsv", window=100000)
getCoaccessible <- function(peak_anno, quant_matrix, window=500000, sample_num = 100){
  checkGeAnno()
  anno <- read.table(peak_anno, head=T, sep="\t", stringsAsFactors=FALSE)
  summit <- anno[,-c(4)]
  colnames(summit) <- c('chr', 'bp1', 'bp2', 'summit')
  category <- anno[,-c(5)]
  colnames(category) <- c('chr', 'bp1', 'bp2', 'category')

  # processing the matrix
  mat1 <- read.table(quant_matrix, row.names=1, head=T, stringsAsFactors=FALSE)
  peak <- data.frame(peak = rownames(mat1))
  peak$chr <- sapply(strsplit(peak$peak,":"), `[`, 1 )
  peak$bed <- sapply(strsplit(peak$peak,":"), `[`, 2 )
  mat2 <- peak %>% tidyr::separate(bed, c("bp1", "bp2"), "-")
  mat2 <- mat2[,c(2:4)]
  mat<- cbind(mat2, mat1)
  rownames(mat) <- NULL
  rownames(mat) <- peaks <- sprintf("%s_%s_%s", mat$chr, mat$bp1, mat$bp2)

  ## regularized correlation
  genomic_coords <- data.frame(CATAnno$genome)
  colnames(genomic_coords) <- c("V1","V2")

  ## Split the chrom
  x <- unique(mat$chr)
  chrom_list <- gtools::mixedsort(x)

  coacc <- lapply(chrom_list, function(x){
    genomic_coords <- genomic_coords[genomic_coords$V1 == x, ]
    indata <- mat[mat$chr ==x, ]
    logfile(sprintf("Co-accessibility with %s.", x))
    res <- run_cicero(indata, window = window, genomic_coords, sample_num = sample_num)
    regcon <- res %>% dplyr::group_by(Peak1, Peak2) %>% dplyr::summarize(coaccess = mean(coaccess))
    return(regcon)
  })

  # cl <- parallel::makeCluster(getOption("cl.cores", cores_N))
  # coacc <- parallel::parLapply(cl, chrom_list, function(x){
  #   genomic_coords <- genomic_coords[genomic_coords$V1 == x, ]
  #   indata <- mat[mat$chr ==x, ]
  #   logfile(sprintf("Co-accessibility with %s.", x))
  #   res <- run_cicero(indata, window = window, genomic_coords)
  #   regcon <- res %>% dplyr::group_by(Peak1, Peak2) %>% dplyr::summarize(coaccess = mean(coaccess))
  #   return(regcon)
  # })
  # parallel::stopCluster(cl)



  all_coacc <- dplyr::bind_rows(coacc) %>% as.data.frame()

  logfile("Preparing the result!")
  all_coacc$Peak1 <- sub("_", ":", all_coacc$Peak1, fixed = TRUE)
  all_coacc$Peak2 <- sub("_", ":", all_coacc$Peak2, fixed = TRUE)
  all_coacc$Peak1 <- sub("_", "-", all_coacc$Peak1, fixed = TRUE)
  all_coacc$Peak2 <- sub("_", "-", all_coacc$Peak2, fixed = TRUE)
  colist <- all_coacc
  colist <- colist[summit[as.character(colist$Peak2),'summit'] > summit[as.character(colist$Peak1),'summit'], ]
  colist$dist <- summit[as.character(colist$Peak2),'summit'] - summit[as.character(colist$Peak1),'summit']
  colist$Center1 <- summit[as.character(colist$Peak1),'summit']
  colist$Center2 <- summit[as.character(colist$Peak2),'summit']
  colist$dist10k <- round(colist$dist / 10000)
  colist$dist2k <- round(colist$dist / 2000)*2000
  qlts <- c(-1, quantile(colist$coaccess, probs=1:2/3), 1)
  colist$colevel <- factor(cut(colist$coaccess, qlts), labels=c('Low','Medium','High'))
  colist$group <- ifelse(category[as.character(colist$Peak1),'category']=="Distal" & category[as.character(colist$Peak2),'category']=="Distal", "Enhancer-Enhancer",
                          ifelse(category[as.character(colist$Peak1),'category']=="Proximal" & category[as.character(colist$Peak2),'category']=="Proximal", "Promoter-Promoter", "Promoter-Enhancer"))
  return(colist)
}
