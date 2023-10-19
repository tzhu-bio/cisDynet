#' A function help us look for which genes these SNPs might affect the expression of.
#'
#' @param snp The SNP coordinates file.
#' @param peak2gene Peak2Gene result obtained by **getPeak2Gene**.
#' @param corr_cutoff Pearson's correlation coefficient threshold for determining significance Peak2Gene links. Default: 0.4.
#'
#' @return
#' @export
#'
#' @examples  getSNP2Gene(snp="~/encode/CAT/test.snp", peak2gene="~/all_peak2gene.txt")
getSNP2Gene <- function(snp, peak2gene, corr_cutoff = 0.4){
  snp_bed <- read.table(snp, head = F)
  colnames(snp_bed) <- c("chrom","start")
  snp_bed$end <- snp_bed$start
  snp1 <- GenomicRanges::makeGRangesFromDataFrame(snp_bed)
  snp2 <- valr::gr_to_bed(snp1)
  p2g <- readRDS(peak2gene)
  p2g$chrom <- sapply(strsplit(p2g$Peak,":"), `[`, 1)
  p2g$bed <- sapply(strsplit(p2g$Peak,":"), `[`, 2)
  p2g1 <- p2g %>% tidyr::separate(bed, c("start", "end"), "-")
  p2g1 <- p2g1[,c("chrom", "start", "end", "Peak","Gene", "correlations","p.value","Type","TSS","Summit2TSS")]
  colnames(p2g1)[4:10] <- c("Peak","Gene","Correlation","P_value","Type","TSS","Distance")
  p2g2 <- p2g1[abs(p2g1$Correlation) >= corr_cutoff, c(1,2,3)]
  p2g2$start <- as.integer(p2g2$start) + 1
  b <- GenomicRanges::makeGRangesFromDataFrame(p2g2)
  c <- valr::gr_to_bed(b)
  res <- valr::bed_intersect(c, snp2)
  res1 <- as.data.frame(res[,c(1,2,3,7)])
  res1$Peak <- sprintf("%s:%s-%s",res1$chrom, res1$start.x,res1$end.x)
  final <- merge(res1, p2g1, by="Peak", all.x=T)[,c("start.y","Peak","Gene","Correlation","P_value","Type","TSS","Distance")]
  colnames(final)[1] <- c("SNP")
  final <- final[abs(final$Correlation) >= corr_cutoff,]
  return(final)
}
