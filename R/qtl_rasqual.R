#' Save a list of matrices into a suitable format for RASQUAL
#'
#' Works with expression and covariates matrices.
#' 
#' @param data_list list of matrices
#' @param output_dir relative path to the output dir
#' @param file_suffix suffix added to each file after their name in the list.
#' @return None
#' @author Kaur Alasoo
#' @export 
saveRasqualMatrices <- function(data_list, output_dir, file_suffix = "expression"){
  #Save data for FastQTL to disk
  
  #Save each matrix as a separate  txt file
  for (sn in names(data_list)){
    file_path = file.path(output_dir, paste(sn,file_suffix, "txt", sep = "."))
    file_path_bin = file.path(output_dir, paste(sn,file_suffix, "bin", sep = "."))
    print(file_path)
    write.table(data_list[[sn]], file_path, quote = FALSE, row.names = TRUE, col.names = FALSE, sep = "\t")
    writeBin(as.double(c(t(data_list[[sn]]))), file_path_bin)
  }
}


rasqualSizeFactorsMatrix <- function(expression_list, factor_name){
  factor_vector = expression_list$norm_factor[,factor_name]
  size_matrix = matrix(rep(factor_vector, nrow(expression_list$counts)), nrow = nrow(expression_list$counts), byrow = TRUE)
  rownames(size_matrix) = rownames(expression_list$counts)
  return(size_matrix)
}

#' Count the number of feature SNPs and cis SNPs overlapping a set of peak calls
#'
#' @param peak_metadata Data frame with peak metadata (required columns: gene_id, chr, start, end)
#' @param snp_coords Data frame with SNP coordinates from a VCF file (required columns: chr, pos, snp_id)
#' @param cis_window Size of the cis window from both sides of the peak.
#'
#' @return Data frame with peak coordinates, cis region coordiantes as well as number of cis and feature snps.
#' @export
countSnpsOverlapingPeaks <- function(peak_metadata, snp_coords, cis_window = 500000){
  
  #Construct peak coords data frame
  peak_coords = dplyr::transmute(peak_metadata, gene_id, chromosome_name = chr, strand = 1, 
                                 exon_starts = start, exon_ends = end) %>%
    dplyr::mutate(range_start = pmax(0, exon_starts - cis_window), range_end = exon_ends + cis_window)
  
  #Construct GRanges objects
  peak_granges = GenomicRanges::GRanges(seqnames = peak_coords$chromosome_name, 
                         ranges = IRanges::IRanges(start = peak_coords$exon_starts, end = peak_coords$exon_ends))
  region_granges = GenomicRanges::GRanges(seqnames = peak_coords$chromosome_name, 
                           ranges = IRanges::IRanges(start = peak_coords$range_start, end = peak_coords$range_end))
  snp_granges = GenomicRanges::GRanges(seqnames = snp_coords$chr, ranges = IRanges::IRanges(start = snp_coords$pos, end = snp_coords$pos))
  
  #Count overlaps
  feature_snps = GenomicRanges::countOverlaps(peak_granges, snp_granges)
  cis_snps = GenomicRanges::countOverlaps(region_granges, snp_granges)
  
  new_peak_coords = dplyr::mutate(peak_coords, feature_snp_count = feature_snps, cis_snp_count = cis_snps)
  return(new_peak_coords)
}