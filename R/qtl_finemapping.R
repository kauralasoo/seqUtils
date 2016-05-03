
#' Add ATAC peak annotations to the eQTL credible set
#'
#' @param eqtl_credible_set 
#' @param atac_peak_metadata 
#' @param atac_tabix_path 
#'
#' @return Credible set with ATAC peak annotations
#' @export
annotateCredibleSet <- function(eqtl_credible_set, atac_peak_metadata, atac_tabix_path){
  
  #Construct Granges for the credible set SNPs
  variant_ranges = dplyr::transmute(eqtl_credible_set, seqnames = chr, start = pos, end = pos, strand = "+", snp_id) %>%
    dataFrameToGRanges()
  
  #Add overlapping peaks
  peak_finemapping = addOverlappingPeaks(eqtl_credible_set, atac_peak_metadata)
  
  #Find most associated peaks for each SNP
  associated_peaks = tabixFetchSNPs(variant_ranges, atac_tabix_path)
  if (!is.null(associated_peaks)){
    associated_peaks = associated_peaks %>% 
      dplyr::group_by(snp_id) %>% 
      dplyr::arrange(p_nominal) %>% 
      dplyr::filter(row_number() == 1) %>% 
      dplyr::ungroup() %>% 
      dplyr::transmute(assoc_peak_id = gene_id, snp_id, p_nominal_peak = p_nominal)
    
    #Add associated peaks
    peak_finemapping = dplyr::left_join(peak_finemapping, associated_peaks, by = "snp_id")
  }else{
    peak_finemapping = dplyr::mutate(peak_finemapping, assoc_peak_id = NA)
  }
  
  return(peak_finemapping)
}

#' Add overlapping ATAC peaks to a credible set of QTL variants
#'
#' @param qtl_credible_set Data frame with credible set of QTL variants 
#' (columns: snp_id, chr, pos)
#' @param atac_peak_metadata Data frame of start end coordinates of ATAC peaks 
#' (columns: gene_id, chr, start, end, strand)
#' @param extend Number of bp added to peak start and end coordinates to increase the chance of finding an overlap.
#'
#' @return Modifed qtl_credible_set data frame with added overlap_peak_id column.
#' @export
addOverlappingPeaks <- function(qtl_credible_set, atac_peak_metadata, extend = 0) {
  
  #Some basic assertions
  #Check that atac_peak_metadata has required columns
  assertthat::assert_that(assertthat::has_name(atac_peak_metadata, "gene_id"))
  assertthat::assert_that(assertthat::has_name(atac_peak_metadata, "chr"))
  assertthat::assert_that(assertthat::has_name(atac_peak_metadata, "start"))
  assertthat::assert_that(assertthat::has_name(atac_peak_metadata, "end"))
  assertthat::assert_that(assertthat::has_name(atac_peak_metadata, "strand"))
  
  #Check that qtl_credible_set has required columns
  assertthat::assert_that(assertthat::has_name(qtl_credible_set, "snp_id"))
  assertthat::assert_that(assertthat::has_name(qtl_credible_set, "chr"))
  assertthat::assert_that(assertthat::has_name(qtl_credible_set, "pos"))
  
  assertthat::assert_that(is.numeric(extend))
  
  #Extend peaks by some buffer sequence
  if(extend > 0){
    atac_peak_metadata = dplyr::mutate(atac_peak_metadata, start = start - extend, end = end + extend)
  }
  
  #construct Granges object with atac peak data
  peak_ranges = dataFrameToGRanges(dplyr::rename(atac_peak_metadata, seqnames = chr))
  
  #Construct Granges for the credible set SNPs
  variant_ranges = dplyr::transmute(qtl_credible_set, seqnames = chr, start = pos, end = pos, strand = "+", snp_id) %>%
    dataFrameToGRanges()
  
  #Find overlapping peaks for each variant
  olaps = GenomicRanges::findOverlaps(variant_ranges, peak_ranges, ignore.strand = TRUE)
  variants_overlaps = variant_ranges[GenomicRanges::queryHits(olaps)]
  peak_overlaps = peak_ranges[GenomicRanges::subjectHits(olaps),]
  overlapping_peaks = dplyr::data_frame(snp_id = variants_overlaps$snp_id, overlap_peak_id = peak_overlaps$gene_id)
  
  #Add overlapping peaks
  peak_finemapping = dplyr::left_join(qtl_credible_set, overlapping_peaks, by = "snp_id")
  return(peak_finemapping)
}
