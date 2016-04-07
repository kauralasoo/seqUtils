
#' Add ATAC peak annotations to the eQTL credible set
#'
#' @param eqtl_credible_set 
#' @param atac_peak_metadata 
#' @param atac_tabix_path 
#'
#' @return Credible set with ATAC peak annotations
#' @export
annotateCredibleSet <- function(eqtl_credible_set, atac_peak_metadata, atac_tabix_path){
  
  #construct Granges object with atac peak data
  peak_ranges = dataFrameToGRanges(dplyr::rename(atac_peak_metadata, seqnames = chr))
  
  #Construct Granges for the credible set SNPs
  variant_ranges = dplyr::transmute(eqtl_credible_set, seqnames = chr, start = pos, end = pos, strand = "+", snp_id) %>%
    dataFrameToGRanges()
  
  #Find overlapping peaks for each variant
  olaps = GenomicRanges::findOverlaps(variant_ranges, peak_ranges, ignore.strand = TRUE)
  variants_overlaps = variant_ranges[GenomicRanges::queryHits(olaps)]
  peak_overlaps = peak_ranges[GenomicRanges::subjectHits(olaps),]
  overlapping_peaks = dplyr::data_frame(snp_id = variants_overlaps$snp_id, overlap_peak_id = peak_overlaps$gene_id)
  
  #Add overlapping peaks
  peak_finemapping = dplyr::left_join(eqtl_credible_set, overlapping_peaks, by = "snp_id")
  
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

