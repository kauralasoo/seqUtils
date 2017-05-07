makeQTLCoveragePlot <- function(qtl_df, str1_df, str2_df, genotypes, exons, gene_metadata = NULL, ...){
  
  assertthat::assert_that(assertthat::has_name(qtl_df, "snp_id"))
  assertthat::assert_that(assertthat::has_name(qtl_df, "phenotype_id"))
  assertthat::assert_that(assertthat::has_name(qtl_df, "other_phenotype_id"))
  assertthat::assert_that(assertthat::has_name(qtl_df, "group_id"))
  
  print(qtl_df$phenotype_id)
  #Select transcripts
  if(!is.null(gene_metadata)){
    tx_names = dplyr::filter(gene_metadata, gene_id == qtl_df$group_id)$transcript_id
  } else{
    tx_names = c(qtl_df$phenotype_id, qtl_df$other_phenotype_id)
  }
  selected_exons = exons[tx_names]
  
  #Extract strand
  selected_strand = as.vector(strand(selected_exons[[1]]))[1]
  
  #Construct track data
  if(selected_strand == "-"){
    bigwig_meta = str1_df
  }else{
    bigwig_meta = str2_df
  }
  track_data = wiggleplotrGenotypeColourGroup(bigwig_meta, qtl_df$snp_id, genotypes, 1)
  
  #Make a coverage plot
  plot = wiggleplotr::plotCoverage(selected_exons, track_data = track_data, fill_palette = getGenotypePalette(), ...)
  return(plot)
}
