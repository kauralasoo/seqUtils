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

extractTrQtlDataFromSE <- function(qtl_row, se, genotype_matrix, variant_information, 
                                   condition_levels = c("Ctrl", "AcLDL"), phenotype_name = "tpm_ratios"){
  #Extraxct sample metadata from SE
  sample_metadata = SummarizedExperiment::colData(se) %>% tbl_df2() %>%
    dplyr::mutate(condition_name = factor(condition_name, levels = condition_levels))
  
  #Extract phenotype matrix
  pheno_matrix = SummarizedExperiment::assays(se)[[phenotype_name]]
  
  #Reformat gene data
  transcript_metadata = tbl_df2(SummarizedExperiment::rowData(se)) %>% 
    dplyr::mutate(gene_id = transcript_id, gene_name = transcript_id)
  
  #Construct data
  data = constructQtlPlotDataFrame(qtl_row$phenotype_id, qtl_row$snp_id, 
                                   pheno_matrix, genotype_matrix, 
                                   sample_metadata, 
                                   transcript_metadata) %>%
    dplyr::left_join(constructGenotypeText(qtl_row$snp_id, variant_information), by = "genotype_value")
  
  return(data)
}


