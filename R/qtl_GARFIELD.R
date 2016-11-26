rasqualSummariesToGarfieldByChr <- function(selected_chr, rasqual_min_pvalues, summary_paths, GRCh37_variants, 
                                            gene_metadata, garfield_coords, cis_dist = 5e5, p_threshold = 1e-5){
  
  #Extract snp positions from variant data
  snp_pos_df = dplyr::select(GRCh37_variants, snp_id, pos)
  #Extract gene chromosomes from metadata
  gene_chr_map = dplyr::select(gene_metadata, gene_id, chr)
  
  #Extract QTLs per condition
  chr_qtl_list = purrr::map(rasqual_min_pvalues, ~dplyr::filter(., p_nominal < p_threshold) %>%
                              dplyr::left_join(gene_chr_map, by = "gene_id") %>%
                              dplyr::filter(chr == selected_chr))
  
  #Construct gene ranges for each condition
  gene_ranges_list = purrr::map(chr_qtl_list, 
                                ~rasqualTools::constructGeneRanges(.,gene_metadata, cis_dist))
  
  #Extract summary stats for each condition
  summaries = purrr::map2(gene_ranges_list, summary_paths, ~rasqualTools::tabixFetchGenes(.x, .y))
  
  #filter QTLs by p-value
  pvalue_hits = purrr::map(summaries, ~purrr::map_df(., ~dplyr::filter(., p_nominal < 1e-5)))
  
  #Identify unique variants and their positions in GRCh37
  unique_snps = purrr::map_df(pvalue_hits, identity) %>% dplyr::select(snp_id) %>% unique()
  subset_positions = dplyr::filter(snp_pos_df, snp_id %in% unique_snps$snp_id)
  
  #Add snp posiitions to pvalue hits
  esnps = purrr::map(pvalue_hits, ~dplyr::select(., snp_id) %>% 
                       unique() %>% 
                       dplyr::left_join(subset_positions, by = "snp_id") %>%
                       dplyr::mutate(esnp = 1) %>%
                       dplyr::select(pos, esnp))
  
  #Add eSNPs to GARFIELD coordinates
  esnp_df_list = purrr::map(esnps, ~dplyr::left_join(garfield_coords, ., by = "pos") %>% 
                              dplyr::mutate(esnp = ifelse(is.na(esnp), 0, esnp)) %>%
                              dplyr::select(esnp))
  
  #Merge GARFIELD fields
  eqtl_df = dplyr::bind_cols(esnp_df_list)
  result = dplyr::bind_cols(garfield_coords, eqtl_df)
  result$pos = paste0(result$pos, " ") #Add space to pos field
  
  return(result)
}
