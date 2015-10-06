
filterEQTLs <- function(data_frame, gene_id_name_map, fdr_cutoff = 0.1){
  dat = dplyr::filter(data_frame, FDR < fdr_cutoff) %>% 
    dplyr::rename(gene_id = gene, snp_id = snps) %>%
    dplyr::group_by(gene_id) %>% 
    dplyr::arrange(pvalue) %>% 
    dplyr::filter(row_number() == 1) %>%
    dplyr::ungroup() %>% 
    dplyr::arrange(pvalue) %>%
    dplyr::left_join(gene_id_name_map, by = "gene_id")
  return(dat)
}