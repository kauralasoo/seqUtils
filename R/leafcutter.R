#Apply cluster-wise bonferroni correction to the fastqtl output
leafcutterBonferroniCorrection <- function(fastqtl_df, cluster_meta){
  result = dplyr::left_join(fastqtl_df, cluster_meta, by = "gene_id") %>% 
    dplyr::mutate(p_bonferroni = p_beta * cluster_size) %>%
    dplyr::mutate(p_bonferroni = pmin(p_bonferroni, 1)) %>%
    dplyr::mutate(p_fdr = p.adjust(p_bonferroni, method = "fdr")) %>%
    dplyr::group_by(cluster_id) %>% dplyr::arrange(cluster_id, p_fdr) %>% 
    dplyr::filter(row_number() == 1) %>% 
    dplyr::filter(p_fdr < 0.1) %>% 
    dplyr::ungroup() %>% 
    dplyr::arrange(p_fdr)
  return(result)
}