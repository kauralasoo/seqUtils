filterDESeqResults <- function(results,gene_metadata, min_padj = 0.01, min_fc = 1, biotype_filter = NULL){
  #Add gene name to the DESeq result and filter out up and downregulated genes.
  
  #Construct a results table
  result_table = results %>% 
    data.frame() %>% 
    dplyr::mutate(gene_id = rownames(results)) %>% 
    tbl_df() %>% 
    dplyr::left_join(gene_metadata, by = "gene_id") %>% 
    dplyr::arrange(padj)
  
  #Find up and down-regulated genes
  up_genes = dplyr::filter(result_table, padj < min_padj, log2FoldChange > min_fc) %>% 
    arrange(-log2FoldChange)
  down_genes = dplyr::filter(result_table, padj < min_padj, log2FoldChange < -min_fc) %>% 
    arrange(log2FoldChange)
  
  #Filter up and down-regulated genes by biotype
  if(!is.null(biotype_filter)){
    up_genes = dplyr::filter(up_genes, gene_biotype == biotype_filter)
    down_genes = dplyr::filter(down_genes, gene_biotype == biotype_filter)
  }
  return(list(up_genes = up_genes, down_genes = down_genes, results_table = result_table))
}

tidyDESeq <- function(result, gene_metadata){
  result_table = result %>% 
    as.data.frame() %>% 
    dplyr::mutate(gene_id = rownames(result)) %>% 
    tbl_df() %>% 
    dplyr::left_join(gene_metadata, by = "gene_id") %>% 
    dplyr::arrange(padj) %>%
    dplyr::select(gene_id, gene_name, everything())
  return(result_table)
}

tidyTopTable <- function(result){
  names = rownames(result)
  result = result %>% dplyr::tbl_df() %>%
    dplyr::mutate(gene_id = names) %>%
    dplyr::select(gene_id, everything())
  return(result)
}
