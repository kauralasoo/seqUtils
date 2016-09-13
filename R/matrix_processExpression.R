regressPrinicpalComponents <- function(data_matrix, n_pcs){
  #Regress out first n principal components from the data matrix
  if(n_pcs > 0){
    pca = prcomp(t(data_matrix))
    pca_explained = (pca$x[,1:n_pcs] %*% t(pca$rotation[,1:n_pcs])) %>%
      scale(center = -1 * pca$center, scale = FALSE) %>% t()
    result = data_matrix - pca_explained
  } else{
    result = data_matrix
  }
  return(result)
}

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


explainPEER <- function(peer_factors, covariates){
  #Identify which covariates are most correlated with hidden factors discovered by PEER
  dataset = cbind(covariates, peer_factors)
  covar_names = colnames(covariates)
  factor_names = colnames(peer_factors)
  results = c()
  for (covar in covar_names){
    for (factor in factor_names){
      formula = as.formula(paste(factor, covar, sep = " ~ "))
      model = lm(formula, data = dataset)
      model_summary = summary(model)
      result = c(covar, factor, model_summary$r.squared, model_summary$coefficients[2,4])
      results = rbind(results, result)
    }
  }
  colnames(results) = c("covariate", "peer_factor", "r_squared", "p_value")
  rownames(results) = c()
  results = as.data.frame(results, stringsAsFactors = FALSE)
  results$r_squared = as.numeric(results$r_squared)
  results$p_value = as.numeric(results$p_value)
  return(results)
}
