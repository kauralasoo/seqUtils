extractConditionFromExpressionList <- function(expression_list, cond_name){
  
  new_sample_meta = dplyr::filter(expression_list$sample_metadata, condition_name == cond_name)
  new_counts = expression_list$counts[,new_sample_meta$sample_id]
  new_tpm = expression_list$tpm[,new_sample_meta$sample_id]
  new_cqn = expression_list$cqn[,new_sample_meta$sample_id]
  new_norm_factors = dplyr::semi_join(expression_list$norm_factors, new_sample_meta, by = "sample_id")
  
  #If covariates are present then filter those as well
  if(!is.null(expression_list$covariates)){
    new_covariates = expression_list$covariates[, new_sample_meta$sample_id]
  } else {
    new_covariates = NULL
  }
  
  #Combine everything into a list
  results_list = list(
    counts = new_counts,
    cqn = new_cqn,
    tpm = new_tpm,
    covariates = new_covariates,
    norm_factors = new_norm_factors,
    sample_metadata = new_sample_meta,
    gene_metadata = expression_list$gene_metadata)
}

extractGenesFromExpressionList <- function(expression_list, gene_ids){
  new_counts = expression_list$counts[gene_ids,]
  new_cqn = expression_list$cqn[gene_ids,]
  new_tpm = expression_list$tpm[gene_ids,]
  new_gene_metadata = dplyr::filter(expression_list$gene_metadata, gene_id %in% gene_ids)
  
  result_list = list(
    counts = new_counts,
    cqn = new_cqn,
    tpm = new_tpm,
    covariates = expression_list$covariates,
    norm_factors = expression_list$norm_factors,
    sample_metadata = expression_list$sample_metadata,
    gene_metadata = new_gene_metadata
  )
}

renameMatrixColumnsInExpressionList <- function(expression_list, old_column_names, new_column_names){
  
  #Change colunn names
  new_counts = extractSubset(expression_list$sample_metadata, expression_list$counts, old_column_names, new_column_names)
  new_tpm = extractSubset(expression_list$sample_metadata, expression_list$tpm, old_column_names, new_column_names)
  new_cqn = extractSubset(expression_list$sample_metadata, expression_list$cqn, old_column_names, new_column_names)
  
  #Update result list
  result_list = expression_list
  result_list$counts = new_counts
  result_list$tpm = new_tpm
  result_list$cqn = new_cqn
  
  #If covariates are present then filter those as well
  if(!is.null(expression_list$covariates)){
    new_covariates = extractSubset(expression_list$sample_metadata, expression_list$covariates, old_column_names, new_column_names)
    result_list$covariates = new_covariates
  }
  
  return(result_list)
}