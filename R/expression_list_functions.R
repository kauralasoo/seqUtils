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

#' Convert named vector into a tidy data_frame.
#' 
#' @param named_vector Named vector.
#' @return data_frame with two columns: value, sample_id.
#' @author Kaur Alasoo
#' @export 
tidyVector <- function(named_vector){
  dplyr::data_frame(value = named_vector, sample_id = names(named_vector))
}

#' For a given gene_id, extract its expression from expression_matrix and join with metadata.
#' 
#' @param gene_id Gene id, corresponds to a row name of expression matrix.
#' @param expression_matrix Matrix of gene exrpression; rows - gene ids, cols - sample_ids.
#' @param sample_metadata Data frame with metadata for each sample.
#' @return data frame that has gene expression in value column and metadata in other columns.
#' @author Kaur Alasoo
#' @export 
constructGeneData <- function(gene_id, expression_matrix, sample_metadata){
  #Construct df of gene expression for lmer analysis
  gene_df = tidyVector(expression_matrix[gene_id,])
  model_data = dplyr::left_join(gene_df, sample_metadata, by = "sample_id") 
  return(model_data)
}