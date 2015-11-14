calculateTPM <- function(counts_matrix, lengths, selected_genes = NULL, fragment_length = 250){
  #Normalize gene expression using TPM method
  
  if (is.null(selected_genes)){
    selected_genes = rownames(counts_matrix)
  }
  
  #Add rownames to lengths df 
  rownames(lengths) = lengths$gene_id
  
  #Calculated scaling factor
  selected_counts = counts_matrix[intersect(rownames(counts_matrix), selected_genes),]
  selected_lengths = lengths[rownames(selected_counts),]$length
  scaling_factor = colSums((selected_counts*fragment_length)/selected_lengths)
  
  #Calculate TPM for each gene
  length_vector = lengths[rownames(counts_matrix),]$length
  tpm = t(t(counts_matrix * fragment_length * 1e6)/scaling_factor)
  tpm = tpm/length_vector
  
  return(tpm)
}

#' Normalize read counts matrix using the cqn method.
#'
#' Quantile normalize read counts while correcting for feature length and GC countent.
#' 
#' @param counts_matrix Matrix of read counts.
#' @param gene_metadata data.frame with at least three columns: gene_id, 
#' percentage_gc_content and length.
#' @return Quantile-normalized and GC-corrected matrix.
#' @author Kaur Alasoo
#' @export 
calculateCQN <- function(counts_matrix, gene_metadata){
  #Normalize read counts using the CQN method.
  expression_cqn = cqn(counts = counts_matrix[gene_metadata$gene_id,], 
                       x = gene_metadata$percentage_gc_content, 
                       lengths = gene_metadata$length, verbose = TRUE)
  expression_norm = expression_cqn$y + expression_cqn$offset
  return(expression_norm)
}

calculateNormFactors <- function(counts_matrix, method = "RLE", output = "rasqual"){
  #Calculate norm factors for a counts matrix
  dge = edgeR::DGEList(counts = counts_matrix)
  dge = edgeR::calcNormFactors(dge, method = method)
  sample_info = dge$samples[,-1]
  if (output == "rasqual"){
    size_matrix = matrix(rep(sample_info$norm.factors, nrow(counts_matrix)), 
                         nrow = nrow(counts_matrix), byrow = TRUE)
    rownames(size_matrix) = rownames(counts_matrix)
    return(size_matrix)
  }else{
    return(sample_info)
  }
}

filterExpressionDataset <- function(dataset, sample_ids = NULL, gene_ids = NULL){
  #Filter expression dataset by sample_ids or gene_ids
  if(!is.null(sample_ids)){
    dataset$design = dataset$design[sample_ids,]
    dataset$exprs_cqn = dataset$exprs_cqn[,sample_ids]
    dataset$exprs_counts = dataset$exprs_counts[,sample_ids]
  }
  if(!is.null(gene_ids)){
    dataset$exprs_cqn = dataset$exprs_cqn[gene_ids,]
    dataset$exprs_counts = dataset$exprs_counts[gene_ids,]
    dataset$gene_metadata = dplyr::filter(dataset$gene_metadata, gene_id %in% gene_ids)
  }
  return(dataset)
}

#' Convert named vector into a tidy data_frame.
#' 
#' @param named_vector Named vector.
#' @return data_frame with two columns: value, sample_id.
#' @author Kaur Alasoo
#' @export 
tidyVector <- function(named_vector){
  data_frame(value = named_vector, sample_id = names(named_vector))
}

#' For a given gene_id, extract its expression from expression_matrix and join with metdata.
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
