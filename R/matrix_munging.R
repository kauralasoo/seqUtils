#Set of functions to work with raw and processed data matrices

#' Normalise counts matrix using transcripts per million (TPM) normalisation
#'
#' @param counts_matrix Matrix of read counts (genes in rows, samples in columns)
#' @param lengths data frame with gene lengths (required columns: gene_id, length)
#' @param selected_genes Subset of genes used to to calculate library size (default: NULL = use all genes)
#' @param fragment_length Mean fragment length
#'
#' @return Read count matrix where entries have been TPM normalized.
#' @export
calculateTPM <- function(counts_matrix, lengths, selected_genes = NULL, fragment_length = 250){
  #Normalize gene expression using TPM method
  
  if (is.null(selected_genes)){
    selected_genes = rownames(counts_matrix)
  }
  
  #Add rownames to lengths df
  lengths = as.data.frame(lengths) #Make sure that lengths is a data frame and not a tbl_db, because tbl_df cannot have row names
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
#' @param return_type If return_type == "normalised" then return normalized expression matrix, otherwise return cqn object.
#' @return Quantile-normalized and GC-corrected matrix.
#' @author Kaur Alasoo
#' @export 
calculateCQN <- function(counts_matrix, gene_metadata, return_type = "normalised"){
  #Normalize read counts using the CQN method.
  expression_cqn = cqn(counts = counts_matrix[gene_metadata$gene_id,], 
                       x = gene_metadata$percentage_gc_content, 
                       lengths = gene_metadata$length, verbose = TRUE)
  #Choose return type
  if(return_type == "normalised"){
    expression_norm = expression_cqn$y + expression_cqn$offset
    return(expression_norm)
  }
  else{
    return(expression_cqn)
  }
}

calculateNormFactors <- function(counts_matrix, method = "RLE"){
  #Calculate norm factors for a counts matrix
  dge = edgeR::DGEList(counts = counts_matrix)
  dge = edgeR::calcNormFactors(dge, method = method)
  sample_info = dge$samples[,-1]
  colnames(sample_info) = c("library_size", "norm_factor")
  result = dplyr::mutate(sample_info, sample_id = rownames(sample_info)) %>%
    dplyr::select(sample_id, everything())
  return(result)
}

calculateMean <- function(matrix, design, factor, sample_id_col = "sample_id"){
  #Calculate the mean value in matrix over all possible factor values.
  
  #If the factor is not a factor then make it a factor.
  if(!is.factor(design[,factor])){
    design[,factor] = factor(design[,factor])
  }
  
  #Set sample_id column as rownames
  rownames(design) = design[,sample_id_col]
  factor = design[,factor]
  levs = levels(factor)
  result = c()
  for (lev in levs){
    filter = factor == lev
    samples = rownames(design[filter,])
    mat = matrix[,samples]
    mat = rowMeans(mat)
    result = cbind(result, mat)
  }
  colnames(result) = levs
  return(data.frame(result))
}

zScoreNormalize <- function(matrix){
  #Normalize expression matrix by z-score
  matrix = matrix - rowMeans(matrix)
  matrix = matrix / apply(matrix, 1, sd)
  return(matrix)
}

replaceNAsWithRowMeans <- function(matrix){
  #replace with row means
  na_pos = which(is.na(matrix), arr.ind = TRUE)
  matrix[na_pos] = rowMeans(matrix, na.rm=TRUE)[na_pos[,1]]
  
  #If there are addional NAs left (whole row NAs) then replace with 0
  matrix[is.na(matrix)] = 0
  return(matrix)
}

performPCA <- function(matrix, design, n_pcs = NULL, feature_id = "sample_id", column_prefix = "", ...){
  #Perform PCA of gene expression matrix add experimental design metadata to the results
  pca = prcomp(t(matrix), ...)
  if(is.null(n_pcs)){
    n_pcs = ncol(matrix)
  }
  pca_mat = as.data.frame(pca$x[,1:n_pcs])
  colnames(pca_mat) = paste0(column_prefix, colnames(pca_mat))
  pca_matrix = pca_mat %>% 
    dplyr::mutate(sample_id = rownames(pca$x)) %>%
    dplyr::rename_(.dots = setNames("sample_id", feature_id)) %>% #Hack to make renaming work
    dplyr::left_join(design, by = feature_id)
  #Calculate variance explained by each component
  var_exp = (pca$sdev^2) / sum(pca$sdev^2)
  return(list(pca_matrix = pca_matrix, pca_object = pca, var_exp = var_exp))
}

#' Force a vector of values into standard normal distribution
#'
#' @param x numeric vector with arbitrary distribution
#'
#' @return Vector with a standard normal distribution
#' @export
quantileNormaliseVector = function(x){
  qnorm(rank(x,ties.method = "random")/(length(x)+1))
}

quantileNormaliseMatrix <- function(matrix){
  quantile_matrix = matrix(0, nrow(matrix), ncol(matrix))
  for (i in seq_along(matrix[1,])){
    quantile_matrix[,i] = quantileNormaliseVector(matrix[,i])
  }
  #Add names
  rownames(quantile_matrix) = rownames(matrix)
  colnames(quantile_matrix) = colnames(matrix)
  return(quantile_matrix)
}

quantileNormaliseCols <- function(matrix,...){
  quantileNormaliseMatrix(matrix, ...)
}

quantileNormaliseRows <- function(matrix,...){
  t(quantileNormaliseMatrix(t(matrix), ...))
}