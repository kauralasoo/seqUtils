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

replaceNAsWithRowMeans <- function(matrix){
  #replace with row means
  na_pos = which(is.na(matrix), arr.ind = TRUE)
  matrix[na_pos] = rowMeans(matrix, na.rm=TRUE)[na_pos[,1]]
  
  #If there are addional NAs left (whole row NAs) then replace with 0
  matrix[is.na(matrix)] = 0
  return(matrix)
}
