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

performPCA <- function(matrix, design, ...){
  #Perform PCA of gene expression matrix add experimental design metadata to the results
  pca = prcomp(t(matrix), ...)
  pca_matrix = as.data.frame(pca$x) %>% 
    dplyr::mutate(sample_id = rownames(pca$x)) %>%
    dplyr::left_join(design, by = "sample_id")
  return(list(pca_matrix = pca_matrix, pca_object = pca))
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
