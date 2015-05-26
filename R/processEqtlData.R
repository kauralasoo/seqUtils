qtlProcessExpression <- function(sample_meta, exprs_cqn, remove_n_pcs = 4, new_column_names = "donor"){
  
  #Keep only sample that are in the metadata column
  condition_data = exprs_cqn[,sample_meta$sample_id]
  
  #Rename the columns with donor id
  colnames(condition_data) = sample_meta[,new_column_names]
  
  #Remove PCs
  if(remove_n_pcs > 0){
    pca = performPCA(condition_data, sample_meta)
    pca_explained = (pca$pca_object$x[,1:remove_n_pcs] %*% t(pca$pca_object$rotation[,1:remove_n_pcs])) %>%
      scale(center = -1 * pca$pca_object$center, scale = FALSE) %>% t()
    condition_data = condition_data - pca_explained
  }
  
  return(condition_data)
}

qtlProcessGenotypes <- function(sample_meta, genotypes, new_column_names = "donor"){
  gt = vcf_file$genotypes[,sample_meta$genotype_id]
  colnames(gt) = sample_meta[,new_column_names]
  return(gt)
}

runMatrixEQTL <- function(exp_data, geno_data, snpspos, genepos){
  #Run matrixeQTL on a prepared data set
  
  #Construct a SlicedData object of the expression data
  expression_sliced = SlicedData$new()
  expression_sliced$CreateFromMatrix(exp_data)
  expression_sliced$ResliceCombined(sliceSize = 2000)
  
  #Create a SlicedData obejct for the genotypes
  snps = SlicedData$new()
  snps$CreateFromMatrix(geno_data)
  snps$ResliceCombined()
  
  #RUN
  me = Matrix_eQTL_main(
    snps = snps,
    gene = expression_sliced,
    cvrt = SlicedData$new(),
    output_file_name = "",
    pvOutputThreshold = 0,  
    output_file_name.cis = NULL,
    pvOutputThreshold.cis = 1e-2,
    snpspos = snpspos,
    genepos = genepos,
    cisDist = 5e5,
    useModel = modelLINEAR, 
    errorCovariance = numeric(), 
    verbose = TRUE,
    pvalue.hist = TRUE,
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = FALSE);
  
  return(me)
}