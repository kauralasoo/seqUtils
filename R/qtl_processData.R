qtlProcessExpression <- function(sample_meta, exprs_cqn, new_column_names = "donor"){
  
  #Keep only sample that are in the metadata column
  condition_data = exprs_cqn[,sample_meta$sample_id]
  
  #Rename the columns with donor id
  sample_meta = as.data.frame(sample_meta) #Ensure that its a data.frame
  colnames(condition_data) = sample_meta[,new_column_names]
  
  return(condition_data)
}

qtlProcessGenotypes <- function(sample_meta, genotypes, new_column_names = "donor"){
  gt = genotypes[,sample_meta$genotype_id]
  sample_meta = as.data.frame(sample_meta) #Ensure that its a data.frame
  colnames(gt) = sample_meta[,new_column_names]
  return(gt)
}

runMatrixEQTL <- function(exp_data, geno_data, snpspos, genepos, covariates = NULL){
  #Run matrixeQTL on a prepared data set
  
  #Construct a SlicedData object of the expression data
  expression_sliced = SlicedData$new()
  expression_sliced$CreateFromMatrix(exp_data)
  expression_sliced$ResliceCombined(sliceSize = 2000)
  
  #Create a SlicedData obejct for the genotypes
  snps = SlicedData$new()
  snps$CreateFromMatrix(geno_data)
  snps$ResliceCombined()
  
  #Add covariates
  cvrt = SlicedData$new()
  if (!is.null(covariates)){
    cvrt$CreateFromMatrix(covariates)
    cvrt$ResliceCombined()
  }
  
  #RUN
  me = Matrix_eQTL_main(
    snps = snps,
    gene = expression_sliced,
    cvrt = cvrt,
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
    pvalue.hist = "qqplot",
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = FALSE);
  
  return(me)
}