extractSubset <- function(sample_meta, data, old_column_names = "sample_id", new_column_names = "donor"){
  #Extract subset of columns from a data.frame and rename column names
  
  sample_meta = as.data.frame(sample_meta)               #Ensure that its a data.frame
  subset_data = data[,sample_meta[,old_column_names]]    #Keep only sample that are in the metadata column
  colnames(subset_data) = sample_meta[,new_column_names] #Rename the columns with donor id
  
  return(subset_data)
}

runMatrixEQTL <- function(exp_data, geno_data, snpspos, genepos, covariates = NULL){
  #Run matrixeQTL on a prepared data set
  
  #Perform some sanity checks
  if(!all(colnames(exp_data) == colnames(geno_data))){
    stop("Column names of expression and genotype data are not equal.")
  }
  
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
    if(!all(colnames(exp_data) == colnames(covariates))){
      stop("Column names of expression and covariates data are not equal.")
    }
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

savePEERData <- function(expression_list, outdir){
  #Save eQTL gene expression data to be used with PEER
  sample_names = names(expression_list)
  for (sample in sample_names){
    matrix = expression_list[[sample]]
    path = file.path(outdir, paste(sample, ".exprs.txt", sep = ""))
    print(path)
    write.table(t(matrix), path, row.names = FALSE, col.names = FALSE, sep = ",")
  }
}

importPEERFactors <- function(file_path, design_matrix, remove_mean = TRUE){
  #Import output from runPEER.py script back into R
  peer_factors = read.table(file_path, sep =",")
  if(remove_mean){ #Remove the first column that contains the mean
    peer_factors = peer_factors[,2:ncol(peer_factors)]
  }
  colnames(peer_factors) = paste("PEER_factor_", c(1:ncol(peer_factors)), sep = "")
  peer_factors = dplyr::mutate(peer_factors, sample_id = design_matrix$sample_id) %>%
    dplyr::select(sample_id, everything())
  return(peer_factors)
}
