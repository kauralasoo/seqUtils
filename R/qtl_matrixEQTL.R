#' A wrapper around MatrixeQTL
#'
#' Prepares SlicedData objects of expression data, genotypes and covariates 
#' and then runs MatrixeQTL on the data.
#' 
#' @param exp_data Matrix of gene expression data (genes in rows, samples in columns).
#' @param geno_data Matrix of genotype data (SNPs in rows, samples in columns).
#' @param snpspos Matrix of SNP coordinates (columns: snpid,chr,pos).
#' @param genepos Matrix of gene coordinates (columns: geneid,chr,left,right).
#' @param covariates Matrix of covariates (samples in columns).
#' @param cisDist cis distance from the gene.
#' @param pvOutputThreshold Maximum p-value to report. Smaller values make the code quicker.
#' @param permute Permute genotype labels before qtl mapping.
#' @param model Specifies which MatrixEQTL model to use. Options: modelLINEAR, modelLINEAR_CROSS and modelANOVA.
#' @return MatrixeQTL result object.
#' @author Kaur Alasoo
#' @export 
runMatrixEQTL <- function(exp_data, geno_data, snpspos, genepos, covariates = NULL, 
                          cisDist = 5e5, pvOutputThreshold = 1e-2, permute = FALSE, model = modelLINEAR){
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
  if(permute == TRUE){
    #Permute column labels of the genotype data
    genotype_labels = colnames(geno_data)
    geno_data = geno_data[,sample(length(genotype_labels), length(genotype_labels))]
    colnames(geno_data) = genotype_labels
  }
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
    pvOutputThreshold.cis = pvOutputThreshold,
    snpspos = snpspos,
    genepos = genepos,
    cisDist = cisDist,
    useModel = model, 
    errorCovariance = numeric(), 
    verbose = TRUE,
    pvalue.hist = "qqplot",
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = FALSE);
  
  return(me)
}

#' Convert gene_metadata table into a gene positions data frame suitable for MatrixEQTL.
#'
#' @param gene_metadata Data frame with gene metadata. Required columns: gene_id, chr, start, end)
#'
#' @return Gene postion data frame for matrixEQTL.
#' @export
constructMatrixEQTLGenePos <- function(gene_metadata){
  res = dplyr::transmute(gene_metadata, geneid = gene_id, chr, left = start, right = end)
  return(res)
}

#' Convert sample_metadata into a covariates matrix suitable for matrixEQTL
#'
#' @param sample_metadata Sample metadata matrix
#' @param covariate_names Vector of columns names that will be used as covariates.
#' @param id_column Name of the column that will be used as column names for the covariates matrix.
#' Default sample id.
#'
#' @return Covariates matrix suitable for matrixEQTL.
#' @export
constructMatrixEQTLCovariates <- function(sample_metadata, covariate_names, id_column = "sample_id"){
  cov_matrix = t(sample_metadata[,covariate_names])
  col_ids = as.data.frame(sample_metadata)[,id_column]
  colnames(cov_matrix) = col_ids
  return(cov_matrix)
}

matrixEQTLExtractCisQTLs <- function(matrixeqtl_object, snp_positions){
  renamed_data = dplyr::transmute(matrixeqtl_object$cis$eqtls, gene_id = as.vector(gene), snp_id = as.vector(snps), 
                                  statistic, p_nominal = pvalue, beta) %>% 
    tbl_df()
  selected_snps = dplyr::filter(snp_positions, snpid %in% renamed_data$snp_id) %>% 
    dplyr::rename(snp_id = snpid)
  joint_result = dplyr::left_join(renamed_data, selected_snps, by = "snp_id")
  return(joint_result)
}
