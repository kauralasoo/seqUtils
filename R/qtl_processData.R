extractSubset <- function(sample_meta, data, old_column_names = "sample_id", new_column_names = "donor"){
  #Extract subset of columns from a data.frame and rename column names
  
  sample_meta = as.data.frame(sample_meta)               #Ensure that its a data.frame
  subset_data = data[,sample_meta[,old_column_names]]    #Keep only sample that are in the metadata column
  colnames(subset_data) = sample_meta[,new_column_names] #Rename the columns with donor id
  
  return(subset_data)
}

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

#' Save list of expression matrices in a format suitable for runPEER.py
#'
#' @param expression_list list of expression matrices
#' @param outdir output directory
#' @param file_suffix suffix of the output file
#'
#' @export
#'
savePEERData <- function(expression_list, outdir, file_suffix = "exprs"){
  #Save eQTL gene expression data to be used with PEER
  sample_names = names(expression_list)
  for (sample in sample_names){
    matrix = expression_list[[sample]]
    path = file.path(outdir, paste(sample, file_suffix ,"txt", sep = "."))
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

#' Calculate Pi1 replicability statistic between two sets of pvalues
#' 
#' Currently expects the following columns: gene_id, qvalue and p_beta.
#' 
#' @param table1 First table maximum p-values per feature.
#' @param table2 Second table of maximum p-values per feature.
#' @param qvalue_thresh qvalue threshold for table1.
#' @return None
#' @author Kaur Alasoo
#' @export 
calculatePi1 <- function(table1, table2, qvalue_thresh = 0.1){
  #Identify significant hits from first table
  table1_hits = dplyr::filter(table1, qvalue < qvalue_thresh)
  #Extract the same genes from the second table
  table2_hits = dplyr::semi_join(table2, table1_hits, by = "gene_id")
  #Estimate the proportion of replicated qtls
  pi1 = 1 - qvalue::qvalue(table2_hits$p_beta)$pi0
  return(pi1)
}

#' Calculate all pairwise Pi1 statistics for a list p-value tables.
#' 
#' Currently expects the following columns: gene_id, qvalue and p_beta.
#' 
#' @param qtl_list List of p-value tables
#' @param qvalue_thresh qvalue threshold for table1.
#' @return None
#' @author Kaur Alasoo
#' @export 
calculatePairwisePi1 <- function(qtl_list, qvalue_thresh = 0.1, tidy = FALSE){
  sample_names = names(qtl_list)
  rep_matrix = matrix(1,length(sample_names),length(sample_names))
  colnames(rep_matrix) = sample_names
  rownames(rep_matrix) = sample_names
  
  #Iterate through all pairs of p-values
  for (sn1 in 1:length(sample_names)){
    for (sn2 in 1:length(sample_names)){
      if (sn1 != sn2){
        rep_matrix[sn1, sn2] = calculatePi1(qtl_list[[sn1]], qtl_list[[sn2]])
      }
    }
  }
  
  #If tidy then return data frme insted of matrix
  if(tidy == TRUE){
    res = as.data.frame(rep_matrix) %>% 
      dplyr::mutate(first = rownames(rep_matrix)) %>% 
      dplyr::select(first, everything()) %>% 
      tidyr::gather("second","pi1",2:(ncol(rep_matrix)+1)) %>% 
      dplyr::arrange(first, second)
    return(res)
  }
  return(rep_matrix)
}

constructMatrixEQTLGenePos <- function(gene_metadata){
  res = dplyr::transmute(gene_metadata, geneid = gene_id, chr, left = start, right = end)
  return(res)
}

