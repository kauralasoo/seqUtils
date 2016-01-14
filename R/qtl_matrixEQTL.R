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

constructMatrixEQTLGenePos <- function(gene_metadata){
  res = dplyr::transmute(gene_metadata, geneid = gene_id, chr, left = start, right = end)
  return(res)
}

#' Use MatrixEQTL to test interaction between genotype and two conditions
#'
#' @param conditions Character vector containing the names of the two conditions
#' @param gene_snp_pairs Data frame with two columns: gene_id, snp_id.
#' @param expression_list List with expression data.
#' @param genotype_list List with genotype data.
#'
#' @return Tidy table with interaction testing results
#' @export
matrixeqtlTestInteraction <- function(conditions, gene_snp_pairs, expression_list, genotype_list){
  #Extract gene and snp ids
  gene_ids = unique(gene_snp_pairs$gene_id)
  snp_ids = unique(gene_snp_pairs$snp_id)
  
  #Filter data
  design_matrix = dplyr::filter(expression_list$sample_metadata, condition_name %in% conditions)
  exp_matrix = expression_list$cqn[gene_ids, design_matrix$sample_id]
  geno_matrix = extractSubset(design_matrix, genotype_list$genotypes[snp_ids,], 
                              old_column_names = "genotype_id", new_column_names = "sample_id")
  snpspos = dplyr::filter(genotype_list$snpspos, snpid %in% snp_ids) %>% as.data.frame()
  genepos = constructMatrixEQTLGenePos(expression_list$gene_metadata) %>% dplyr::filter(geneid %in% gene_ids)
  
  #Construct condition covariate
  cov = dplyr::mutate(design_matrix, condition_cov = ifelse(condition_name == conditions[1], 0, 1))$condition_cov
  cov_matrix = t(as.matrix(cov))
  colnames(cov_matrix) = design_matrix$sample_id
  
  interaction_res = runMatrixEQTL(exp_matrix, geno_matrix, snpspos, genepos, cov_matrix, 
                                  pvOutputThreshold = 1, model = modelLINEAR_CROSS)
  interaction_table = interaction_res$cis$eqtls %>% dplyr::rename(gene_id = gene, snp_id = snps, p_nominal = pvalue) %>% 
    dplyr::mutate(gene_id = as.character(gene_id), snp_id = as.character(snp_id)) %>% 
    dplyr::semi_join(gene_snp_pairs, by = c("gene_id", "snp_id")) %>%
    dplyr::mutate(p_fdr = p.adjust(p_nominal, "fdr")) %>%
    dplyr::arrange(p_fdr) %>%
    dplyr::select(-FDR)
  return(interaction_table)
}

filterEQTLs <- function(data_frame, gene_id_name_map, fdr_cutoff = 0.1){
  dat = dplyr::filter(data_frame, FDR < fdr_cutoff) %>% 
    dplyr::rename(gene_id = gene, snp_id = snps) %>%
    dplyr::group_by(gene_id) %>% 
    dplyr::arrange(pvalue) %>% 
    dplyr::filter(row_number() == 1) %>%
    dplyr::ungroup() %>% 
    dplyr::arrange(pvalue) %>%
    dplyr::left_join(gene_id_name_map, by = "gene_id")
  return(dat)
}

