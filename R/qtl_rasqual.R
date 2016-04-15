#' Import rasqual output table into R
#'
#' Skipped gene-SNP pairs are automatically removed and 
#' chisq statistic is converted into p-value (p_nominal).
#'
#' @param path Bath to rasqual output file.
#'
#' @return data_frame
#' @export
importRasqualTable <- function(path){
  rasqual_results = readr::read_delim(path, delim = "\t", col_types = "ccciddddddddiii",col_names = FALSE)
  colnames(rasqual_results) = c("gene_id", "snp_id", "chr", "pos", "allele_freq", "HWE", "IA", "chisq", "effect_size", "delta", "phi", "overdisp", "n_feature_snps", "n_cis_snps", "converged")
  
  rasqual_pvalues = dplyr::filter(rasqual_results, snp_id != "SKIPPED") %>%
    dplyr::mutate(p_nominal = pchisq(chisq, df = 1, lower = FALSE))
  
  return(rasqual_pvalues)
}

#' Find SNP with minimal p-value per gene.
#'
#' Group QTL matrix by gene, sort by p_nominal, keep SNP with smalles p-value,
#' Correct that using bonferroni correction and then apply FDR correction across genes.
#'
#' @param qtl_df Data frame with QTL mapping results 
#' (required columns: gene_id, p_nominal, n_cis_snps)
#'
#' @return Only gene-SNP pairs with minimal p-values,
#'  added columns: p_bonferroni, p_fdr.
#' @export
findMinimalSnpPvalues <- function(qtl_df){
  result = dplyr::group_by(qtl_df, gene_id) %>%
    dplyr::arrange(p_nominal) %>% 
    dplyr::filter(row_number() == 1) %>%
    dplyr::mutate(p_bonferroni = p.adjust(p_nominal, "bonferroni", n_cis_snps)) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(p_nominal) %>%
    dplyr::mutate(p_fdr = p.adjust(p_bonferroni, "fdr"))
  return(result)
}

#' Write multiple files for RASQUAL onto disk.
#'
#' @param condition_list Named list of expression lists, each expression list needs to
#' contain at least the following elements: 'counts' matrix, 'sample_metadata' df,
#' @param rasqual_input_folder Path to the RASQUAL input folder.
#' @param max_batch_size Maximal number of feaures to be included in a single batch.
#'
#' @export
exportDataForRasqual <- function(condition_list, rasqual_input_folder, max_batch_size = 50){
  if(!is.list(condition_list)){
    stop("Input condiiton_list must be a list.")
  }
  if(length(condition_list) < 1){
    stop("The conditon_list must have at least one element.")
  }
  
  #Extract and save read count matrices
  counts_list = lapply(condition_list, function(x){x$counts})
  saveRasqualMatrices(counts_list, rasqual_input_folder, file_suffix = "expression")
  
  #Extract sample-genotype map for each condition
  sg_map = lapply(condition_list, function(x){ dplyr::select(x$sample_metadata, sample_id, genotype_id) })
  saveFastqtlMatrices(sg_map, rasqual_input_folder, file_suffix = "sg_map", col_names = FALSE)
  
  #Export library size
  library_size_list = lapply(counts_list, rasqualCalculateSampleOffsets, condition_list[[1]]$gene_metadata, gc_correct = FALSE)
  saveRasqualMatrices(library_size_list, rasqual_input_folder, file_suffix = "library_size")
  
  #Export GC-corrected library sizes
  gc_library_size_list = lapply(counts_list, rasqualCalculateSampleOffsets, condition_list[[1]]$gene_metadata)
  saveRasqualMatrices(gc_library_size_list, rasqual_input_folder, file_suffix = "gc_library_size")
  
  #Calculate covariates using Natsuhiko's SVD code
  covariates_list = lapply(condition_list, function(x){
    sf = rasqualCalculateSampleOffsets(x$counts, x$gene_metadata)
    covariates = rasqualMakeCovariates(x$counts, sf)
    return(covariates)
  })
  saveRasqualMatrices(covariates_list, rasqual_input_folder, file_suffix = "svd_covariates")
  
  #Extract covariates from sample metadata
  meta_cov_list = lapply(condition_list, function(x){
    meta_matrix = dplyr::select(x$sample_metadata, sample_id, sex_binary, PEER_factor_1:PEER_factor_10)
    cov_matrix = rasqualMetadataToCovariates(meta_matrix)[,1:5]
    return(cov_matrix)
  })
  saveRasqualMatrices(meta_cov_list, rasqual_input_folder, file_suffix = "PEER_covariates")
  
  #Extract covariates from sample metadata
  meta_cov_list = lapply(condition_list, function(x){
    meta_matrix = dplyr::select(x$sample_metadata, sample_id, sex_binary, PEER_factor_1:PEER_factor_10)
    cov_matrix = rasqualMetadataToCovariates(meta_matrix)[,1:3]
    return(cov_matrix)
  })
  saveRasqualMatrices(meta_cov_list, rasqual_input_folder, file_suffix = "PEER_covariates_n3")
  
  #Save feature names to disk
  feature_names = rownames(counts_list[[1]])
  write.table(feature_names, file.path(rasqual_input_folder, "feature_names.txt"), 
              row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  #Construct a list of batches for rasqual
  feature_batches = rasqualConstructGeneBatches(condition_list[[1]]$gene_metadata, max_batch_size)
  write.table(feature_batches, file.path(rasqual_input_folder, "feature_batches.txt"), 
              row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
}


