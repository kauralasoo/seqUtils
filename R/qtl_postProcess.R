
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

#' Calculate emiprical permutation p-values
#'
#' For each element in the te vector of original p-values, count the number of permutation p-values that are smaller than that.
#' 
#' @param raw_pvalues Vector of raw p-values -log10 form the original data.
#' @param output_dir Vector of maximum -log10 pvalues from the permuted data.
#' @return None
#' @author Kaur Alasoo
#' @export 
calculatePermutationPvalues <- function(raw_pvalues, max_permutation_pvalues){
  
  #Count the number of permutations
  n_perm = length(max_permutation_pvalues)
  
  #Count the number of pvalues
  n_pvalues = length(raw_pvalues)
  
  #Construct matrix of permutation pvalues
  perm_mat = matrix(rep(sort(max_permutation_pvalues), n_pvalues), nrow = n_pvalues, byrow = TRUE)
  comparison = perm_mat - raw_pvalues
  comparison[comparison > 0] = 1
  comparison[comparison <= 0] = 0
  permutation_pvalues = (rowSums(comparison) + 1) / (n_perm + 1)
  
  return(permutation_pvalues)
}
