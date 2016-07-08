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


#' Filter SNPs based on R2
#'
#' Identify genes that have two different lead SNPs and remove one of them if 
#' the R2 between the two SNPs is greater than R2_thresh.
#'
#' @param feature_snp_pairs Data frame with two columns: gene_id, snp_id.
#' @param genotypes Genotype matrix (SNPs in rows, samples in columns)
#' @param R2_thresh R2 threshold, default 0.8.
#'
#' @return Filtered feature_snp_pairs that does not cotain highly correlated SNPs for single feature.
#' @export
filterHitsR2 <- function(feature_snp_pairs, genotypes, R2_thresh = 0.8){
  
  #Filter genotypes matrix to only include the relevant variants
  genotypes = genotypes[unique(feature_snp_pairs$snp_id),]
  
  #Count SNPs per gene
  snp_count = dplyr::group_by(feature_snp_pairs, gene_id) %>% dplyr::summarise(snp_count = length(snp_id))
  single_snps = dplyr::semi_join(feature_snp_pairs, dplyr::filter(snp_count, snp_count == 1), by = "gene_id")
  multi_snps = dplyr::anti_join(feature_snp_pairs, dplyr::filter(snp_count, snp_count == 1), by = "gene_id")
  
  #Filter genes with multiple genes
  multi_snp_list = plyr::dlply(multi_snps, "gene_id")
  multi_snp_filtered_list = lapply(multi_snp_list, filterGeneR2, genotypes, R2_thresh)
  multi_snp_filtered = plyr::ldply(multi_snp_filtered_list, .id = NULL)
  result = rbind(single_snps, multi_snp_filtered)
  return(result)
}

#Helper function for filterHitsR2
filterGeneR2 <- function(gene_df, genotypes, r2_thresh){
  genotype_matrix = t(genotypes[gene_df$snp_id,])
  r2 = cor(genotype_matrix, use = "pairwise.complete.obs")^2
  r2 = r2 - diag(1, nrow(r2))
  r2[lower.tri(r2)] = 0
  gene_df[colSums(r2 > r2_thresh) == 0,]
}

#' Calculate R2 between the lead SNP and all other SNPs in the data frame and adds corresponding column.
#'
#' Assumes that lead SNP is in the first row of the gene_df data frame 
#' (data frame is sorted by p-value).
#'
#' @param gene_df Data frame of QTL snps (required column: gene_id)
#' @param genotypes Genotype matrix
#'
#' @return gene_df data frame with R2 column added.
#' @export
addR2FromLead <- function(gene_df, genotypes){
  
  assertthat::assert_that(!is.null(gene_df))
  assertthat::assert_that(assertthat::has_name(gene_df, "snp_id"))
  
  #Extract genotype matrix
  genotype_matrix = t(genotypes[gene_df$snp_id,])
  
  #If more than one SNP then calculate R2
  if( nrow(gene_df) > 1 ){ 
    #Calculate R2 between the first SNP in the genotype matrix and all other SNPs
    r2 = apply(genotype_matrix, 2, function(x, y){ cor(x,y,use = "pairwise.complete.obs") }, genotype_matrix[,1])^2
    gene_df = dplyr::mutate(gene_df, R2 = r2)
  } else { #If only one SNP then set R2 to 1
    gene_df = dplyr::mutate(gene_df, R2 = 1)
  }
  return(gene_df)
}

calculateR2FromLead <- function(snp_ids, genotypes){
  
  #Extract genotype matrix
  genotype_matrix = t(genotypes[snp_ids,])
  
  #If more than one SNP then calculate R2
  if( length(snp_ids) > 1 ){ 
    #Calculate R2 between the first SNP in the genotype matrix and all other SNPs
    r2 = apply(genotype_matrix, 2, function(x, y){ cor(x,y,use = "pairwise.complete.obs") }, genotype_matrix[,1])^2
  } else { #If only one SNP then set R2 to 1
    r2 = 1
  }
  return(r2)
}

#' Add expected p-value for Q-Q plots
#'
#' @param pvalue_df 
#'
#' @return pvalue_df with p_expected column
#' @export
addExpectedPvalue <- function(pvalue_df){
  dplyr::arrange(pvalue_df, p_eigen) %>%
    dplyr::mutate(p_expected = c(1:length(p_eigen))/length(p_eigen))
}

#' Extract all QTLs at a specific FDR level from a list of min pvalues by condition
#'
#' Multiple variants per gene are sorted by p-value
#'
#' @param min_pvalue_list List of QTLs per condition
#' @param fdr_cutoff 
#'
#' @return Data frame of QTLs
#' @export
extractQTLsFromList <- function(min_pvalue_list, fdr_cutoff = 0.1){
  min_hits = purrr::map(min_pvalue_list, ~dplyr::filter(.,p_fdr < fdr_cutoff))
  qtl_df = purrr::map_df(min_hits, identity, .id = "condition_name") %>% 
    dplyr::group_by(gene_id) %>% 
    dplyr::arrange(p_nominal) %>%
    dplyr::ungroup()
  return(qtl_df)
}
