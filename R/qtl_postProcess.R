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

#' Filter QTL interaction test results based on RASQUAL effect size
#'
#' @param conditions Character vector of condition names.
#' @param interaction_result Output from the matrixeqtlTestInteraction function.
#' @param pvalue_list List of RASQUAL gene-SNP results by condition.
#' @param fdr_thresh Interaction test FDR threshold. 
#' 
#' @return Data frame augmented with Rasqual effect sizes and p-values in each condition.
#' @export
filterInteractionResults <- function(conditions, interaction_result, pvalue_list, fdr_thresh = 0.1){
  
  #Filter interaction tests by p-value
  interaction_candidates = dplyr::filter(interaction_result, p_fdr <= fdr_thresh)
  
  #Extract RASQUAL results for significant gene-SNP pairs
  res_A = dplyr::semi_join(pvalue_list[[ conditions[1] ]], interaction_candidates, by = c("gene_id", "snp_id")) %>%
    dplyr::select(gene_id, snp_id, beta, p_nominal) %>%
    dplyr::rename(beta_A = beta, p_nominal_A = p_nominal)
  res_B = dplyr::semi_join(pvalue_list[[ conditions[2] ]], interaction_candidates, by = c("gene_id", "snp_id")) %>%
    dplyr::select(gene_id, snp_id, beta, p_nominal) %>%
    dplyr::rename(beta_B = beta, p_nominal_B = p_nominal)
  
  #Join all of the data together and filter
  res = dplyr::left_join(interaction_candidates, res_A, by = c("gene_id", "snp_id")) %>% 
    tbl_df() %>% 
    dplyr::left_join(res_B, by =c("gene_id","snp_id")) %>% 
    dplyr::mutate(beta_diff = beta_B - beta_A) %>%
    dplyr::mutate(abs_beta_min = pmin(abs(beta_A), abs(beta_B))) %>%
    dplyr::group_by(gene_id, snp_id) %>%
    dplyr::mutate(min_condition = which.min(c(abs(beta_A), abs(beta_B)))) %>%
    dplyr::mutate(min_condition_name = ifelse(min_condition == 1, conditions[1], conditions[2]))
  return(res)
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

testInterctionsBetweenPairs <- function(condition_pair, rasqual_min_hits, combined_expression_data, covariate_names, vcf_file, fdr_thresh = 0.1){
  #Extraxt gene name map
  gene_name_map = dplyr::select(combined_expression_data$gene_metadata, gene_id, gene_name)
  
  #Find independent pairs of SNPs
  min_pvalue_df = ldply(rasqual_min_hits[condition_pair], .id = "condition_name")
  joint_pairs = dplyr::select(min_pvalue_df, gene_id, snp_id) %>% unique()
  joint_pairs_filtered = filterHitsR2(joint_pairs, vcf_file$genotypes, 0.2)

  #Test pairwise interactions using matrixeQTL
  interaction_df = matrixeqtlTestInteraction(condition_pair, joint_pairs_filtered, combined_expression_data, vcf_file, covariate_names)
  interaction_effects = filterInteractionResults(condition_pair, interaction_df, rasqual_selected_pvalues, fdr_thresh = fdr_thresh) %>%
    dplyr::left_join(gene_name_map, by = "gene_id")
  return(interaction_effects)
}

#' Test for interaction between genotype and condition using ANOVA model
#'
#' @param gene_id Tested gene id
#' @param snp_id Tested SNP id
#' @param eqtl_data_list Exression dataset
#' @param vcf_file VCF file from gdsToMatrix() function
#' @param qtl_formula Formula for the model with just genotype and condition terms
#' @param interaction_formula Formula for the model with interaction term betweeb genotype and condition
#'
#' @return Either a pvalue (if return_value == "ponly") or the full linear model object.
#' @export
testInteraction <- function(gene_id, snp_id, trait_matrix, sample_metadata, vcf_file, qtl_formula, interaction_formula, return_value = "ponly"){

  #Extract data
  exp_data = data_frame(sample_id = colnames(trait_matrix), expression = trait_matrix[gene_id,])
  geno_data = data_frame(genotype_id = colnames(vcf_file$genotypes), genotype = vcf_file$genotypes[snp_id,])
  
  sample_data = dplyr::left_join(sample_metadata, exp_data, by = "sample_id") %>%
    dplyr::left_join(geno_data, by = "genotype_id")
  
  #apply two models to the data and compare them using anova
  no_interaction = lm(qtl_formula, as.data.frame(sample_data))
  interaction = lm(interaction_formula, as.data.frame(sample_data))
  res = anova(no_interaction, interaction)
  
  #Return value
  if(return_value == "ponly"){
    return(res[[6]])
  }
  else{
    return(list(anova = res, interaction_model = interaction))
  }
}

#Run testInteraction on a data.frame of gene-SNP pairs
testMultipleInteractions <- function(snps_df, trait_matrix, sample_metadata, vcf_file, qtl_formula, interaction_formula, return_value = "ponly"){
  #Plot eQTL results for a list of gene and SNP pairs.
  result = list()
  for(i in 1:nrow(snps_df)){
    gene_id = snps_df[i,]$gene_id
    snp_id = snps_df[i,]$snp_id
    print(gene_id)
    test = testInteraction(gene_id, snp_id, trait_matrix, sample_metadata, vcf_file, qtl_formula, interaction_formula, return_value)
    result[[paste(gene_id,snp_id, sep = ":")]] = test
  }
  return(result)
}

#Post-process results from testMultipleInteractions
postProcessInteractionPvalues <- function(pvalue_list){
  res = plyr::ldply(pvalue_list, .id = "id") %>% 
    tidyr::separate(id, into = c("gene_id", "snp_id"), sep = ":") %>%
    dplyr::rename(p_nominal = V2) %>% tbl_df() %>%
    dplyr::select(gene_id, snp_id, p_nominal) %>%
    dplyr::mutate(p_fdr = qvalue::qvalue(p_nominal)$qvalues) %>%
    dplyr::arrange(p_nominal)
  return(res)
}

