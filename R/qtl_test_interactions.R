#' Use MatrixEQTL to test interaction between genotype and two conditions
#'
#' @param conditions Character vector containing the names of the two conditions
#' @param gene_snp_pairs Data frame with two columns: gene_id, snp_id.
#' @param expression_list List with expression data.
#' @param genotype_list List with genotype data.
#' @param covariates_names Names of the columns in the expression_list$sample_metadata table 
#' that will be used as covariates for the matrixEQTL model.
#' 
#' @return Tidy table with interaction testing results
#' @export
matrixeqtlTestInteraction <- function(conditions, gene_snp_pairs, expression_list, genotype_list, covariate_names = NULL){
  #Extract gene and snp ids
  gene_ids = unique(gene_snp_pairs$gene_id)
  snp_ids = unique(gene_snp_pairs$snp_id)
  
  #Filter data
  design_matrix = dplyr::filter(expression_list$sample_metadata, condition_name %in% conditions)
  exp_matrix = expression_list$cqn[gene_ids, design_matrix$sample_id]
  geno_matrix = extractSubset(design_matrix, genotype_list$genotypes[snp_ids,], 
                              old_column_names = "genotype_id", new_column_names = "sample_id")
  snpspos = dplyr::filter(genotype_list$snpspos, snpid %in% snp_ids) %>% as.data.frame()
  genepos = constructMatrixEQTLGenePos(expression_list$gene_metadata) %>% 
    dplyr::filter(geneid %in% gene_ids) %>%
    as.data.frame()
  
  #Construct condition covariate
  design_matrix = dplyr::mutate(design_matrix, condition_cov = ifelse(condition_name == conditions[1], 0, 1))
  joint_covariate_names = c(covariate_names, "condition_cov")
  cov_matrix = constructMatrixEQTLCovariates(design_matrix, joint_covariate_names)
  
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



#' Test for interaction between genotype and condition using ANOVA model
#'
#' @param gene_id Tested gene id
#' @param snp_id Tested SNP id
#' @param trait_matrix expression or chromatin accessibility matrix
#' @param sample_metadata data frame with sample metadata and covariates
#' @param vcf_file VCF file from gdsToMatrix() function
#' @param qtl_formula Formula for the model with just genotype and condition terms
#' @param interaction_formula Formula for the model with interaction term betweeb genotype and condition
#'
#' @return Either a pvalue (if return_value == "ponly") or the full linear model object.
#' @export
testInteraction <- function(gene_id, snp_id, trait_matrix, sample_metadata, vcf_file, qtl_formula,
                            interaction_formula, return_value = "ponly"){
  
  #Some basic assertions
  assertthat::assert_that(is.matrix(trait_matrix))
  assertthat::assert_that(assertthat::has_name(sample_metadata, "sample_id"))
  assertthat::assert_that(assertthat::has_name(sample_metadata, "genotype_id"))
  
  #Extract data
  exp_data = data_frame(sample_id = colnames(trait_matrix), expression = trait_matrix[gene_id,])
  geno_data = data_frame(genotype_id = colnames(vcf_file$genotypes), genotype = vcf_file$genotypes[snp_id,])
  
  if(length(intersect(exp_data$sample_id, sample_metadata$sample_id)) != length(exp_data$sample_id)){
    stop("sample_id columns for trait_matrix and sample_metadata do not match.")
  } 

  sample_data = dplyr::left_join(sample_metadata, exp_data, by = "sample_id") %>%
    dplyr::left_join(geno_data, by = "genotype_id") %>%
    dplyr::filter(!is.na(genotype)) #Remove NA genotypes, because lm does not like them

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

#' Test for interaction between genotype and condition using ANOVA and lme4
#'
#' @param gene_id Tested gene id
#' @param snp_id Tested SNP id
#' @param trait_matrix expression or chromatin accessibility matrix
#' @param sample_metadata data frame with sample metadata and covariates
#' @param vcf_file VCF file from gdsToMatrix() function
#' @param qtl_formula Formula for the model with just genotype and condition terms
#' @param interaction_formula Formula for the model with interaction term betweeb genotype and condition
#'
#' @return Either a pvalue (if return_value == "ponly") or the full linear model object.
#' @export
testInteractionLme4 <- function(gene_id, snp_id, trait_matrix, sample_metadata, vcf_file, qtl_formula,
                            interaction_formula, return_value = "ponly"){
  
  #Some basic assertions
  assertthat::assert_that(is.matrix(trait_matrix))
  assertthat::assert_that(assertthat::has_name(sample_metadata, "sample_id"))
  assertthat::assert_that(assertthat::has_name(sample_metadata, "genotype_id"))
  
  #Extract data
  exp_data = data_frame(sample_id = colnames(trait_matrix), expression = trait_matrix[gene_id,])
  geno_data = data_frame(genotype_id = colnames(vcf_file$genotypes), genotype = vcf_file$genotypes[snp_id,])
  
  sample_data = dplyr::left_join(sample_metadata, exp_data, by = "sample_id") %>%
    dplyr::left_join(geno_data, by = "genotype_id")
  
  #apply two models to the data and compare them using anova
  no_interaction = lme4::lmer(qtl_formula, sample_data, REML = FALSE)
  interaction = lme4::lmer(interaction_formula, sample_data, REML = FALSE)
  res = anova(no_interaction, interaction)
  
  #Return value
  if(return_value == "ponly"){
    return(res[[8]])
  }
  else{
    return(list(anova = res, interaction_model = interaction))
  }
}


#Run testInteraction on a data.frame of gene-SNP pairs
testMultipleInteractions <- function(snps_df, trait_matrix, sample_metadata, vcf_file, qtl_formula, 
                                     interaction_formula, return_value = "ponly", id_field_separator = ";",
                                     lme4 = FALSE){
  
  #Filter the VCF file for quicker testing
  genotypes = vcf_file$genotypes[unique(snps_df$snp_id),]
  snps_pos = dplyr::filter(vcf_file$snpspos, snpid %in% rownames(genotypes))
  filtered_vcf = list(snpspos = snps_pos, genotypes = genotypes)
  
  #Plot eQTL results for a list of gene and SNP pairs.
  result = list()
  for(i in 1:nrow(snps_df)){
    gene_id = snps_df[i,]$gene_id
    snp_id = snps_df[i,]$snp_id
    print(gene_id)
    if (lme4 == TRUE){
      test = testInteractionLme4(gene_id, snp_id, trait_matrix, sample_metadata, filtered_vcf, qtl_formula, interaction_formula, return_value)
    } else{
      test = testInteraction(gene_id, snp_id, trait_matrix, sample_metadata, filtered_vcf, qtl_formula, interaction_formula, return_value)
    }
    result[[paste(gene_id,snp_id, sep = id_field_separator)]] = test
  }
  return(result)
}

#Post-process results from testMultipleInteractions
postProcessInteractionPvalues <- function(pvalue_list, id_field_separator = ";"){
  res = plyr::ldply(pvalue_list, .id = "id") %>% 
    tidyr::separate(id, into = c("gene_id", "snp_id"), sep = id_field_separator) %>%
    dplyr::rename(p_nominal = V2) %>% tbl_df() %>%
    dplyr::select(gene_id, snp_id, p_nominal) %>%
    dplyr::mutate(p_fdr = p.adjust(p_nominal,"fdr")) %>%
    dplyr::mutate(qvalue = qvalue::qvalue(p_nominal)$qvalues) %>%
    dplyr::arrange(p_nominal)
  return(res)
}

#' Test for threeway interactions between two features, genotype and environment.
#' 
#' The two features (master and dependent) are assumed to be regulated by the same 
#' causal genetic variant.
#'
#' @param feature_pairs data frame with master_id, dependent_id and snp_id
#' @param trait_matrix expression or chromatin accessibility matrix
#' @param sample_metadata data frame with sample metadata and covariates
#' @param vcf_file VCF file from gdsToMatrix() function
#' @param model0 Null model
#' @param model1 Model of interest that is compared to the null model
#' @param p_only Specifies return type
#'
#' @return If TRUE, then only p-value is returned, otherwise returns a list of model data and fits.
#' @export
testThreewayInteraction <- function(feature_pairs, trait_matrix, sample_metadata, vcf_file, model0, model1, p_only = TRUE){
  #Make sure that feature pairs is a data frame with only one row
  assertthat::assert_that(nrow(feature_pairs) == 1)
  assertthat::assert_that(assertthat::has_name(feature_pairs, "master_id"))
  assertthat::assert_that(assertthat::has_name(feature_pairs, "dependent_id"))
  assertthat::assert_that(assertthat::has_name(feature_pairs, "snp_id"))
  
  #Extract trait data
  master_cqn = data_frame(sample_id = colnames(trait_matrix), 
                          cqn = trait_matrix[feature_pairs$master_id,], 
                          peak_type = "master")
  dependent_cqn = data_frame(sample_id = colnames(trait_matrix), 
                             cqn = trait_matrix[feature_pairs$dependent_id,], 
                             peak_type = "dependent")
  cqn_data = rbind(master_cqn, dependent_cqn)
  
  #genotypes
  geno_data = data_frame(genotype_id = colnames(vcf_file$genotypes), 
                         genotype = vcf_file$genotypes[feature_pairs$snp_id,])
  
  #Everything
  joint_data = dplyr::left_join(cqn_data, sample_metadata, by = "sample_id") %>% 
    dplyr::left_join(geno_data, by = "genotype_id")
  
  #Compare models
  model0_fit = lm(model0, joint_data)
  model1_fit = lm(model1, joint_data)
  result = anova(model0_fit, model1_fit)
  
  #Return value
  if(p_only){
    return(result[[6]][2])
  }
  else{
    return(list(anova = result, interaction_model = interaction, data = joint_data))
  }
}

testInteractionByChromosome <- function(chr, joint_pairs, ...){
  
  #Import genotypes
  vcf_file = seqUtils::gdsToMatrix(paste0("processed/Fairfax/geno_by_chr/",chr,".gds"))
  
  #Filter unique snps per probe
  chr_pairs = dplyr::filter(joint_pairs, pheno_chr == chr) %>%
    dplyr::select(gene_id, snp_id)
  filtered_pairs = filterHitsR2(chr_pairs, vcf_file$genotypes, .8) %>% dplyr::tbl_df()
  
  #Test interactions
  interaction_results = testMultipleInteractions(snps_df = filtered_pairs,
                                                 vcf_file = vcf_file, ...)
  return(interaction_results)
}

