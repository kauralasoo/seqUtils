extractBetasFromList <- function(gene_snp_pairs, rasqual_result_list){
  betas_list = lapply(rasqual_result_list, function(x) dplyr::select(x, gene_id, snp_id, beta))
  res = gene_snp_pairs
  for (i in seq_along(names(rasqual_result_list))){
    res = dplyr::left_join(res, betas_list[[i]], by = c("gene_id", "snp_id"))
  }
  #Rename newly added columns with the names of the list
  colnames(res)[(ncol(gene_snp_pairs)+1):ncol(res)] = names(rasqual_result_list)
  return(res)
}

calculateBetaDiffMatrix <- function(beta_matrix, baseline_column = "naive"){
  beta_diff_matrix = dplyr::select(beta_matrix, -gene_id, -snp_id) %>% as.data.frame() #Remvoe gene id and SNP ids
  beta_diff_matrix = beta_diff_matrix - beta_diff_matrix[,baseline_column] #Calculate diff compared to baseline
  colnames(beta_diff_matrix) = paste(colnames(beta_diff_matrix), "diff", sep = "_") #Rename columns
  beta_diff_matrix$max_abs_diff = apply(abs(beta_diff_matrix), 1, max) #Max diff
  beta_diff_matrix = cbind(beta_matrix[,c("gene_id", "snp_id")], beta_diff_matrix)
  
  return(beta_diff_matrix)
}

clusterBetasKmeans <- function(beta_df, k){
  beta_matrix = dplyr::select(beta_df, -gene_id, -snp_id)
  beta_matrix = beta_matrix/apply(beta_matrix, 1, max) #Scale by maximum beta for each gene-snp pair
  clustering = kmeans(beta_matrix, k, nstart = 100)
  beta_df$cluster_id = clustering$cluster
  return(beta_df)
}

calculateClusterSizes <- function(clusters_df, selected_condition = "naive"){
  cluster_sizes = dplyr::select(clusters_df, cluster_id, condition_name) %>% 
    dplyr::filter(condition_name == selected_condition) %>% 
    dplyr::group_by(cluster_id) %>% 
    dplyr::summarise(count = length(cluster_id), condition_name = condition_name[1])
  return(cluster_sizes)
}

calculateClusterMeans <- function(clusters_df){
  cluster_means = dplyr::group_by(clusters_df, cluster_id, condition_name) %>% 
    dplyr::summarise(beta_mean = mean(beta), beta_sd = sd(beta))
  return(cluster_means)
}

#Flip the sign of the beta so that the maximal beta is positive
betaCorrectSign <- function(beta_df){
  #Calculate the sign for the maximal effect size
  max_sign = dplyr::group_by(beta_df, gene_id, snp_id) %>% 
    dplyr::arrange(-abs(beta)) %>% 
    dplyr::filter(row_number() == 1) %>%
    dplyr::mutate(max_sign = sign(beta)) %>% 
    dplyr::select(gene_id, snp_id, max_sign)
  
  #Flip the sign of all of the effect sizes
  beta_correct_sign = dplyr::left_join(beta_df, max_sign, by = c("gene_id", "snp_id")) %>%
    dplyr::arrange(gene_id, snp_id) %>% 
    dplyr::transmute(gene_id, snp_id, condition_name, beta = beta*max_sign)
  return(beta_correct_sign)
}

extractAndProcessBetas <- function(gene_snp_pairs, rasqual_results_list, baseline_column = "naive"){
  #Extract effect sizes for all gene-snp pairs from RASQUAL data
  beta_matrix = extractBetasFromList(gene_snp_pairs, rasqual_results_list) %>% dplyr::ungroup()
  
  #Convert Beta matrix into a df and revert signs to the sign of the maximal effect size
  beta_df = tidyr::gather(beta_matrix, condition_name, beta, naive:IFNg_SL1344)
  beta_correct_sign = betaCorrectSign(beta_df)
  beta_correct_sign_matrix = tidyr::spread(beta_correct_sign, condition_name, beta)
  
  #Calculate maximum absolute diff in betas between conditions
  beta_diff_matrix = calculateBetaDiffMatrix(beta_correct_sign_matrix, baseline_column) %>% dplyr::tbl_df()
  
  #Merge summary stats
  beta_summaries = dplyr::group_by(beta_correct_sign, gene_id, snp_id) %>% 
    dplyr::summarise(max_abs_beta = max(abs(beta)), min_abs_beta = min(abs(beta))) %>%
    dplyr::ungroup() %>%
    dplyr::left_join(beta_diff_matrix, by = c("gene_id","snp_id"))
  beta_summaries_matrix = dplyr::left_join(beta_correct_sign_matrix, beta_summaries, by = c("gene_id","snp_id"))
  
  return(list(beta_df = beta_correct_sign, beta_summaries = beta_summaries_matrix))
}

findMostAssociatedPeakPerSNP <- function(qtl_hits, rasqual_results){
  results = dplyr::semi_join(rasqual_results, qtl_hits, by = "snp_id") %>%
    dplyr::group_by(snp_id) %>% dplyr::arrange(desc(chisq)) %>%
    dplyr::mutate(n_tests = length(gene_id)) %>%
    dplyr::filter(row_number() == 1) %>% 
    dplyr::select(gene_id, snp_id, p_nominal, beta, n_tests) %>%
    dplyr::mutate(p_bonferroni = p.adjust(p_nominal, "bonferroni", n_tests)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(p_fdr = p.adjust(p_bonferroni, "fdr"))
  return(results)
}

mergeATACandRNAEffects <- function(atac_beta_df, rna_qtls, rna_beta_df){
  #Rename gene_id to peak_id
  atac_beta_df_renamed = dplyr::rename(atac_beta_df, peak_id = gene_id) %>% dplyr::mutate(phenotype = "ATAC")
  #Filter out RNA effect sizes
  rna_effects = dplyr::semi_join(rna_beta_df, rna_qtls, by = c("gene_id", "snp_id")) %>% 
    dplyr::semi_join(atac_beta_df_renamed, by = "snp_id") %>% 
    dplyr::mutate(phenotype = "RNA") %>% 
    dplyr::left_join(dplyr::select(atac_beta_df_renamed, peak_id, snp_id) %>% unique(), by = "snp_id") %>% 
    dplyr::select(gene_id, snp_id, peak_id, everything())
  #Filter out ATAC effect sizes
  atac_effects = dplyr::left_join(atac_beta_df_renamed, dplyr::select(rna_effects, gene_id, snp_id, peak_id) %>% unique(), 
                                  by = c("peak_id", "snp_id"))
  joint_effects = dplyr::bind_rows(rna_effects, atac_effects)
  return(joint_effects)
}

#' Add gene_id and gene_name to atac hits and rank them by "rank_by" column
rankAtacSummaries <- function(atac_beta_summaries, joint_effects, rank_by){
  atac_beta_summaries = dplyr::rename(atac_beta_summaries, peak_id = gene_id) #Rename gene_id to peak_id
  id_map = dplyr::select(joint_effects, peak_id, gene_id, gene_name, snp_id) %>% unique() #Select columns
  peak_gene_match = dplyr::left_join(atac_beta_summaries, id_map, by = c("peak_id", "snp_id")) %>% 
    dplyr::arrange_(rank_by)
  return(peak_gene_match)
}

prepareBetasDf <- function(rna_qtls_df, rna_betas, atac_rasqual_list, 
                           gene_name_map, appear_condition = "IFNg", rank_by = "IFNg_dff", baseline_column = "naive"){
  
  #Find msot associatied peak for each SNP
  matched_atac_peaks = findMostAssociatedPeakPerSNP(rna_qtls_df, atac_rasqual_list[[appear_condition]]) %>% dplyr::filter(p_fdr < 0.1)
  
  #Extract betas for ATAC peaks
  atac_betas = extractAndProcessBetas(dplyr::select(matched_atac_peaks, gene_id, snp_id), atac_rasqual_list, baseline_column)
  
  #Merge ATAC and RNA betas into a single data frame
  joint_betas =  mergeATACandRNAEffects(atac_betas$beta_df, rna_qtls_df, rna_betas$beta_df) %>%
    dplyr::left_join(gene_name_map, by = "gene_id")
  
  #Rank genes based on difference in ATAC effect size
  atac_ranked_summaries = rankAtacSummaries(atac_betas$beta_summaries, joint_betas, rank_by)
  joint_betas = dplyr::mutate(joint_betas, gene_name = factor(gene_name, levels = atac_ranked_summaries$gene_name))
  
  return(joint_betas)
}


