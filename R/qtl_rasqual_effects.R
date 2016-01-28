extractBetasFromList <- function(gene_snp_pairs, rasqual_result_list){
  betas_list = lapply(rasqual_result_list, function(x) dplyr::select(x, gene_id, snp_id, beta))
  res = gene_snp_pairs
  for (i in seq_along(names(rasqual_result_list))){
    res = dplyr::left_join(res, betas_list[[i]], by = c("gene_id", "snp_id"))
  }
  colnames(res)[3:ncol(res)] = names(rasqual_result_list)
  return(res)
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