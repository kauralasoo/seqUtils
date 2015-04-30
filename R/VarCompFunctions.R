varianceExplained <- function(lmer_model){
  #Calculate the percentage of variance explained by different factors
  variance = as.data.frame(VarCorr(lmer_model))
  var_percent = dplyr::mutate(variance, percent_variance = vcov/sum(vcov)) %>% 
    dplyr::select(grp, percent_variance) %>% 
    dplyr::mutate(type = "gene")
  var_row = tidyr::spread(var_percent, grp, percent_variance)
  return(var_row)  
}

constructGeneData <- function(gene_id, expression_matrix, design_matrix, metadata){
  #Construct df of gene expression for lmer analysis
  gene_exp = exprs_cqn[gene_id,]
  gene_df = data_frame(sample_id = names(gene_exp), exp_value = gene_exp)
  model_data = dplyr::left_join(gene_df, design_matrix, by = "sample_id") %>%
    dplyr::left_join(metadata, by = c("donor", "replicate")) %>%
    dplyr::mutate(ng_ul_standard = log(ng_ul_mean,2) - mean(log(ng_ul_mean,2))) %>% #Standardise RNA concentration measure
    dplyr::mutate(ng_ul_standard = ng_ul_standard/sd(ng_ul_standard))
  return(model_data)
}

