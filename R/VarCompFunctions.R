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

estimateVarianceExplained <- function(model_data, model_function){
  #Apply an lmer4 model (model_function) to a dataset and report the proportion of variance explained by each component.
  #Set a flag to idicate if the model converged or not.
  tryCatch({
    var_exp = model_function(model_data) %>% varianceExplained()
    var_exp$converged = TRUE
    return(var_exp)
  },
  warning = function(c){
    var_exp = model_function(model_data) %>% varianceExplained() %>% suppressWarnings()
    var_exp$converged = FALSE
    return(var_exp)
  })
}

binGenesByResidual <- function(variance_table, n_bins = 20){
  #Bin genes by the proportion of variance explained by the residual
  var_table = dplyr::arrange(variance_table, residual) %>% 
    dplyr::mutate(residual_bin = n_bins-floor(residual*n_bins)) %>%
    dplyr::mutate(residual_bin = as.character(105-residual_bin*5)) %>%
    dplyr::mutate(residual_bin = factor(residual_bin, levels = rev(unique(residual_bin))))
  return(var_table)
}

meanVarianceWithinBins <- function(binned_table, binning_variable = "residual_bin"){
  #Caluclate mean variance within bins and add another bin for total variance
  factors = setdiff(colnames(binned_table), c("gene_id", binning_variable))
  var_gathered = tidyr::gather_(binned_table, "component", "var_explained", factors)
  
  #Calculate mean variance within bins
  var_summarised = group_by_(var_gathered, binning_variable, "component") %>% 
    dplyr::summarise(var_explained = mean(var_explained))
  
  #Calcualate total variance explaiend by each factor
  var_total = group_by_(var_gathered, "component") %>%
    dplyr::summarise(var_explained = mean(var_explained)) %>%
    dplyr::transmute(bin = "Total", component, var_explained)
  var_total = dplyr::rename_(var_total, .dots = setNames(list(quote(bin)), binning_variable))
            
  #Bind the two together
  var_summarised = rbind(var_summarised, var_total)
  return(var_summarised)
}

maximumFactorPerGene <- function(variance_table){
  #For each gene, identify the component that explains the most variance for that gene
  factors = setdiff(colnames(variance_table), "gene_id")
  
  maximum_factor = tidyr::gather_(variance_table, "component_max", "var_explained_max", factors) %>% 
    group_by(gene_id) %>% 
    arrange(-var_explained_max) %>% 
    filter(row_number() == 1)
  
  return(maximum_factor)
}


plotBinnedVariance <- function(var_summarised){
  #Make stacked barplot of the variance explained
  plot = ggplot(var_summarised, aes(x = residual_bin, y = var_explained, fill = component)) + 
    geom_bar(stat="identity") +
    ylab("% variance explained") +
    xlab("Residual variance bin")
  return(plot)
}