#' Calculate the proportion of variance explaned by different factors in a lme4 model
varianceExplained <- function(lmer_model){
  variance = as.data.frame(lme4::VarCorr(lmer_model))
  var_percent = dplyr::mutate(variance, percent_variance = vcov/sum(vcov)) %>% 
    dplyr::select(grp, percent_variance) %>% 
    dplyr::mutate(type = "gene")
  var_row = tidyr::spread(var_percent, grp, percent_variance)
  return(var_row)  
}

#' Estimate propotion of variance explained by each compoent in a model.
#'
#' #Apply an lme4 model (model_function) to a data frame and 
#' report the proportion of variance explained by each component.
#' 
#' @param model_data data frame with all of the neccessary columns for the model.
#' @param model_function Function that when applied to a model_data data frame returns an
#' lme4 object.
#' @return data frame containing the proportion of variance explained by each component.
#' @author Kaur Alasoo
#' @export 
estimateVarianceExplained <- function(model_data, model_function){
 
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