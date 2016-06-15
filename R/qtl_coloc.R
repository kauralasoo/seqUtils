#' Wrapper around coloc.abf to directly work with tidy data frames
#'
#' Group QTL database by gene, sort by p_nominal, keep SNP with smallest p-value,
#' correct that using bonferroni correction and then apply FDR correction across genes.
#'
#' @param df1 Summary statistics from the first dataset
#' (required columns: snp_id, p_nominal, beta, MAF).
#' @param df1 Summary statistics from the second dataset
#' (required columns: snp_id, p_nominal, beta, MAF).
#' @param n1 Sample size of the first dataset.
#' @param n2 Sample size of the second dataset.
#' @param p1 prior probability a SNP is associated with trait 1, default 1e-4
#' @param p2 prior probability a SNP is associated with trait 2, default 1e-4
#' @param p12 prior probability a SNP is associated with both traits, default 1e-5
#'
#' @return List of colocalisation results (Same as coloc.abf).
#' @export
testColoc <- function(df1, df2, n1, n2, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5){
  
  #Test for colocalisation between two sets of p-values
  result = coloc::coloc.abf(dataset1 = list(pvalues = df1$p_nominal, N = n1, 
                                            beta = df1$beta, MAF = df1$MAF, 
                                            snp = df1$snp_id, type = "quant"), 
                            dataset2 = list(pvalues = df2$p_nominal, N = n2, 
                                            beta = df2$beta, MAF = df2$MAF,
                                            snp = df2$snp_id, type = "quant"),
                            p1 = p1, p2 = p2, p12 = p12)
  return(result)
}

#' Calculate posterior probabilities from association summary stats
#'
#' @param dataset Association summary statistics
#' (required columns: snp_id, p_nominal, effect_size, MAF).
#' @param n Sample size of the dataset
#'
#' @return Original data frame with one additional column containing posterior probabilities.
#' @export
addAssociationPosterior <- function(dataset, n){
  
  #Use the coloc package to calculate ABFs
  coloc_res = testColoc(dataset, dataset, n, n)
  post_df = dplyr::transmute(coloc_res$results, snp_id = snp, lABF = lABF.df1) %>% 
    dplyr::mutate(posterior = exp(lABF)/sum(exp(lABF))) %>% 
    dplyr::select(snp_id, lABF, posterior) %>%
    dplyr::mutate(snp_id = as.character(snp_id))
  result = dplyr::left_join(dataset, post_df, by = "snp_id")
  return(result)
}


#' Construct the credible set of associated SNPs for each gene.
#'
#' @param dataset Output from addAssociationPosterior
#' @param threshold Cumulative posterior porbability of most associated SNPs.
#'
#' @return Dataset that only contains the credible set of associated SNPs.
#' @export
constructCredibleSet <- function(dataset, threshold = 0.95){
  result = dplyr::arrange(dataset, -posterior) %>% 
    dplyr::mutate(cum_posterior = cumsum(posterior))
  #Identify the number of top associated SNPs whose cumulative posterior is above the threshold
  n_top_rows = nrow(result) - (length(which(result$cum_posterior >= threshold))-1)
  result = dplyr::filter(result, row_number() <= n_top_rows)
  return(result)
}


