#' Convert FastQTL p-values into a data frame of z-scores
#'
#' @param pvalue_df data.fram of p-values (columns: p_nominal, beta)
#' @param n_top_snps Number of top associated SNPs to be included.
#'
#' @return data frame of z-scores
#' @export
pToZ <- function(pvalue_df, n_top_snps = 200){
  result = pvalue_df %>% 
    dplyr::arrange(p_nominal) %>%
    head(n_top_snps) %>% 
    dplyr::arrange(pos) %>%
    dplyr::mutate(z = qnorm(p_nominal)*sign(beta)) %>%
    dplyr::select(snp_id, z)
  return(result)
}

#' Calculate LD matrix based on a data frame of z-scores and genotypes
#'
#' @param z_df data frame of z-scores
#' @param genotypes matrix of genotypes
#'
#' @return matrix of pearson correlaton between genotypes
#' @export
zToLD <- function(z_df, genotypes){
  cor_mat = cor(t(genotypes[z_df$snp_id,]), method = "pearson", use = "pairwise.complete.obs")
  return(cor_mat)
}

saveFinemapMatrices = function(data_list, output_dir, file_suffix = "z"){
  #Save data for FastQTL to disk
  
  #Save each matrix as a separate txt file
  for (sn in names(data_list)){
    file_path = file.path(output_dir, paste(sn, file_suffix, sep = "."))
    print(file_path)
    write.table(data_list[[sn]], file_path, sep = " ", quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
}


finemapConstructMeta <- function(gene_ids, path, n_inds){
  meta = data_frame(z = paste0(path, gene_ids, ".z"),
                        ld = paste0(path, gene_ids, ".ld"),
                        snp = paste0(path, gene_ids, ".snp"),
                        config = paste0(path, gene_ids, ".config"),
                        k = paste0(path, "prior.k"),
                        log = paste0(path, gene_ids, ".log"), 
                        "n-ind" = n_inds)
  return(meta)
}

