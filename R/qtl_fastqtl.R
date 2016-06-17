#' Save a list of matrices into a suitable format for FastQTL
#'
#' Works with expression and covariates matrices.
#' 
#' @param data_list list of matrices
#' @param output_dir relative path to the output dir
#' @param file_suffix suffix added to each file after their name in the list.
#' @return None
#' @author Kaur Alasoo
#' @export 
saveFastqtlMatrices <- function(data_list, output_dir, file_suffix = "expression", col_names = TRUE){
  #Save data for FastQTL to disk

  #Save each matrix as a separate  txt file
  for (sn in names(data_list)){
    file_path = file.path(output_dir, paste(sn,file_suffix, "txt", sep = "."))
    print(file_path)
    write.table(data_list[[sn]], file_path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = col_names)
  }
}

#' Import fastQTL output table into R.
#'
#' Detect if the table is from nominal run or permutation run and add proper column names.
#' 
#' @param file_path Path to the fastQTL output file
#' @return data_frame containing gene_ids, snp ids and p-values.
#' @author Kaur Alasoo
#' @export 
importFastQTLTable <- function(file_path){
  table = read.table(file_path, stringsAsFactors = FALSE)
  if(ncol(table) == 11){
    colnames(table) = c("gene_id", "n_cis_snps", "beta1", "beta2", "dummy", "snp_id", "distance","p_nominal", "slope","p_perm","p_beta")
    table = table %>% tbl_df() %>% 
      dplyr::filter(!is.na(p_beta)) %>%
      dplyr::mutate(p_fdr = p.adjust(p_beta, method = "fdr")) %>%
      dplyr::mutate(qvalue = qvalue::qvalue(p_beta)$qvalues) %>%
      dplyr::arrange(p_fdr)
  }
  return(table)
}

#' Construct gene position data.frame suitable for FastQTL
#' 
#' @param gene_metadata Gene metadata data frame from expression list (columns: chr, start, end, gene_id)
#' @return data_frame suitable for FastQTL
#' @author Kaur Alasoo
#' @export 
constructFastQTLGenePos <- function(gene_metadata){
	genepos = dplyr::select(gene_metadata, chr, start, end, gene_id) %>% 
	  dplyr::arrange(chr, start) %>%
    dplyr::rename_("left" = "start", "right" = "end", "geneid" = "gene_id", "#chr" = "chr")
  return(genepos)
}

#' Add gene position to the FastQTL expression matrix
#' 
#' @param matrix expression matrix
#' @param genepos data frame form constructFastQTLGenePos
#' @return data_frame suitable for FastQTL
#' @author Kaur Alasoo
#' @export 
prepareFastqtlMatrix <- function(matrix, genepos){
  res = dplyr::mutate(as.data.frame(matrix), geneid = rownames(matrix)) %>%
    dplyr::select(geneid, everything()) %>%
    dplyr::left_join(genepos, ., by = "geneid") %>%
    dplyr::arrange()
}

fastqtlMetadataToCovariates <- function(metadata){
  genotype_ids = metadata$genotype_id
  genotype_position = which(colnames(metadata) == "genotype_id")
  covariate_ids = colnames(metadata)[-genotype_position]
  cov_matrix = t(metadata[-genotype_position]) %>% as.data.frame()
  colnames(cov_matrix) = genotype_ids
  cov_matrix = dplyr::mutate(cov_matrix, id = covariate_ids) %>% dplyr::select(id, everything())
  return(cov_matrix)
}

exportDataForFastQTL <- function(condition_list, fastqtl_input_folder, n_chunks = 25){
  
  #Save expression datasets to disk
  fastqtl_genepos = constructFastQTLGenePos(condition_list[[1]]$gene_metadata)
  cqn_list = lapply(condition_list, function(x){x$cqn})
  fastql_cqn_list = lapply(cqn_list, prepareFastqtlMatrix, fastqtl_genepos)
  saveFastqtlMatrices(fastql_cqn_list, fastqtl_input_folder, file_suffix = "cqn_expression")
  
  #Extract covariates from sample metadata
  meta_cov_list = lapply(condition_list, function(x){
    meta_matrix = dplyr::select(x$sample_metadata, genotype_id, sex_binary, PEER_factor_1:PEER_factor_10)
    cov_matrix = fastqtlMetadataToCovariates(meta_matrix)[1:5,]
    return(cov_matrix)
  })
  saveFastqtlMatrices(meta_cov_list, fastqtl_input_folder, file_suffix = "covariates")
  
  #Construct chunks table
  chunks_matrix = data.frame(chunk = seq(1:n_chunks), n = n_chunks)
  write.table(chunks_matrix, file.path(fastqtl_input_folder, "chunk_table.txt"), 
              row.names = FALSE, quote = FALSE, col.names = FALSE, sep = " ")
}

fastQTLCorrectEigenMT <- function(fastqtl_results, n_tests){
  res = dplyr::left_join(fastqtl_results, n_tests, by = "gene_id") %>% 
    dplyr::mutate(n_tests = ifelse(is.na(n_tests), n_cis_snps, n_tests)) %>% 
    dplyr::group_by(gene_id) %>% 
    dplyr::mutate(p_eigen = p.adjust(p_nominal, "bonferroni", n = n_tests)) %>%
    dplyr::ungroup() %>% 
    dplyr::mutate(p_eigen_fdr = p.adjust(p_eigen, method = "fdr"))
  return(res)
}


#' Fetch particular genes from tabix indexed FastQTL output file.
#'
#' @param gene_ranges GRanges object with coordinates of the cis regions around genes.
#' @param tabix_file Tabix-indexed fastqtl output file.
#'
#' @return List of data frames containing Rasqual results for each gene.
#' @export
fastqtlTabixFetchGenes <- function(gene_ranges, tabix_file){
  #Set column names for rasqual
  fastqtl_columns = c("gene_id","chr","pos","snp_id","distance","p_nominal","beta")
  
  result = list()
  for (i in seq_along(gene_ranges)){
    selected_gene_id = gene_ranges[i]$gene_id
    print(i)
    tabix_table = scanTabixDataFrame(tabix_file, gene_ranges[i], col_names = fastqtl_columns)[[1]] %>%
      dplyr::filter(gene_id == selected_gene_id)
    
    #Add additional columns
    result[[selected_gene_id]] = tabix_table
  }
  return(result)
}